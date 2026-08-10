#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Interpolation des NoData et lissage final."""
import os
from multiprocessing import Pool, cpu_count

import numpy as np
import rasterio
from loguru import logger
from scipy.ndimage import uniform_filter, label, binary_dilation
from scipy.spatial import cKDTree
from tqdm import tqdm

from mqga.tiling import init_worker


def _reapply_protect_mask(chem_raster, protect_mask_path, no_data=-9999):
	"""Remet nodata sur les pixels du masque protégé (ex. zones de décrochage)."""
	if not protect_mask_path:
		return
	from mqga.decrochage import apply_mask_as_nodata
	apply_mask_as_nodata(
		chem_raster, protect_mask_path, chem_out=chem_raster, no_data=no_data
	)
	logger.info("Masque protégé réappliqué (nodata): {}", protect_mask_path)


def _load_protect_bool(protect_mask_path, shape):
	"""True = pixel protégé (décrochage) : ne pas interpoler. False partout si pas de masque."""
	if not protect_mask_path:
		return np.zeros(shape, dtype=bool)
	with rasterio.open(protect_mask_path, "r") as src:
		mask = src.read(1)
	if mask.shape != shape:
		raise ValueError(
			f"Masque protégé {mask.shape} incompatible avec raster {shape}"
		)
	return mask != 0


def _process_block_idw(args):
	"""
	Fonction helper pour le traitement parallèle des blocs IDW.
	Doit être au niveau du module pour être picklable par multiprocessing.
	"""
	from scipy.spatial import cKDTree
	
	(block_y, block_x, chem_in_local, height_local, width_local, 
	 block_size, search_radius, no_data, power, protect_mask_path) = args
	
	# Ouvrir le fichier dans le worker
	with rasterio.open(chem_in_local, 'r') as src_local:
		# Limites du bloc avec padding
		row_start = max(0, block_y * block_size - search_radius)
		row_end = min(height_local, (block_y + 1) * block_size + search_radius)
		col_start = max(0, block_x * block_size - search_radius)
		col_end = min(width_local, (block_x + 1) * block_size + search_radius)
		
		# Zone de traitement
		process_row_start = block_y * block_size
		process_row_end = min(height_local, (block_y + 1) * block_size)
		process_col_start = block_x * block_size
		process_col_end = min(width_local, (block_x + 1) * block_size)
		
		# Lire le bloc
		window = rasterio.windows.Window(col_start, row_start, 
										  col_end - col_start, 
										  row_end - row_start)
		block_data = src_local.read(1, window=window).astype(np.float32)
		block_protected = np.zeros(block_data.shape, dtype=bool)
		if protect_mask_path:
			with rasterio.open(protect_mask_path, 'r') as src_prot:
				block_protected = src_prot.read(1, window=window) != 0
	
	# Masques
	mask_valid = (block_data != no_data) & ~np.isnan(block_data)
	mask_nodata = ~mask_valid
	
	# Offsets
	process_row_off = process_row_start - row_start
	process_col_off = process_col_start - col_start
	process_height = process_row_end - process_row_start
	process_width = process_col_end - process_col_start
	
	# Masque zone de traitement
	process_mask = np.zeros_like(mask_nodata, dtype=bool)
	process_mask[process_row_off:process_row_off + process_height,
				 process_col_off:process_col_off + process_width] = True
	
	# Ne pas interpoler les pixels protégés (décrochage) : resteront nodata
	process_nodata = mask_nodata & process_mask & ~block_protected
	
	# Créer le bloc de résultat
	result_block = block_data[process_row_off:process_row_off + process_height,
							  process_col_off:process_col_off + process_width].copy()
	
	if np.any(process_nodata) and np.any(mask_valid):
		# Coordonnées
		block_rows, block_cols = np.meshgrid(
			np.arange(block_data.shape[0]), 
			np.arange(block_data.shape[1]), 
			indexing='ij'
		)
		
		# Points valides et leurs valeurs
		points_valid = np.column_stack([block_rows[mask_valid], block_cols[mask_valid]])
		values_valid = block_data[mask_valid]
		
		# Points nodata à interpoler
		points_nodata = np.column_stack([block_rows[process_nodata], block_cols[process_nodata]])
		
		# Construire un arbre KD pour recherche rapide des voisins
		tree = cKDTree(points_valid)
		
		# Trouver tous les voisins dans le rayon (vectorisé)
		distances_list, indices_list = tree.query(points_nodata, k=min(10, len(points_valid)), 
												  distance_upper_bound=search_radius)
		
		# Traiter chaque point nodata (vectorisé par batch)
		for i, (distances, indices) in enumerate(zip(distances_list, indices_list)):
			# Filtrer les distances infinies
			valid_mask = np.isfinite(distances) & (distances > 0)
			
			if np.any(valid_mask):
				valid_distances = distances[valid_mask]
				valid_indices = indices[valid_mask]
				valid_values = values_valid[valid_indices]
				
				# Éviter division par zéro
				valid_distances = np.maximum(valid_distances, 0.1)
				
				# Poids = 1 / distance^power (vectorisé)
				weights = 1.0 / (valid_distances ** power)
				
				# IDW (vectorisé)
				interp_value = np.sum(weights * valid_values) / np.sum(weights)
				
				# Mettre à jour le résultat
				row_idx, col_idx = points_nodata[i]
				rel_row = int(row_idx - process_row_off)
				rel_col = int(col_idx - process_col_off)
				if 0 <= rel_row < result_block.shape[0] and 0 <= rel_col < result_block.shape[1]:
					result_block[rel_row, rel_col] = interp_value
	
	return (process_col_start, process_row_start, process_width, process_height, result_block)

#############################################################################################################################	
def interpolate_nodata_idw_vectorized(
	chem_in, chem_out, no_data=-9999, search_radius=50, power=2, block_size=2000, n_jobs=None,
	protect_mask_path=None,
):
	"""
	Version optimisée et vectorisée de l'interpolation IDW.
	Beaucoup plus rapide que la version originale grâce à la vectorisation numpy.
	
	Args:
		chem_in: Chemin vers l'image d'entrée
		chem_out: Chemin vers l'image de sortie
		no_data: Valeur nodata
		search_radius: Rayon de recherche pour les pixels valides (par défaut 50)
		power: Puissance pour la pondération (par défaut 2)
		block_size: Taille des blocs pour le traitement (augmenté à 2000 pour meilleure performance)
		n_jobs: Nombre de processus parallèles (None = auto)
		protect_mask_path: masque (non nul) exclus du fill et remis en nodata après
	"""
	if n_jobs is None:
		n_jobs = max(1, cpu_count() - 1)
	
	# Supprimer le fichier de sortie
	if os.path.exists(chem_out):
		try:
			os.remove(chem_out)
		except OSError:
			pass
	
	chem_out_aux = chem_out + '.aux.xml'
	if os.path.exists(chem_out_aux):
		try:
			os.remove(chem_out_aux)
		except OSError:
			pass
	
	# Ouvrir l'image source
	with rasterio.open(chem_in, 'r') as src:
		metadata = src.meta.copy()
		height, width = src.height, src.width
		
		# Créer l'image de sortie
		with rasterio.open(chem_out, 'w', **metadata) as dst:
			# Traiter par blocs
			n_blocks_y = (height + block_size - 1) // block_size
			n_blocks_x = (width + block_size - 1) // block_size
			total_blocks = n_blocks_y * n_blocks_x
			
			logger.info("Traitement par blocs (IDW vectorisé): {}x{} blocs de {}x{} pixels ({} processus)", n_blocks_x, n_blocks_y, block_size, block_size, n_jobs)
			
			# Traiter les blocs en parallèle
			# Passer tous les paramètres nécessaires à la fonction globale
			block_args = [(block_y, block_x, chem_in, height, width, block_size, search_radius, no_data, power, protect_mask_path) 
						  for block_y in range(n_blocks_y) for block_x in range(n_blocks_x)]
			
			total_nodata_interpolated = 0
			with Pool(processes=n_jobs, initializer=init_worker) as pool:
				results = list(tqdm(pool.imap(_process_block_idw, block_args), total=total_blocks, desc="Interpolation IDW"))
			
			# Écrire les résultats
			for col_start, row_start, width, height, result_block in results:
				write_window = rasterio.windows.Window(col_start, row_start, width, height)
				dst.write(result_block.astype(metadata['dtype']), 1, window=write_window)
				total_nodata_interpolated += np.sum((result_block != no_data) & ~np.isnan(result_block))
	
	logger.info("Interpolation terminée (IDW vectorisé): {} pixels nodata interpolés.", total_nodata_interpolated)
	_reapply_protect_mask(chem_out, protect_mask_path, no_data=no_data)

def interpolate_nodata_hybrid(
	chem_in, chem_out, no_data=-9999,
	connectivity=4, seuil_percent=50,
	poids=1, rayon=50, n=1, block_size=2000,
	protect_mask_path=None,
	vcalc_mode="p90",
):
	"""
	Interpole les pixels nodata avec la méthode hybride de xingng.
	Combine interpolation locale sur les pixels de bord et constante statistique.
	
	Args:
		chem_in: Chemin vers l'image d'entrée
		chem_out: Chemin vers l'image de sortie
		no_data: Valeur nodata
		connectivity: Connexité (4 ou 8, par défaut 4)
		seuil_percent: Pourcentage de valeurs minimales à exclure si vcalc_mode="min" (défaut 50)
		poids: Puissance pour la pondération IDW (défaut 1 = linéaire)
		rayon: Rayon de recherche pour l'interpolation locale (défaut 50)
		n: Facteur de pondération entre interpolation et constante (défaut 1)
		block_size: Taille des blocs pour le traitement (défaut 2000, non utilisé actuellement)
		protect_mask_path: masque (non nul) exclus du fill et remis en nodata après
		vcalc_mode: Constante de trou — "p90" = P90(ε_bord) (défaut), "min" = historique
	"""
	# Supprimer le fichier de sortie
	if os.path.exists(chem_out):
		try:
			os.remove(chem_out)
		except OSError:
			pass
	
	chem_out_aux = chem_out + '.aux.xml'
	if os.path.exists(chem_out_aux):
		try:
			os.remove(chem_out_aux)
		except OSError:
			pass
	
	# Définir la structure de connexité
	if connectivity == 4:
		structure = np.array([[0, 1, 0],
							  [1, 1, 1],
							  [0, 1, 0]], dtype=bool)
	else:  # connexité 8
		structure = np.ones((3, 3), dtype=bool)
	
	# Ouvrir l'image source
	with rasterio.open(chem_in, 'r') as src:
		metadata = src.meta.copy()
		height, width = src.height, src.width
		
		# Lire toute l'image (nécessaire pour identifier les trous connexes)
		data = src.read(1).astype(np.float32)
	
	# Identifier les pixels nodata et valides
	mask_nodata = (data == no_data) | np.isnan(data)
	mask_valid = ~mask_nodata
	# Décrochage / masque protégé : rester nodata, ne pas entrer dans les trous à boucher
	mask_protected = _load_protect_bool(protect_mask_path, data.shape)
	mask_to_fill = mask_nodata & ~mask_protected
	n_nodata = int(np.sum(mask_nodata))
	n_protected = int(np.sum(mask_protected))
	n_protected_skip = int(np.sum(mask_nodata & mask_protected))
	n_to_fill = int(np.sum(mask_to_fill))
	logger.info(
		"Interpolation hybride: protect_mask={}, nodata={}, protégés={}, "
		"exclus du fill={}, à interpoler={}",
		protect_mask_path or "(aucun)",
		n_nodata,
		n_protected,
		n_protected_skip,
		n_to_fill,
	)
	
	if not np.any(mask_to_fill):
		with rasterio.open(chem_out, 'w', **metadata) as dst:
			dst.write(data, 1)
		logger.info("Aucun pixel nodata à interpoler (hors zones protégées).")
		_reapply_protect_mask(chem_out, protect_mask_path, no_data=no_data)
		return
	
	# Identifier les trous connexes (composantes connexes de nodata non protégés)
	labeled_holes, num_holes = label(mask_to_fill, structure=structure)
	logger.info("Interpolation hybride: {} trou(s) connexe(s) à boucher.", num_holes)
	
	# print(f"Traitement de {num_holes} trou(s) connexe(s)...")  # Désactivé pour ne garder que la barre de progression
	
	# Créer une copie pour le résultat
	result = data.copy()
	
	# Traiter chaque trou
	for hole_id in tqdm(range(1, num_holes + 1), desc="Bouchage des trous"):
		# Masque du trou actuel
		mask_hole = (labeled_holes == hole_id)
		
		# Pixels de bord de ce trou spécifique
		# Un pixel de bord est un pixel valide adjacent à ce trou
		hole_dilated = binary_dilation(mask_hole, structure=structure)
		border_mask = hole_dilated & mask_valid & ~mask_hole
		
		# Collecter les valeurs des pixels de bord
		border_values = data[border_mask].tolist()
		
		if len(border_values) == 0:
			# Pas de pixels de bord, on ne peut pas boucher ce trou
			# print(f"  Attention: Trou {hole_id} n'a pas de pixels de bord, ignoré.")  # Désactivé pour ne garder que la barre de progression
			continue
		
		# Calculer V_calc pour ce trou (ancre du centre)
		# p90 (défaut): P90(ε_bord) — moins optimiste qu'un min, moins sensible qu'un max
		# min: historique — exclure seuil% des plus faibles, puis minimum des restantes
		if vcalc_mode == "min":
			values_sorted = sorted(border_values)
			n_exclude = int(len(values_sorted) * seuil_percent / 100)
			if n_exclude >= len(values_sorted):
				n_exclude = len(values_sorted) - 1
			if n_exclude < 0:
				n_exclude = 0
			values_filtered = values_sorted[n_exclude:]
			V_calc = min(values_filtered) if len(values_filtered) > 0 else values_sorted[-1]
		else:
			V_calc = float(np.percentile(np.asarray(border_values, dtype=np.float32), 90))
		
		# Coordonnées des pixels de bord pour ce trou
		border_coords = np.column_stack(np.where(border_mask))
		if len(border_coords) == 0:
			continue
		
		border_coords_float = border_coords.astype(np.float32)
		border_values_array = np.array(border_values)
		
		# Construire un arbre KD pour recherche rapide des distances
		tree = cKDTree(border_coords_float)
		
		# Coordonnées des pixels nodata de ce trou
		hole_coords = np.column_stack(np.where(mask_hole))
		hole_coords_float = hole_coords.astype(np.float32)
		
		# Traiter par batch pour réduire la mémoire
		batch_size = 5000
		for i in range(0, len(hole_coords), batch_size):
			batch_coords = hole_coords_float[i:i+batch_size]
			
			# Calculer les distances au bord le plus proche pour chaque pixel nodata
			distances_to_border, _ = tree.query(batch_coords, k=1)
			
			# Calculer K pour chaque pixel : K = min(1, d/rayon)
			K = np.minimum(1.0, distances_to_border / rayon)
			
			# Initialiser V_interpole avec V_calc (fallback)
			V_interpole = np.full(len(batch_coords), V_calc, dtype=np.float32)
			
			# Trouver les pixels dans le rayon (d <= rayon)
			mask_in_radius = distances_to_border <= rayon
			
			if np.any(mask_in_radius):
				# Pour ces pixels, calculer l'interpolation IDW sur les pixels de bord
				coords_in_radius = batch_coords[mask_in_radius]
				
				for j, coord in enumerate(coords_in_radius):
					# Trouver les pixels de bord dans le rayon
					distances_to_border_points, indices = tree.query(
						coord.reshape(1, -1), 
						k=min(len(border_coords), 20),
						distance_upper_bound=rayon
					)
					
					# Filtrer les distances infinies et zéro
					valid_mask = np.isfinite(distances_to_border_points[0]) & (distances_to_border_points[0] > 0)
					
					if np.any(valid_mask):
						valid_distances = distances_to_border_points[0][valid_mask]
						valid_indices = indices[0][valid_mask]
						valid_values = border_values_array[valid_indices]
						
						# Éviter division par zéro
						valid_distances = np.maximum(valid_distances, 0.1)
						
						# Poids = 1 / distance^poids (IDW)
						weights = 1.0 / (valid_distances ** poids)
						
						# Interpolation IDW
						V_interpole[np.where(mask_in_radius)[0][j]] = np.sum(weights * valid_values) / np.sum(weights)
			
			# Calculer V final pour ce batch
			# V = (1 - K^n) * V_interpole + K^n * V_calc
			if n == 1:
				V = (1 - K) * V_interpole + K * V_calc
			elif n > 0:
				K_power = K ** n
				V = (1 - K_power) * V_interpole + K_power * V_calc
			elif n == -1:
				# Cas spécial n=-1 équivaut à n=1
				V = (1 - K) * V_interpole + K * V_calc
			else:  # n < 0 et n != -1
				K_power = K ** abs(n)
				denom = ((1 - K) ** abs(n)) + K_power
				V = ((1 - K) * V_interpole + K_power * V_calc) / denom
			
			# Mettre à jour le résultat
			for k, coord in enumerate(batch_coords):
				row, col = int(coord[0]), int(coord[1])
				if 0 <= row < height and 0 <= col < width:
					result[row, col] = V[k]
		
		# print(f"  Trou {hole_id}: {np.sum(mask_hole)} pixels bouchés, V_calc={V_calc:.2f}")  # Désactivé pour ne garder que la barre de progression
	
	# Sauvegarder le résultat
	metadata['dtype'] = result.dtype
	with rasterio.open(chem_out, 'w', **metadata) as dst:
		dst.write(result, 1)
	
	nodata_filled = int(np.sum(mask_to_fill & (result != no_data)))
	logger.info("Interpolation terminée (hybride): {} pixels nodata interpolés.", nodata_filled)
	_reapply_protect_mask(chem_out, protect_mask_path, no_data=no_data)

def apply_moving_average(chem_in, chem_out, window_size=50, no_data=-9999, protect_mask_path=None):
	"""
	Applique une moyenne sur une fenêtre glissante à l'image.
	
	Args:
		chem_in: Chemin vers l'image d'entrée
		chem_out: Chemin vers l'image de sortie
		window_size: Taille de la fenêtre glissante (par défaut 50x50)
		no_data: Valeur nodata (par défaut -9999)
	"""
	# Supprimer le fichier de sortie s'il existe déjà pour garantir l'écrasement
	if os.path.exists(chem_out):
		try:
			os.remove(chem_out)
		except OSError:
			pass
	
	# Supprimer aussi les fichiers auxiliaires (.aux.xml) s'ils existent
	chem_out_aux = chem_out + '.aux.xml'
	if os.path.exists(chem_out_aux):
		try:
			os.remove(chem_out_aux)
		except OSError:
			pass
	
	# Lire l'image
	with rasterio.open(chem_in, 'r') as src:
		data = src.read(1).astype(np.float32)
		metadata = src.meta.copy()
	
	# Créer un masque pour les pixels nodata
	mask_nodata = (data == no_data) | np.isnan(data)
	
	# Créer une copie des données en float64 pour les calculs
	data_float = data.astype(np.float64)
	
	# Remplacer les nodata par 0 pour le calcul de la somme
	data_sum = data_float.copy()
	data_sum[mask_nodata] = 0.0
	
	# Créer un masque de poids (1 pour valide, 0 pour nodata)
	weights = (~mask_nodata).astype(np.float64)
	
	# Calculer la somme et le nombre de pixels valides dans chaque fenêtre
	sum_window = uniform_filter(data_sum, size=window_size, mode='constant', cval=0.0)
	count_window = uniform_filter(weights, size=window_size, mode='constant', cval=0.0)
	
	# Calculer la moyenne uniquement là où il y a des pixels valides
	# Éviter la division par zéro et supprimer les warnings
	with np.errstate(divide='ignore', invalid='ignore'):
		result = np.where(count_window > 0, sum_window / count_window, no_data)
	
	# Conserver les nodata originaux (si un pixel était nodata, il reste nodata)
	result[mask_nodata] = no_data
	
	# Convertir en float32 pour la sauvegarde
	result = result.astype(np.float32)
	
	# Sauvegarder l'image résultante
	metadata['dtype'] = result.dtype
	with rasterio.open(chem_out, 'w', **metadata) as dst:
		dst.write(result, 1)

	_reapply_protect_mask(chem_out, protect_mask_path, no_data=no_data)

