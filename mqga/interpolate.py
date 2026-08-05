#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Interpolation des NoData et lissage final."""
import os
from multiprocessing import Pool, cpu_count

import numpy as np
import rasterio
from scipy.ndimage import generic_filter, uniform_filter, label, binary_dilation
from scipy.interpolate import LinearNDInterpolator, griddata
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from tqdm import tqdm
from loguru import logger

from mqga.tiling import init_worker


def interpolate_nodata_with_linearnd(chem_in, chem_out, no_data=-9999, block_size=1000):
	"""
	Interpole les pixels nodata (valeur -9999) en utilisant LinearNDInterpolator de scipy.
	Version optimisée qui traite l'image par blocs pour réduire la consommation mémoire.
	
	Args:
		chem_in: Chemin vers l'image d'entrée avec des pixels nodata
		chem_out: Chemin vers l'image de sortie avec les pixels nodata interpolés
		no_data: Valeur nodata (par défaut -9999)
		block_size: Taille des blocs pour le traitement (par défaut 1000x1000)
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
	
	# Ouvrir l'image source
	with rasterio.open(chem_in, 'r') as src:
		metadata = src.meta.copy()
		height, width = src.height, src.width
		
		# Créer l'image de sortie
		with rasterio.open(chem_out, 'w', **metadata) as dst:
			# Traiter l'image par blocs pour réduire la consommation mémoire
			n_blocks_y = (height + block_size - 1) // block_size
			n_blocks_x = (width + block_size - 1) // block_size
			total_blocks = n_blocks_y * n_blocks_x
			
			# print(f"Traitement par blocs: {n_blocks_x}x{n_blocks_y} blocs de {block_size}x{block_size} pixels")
			
			block_count = 0
			total_nodata_interpolated = 0
			
			for block_y in range(n_blocks_y):
				for block_x in range(n_blocks_x):
					block_count += 1
					
					# Calculer les limites du bloc avec padding pour avoir des pixels valides autour
					row_start = max(0, block_y * block_size - block_size // 4)
					row_end = min(height, (block_y + 1) * block_size + block_size // 4)
					col_start = max(0, block_x * block_size - block_size // 4)
					col_end = min(width, (block_x + 1) * block_size + block_size // 4)
					
					# Zone de traitement (sans le padding)
					process_row_start = block_y * block_size
					process_row_end = min(height, (block_y + 1) * block_size)
					process_col_start = block_x * block_size
					process_col_end = min(width, (block_x + 1) * block_size)
					
					# Lire le bloc avec padding
					window = rasterio.windows.Window(col_start, row_start, 
													  col_end - col_start, 
													  row_end - row_start)
					block_data = src.read(1, window=window).astype(np.float32)
					
					# Identifier les pixels valides et nodata dans le bloc
					mask_valid = (block_data != no_data) & ~np.isnan(block_data)
					mask_nodata = ~mask_valid
					
					# Calculer les offsets dans le bloc
					process_row_off = process_row_start - row_start
					process_col_off = process_col_start - col_start
					process_height = process_row_end - process_row_start
					process_width = process_col_end - process_col_start
					
					# Masque pour la zone de traitement (sans le padding)
					process_mask = np.zeros_like(mask_nodata, dtype=bool)
					process_mask[process_row_off:process_row_off + process_height,
								 process_col_off:process_col_off + process_width] = True
					
					# Pixels nodata uniquement dans la zone de traitement
					process_nodata = mask_nodata & process_mask
					
					if not np.any(process_nodata):
						# Pas de nodata à traiter dans ce bloc, copier directement
						result_block = block_data[process_row_off:process_row_off + process_height,
												  process_col_off:process_col_off + process_width].copy()
					else:
						# Créer les coordonnées relatives dans le bloc
						block_rows, block_cols = np.meshgrid(
							np.arange(block_data.shape[0]), 
							np.arange(block_data.shape[1]), 
							indexing='ij'
						)
						
						# Coordonnées des pixels valides dans le bloc
						points_valid = np.column_stack([block_rows[mask_valid], block_cols[mask_valid]])
						values_valid = block_data[mask_valid]
						
						# Coordonnées des pixels nodata à interpoler (seulement dans la zone de traitement)
						points_nodata = np.column_stack([block_rows[process_nodata], block_cols[process_nodata]])
						
						if len(points_valid) < 3:
							# Pas assez de points valides pour interpoler, garder nodata
							result_block = block_data[process_row_off:process_row_off + (process_row_end - process_row_start),
													  process_col_off:process_col_off + (process_col_end - process_col_start)].copy()
						else:
							# Créer l'interpolateur pour ce bloc
							interpolator = LinearNDInterpolator(points_valid, values_valid)
							
							# Interpoler les valeurs pour les pixels nodata
							interpolated_values = interpolator(points_nodata)
							
							# Créer une copie du bloc de résultat (zone de traitement uniquement)
							result_block = block_data[process_row_off:process_row_off + process_height,
													  process_col_off:process_col_off + process_width].copy()
							
							# Créer un masque pour les pixels nodata dans la zone de traitement
							process_nodata_local = process_nodata[process_row_off:process_row_off + process_height,
																   process_col_off:process_col_off + process_width]
							
							# Remplacer les pixels nodata par les valeurs interpolées
							result_block[process_nodata_local] = interpolated_values
							
							# Gérer les NaN (hors du domaine convexe)
							mask_nan = np.isnan(interpolated_values)
							if np.any(mask_nan):
								# Remplacer les NaN par nodata
								result_block[process_nodata_local][mask_nan] = no_data
							
							total_nodata_interpolated += np.sum(~np.isnan(interpolated_values))
					
					# Écrire le bloc de résultat
					write_window = rasterio.windows.Window(process_col_start, process_row_start,
															process_col_end - process_col_start,
															process_row_end - process_row_start)
					dst.write(result_block.astype(metadata['dtype']), 1, window=write_window)
					
					# Afficher la progression
					if block_count % 10 == 0 or block_count == total_blocks:
						progress = (block_count / total_blocks) * 100
						print(f"  Progression: {block_count}/{total_blocks} blocs ({progress:.1f}%) - {total_nodata_interpolated} pixels interpolés", end='\r')
			
			print()  # Nouvelle ligne après la progression
	
	logger.info("Interpolation terminée (LinearNDInterpolator): {} pixels nodata interpolés.", total_nodata_interpolated)

#############################################################################################################################	
def interpolate_nodata_griddata(chem_in, chem_out, no_data=-9999, block_size=1000):
	"""
	Interpole les pixels nodata avec griddata (méthode 'linear').
	Plus rapide que LinearNDInterpolator.
	
	Args:
		chem_in: Chemin vers l'image d'entrée avec des pixels nodata
		chem_out: Chemin vers l'image de sortie avec les pixels nodata interpolés
		no_data: Valeur nodata (par défaut -9999)
		block_size: Taille des blocs pour le traitement (par défaut 1000x1000)
	"""
	# Supprimer le fichier de sortie s'il existe déjà
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
			# Traiter l'image par blocs
			n_blocks_y = (height + block_size - 1) // block_size
			n_blocks_x = (width + block_size - 1) // block_size
			total_blocks = n_blocks_y * n_blocks_x
			
			logger.info("Traitement par blocs (griddata): {}x{} blocs de {}x{} pixels", n_blocks_x, n_blocks_y, block_size, block_size)
			
			block_count = 0
			total_nodata_interpolated = 0
			
			for block_y in range(n_blocks_y):
				for block_x in range(n_blocks_x):
					block_count += 1
					
					# Calculer les limites du bloc avec padding
					row_start = max(0, block_y * block_size - block_size // 4)
					row_end = min(height, (block_y + 1) * block_size + block_size // 4)
					col_start = max(0, block_x * block_size - block_size // 4)
					col_end = min(width, (block_x + 1) * block_size + block_size // 4)
					
					# Zone de traitement (sans le padding)
					process_row_start = block_y * block_size
					process_row_end = min(height, (block_y + 1) * block_size)
					process_col_start = block_x * block_size
					process_col_end = min(width, (block_x + 1) * block_size)
					
					# Lire le bloc avec padding
					window = rasterio.windows.Window(col_start, row_start, 
													  col_end - col_start, 
													  row_end - row_start)
					block_data = src.read(1, window=window).astype(np.float32)
					
					# Identifier les pixels valides et nodata
					mask_valid = (block_data != no_data) & ~np.isnan(block_data)
					mask_nodata = ~mask_valid
					
					# Calculer les offsets
					process_row_off = process_row_start - row_start
					process_col_off = process_col_start - col_start
					process_height = process_row_end - process_row_start
					process_width = process_col_end - process_col_start
					
					# Masque pour la zone de traitement
					process_mask = np.zeros_like(mask_nodata, dtype=bool)
					process_mask[process_row_off:process_row_off + process_height,
								 process_col_off:process_col_off + process_width] = True
					
					process_nodata = mask_nodata & process_mask
					
					if not np.any(process_nodata):
						# Pas de nodata, copier directement
						result_block = block_data[process_row_off:process_row_off + process_height,
												  process_col_off:process_col_off + process_width].copy()
					else:
						# Créer les coordonnées
						block_rows, block_cols = np.meshgrid(
							np.arange(block_data.shape[0]), 
							np.arange(block_data.shape[1]), 
							indexing='ij'
						)
						
						# Coordonnées des pixels valides
						points_valid = np.column_stack([block_rows[mask_valid], block_cols[mask_valid]])
						values_valid = block_data[mask_valid]
						
						# Coordonnées des pixels nodata
						points_nodata = np.column_stack([block_rows[process_nodata], block_cols[process_nodata]])
						
						if len(points_valid) < 3:
							# Pas assez de points valides
							result_block = block_data[process_row_off:process_row_off + process_height,
													  process_col_off:process_col_off + process_width].copy()
						else:
							# Utiliser griddata avec méthode 'linear'
							interpolated_values = griddata(
								points_valid, values_valid, points_nodata,
								method='linear', fill_value=np.nan
							)
							
							# Créer le bloc de résultat
							result_block = block_data[process_row_off:process_row_off + process_height,
													  process_col_off:process_col_off + process_width].copy()
							
							# Masque local pour les nodata
							process_nodata_local = process_nodata[process_row_off:process_row_off + process_height,
																   process_col_off:process_col_off + process_width]
							
							# Remplacer les pixels nodata
							result_block[process_nodata_local] = interpolated_values
							
							# Gérer les NaN
							mask_nan = np.isnan(interpolated_values)
							if np.any(mask_nan):
								result_block[process_nodata_local][mask_nan] = no_data
							
							total_nodata_interpolated += np.sum(~np.isnan(interpolated_values))
					
					# Écrire le bloc
					write_window = rasterio.windows.Window(process_col_start, process_row_start,
															process_col_end - process_col_start,
															process_row_end - process_row_start)
					dst.write(result_block.astype(metadata['dtype']), 1, window=write_window)
					
					# Progression
					if block_count % 10 == 0 or block_count == total_blocks:
						progress = (block_count / total_blocks) * 100
						print(f"  Progression: {block_count}/{total_blocks} blocs ({progress:.1f}%) - {total_nodata_interpolated} pixels interpolés", end='\r')
			
			print()
	
	logger.info("Interpolation terminée (griddata): {} pixels nodata interpolés.", total_nodata_interpolated)

#############################################################################################################################	
def interpolate_nodata_idw(chem_in, chem_out, no_data=-9999, search_radius=50, power=2, block_size=1000):
	"""
	Interpole les pixels nodata avec IDW (Inverse Distance Weighting).
	Rapide et efficace.
	
	Args:
		chem_in: Chemin vers l'image d'entrée
		chem_out: Chemin vers l'image de sortie
		no_data: Valeur nodata
		search_radius: Rayon de recherche pour les pixels valides (par défaut 50)
		power: Puissance pour la pondération (par défaut 2)
		block_size: Taille des blocs pour le traitement
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
			
			logger.info("Traitement par blocs (IDW): {}x{} blocs de {}x{} pixels", n_blocks_x, n_blocks_y, block_size, block_size)
			
			block_count = 0
			total_nodata_interpolated = 0
			
			for block_y in range(n_blocks_y):
				for block_x in range(n_blocks_x):
					block_count += 1
					
					# Limites du bloc avec padding
					row_start = max(0, block_y * block_size - block_size // 4)
					row_end = min(height, (block_y + 1) * block_size + block_size // 4)
					col_start = max(0, block_x * block_size - block_size // 4)
					col_end = min(width, (block_x + 1) * block_size + block_size // 4)
					
					# Zone de traitement
					process_row_start = block_y * block_size
					process_row_end = min(height, (block_y + 1) * block_size)
					process_col_start = block_x * block_size
					process_col_end = min(width, (block_x + 1) * block_size)
					
					# Lire le bloc
					window = rasterio.windows.Window(col_start, row_start, 
													  col_end - col_start, 
													  row_end - row_start)
					block_data = src.read(1, window=window).astype(np.float32)
					
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
					
					process_nodata = mask_nodata & process_mask
					
					# Créer le bloc de résultat
					result_block = block_data[process_row_off:process_row_off + process_height,
											  process_col_off:process_col_off + process_width].copy()
					
					if np.any(process_nodata):
						# Coordonnées
						block_rows, block_cols = np.meshgrid(
							np.arange(block_data.shape[0]), 
							np.arange(block_data.shape[1]), 
							indexing='ij'
						)
						
						points_valid = np.column_stack([block_rows[mask_valid], block_cols[mask_valid]])
						values_valid = block_data[mask_valid]
						points_nodata = np.column_stack([block_rows[process_nodata], block_cols[process_nodata]])
						
						if len(points_valid) > 0:
							# Traiter par batch pour réduire la mémoire
							batch_size = 5000
							for i in range(0, len(points_nodata), batch_size):
								batch_points = points_nodata[i:i+batch_size]
								
								# Calculer les distances
								distances = cdist(batch_points, points_valid)
								
								# Trouver les voisins dans le rayon
								mask_near = distances <= search_radius
								
								# IDW pour chaque pixel
								for j in range(len(batch_points)):
									near_mask = mask_near[j]
									if np.any(near_mask):
										near_distances = distances[j, near_mask]
										near_values = values_valid[near_mask]
										
										# Éviter division par zéro
										near_distances = np.maximum(near_distances, 0.1)
										
										# Poids = 1 / distance^power
										weights = 1.0 / (near_distances ** power)
										interp_value = np.sum(weights * near_values) / np.sum(weights)
										
										# Mettre à jour le résultat
										row_idx, col_idx = batch_points[j]
										rel_row = int(row_idx - process_row_off)
										rel_col = int(col_idx - process_col_off)
										if 0 <= rel_row < result_block.shape[0] and 0 <= rel_col < result_block.shape[1]:
											result_block[rel_row, rel_col] = interp_value
											total_nodata_interpolated += 1
					
					# Écrire le bloc
					write_window = rasterio.windows.Window(process_col_start, process_row_start,
															process_col_end - process_col_start,
															process_row_end - process_row_start)
					dst.write(result_block.astype(metadata['dtype']), 1, window=write_window)
					
					# Progression
					if block_count % 10 == 0 or block_count == total_blocks:
						progress = (block_count / total_blocks) * 100
						print(f"  Progression: {block_count}/{total_blocks} blocs ({progress:.1f}%) - {total_nodata_interpolated} pixels interpolés", end='\r')
			
			print()
	
	logger.info("Interpolation terminée (IDW): {} pixels nodata interpolés.", total_nodata_interpolated)

#############################################################################################################################	
def _process_block_idw(args):
	"""
	Fonction helper pour le traitement parallèle des blocs IDW.
	Doit être au niveau du module pour être picklable par multiprocessing.
	"""
	from scipy.spatial import cKDTree
	
	(block_y, block_x, chem_in_local, height_local, width_local, 
	 block_size, search_radius, no_data, power) = args
	
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
	
	process_nodata = mask_nodata & process_mask
	
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
def interpolate_nodata_idw_vectorized(chem_in, chem_out, no_data=-9999, search_radius=50, power=2, block_size=2000, n_jobs=None):
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
			block_args = [(block_y, block_x, chem_in, height, width, block_size, search_radius, no_data, power) 
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

#############################################################################################################################	
def interpolate_nodata_fast(chem_in, chem_out, no_data=-9999, max_iterations=5):
	"""
	Interpolation rapide utilisant scipy.ndimage pour remplir les nodata.
	Très rapide mais moins précise que les méthodes d'interpolation spatiale.
	Utilise une approche de propagation itérative.
	
	Args:
		chem_in: Chemin vers l'image d'entrée
		chem_out: Chemin vers l'image de sortie
		no_data: Valeur nodata
		max_iterations: Nombre maximum d'itérations (par défaut 5)
	"""
	from scipy.ndimage import binary_dilation, uniform_filter
	
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
	
	# Lire l'image
	with rasterio.open(chem_in, 'r') as src:
		data = src.read(1).astype(np.float32)
		metadata = src.meta.copy()
	
	mask_nodata = (data == no_data) | np.isnan(data)
	
	if not np.any(mask_nodata):
		with rasterio.open(chem_out, 'w', **metadata) as dst:
			dst.write(data, 1)
		logger.info("Aucun pixel nodata à interpoler.")
		return
	
	# Copier les données
	result = data.copy()
	
	# Itérer pour remplir progressivement les nodata
	for iteration in range(max_iterations):
		# Masque des nodata restants
		current_nodata = (result == no_data) | np.isnan(result)
		
		if not np.any(current_nodata):
			break
		
		# Appliquer un filtre uniforme (moyenne) pour propager les valeurs
		# Utiliser un masque pour ne traiter que les zones nodata
		filtered = uniform_filter(result, size=3, mode='constant', cval=no_data)
		
		# Remplacer seulement les nodata par les valeurs filtrées
		result[current_nodata] = filtered[current_nodata]
		
		# Remettre no_data où il n'y a toujours pas de valeur valide
		still_nodata = (result == no_data) | np.isnan(result)
		result[still_nodata] = no_data
	
	# Sauvegarder
	metadata['dtype'] = result.dtype
	with rasterio.open(chem_out, 'w', **metadata) as dst:
		dst.write(result, 1)
	
	nodata_filled = np.sum(mask_nodata & (result != no_data))
	logger.info("Interpolation terminée (rapide): {} pixels nodata interpolés en {} itérations.", nodata_filled, iteration+1)

#############################################################################################################################	
def interpolate_nodata_window(chem_in, chem_out, no_data=-9999, window_size=20):
	"""
	Interpole les nodata avec une moyenne pondérée dans une fenêtre glissante.
	Très rapide mais moins précise.
	
	Args:
		chem_in: Chemin vers l'image d'entrée
		chem_out: Chemin vers l'image de sortie
		no_data: Valeur nodata
		window_size: Taille de la fenêtre de recherche (par défaut 20)
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
	
	# Lire l'image
	with rasterio.open(chem_in, 'r') as src:
		data = src.read(1).astype(np.float32)
		metadata = src.meta.copy()
	
	mask_nodata = (data == no_data) | np.isnan(data)
	
	if not np.any(mask_nodata):
		with rasterio.open(chem_out, 'w', **metadata) as dst:
			dst.write(data, 1)
		logger.info("Aucun pixel nodata à interpoler.")
		return
	
	# Fonction pour remplir les nodata dans une fenêtre
	def fill_nodata(window):
		center = window[window_size//2, window_size//2]
		if center != no_data and not np.isnan(center):
			return center
		
		valid = (window != no_data) & ~np.isnan(window)
		if np.any(valid):
			return np.mean(window[valid])
		return no_data
	
	logger.info("Interpolation par fenêtre glissante ({}x{})...", window_size, window_size)
	result = generic_filter(data, fill_nodata, size=window_size, mode='constant', cval=no_data)
	
	# Garder les nodata qui n'ont pas pu être interpolés
	result[mask_nodata & (result == no_data)] = no_data
	
	# Sauvegarder
	metadata['dtype'] = result.dtype
	with rasterio.open(chem_out, 'w', **metadata) as dst:
		dst.write(result, 1)
	
	nodata_filled = np.sum((mask_nodata) & (result != no_data))
	logger.info("Interpolation terminée (fenêtre): {} pixels nodata interpolés.", nodata_filled)

#############################################################################################################################	
def interpolate_nodata_hybrid(chem_in, chem_out, no_data=-9999, 
                               connectivity=4, seuil_percent=50, 
                               poids=1, rayon=50, n=1, block_size=2000):
	"""
	Interpole les pixels nodata avec la méthode hybride de xingng.
	Combine interpolation locale sur les pixels de bord et constante statistique.
	Équivalent à xingng -FB:2:C:50,1:1:50:1
	
	Args:
		chem_in: Chemin vers l'image d'entrée
		chem_out: Chemin vers l'image de sortie
		no_data: Valeur nodata
		connectivity: Connexité (4 ou 8, par défaut 4)
		seuil_percent: Pourcentage de valeurs minimales à exclure pour V_calc (défaut 50)
		poids: Puissance pour la pondération IDW (défaut 1 = linéaire)
		rayon: Rayon de recherche pour l'interpolation locale (défaut 50)
		n: Facteur de pondération entre interpolation et constante (défaut 1)
		block_size: Taille des blocs pour le traitement (défaut 2000, non utilisé actuellement)
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
	
	if not np.any(mask_nodata):
		with rasterio.open(chem_out, 'w', **metadata) as dst:
			dst.write(data, 1)
		logger.info("Aucun pixel nodata à interpoler.")
		return
	
	# Identifier les trous connexes (composantes connexes de nodata)
	labeled_holes, num_holes = label(mask_nodata, structure=structure)
	
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
		
		# Calculer V_calc pour ce trou
		# Exclure seuil% des plus faibles valeurs, puis prendre le minimum
		values_sorted = sorted(border_values)
		n_exclude = int(len(values_sorted) * seuil_percent / 100)
		if n_exclude >= len(values_sorted):
			n_exclude = len(values_sorted) - 1
		if n_exclude < 0:
			n_exclude = 0
		values_filtered = values_sorted[n_exclude:]
		V_calc = min(values_filtered) if len(values_filtered) > 0 else values_sorted[-1]
		
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
	
	nodata_filled = np.sum(mask_nodata & (result != no_data))
	logger.info("Interpolation terminée (hybride): {} pixels nodata interpolés.", nodata_filled)

#############################################################################################################################	
def apply_moving_average(chem_in, chem_out, window_size=50, no_data=-9999):
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

