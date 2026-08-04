#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Assemblage des tuiles et fondu des recouvrements."""
import os

import numpy as np
import rasterio
from rasterio.merge import merge
from rasterio.windows import from_bounds
from tqdm import tqdm


def calculate_overlap_bounds(src1_path, src2_path):
	"""Calcule les bounds de recouvrement entre deux dalles (intersection)."""
	with rasterio.open(src1_path) as src1, rasterio.open(src2_path) as src2:
		bounds1 = src1.bounds
		bounds2 = src2.bounds
		
		# Intersection des bounds
		left = max(bounds1.left, bounds2.left)
		right = min(bounds1.right, bounds2.right)
		bottom = max(bounds1.bottom, bounds2.bottom)
		top = min(bounds1.top, bounds2.top)
		
		# Vérifier qu'il y a bien un recouvrement
		if left >= right or bottom >= top:
			return None
			
		return (left, bottom, right, top)  # X_0, Y_0, X_1, Y_1

#################################################################################################### 
def create_weight_image_horizontal(overlap_bounds, target_shape):
	"""Crée une image de poids basée sur C/NC (colonne/nombre de colonnes).
	Retourne un array numpy avec des valeurs entre 0 (gauche) et 1 (droite).
	
	Args:
		overlap_bounds: (left, bottom, right, top) - bounds de recouvrement
		target_shape: (height, width) - forme de l'array cible
	"""
	height, width = target_shape
	
	# Créer un array de poids basé sur la position en colonne
	# C/NC : numéro de colonne / nombre de colonnes
	col_indices = np.arange(width, dtype=np.float32)
	weight_array = col_indices / (width - 1) if width > 1 else np.ones(width)
	
	# Étendre sur toutes les lignes
	weights = np.tile(weight_array, (height, 1))
	
	return weights

#################################################################################################### 
def create_weight_image_vertical(overlap_bounds, target_shape):
	"""Crée une image de poids basée sur L/NL (ligne/nombre de lignes).
	Retourne un array numpy avec des valeurs entre 0 (haut) et 1 (bas).
	
	Args:
		overlap_bounds: (left, bottom, right, top) - bounds de recouvrement
		target_shape: (height, width) - forme de l'array cible
	"""
	height, width = target_shape
	
	# Créer un array de poids basé sur la position en ligne
	# L/NL : numéro de ligne / nombre de lignes
	row_indices = np.arange(height, dtype=np.float32)
	weight_array = row_indices / (height - 1) if height > 1 else np.ones(height)
	
	# Étendre sur toutes les colonnes
	weights = np.tile(weight_array.reshape(-1, 1), (1, width))
	
	return weights

#################################################################################################### 
def weighted_blend_overlap(dalle1_path, dalle2_path, weight1, weight2, overlap_bounds, output_path):
	"""Fait la moyenne pondérée de deux dalles dans la zone de recouvrement.
	Formule: I1*I2 + I3*I4 où I1=dalle1, I2=poids1, I3=dalle2, I4=poids2
	L'image de sortie a les bounds définis par overlap_bounds (équivalent à -cg:)."""
	with rasterio.open(dalle1_path) as src1, rasterio.open(dalle2_path) as src2:
		# Calculer les fenêtres de recouvrement pour chaque dalle
		window1 = src1.window(*overlap_bounds)
		window2 = src2.window(*overlap_bounds)
		
		# Arrondir les fenêtres pour éviter les problèmes de taille dus aux arrondis flottants
		window1 = window1.round_lengths().round_offsets()
		window2 = window2.round_lengths().round_offsets()
		
		# Lire les données dans la zone de recouvrement
		data1 = src1.read(1, window=window1).astype(np.float32)
		data2 = src2.read(1, window=window2).astype(np.float32)
		
		# Vérifier que les deux images ont la même taille
		if data1.shape != data2.shape:
			raise ValueError(f"Tailles incompatibles: data1={data1.shape}, data2={data2.shape}")
		
		# Redimensionner les poids si nécessaire pour correspondre à la taille réelle des données
		if weight1.shape != data1.shape:
			# Utiliser interpolation pour redimensionner les poids
			from scipy.ndimage import zoom
			zoom_factors = (data1.shape[0] / weight1.shape[0], data1.shape[1] / weight1.shape[1])
			weight1 = zoom(weight1, zoom_factors, order=1, mode='nearest')
			weight2 = zoom(weight2, zoom_factors, order=1, mode='nearest')
		
		# Vérifier que les poids ont maintenant la bonne taille
		if weight1.shape != data1.shape or weight2.shape != data1.shape:
			raise ValueError(f"Taille des poids incompatible après redimensionnement: poids1={weight1.shape}, poids2={weight2.shape}, données={data1.shape}")
		
		# Gérer les no-data
		no_data1 = src1.nodata if src1.nodata is not None else -9999
		no_data2 = src2.nodata if src2.nodata is not None else -9999
		
		# Masques pour les valeurs valides
		valid1 = (data1 != no_data1) & ~np.isnan(data1)
		valid2 = (data2 != no_data2) & ~np.isnan(data2)
		
		# Calculer la moyenne pondérée
		result = np.full_like(data1, no_data1, dtype=np.float32)
		
		# Cas où les deux valeurs sont valides
		both_valid = valid1 & valid2
		result[both_valid] = data1[both_valid] * weight1[both_valid] + data2[both_valid] * weight2[both_valid]
		
		# Cas où seule la dalle1 est valide
		only1 = valid1 & ~valid2
		result[only1] = data1[only1]
		
		# Cas où seule la dalle2 est valide
		only2 = valid2 & ~valid1
		result[only2] = data2[only2]
		
		# Créer le transform pour l'image de sortie avec les bounds de recouvrement
		# overlap_bounds = (left, bottom, right, top)
		left, bottom, right, top = overlap_bounds
		height, width = result.shape
		
		# Calculer la résolution
		pixel_size_x = (right - left) / width
		pixel_size_y = (top - bottom) / height
		
		# Créer le transform (coin haut-gauche)
		from rasterio.transform import Affine
		transform = Affine(pixel_size_x, 0.0, left,
						   0.0, -pixel_size_y, top)
		
		# Métadonnées pour l'image de sortie
		metadata = src1.meta.copy()
		metadata.update({
			'height': height,
			'width': width,
			'transform': transform,
			'dtype': result.dtype,
			'nodata': no_data1,
			'compress': 'lzw'
		})
		
		# Écrire l'image de sortie
		with rasterio.open(output_path, 'w', **metadata) as dst:
			dst.write(result, 1)


def assemble_tiles_and_overlaps(masks, overlaps, output_path):

    # --- 1) MOSAÏQUE DES TUILES PRINCIPALES ---
    src_masks = [rasterio.open(p) for p in masks]
    nodata = src_masks[0].nodata if src_masks[0].nodata is not None else -9999

    mosaic, mosaic_transform = merge(
        src_masks,
        nodata=nodata,
        method="last"
    )

    meta = src_masks[0].meta.copy()
    meta.update({
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": mosaic_transform,
        "nodata": nodata,
        "compress": "lzw"
    })

    # Convertir la mosaïque en édition locale
    final = mosaic.astype(np.float32)

    # --- 2) APPLICATION DES PATCHS DE RECOUVREMENT ---
    for patch_path in overlaps:
        with rasterio.open(patch_path) as patch:

            # Fenêtre d'insertion (patch → mosaïque)
            win = from_bounds(
                patch.bounds.left,
                patch.bounds.bottom,
                patch.bounds.right,
                patch.bounds.top,
                transform=mosaic_transform
            )
            # Rasterio peut retourner des offsets/taille flottants : on les arrondit pour indexer numpy
            win = win.round_offsets().round_lengths()
            row_off, col_off = int(win.row_off), int(win.col_off)
            height, width = int(win.height), int(win.width)

            # Lecture du patch à la taille exacte de la fenêtre
            patch_arr = patch.read(1, out_shape=(height, width))

            # Extraction du bloc correspondant dans la mosaïque
            final_block = final[0, row_off:row_off+height,
                                   col_off:col_off+width]

            # Pixels valides = non nodata
            valid = patch_arr != nodata

            # Remplacement
            final_block[valid] = patch_arr[valid]

            # Réécriture dans la mosaïque
            final[0, row_off:row_off+height,
                     col_off:col_off+width] = final_block

    # --- 3) ENREGISTREMENT ---
    with rasterio.open(output_path, "w", **meta) as dst:
        dst.write(final)

    # Fermeture des tuiles principales
    for s in src_masks:
        s.close()


def assemble_lines_and_overlaps(lines, overlaps, output_path):
    """
    Variante verticale : assemble les lignes mosaïquées, puis applique les patches verticaux.
    Comportement nodata : les pixels nodata des patches n'écrasent pas les pixels valides existants.
    """
    src_lines = [rasterio.open(p) for p in lines]
    nodata = src_lines[0].nodata if src_lines[0].nodata is not None else -9999

    mosaic, mosaic_transform = merge(
        src_lines,
        nodata=nodata,
        method="last"
    )

    meta = src_lines[0].meta.copy()
    meta.update({
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": mosaic_transform,
        "nodata": nodata,
        "compress": "lzw"
    })

    final = mosaic.astype(np.float32)

    for patch_path in overlaps:
        with rasterio.open(patch_path) as patch:
            win = from_bounds(
                patch.bounds.left,
                patch.bounds.bottom,
                patch.bounds.right,
                patch.bounds.top,
                transform=mosaic_transform
            )
            win = win.round_offsets().round_lengths()
            row_off, col_off = int(win.row_off), int(win.col_off)
            height, width = int(win.height), int(win.width)

            patch_arr = patch.read(1, out_shape=(height, width))
            final_block = final[0, row_off:row_off+height,
                                   col_off:col_off+width]
            valid = patch_arr != nodata
            final_block[valid] = patch_arr[valid]
            final[0, row_off:row_off+height,
                     col_off:col_off+width] = final_block

    with rasterio.open(output_path, "w", **meta) as dst:
        dst.write(final)

    for s in src_lines:
        s.close()


def Make_Assemblage_FINAL(chem_out, NbreDalleX, NbreDalleY, RepTra):
	"""
	Assemble les dalles en utilisant rasterio (version open source).
	"""
	####################################################################################################################################################
	## on raboute tout d'abord 2 dalles côte à côte (en X)	   #########################################################################################
	####################################################################################################################################################
	
	# Boucle avec barre de progression pour le raboutage des dalles côte à côte
	for y in tqdm(range(NbreDalleY), desc="Raboutage en colonnes"):
		for x in tqdm(range(NbreDalleX-1), desc="Progression en Colonne", leave=False):
			
			## Nom de la dalle courante	
			chem_MASK_QUALITY_dalle_xy = os.path.join(RepTra, "Dalle_%s_%s" % (x, y), "MASK_%s_%s.tif" % (x, y))
			chem_MASK_QUALITY_dalle_xy_droite = os.path.join(RepTra, "Dalle_%s_%s" % (x+1, y), "MASK_%s_%s.tif" % (x+1, y))
			
			# Calculer les bounds de recouvrement
			overlap_bounds = calculate_overlap_bounds(chem_MASK_QUALITY_dalle_xy, chem_MASK_QUALITY_dalle_xy_droite)
			
			if overlap_bounds is None:
				print(f"Attention: Pas de recouvrement entre dalle ({x},{y}) et ({x+1},{y})")
				continue
			
			# Lire la shape des données pour créer les poids de la bonne taille
			with rasterio.open(chem_MASK_QUALITY_dalle_xy_droite) as src:
				window = src.window(*overlap_bounds)
				target_shape = (int(window.height), int(window.width))
			
			# Créer les images de poids (horizontal: C/NC)
			weight_droite = create_weight_image_horizontal(overlap_bounds, target_shape)
			weight_gauche = 1.0 - weight_droite
			
			# Faire la moyenne pondérée dans la zone de recouvrement
			chem_reconstruction = os.path.join(RepTra, 'reconstruction_dalle_%s_%s_%s.tif' % (x, x+1, y))
			weighted_blend_overlap(
				chem_MASK_QUALITY_dalle_xy,
				chem_MASK_QUALITY_dalle_xy_droite,
				weight_gauche,
				weight_droite,
				overlap_bounds,
				chem_reconstruction
			)
			
	####################################################################################################################################################		
	## on raboute toutes les dalles sur une même rangée		#########################################################################################
	####################################################################################################################################################

	## on réassemble tout	
	for y in tqdm(range(NbreDalleY), desc="Raboutage en ligne"):
		
		# Construire la liste des images à assembler: dalles originales + images de transition
		image_paths = []
		overlaps = []
		# Ajouter les dalles originales
		for x in range(NbreDalleX):
			image_paths.append(os.path.join(RepTra, "Dalle_%s_%s" % (x, y), "MASK_%s_%s.tif" % (x, y)))
		
		# Ajouter les images de transition
		for x in range(NbreDalleX-1):
			overlaps.append(os.path.join(RepTra, 'reconstruction_dalle_%s_%s_%s.tif' % (x, x+1, y)))
		
		# Assembler de gauche à droite
		chem_final_tmp = os.path.join(RepTra, 'reconstruction_dalle_%s.tif' % y)

		# print("image_paths = ", image_paths)
		# print("chem_final_tmp = ", chem_final_tmp)
		# print("overlaps = ", overlaps)

		try:
			assemble_tiles_and_overlaps(image_paths, overlaps, chem_final_tmp)
		except Exception as e:
			print("ERREUR dans assemble_horizontal:", type(e), e)
			import traceback
			traceback.print_exc()
			print("Vérifiez que la fonction assemble_horizontal est bien importée/définie, qu'il n'y a pas d'erreur de nom de variable ou d'accès aux fichiers ci-dessus.")
			print("Voici la liste des images à assembler, pour vérification des accès fichiers :")
			for ip in image_paths:
				print("  - ", ip, "-->", os.path.exists(ip))
			print("chem_final_tmp =", chem_final_tmp, "--> dossier existe ?", os.path.exists(os.path.dirname(chem_final_tmp)))
			raise  # relancer l'exception pour arrêt si debug

	####################################################################################################################################################
	## on assemble les rangées entre elles		 #####################################################################################################
	####################################################################################################################################################
	
	## on réassemble tout	
	for y in tqdm(range(NbreDalleY-1), desc="Assemblage final"):
				
		chem_dalle_y = os.path.join(RepTra, 'reconstruction_dalle_%s.tif' % y)
		chem_dalle_y_dessous = os.path.join(RepTra, 'reconstruction_dalle_%s.tif' % (y+1))
		
		# Calculer les bounds de recouvrement
		overlap_bounds = calculate_overlap_bounds(chem_dalle_y, chem_dalle_y_dessous)
		
		if overlap_bounds is None:
			print(f"Attention: Pas de recouvrement entre ligne {y} et {y+1}")
			continue
		
		# Lire la shape des données pour créer les poids de la bonne taille
		with rasterio.open(chem_dalle_y_dessous) as src:
			window = src.window(*overlap_bounds)
			target_shape = (int(window.height), int(window.width))
		
		# Créer les images de poids (vertical: L/NL)
		weight_bas = create_weight_image_vertical(overlap_bounds, target_shape)
		weight_haut = 1.0 - weight_bas
		
		# Faire la moyenne pondérée dans la zone de recouvrement
		chem_reconstruction = os.path.join(RepTra, 'reconstruction_dalle_%s_%s.tif' % (y, y+1))
		weighted_blend_overlap(
			chem_dalle_y,
			chem_dalle_y_dessous,
			weight_haut,
			weight_bas,
			overlap_bounds,
			chem_reconstruction
		)
			
	####################################################################################################################################################
	## Assemblage final de toutes les lignes		 #####################################################################################################
	####################################################################################################################################################
	
	# Assemblage final vertical : lignes mosaïquées + patches verticaux
	image_paths = []
	overlaps = []
	
	# Ajouter les lignes assemblées
	for y in range(NbreDalleY):
		chem_dalle_y = os.path.join(RepTra, 'reconstruction_dalle_%s.tif' % y)
		image_paths.append(chem_dalle_y)
		
	# Ajouter les images de transition verticales
	for y in range(NbreDalleY-1):
		overlaps.append(os.path.join(RepTra, 'reconstruction_dalle_%s_%s.tif' % (y, y+1)))
	
	# Assemblage final (vertical)
	assemble_lines_and_overlaps(image_paths, overlaps, chem_out)

