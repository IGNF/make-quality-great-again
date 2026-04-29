#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, sys, glob, shutil, datetime, time, math, re
#from lxml import etree
#from datetime import date
#from osgeo import ogr
#from shapely.geometry import Point, LineString, Polygon, MultiPolygon
#from shapefile import Shapefile
#from shapefile.const import SHPT_POINT, SHPT_POLYGON, SHPT_POLYLINE
import shutil
import argparse
from osgeo import ogr, gdal
#from geojson import Point, Polygon, MultiPolygon, Feature, FeatureCollection, dump
#from shapely.geometry import Polygon
# from rasterstats import zonal_stats  # Non utilisé actuellement
#
import numpy as np
import rasterio
from rasterio.merge import merge
from rasterio.windows import from_bounds
import imageio
import scipy as sp
#import scipy.misc as sm
import scipy.sparse.linalg as ssl
from scipy.ndimage import generic_filter, uniform_filter, label, binary_dilation
from scipy.interpolate import LinearNDInterpolator, griddata
from scipy.spatial.distance import cdist
from scipy.spatial import cKDTree

import time

from multiprocessing import Process,Pool,cpu_count
import concurrent.futures

from tqdm import tqdm

import signal

# Outils OTB (non utilisés actuellement dans le code, commentés pour référence)
# otbcli_ConnectedComponentSegmentation='/home/OTB/OTB-7.4.0-Linux64/bin/otbcli_ConnectedComponentSegmentation'
# otbcli_BinaryMorphologicalOperation='/home/OTB/OTB-7.4.0-Linux64/bin/otbcli_BinaryMorphologicalOperation'
 
############################################################################################################################		 
def GetValue(listInfo,chaine):
	
	for l in listInfo:
		lig = l.strip()
		if lig.find(chaine) != -1:
			return lig.split()[-1]
			
############################################################################################################################		 
############################################################################################################################		 
def GetInfo(cheminTIF):
	"""
	Récupère les métadonnées d'une image GeoTIFF.
	Version open source utilisant rasterio.
	
	Returns:
		[PasX, PasY, Projection, X_0, X_1, Y_0, Y_1, phasage, NbreCol, NbreLig, GModel, GRaster]
	"""
	# Ouvrir l'image avec rasterio
	with rasterio.open(cheminTIF, 'r') as src:
		# Pas en X et Y (résolution)
		transform = src.transform
		PasX = abs(transform[0])  # pixel width
		PasY = abs(transform[4])  # pixel height (généralement négatif, on prend la valeur absolue)
		
		# Bounds (positions)
		bounds = src.bounds
		X_0 = bounds.left   # GAUCHE
		X_1 = bounds.right   # DROITE
		Y_0 = bounds.bottom  # BAS
		Y_1 = bounds.top     # HAUT
		
		# Nombre de colonnes et lignes
		NbreCol = float(src.width)
		NbreLig = float(src.height)
		
		# Code EPSG / Projection
		if src.crs is not None:
			Projection = int(src.crs.to_epsg()) if src.crs.to_epsg() is not None else -1
		else:
			Projection = -1
		
		# GTModelTypeGeoKey et GTRasterTypeGeoKey
		# Ces clés sont dans les tags GeoTIFF, généralement dans les tags de la bande
		GModel = -1
		GRaster = -1
		# Essayer d'abord les tags du dataset
		if hasattr(src, 'tags') and src.tags():
			tags = src.tags()
			if 'GTModelTypeGeoKey' in tags:
				try:
					GModel = int(tags['GTModelTypeGeoKey'])
				except (ValueError, TypeError):
					pass
			if 'GTRasterTypeGeoKey' in tags:
				try:
					GRaster = int(tags['GTRasterTypeGeoKey'])
				except (ValueError, TypeError):
					pass
		# Si pas trouvé, essayer les tags de la première bande
		if (GModel == -1 or GRaster == -1) and hasattr(src, 'tags'):
			try:
				band_tags = src.tags(1)
				if band_tags:
					if 'GTModelTypeGeoKey' in band_tags and GModel == -1:
						try:
							GModel = int(band_tags['GTModelTypeGeoKey'])
						except (ValueError, TypeError):
							pass
					if 'GTRasterTypeGeoKey' in band_tags and GRaster == -1:
						try:
							GRaster = int(band_tags['GTRasterTypeGeoKey'])
						except (ValueError, TypeError):
							pass
			except:
				pass
		
		# Phasage : déterminer si c'est "HG" (haut-gauche) ou "CP" (centre pixel)
		phasage = "HG"  # Par défaut
		if hasattr(src, 'tags') and src.tags():
			# Vérifier les tags pour le phasage
			tags = src.tags()
			if 'AREA_OR_POINT' in tags and tags['AREA_OR_POINT'] == 'Point':
				phasage = "CP"
	
	return [PasX, PasY, Projection, X_0, X_1, Y_0, Y_1, phasage, NbreCol, NbreLig, GModel, GRaster] 
	
#############################################################################################################################	
def read_as_2D_float(filename,no_data):
	
	with rasterio.open(filename, 'r') as dataset:
		data= dataset.read(1).astype(np.float32)
		data[data==no_data] = np.nan
		
	return data

#############################################################################################################################	
def get_pied_histo(histo, seuil_pied):
	for i, elt in enumerate(histo):
		if elt > seuil_pied:
			return i

#############################################################################################################################	
def get_haut_histo(histo, seuil_haut):
	for i, elt in enumerate(histo):
		if elt > seuil_haut:
			return i

#############################################################################################################################	
def calculate_cdf_percent(pixel_array, percentile):
	# Exclude the no-data value and calculate the 5% value of the CDF
	pixel_array = pixel_array[pixel_array != -9999]
	if len(pixel_array) == 0:
		return -9999  # Return no-data value if the window only contains no-data values
	sorted_pixels = np.sort(pixel_array)
	index_5_percent = int(np.ceil(percentile * len(sorted_pixels))) - 1
	return sorted_pixels[max(0, index_5_percent)]
	
#############################################################################################################################	
# def process_image(image, percentile):
# 	# Pad image to handle the borders
# 	padded_image = np.pad(image, 50, mode='constant', constant_values=-9999)
	
# 	# Use generic_filter from scipy.ndimage to apply the function over a 101x101 window
# 	result = generic_filter(padded_image, lambda x: calculate_cdf_percent(x, percentile), size=(101, 101), mode='constant', cval=-9999)
	
# 	# Crop the padded area off the result
# 	return result[50:-50, 50:-50]
	
#############################################################################################################################	
def process_image(image,dl,no_data,percentile):
	# Pad image to handle the borders
	padded_image = np.pad(image, dl, mode='constant', constant_values=no_data)
	# Use generic_filter from scipy.ndimage to apply the function over a 101x101 window
	result = generic_filter(padded_image, lambda x: calculate_cdf_percent(x, percentile), size=(2*dl+1, 2*dl+1), mode='constant', cval=no_data)
	# Crop the padded area off the result
	return result[dl:-dl, dl:-dl]
	
#############################################################################################################################	
def save_ABSOLUTE_image_with_same_geometry(image, output_filename, src_filename):
    # Calculer la valeur absolue de l'image
    abs_image = np.abs(image)
    
    # Ouvrez l'image source pour lire sa géométrie
    with rasterio.open(src_filename) as src:
        metadata = src.meta.copy()  # Copiez les métadonnées de l'image source
        
    # Mettez à jour les métadonnées avec les nouvelles dimensions si nécessaire
    metadata['height'], metadata['width'] = abs_image.shape
    metadata['dtype'] = abs_image.dtype  # Assurez-vous que le type de données correspond à l'image de sortie
    
    # Utilisez les métadonnées copiées pour écrire l'image dans un fichier .tif
    with rasterio.open(output_filename, 'w', **metadata) as dst:
        dst.write(abs_image, 1)  # Écrit l'image dans la première bande en assumant qu'il s'agit d'une image à une seule bande
        
#############################################################################################################################	
def save_image_with_same_geometry(image, output_filename, src_filename):
	# Ouvrez l'image source pour lire sa géométrie
	with rasterio.open(src_filename) as src:
		metadata = src.meta.copy()  # Copiez les métadonnées de l'image source
		
	# Mettez à jour les métadonnées avec les nouvelles dimensions si nécessaire
	metadata['height'], metadata['width'] = image.shape
	metadata['dtype'] = image.dtype  # Assurez-vous que le type de données correspond à l'image de sortie
	
	# Utilisez les métadonnées copiées pour écrire l'image dans un fichier .tif
	with rasterio.open(output_filename, 'w', **metadata) as dst:
		dst.write(image, 1)  # Écrit l'image dans la première bande en assumant qu'il s'agit d'une image à une seule bande

#############################################################################################################################	
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
	
	print(f"Interpolation terminée (LinearNDInterpolator): {total_nodata_interpolated} pixels nodata interpolés.")

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
			
			print(f"Traitement par blocs (griddata): {n_blocks_x}x{n_blocks_y} blocs de {block_size}x{block_size} pixels")
			
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
	
	print(f"Interpolation terminée (griddata): {total_nodata_interpolated} pixels nodata interpolés.")

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
			
			print(f"Traitement par blocs (IDW): {n_blocks_x}x{n_blocks_y} blocs de {block_size}x{block_size} pixels")
			
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
	
	print(f"Interpolation terminée (IDW): {total_nodata_interpolated} pixels nodata interpolés.")

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
			
			print(f"Traitement par blocs (IDW vectorisé): {n_blocks_x}x{n_blocks_y} blocs de {block_size}x{block_size} pixels ({n_jobs} processus)")
			
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
	
	print(f"Interpolation terminée (IDW vectorisé): {total_nodata_interpolated} pixels nodata interpolés.")

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
		print("Aucun pixel nodata à interpoler.")
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
	print(f"Interpolation terminée (rapide): {nodata_filled} pixels nodata interpolés en {iteration+1} itérations.")

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
		print("Aucun pixel nodata à interpoler.")
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
	
	print(f"Interpolation par fenêtre glissante ({window_size}x{window_size})...")
	result = generic_filter(data, fill_nodata, size=window_size, mode='constant', cval=no_data)
	
	# Garder les nodata qui n'ont pas pu être interpolés
	result[mask_nodata & (result == no_data)] = no_data
	
	# Sauvegarder
	metadata['dtype'] = result.dtype
	with rasterio.open(chem_out, 'w', **metadata) as dst:
		dst.write(result, 1)
	
	nodata_filled = np.sum((mask_nodata) & (result != no_data))
	print(f"Interpolation terminée (fenêtre): {nodata_filled} pixels nodata interpolés.")

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
		print("Aucun pixel nodata à interpoler.")
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
	print(f"Interpolation terminée (hybride): {nodata_filled} pixels nodata interpolés.")

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

#################################################################################################### 
### calcule le nombre de dalles en X et en Y en fonction des paramètres de chantier
def CalculNombreDallesXY_NEW(NbColonnes,NbLignes,dalleX,dalleY):
				
	### Calcul du nombre de dalles en X		
	if (NbColonnes%dalleX == 0):
		NbreDalleX=NbColonnes/dalleX
	else: 
		NbreDalleX=int(NbColonnes/dalleX)+1
		
	### Calcul du nombre de dalles en Y		
	if (NbLignes%dalleY == 0):
		NbreDalleY=NbLignes/dalleY
	else: 
		NbreDalleY=int(NbLignes/dalleY)+1
	
	return (NbreDalleX,NbreDalleY)  
	
#################################################################################################### 
### calcule le nombre de dalles en X et en Y en fonction des paramètres de chantier
def CalculNombreDallesXY(NbColonnes,NbLignes,Taille_dalle,Recouv_entre_dalles):
		
	### pour calculer le nombre de dalles en X et en Y		
	NumX=NbColonnes-Taille_dalle
	NumY=NbLignes-Taille_dalle
	Denom=Taille_dalle-Recouv_entre_dalles
		
	### Calcul du nombre de dalles en X		
	if (NumX%Denom == 0):
		NbreDalleX=int(NumX/Denom+1)
	else: 
		NbreDalleX=int((NumX/Denom)+1)+1
		
	### Calcul du nombre de dalles en Y		
	if (NumY%Denom == 0):
		NbreDalleY=int(NumY/Denom+1)
	else: 
		NbreDalleY=int((NumY/Denom)+1)+1
		
	return (NbreDalleX,NbreDalleY) 
	
#################################################################################################### 
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


#################################################################################################### 
def assemble_horizontal_OLD(image_paths, output_path):
	"""Assemble des images de gauche à droite. Les pixels suivants écrasent les précédents en cas de recouvrement.
	Les pixels NoData ne remplacent jamais les pixels valides grâce à skip_empty=True."""
	from rasterio.merge import merge
	
	# Trier les images par leur position X (left) pour garantir l'ordre de gauche à droite
	image_bounds = []
	for path in image_paths:
		with rasterio.open(path) as src:
			image_bounds.append((src.bounds.left, path))
	
	# Trier par position X (left)
	image_bounds.sort(key=lambda x: x[0])
	sorted_paths = [path for _, path in image_bounds]
	
	# Récupérer la valeur nodata de la première image (on suppose qu'elles sont identiques)
	with rasterio.open(sorted_paths[0]) as first:
		nodata_val = first.nodata if first.nodata is not None else -9999
	
	# Ouvrir toutes les images dans l'ordre trié
	srcs = [rasterio.open(path) for path in sorted_paths]
	
	try:
		# Utiliser merge pour assembler (method='last' pour que les pixels suivants écrasent les précédents)
		# skip_empty=True empêche les pixels NoData d'écraser les pixels valides
		# L'ordre dans la liste détermine la priorité : les dernières images écrasent les premières
		mosaic, out_trans = merge(
			srcs,
			method='last',
			nodata=nodata_val,
			skip_empty=True  # Ne pas remplacer les pixels valides par du NoData
		)
		
		# Métadonnées de sortie
		out_meta = srcs[0].meta.copy()
		out_meta.update({
			'driver': 'GTiff',
			'height': mosaic.shape[1],
			'width': mosaic.shape[2],
			'transform': out_trans,
			'compress': 'lzw',
			'nodata': nodata_val
		})
		
		# Écrire l'image assemblée
		with rasterio.open(output_path, 'w', **out_meta) as dst:
			dst.write(mosaic)
	finally:
		# Fermer toutes les sources
		for src in srcs:
			src.close()

#################################################################################################### 	
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

	
#################################################################################################### 
def crop_tile(args):
	"""
	Fonction pour découper une dalle d'une image (utilisée en parallèle).
	
	Args:
		args: tuple (chem_in, Chem_decoup, ligmin, ligmax, colmin, colmax)
	"""
	chem_in, Chem_decoup, ligmin, ligmax, colmin, colmax = args
	
	# Ouvrir l'image source
	with rasterio.open(chem_in, 'r') as src:
		# Calculer les dimensions de la fenêtre
		# ligmin inclus, ligmax exclusif -> height = ligmax - ligmin
		# colmin inclus, colmax exclusif -> width = colmax - colmin
		height = ligmax - ligmin
		width = colmax - colmin
		
		# Tronquer si la fenêtre dépasse les limites de l'image
		max_row = src.height
		max_col = src.width
		
		# Si la fenêtre est complètement en dehors de l'image, on ne fait rien
		if ligmin >= max_row or colmin >= max_col or ligmax <= 0 or colmax <= 0:
			print(f"Attention: Fenêtre complètement hors limites pour {Chem_decoup}")
			return
		
		# Ajuster les offsets pour qu'ils soient >= 0
		row_off = max(0, ligmin)
		col_off = max(0, colmin)
		
		# Ajuster la taille si la fenêtre dépasse les limites
		if row_off + height > max_row:
			height = max_row - row_off
		if col_off + width > max_col:
			width = max_col - col_off
		
		# Vérifier que la fenêtre est valide après ajustement
		if height <= 0 or width <= 0:
			print(f"Attention: Fenêtre invalide pour {Chem_decoup} (height={height}, width={width})")
			return
		
		# Créer la fenêtre de découpage
		window = rasterio.windows.Window(col_off=col_off, row_off=row_off, width=width, height=height)
		
		# Lire les données dans la fenêtre
		data = src.read(1, window=window)
		
		# Calculer le nouveau transform pour cette fenêtre
		transform = rasterio.windows.transform(window, src.transform)
		
		# Copier les métadonnées et les mettre à jour
		metadata = src.meta.copy()
		metadata.update({
			'height': height,
			'width': width,
			'transform': transform,
			'compress': 'lzw'
		})
		
		# Écrire la dalle découpée
		with rasterio.open(Chem_decoup, 'w', **metadata) as dst:
			dst.write(data, 1)

#################################################################################################### 
def init_rep_tra(RepTra, clean=False):
	"""
	Initialise le répertoire de travail temporaire.
	
	Args:
		RepTra: Chemin vers le répertoire de travail
		clean: Si True, supprime le contenu du répertoire s'il existe déjà
	"""
	# Créer le répertoire parent si nécessaire
	parent_dir = os.path.dirname(RepTra) if os.path.dirname(RepTra) else '.'
	if parent_dir and not os.path.exists(parent_dir):
		os.makedirs(parent_dir, exist_ok=True)
	
	# Si le répertoire existe déjà
	if os.path.exists(RepTra):
		if clean:
			print(f"Nettoyage du répertoire temporaire: {RepTra}")
			# Supprimer tout le contenu
			for item in os.listdir(RepTra):
				item_path = os.path.join(RepTra, item)
				if os.path.isdir(item_path):
					shutil.rmtree(item_path)
				else:
					os.remove(item_path)
			print(f"  ✓ Répertoire nettoyé")
		else:
			print(f"Répertoire temporaire existant: {RepTra} (contenu conservé)")
	else:
		# Créer le répertoire
		os.makedirs(RepTra, exist_ok=True)
		print(f"Répertoire temporaire créé: {RepTra}")

#################################################################################################### 
def MakeDecoupage(chem_in, RepTra, NbreDalleX, NbreDalleY, iTailleparcelle, iTailleRecouvrement, iNbreCPU):
	"""
	Découpe l'image en dalles en utilisant rasterio (version open source).
	"""
	tasks = []
	
	for x in range(NbreDalleX):
		for y in range(NbreDalleY):
			#créer le nom du répertoire
			RepDalleXY=os.path.join(RepTra,"Dalle_%s_%s"%(x,y))
			#créer le répertoire, s'il n'existe pas déjà (avec création récursive)
			os.makedirs(RepDalleXY, exist_ok=True)
								
			#Détermination de col_min, col_max, lig_min, lig_max pour la dalle XY
			colminDalleXY=x*(iTailleparcelle-iTailleRecouvrement)
			colmaxDalleXY=x*(iTailleparcelle-iTailleRecouvrement)+iTailleparcelle
			ligminDalleXY=y*(iTailleparcelle-iTailleRecouvrement)
			ligmaxDalleXY=y*(iTailleparcelle-iTailleRecouvrement)+iTailleparcelle
				
			#fichier out mns
			Chem_decoup=os.path.join(RepDalleXY,"IN_%s_%s.tif"%(x,y))
			
			# Préparer les arguments pour la fonction de découpage
			args = (chem_in, Chem_decoup, ligminDalleXY, ligmaxDalleXY, colminDalleXY, colmaxDalleXY)
			tasks.append(args)
			
	# Initialize the pool
	with Pool(processes=iNbreCPU, initializer=init_worker) as pool:
		results = list(tqdm(pool.imap_unordered(crop_tile, tasks), total=len(tasks), desc="Découpage en parallèle des dalles"))
		
	return
			
#################################################################################################### 
def diff_2_mask_quality(args):
	chem_in, chem_out, dl, no_data, percentile = args
	#print("chem_in    >>> ",chem_in)
	#print("chem_out   >>> ",chem_out)
	#print("dl         >>> ",dl)
	#print("no_data    >>> ",no_data)
	#print("percentile >>> ",percentile)
	
	data_in = read_as_2D_float(chem_in, no_data)
	result = process_image(data_in, dl, no_data, percentile)
	save_ABSOLUTE_image_with_same_geometry(result, chem_out, chem_in)
	return

#################################################################################################### 
# def diff_2_mask_quality_BIS(chem_in, chem_out, dl, no_data, percentile):

# 	data_in = read_as_2D_float(chem_in,no_data)
# 	result = process_image(data_in,dl,no_data)
# 	save_ABSOLUTE_image_with_same_geometry(result, chem_out, chem_in)
# 	return
	
#################################################################################################### 
def create_negative_image(chem_in, chem_out, no_data=-9999):
	"""
	Crée une version "négative" de l'image : remplace les pixels > 0 par no_data.
	
	Args:
		chem_in: Chemin vers l'image d'entrée
		chem_out: Chemin vers l'image de sortie
		no_data: Valeur nodata (par défaut -9999)
	"""
	# Lire l'image
	with rasterio.open(chem_in, 'r') as src:
		data = src.read(1).astype(np.float32)
		metadata = src.meta.copy()
	
	# Créer une copie des données
	result = data.copy()
	
	# Remplacer les pixels > 0 par no_data
	# Les pixels <= 0 sont conservés, ainsi que les pixels nodata existants
	mask_positive = (data > 0) & (data != no_data) & ~np.isnan(data)
	result[mask_positive] = no_data
	
	# Sauvegarder l'image résultante
	metadata['dtype'] = result.dtype
	with rasterio.open(chem_out, 'w', **metadata) as dst:
		dst.write(result, 1)
	
#################################################################################################### 
def init_worker():
	signal.signal(signal.SIGINT, signal.SIG_IGN)
	
#################################################################################################### 
def DoParallel(RepTra, NbreDalleX, NbreDalleY, dl, no_data, percentile, iNbreCPU):
		
	tasks = []
	
	## lancement sur chaque dalle
	for x in range(NbreDalleX):
		for y in range(NbreDalleY):
			#créer le nom du répertoire
			RepDalleXY=os.path.join(RepTra,"Dalle_%s_%s"%(x,y))
			#fichier out mns
			chem_in_dalle=os.path.join(RepDalleXY,"IN_%s_%s.tif"%(x,y))
			chem_in_dalle_NEG=os.path.join(RepDalleXY,"IN_%s_%s_NEG.tif"%(x,y))
			chem_out_dalle=os.path.join(RepDalleXY,"MASK_%s_%s.tif"%(x,y))
			#
			# Créer la version négative 
			create_negative_image(chem_in_dalle, chem_in_dalle_NEG, no_data)
			#
			args=(chem_in_dalle_NEG, chem_out_dalle, dl, no_data, percentile)
			tasks.append(args)

	# Initialize the pool
	with Pool(processes=iNbreCPU, initializer=init_worker) as pool:
		# Use tqdm to show the progress bar
		results = list(tqdm(pool.imap_unordered(diff_2_mask_quality, tasks), total=len(tasks),desc="Calcul des masques de qualité par dalle en //"))
				

#################################################################################################### 
def DoParallel_OLD(RepTra, NbreDalleX, NbreDalleY, dl, no_data, percentile, iNbreCPU):
		
	tasks = []
	
	## lancement sur chaque dalle
	for x in range(NbreDalleX):
		for y in range(NbreDalleY):
			#créer le nom du répertoire
			RepDalleXY=os.path.join(RepTra,"Dalle_%s_%s"%(x,y))
			#fichier out mns
			chem_in_dalle=os.path.join(RepDalleXY,"IN_%s_%s.tif"%(x,y))
			chem_out_dalle=os.path.join(RepDalleXY,"MASK_%s_%s.tif"%(x,y))
			#
			args=(chem_in_dalle, chem_out_dalle, dl, no_data, percentile)
			tasks.append(args)
		
	# Initialize the pool
	with Pool(processes=iNbreCPU, initializer=init_worker) as pool:
		# Use tqdm to show the progress bar
		results = list(tqdm(pool.imap_unordered(diff_2_mask_quality, tasks), total=len(tasks),desc="Calcul des masques de qualité par dalle en //"))
				
#############################################################################################################################	
#############################################################################################################################	
#############################################################################################################################					
if __name__ == '__main__':

	try:
		#
		start_time = time.time()
		#
		parser = argparse.ArgumentParser(description='MAKE QUALITY GREAT AGAIN ALL ZONE - Version Spatialisée & Parallélisée')
		parser.add_argument("-diff", type=str, help="Différence DSM/DTM en entrée")
		parser.add_argument("-out", type=str, help="Masque de Qualité en sortie")
		parser.add_argument("-no", type=int, default=-9999, help="Valeur de No Data")
		#parser.add_argument("-bin", type=float, default=0.1, help="Valeur de bin")
		parser.add_argument("-per", type=float, default=0.05, help="Valeur de percentile")
		#parser.add_argument("-count", type=int, default=100, help="Nombre de points minimal pour que la statistique soit valide")
		parser.add_argument("-demiwinl", type=int, default=50, help="Demie-taille en ligne de la fenêtre d'analyse")
		parser.add_argument("-demiwinc", type=int, default=50, help="Demie-taille en colonne de la fenêtre d'analyse")
		#parser.add_argument("-offl", type=int, default=25, help="Décalage en ligne entre 2 analyses / Facteur de sous-échantillonnage en ligne")
		#parser.add_argument("-offc", type=int, default=25, help="Décalage en colonne entre 2 analyses / Facteur de sous-échantillonnage en colonne")
		parser.add_argument("-tile", type=int, default=500, help="Tile / Taille de la tuile")
		parser.add_argument("-pad", type=int, default=50, help="Pad / Recouvrement entre tuiles")	
		parser.add_argument("-RepTra", type=str, help="Répertoire de Travail")
		parser.add_argument("-cpu", type=int, help="Nombre de CPU diponibles")
		parser.add_argument("-winavg", type=int, default=50, help="Taille de la fenêtre glissante pour la moyenne (par défaut 50x50)")
		parser.add_argument("-interp", type=str, default="idw", 
							choices=["griddata", "idw", "idw_old", "window", "linearnd", "fast", "hybrid"],
							help="Méthode d'interpolation pour les pixels nodata (défaut: idw - optimisé et parallélisé, hybrid = méthode hybride xingng)")
		parser.add_argument("-clean", action='store_true', 
							help="Supprimer le contenu du répertoire temporaire s'il existe déjà")
		
		args = parser.parse_args(sys.argv[1:])
		#
		chem_in=args.diff
		chem_out=args.out
		no_data=args.no
		#bin_step=args.bin
		percentile=args.per
		#count_valid=args.count
		demiwinsize_lig=args.demiwinl
		demiwinsize_col=args.demiwinc
		#offset_lig=args.offl
		#offset_col=args.offc
		RepTra=args.RepTra
		iNbreCPU=args.cpu
		#print(iNbreCPU)
		iTailleparcelle=args.tile
		iTailleRecouvrement=args.pad
		winavg_size=args.winavg
		dl=demiwinsize_lig
		dc=demiwinsize_col
		
		### Initialisation du répertoire temporaire
		clean_rep = args.clean if hasattr(args, 'clean') else False
		init_rep_tra(RepTra, clean=clean_rep)
		
		### Découpage
		infos=GetInfo(chem_in)
		
		NbreCol=infos[8]
		NbreLig=infos[9]
				
		NombreDallesXY=CalculNombreDallesXY(NbreCol,NbreLig,iTailleparcelle,iTailleRecouvrement)
		NbreDalleX=NombreDallesXY[0]
		NbreDalleY=NombreDallesXY[1]
		
		#Decoupage en parallèle
		MakeDecoupage(chem_in, RepTra, NbreDalleX, NbreDalleY, iTailleparcelle, iTailleRecouvrement, iNbreCPU)
		
		#Traitement/Calcul en parallèle
		DoParallel(RepTra, NbreDalleX, NbreDalleY, dl, no_data, percentile, iNbreCPU)
		
		#### Assemblage final (version open source)
		# Le fichier d'assemblage est temporaire (fini par _tmp.tif)
		chem_out_tmp = chem_out.replace('.tif', '_tmp.tif')
		Make_Assemblage_FINAL(chem_out_tmp, NbreDalleX, NbreDalleY, RepTra)
		
		#### Post-traitement 2: Bouchage des zones avec no data (interpolation)
		chem_out_clean_bouchage = chem_out.replace('.tif', '_clean_bouchage.tif')
		interp_method = args.interp
		print(f"Post-traitement: Interpolation des pixels nodata (méthode: {interp_method})...")
		
		if interp_method == "griddata":
			interpolate_nodata_griddata(chem_out_tmp, chem_out_clean_bouchage, no_data)
		elif interp_method == "idw":
			# Utiliser la version vectorisée et parallélisée (beaucoup plus rapide)
			interpolate_nodata_idw_vectorized(chem_out_tmp, chem_out_clean_bouchage, no_data, n_jobs=iNbreCPU)
		elif interp_method == "idw_old":
			# Ancienne version (lente, conservée pour compatibilité)
			interpolate_nodata_idw(chem_out_tmp, chem_out_clean_bouchage, no_data)
		elif interp_method == "window":
			interpolate_nodata_window(chem_out_tmp, chem_out_clean_bouchage, no_data)
		elif interp_method == "linearnd":
			interpolate_nodata_with_linearnd(chem_out_tmp, chem_out_clean_bouchage, no_data)
		elif interp_method == "fast":
			# Méthode très rapide (moins précise mais beaucoup plus rapide)
			interpolate_nodata_fast(chem_out_tmp, chem_out_clean_bouchage, no_data)
		elif interp_method == "hybrid":
			# Méthode hybride de xingng (équivalent à -FB:2:C:50,1:1:50:1)
			# Combine interpolation locale sur pixels de bord et constante statistique
			interpolate_nodata_hybrid(
				chem_out_tmp, chem_out_clean_bouchage, no_data,
				connectivity=4, seuil_percent=50, poids=1, rayon=50, n=1
			)
		else:
			print(f"Attention: Méthode d'interpolation '{interp_method}' non reconnue. Utilisation de 'fast' par défaut.")
			interpolate_nodata_fast(chem_out_tmp, chem_out_clean_bouchage, no_data)

		#### Post-traitement 3: Moyenne sur fenêtre glissante (open source)
		# Le fichier final utilise le nom spécifié dans --out
		apply_moving_average(chem_out_clean_bouchage, chem_out, window_size=winavg_size, no_data=no_data)
		
		print("Fichier final généré: %s" % chem_out)
		
		# Afficher le temps total d'exécution arrondi à la seconde
		elapsed_time = time.time() - start_time
		elapsed_seconds = int(round(elapsed_time))
		elapsed_minutes = elapsed_seconds // 60
		elapsed_secs = elapsed_seconds % 60
		if elapsed_minutes > 0:
			print(f"Temps total d'exécution: {elapsed_minutes} minute(s) et {elapsed_secs} seconde(s)")
		else:
			print(f"Temps total d'exécution: {elapsed_seconds} seconde(s)")
				  
	except (RuntimeError, TypeError, NameError):
		print ("ERREUR: ", NameError)


