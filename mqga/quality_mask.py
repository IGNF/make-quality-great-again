#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calcul du masque de qualité (percentile local)."""
import numpy as np
import rasterio
from scipy.ndimage import generic_filter

from mqga.io_raster import (
    read_as_2D_float,
    save_ABSOLUTE_image_with_same_geometry,
)


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
	
