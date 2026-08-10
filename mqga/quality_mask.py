#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calcul du masque de qualité (percentile local ou MAD)."""
import numpy as np
import rasterio
from scipy.ndimage import generic_filter

from mqga.io_raster import (
    read_as_2D_float,
    save_ABSOLUTE_image_with_same_geometry,
)

# σ ≈ 1.4826·MAD ; LE90 ≈ 1.645·σ → ε ≈ 2.44·MAD
DEFAULT_MAD_K = 2.44

# Effectif minimal de pixels négatifs valides dans la fenêtre locale
# (STANAG / historique gs2_mnt_mq_local.filter_stats size=167).
# En dessous → NoData (la statistique n'est pas considérée fiable).
DEFAULT_MIN_VALID = 167


def calculate_cdf_percent(pixel_array, percentile, no_data=-9999, min_valid=DEFAULT_MIN_VALID):
	# Exclude the no-data value and calculate the percentile of the CDF
	pixel_array = np.asarray(pixel_array, dtype=np.float64)
	pixel_array = pixel_array[np.isfinite(pixel_array) & (pixel_array != no_data)]
	if len(pixel_array) < min_valid:
		return no_data
	sorted_pixels = np.sort(pixel_array)
	index_5_percent = int(np.ceil(percentile * len(sorted_pixels))) - 1
	return sorted_pixels[max(0, index_5_percent)]


def calculate_mad(pixel_array, no_data=-9999, min_valid=DEFAULT_MIN_VALID):
	"""
	MAD sur les échantillons valides de la fenêtre (partie négative déjà filtrée en amont).
	MAD = median(|x − median(x)|)
	"""
	arr = np.asarray(pixel_array, dtype=np.float64)
	v = arr[np.isfinite(arr) & (arr != no_data)]
	if v.size < min_valid:
		return no_data
	med = np.median(v)
	return float(np.median(np.abs(v - med)))


def calculate_mad_le90(
	pixel_array, mad_k=DEFAULT_MAD_K, no_data=-9999, min_valid=DEFAULT_MIN_VALID,
):
	"""k·MAD (défaut mad_k=2.44 ≈ LE90 sous hypothèse gaussienne). Le biais |b| est ajouté ensuite."""
	mad = calculate_mad(pixel_array, no_data=no_data, min_valid=min_valid)
	if mad == no_data:
		return no_data
	return mad_k * mad


def _apply_additive_bias(result, bias, no_data, stat):
	"""ε ← |b| + ε_stat (MAD déjà ≥0 ; percentile ≤0 → on décale avant abs)."""
	b = abs(float(bias))
	if b == 0:
		return result
	out = np.array(result, dtype=np.float64, copy=True)
	valid = np.isfinite(out) & (out != no_data)
	if stat == "mad":
		out[valid] = out[valid] + b
	else:
		out[valid] = -(np.abs(out[valid]) + b)
	return out


def process_image(
	image, dl, no_data, percentile, stat="mad", mad_k=DEFAULT_MAD_K, bias=0.0,
	min_valid=DEFAULT_MIN_VALID,
):
	# Pad image to handle the borders
	padded_image = np.pad(image, dl, mode='constant', constant_values=no_data)
	# Use generic_filter from scipy.ndimage to apply the function over a local window
	if stat == "mad":
		fn = lambda x: calculate_mad_le90(
			x, mad_k=mad_k, no_data=no_data, min_valid=min_valid,
		)
	else:
		fn = lambda x: calculate_cdf_percent(
			x, percentile, no_data=no_data, min_valid=min_valid,
		)
	result = generic_filter(
		padded_image,
		fn,
		size=(2 * dl + 1, 2 * dl + 1),
		mode='constant',
		cval=no_data,
	)
	# Crop the padded area off the result
	result = result[dl:-dl, dl:-dl]
	return _apply_additive_bias(result, bias, no_data, stat)


def diff_2_mask_quality(args):
	# compat: 5 (ancien) / 7 (stat+mad_k) / 8 (+bias) / 9 (+min_valid)
	if len(args) == 5:
		chem_in, chem_out, dl, no_data, percentile = args
		stat, mad_k, bias = "percentile", DEFAULT_MAD_K, 0.0
		min_valid = DEFAULT_MIN_VALID
	elif len(args) == 7:
		chem_in, chem_out, dl, no_data, percentile, stat, mad_k = args
		bias = 0.0
		min_valid = DEFAULT_MIN_VALID
	elif len(args) == 8:
		chem_in, chem_out, dl, no_data, percentile, stat, mad_k, bias = args
		min_valid = DEFAULT_MIN_VALID
	else:
		chem_in, chem_out, dl, no_data, percentile, stat, mad_k, bias, min_valid = args
	data_in = read_as_2D_float(chem_in, no_data)
	result = process_image(
		data_in, dl, no_data, percentile, stat=stat, mad_k=mad_k, bias=bias,
		min_valid=min_valid,
	)
	save_ABSOLUTE_image_with_same_geometry(
		result, chem_out, chem_in, no_data=no_data
	)
	return


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
