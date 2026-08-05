#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Lecture / écriture rasters et utilitaires I/O."""
import os
import shutil
from dataclasses import dataclass

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject


@dataclass(frozen=True)
class RasterInfo:
	"""Métadonnées utiles d'un GeoTIFF."""
	pas_x: float
	pas_y: float
	projection: int
	x_0: float  # gauche
	x_1: float  # droite
	y_0: float  # bas
	y_1: float  # haut
	phasage: str
	nbre_col: float
	nbre_lig: float
	g_model: int
	g_raster: int


def GetInfo(cheminTIF):
	"""
	Récupère les métadonnées d'une image GeoTIFF.
	Version open source utilisant rasterio.

	Returns:
		RasterInfo
	"""
	# Ouvrir l'image avec rasterio
	with rasterio.open(cheminTIF, 'r') as src:
		# Pas en X et Y (résolution)
		transform = src.transform
		pas_x = abs(transform[0])  # pixel width
		pas_y = abs(transform[4])  # pixel height (généralement négatif, on prend la valeur absolue)

		# Bounds (positions)
		bounds = src.bounds
		x_0 = bounds.left   # GAUCHE
		x_1 = bounds.right  # DROITE
		y_0 = bounds.bottom  # BAS
		y_1 = bounds.top     # HAUT

		# Nombre de colonnes et lignes
		nbre_col = float(src.width)
		nbre_lig = float(src.height)

		# Code EPSG / Projection
		if src.crs is not None:
			projection = int(src.crs.to_epsg()) if src.crs.to_epsg() is not None else -1
		else:
			projection = -1

		# GTModelTypeGeoKey et GTRasterTypeGeoKey
		# Ces clés sont dans les tags GeoTIFF, généralement dans les tags de la bande
		g_model = -1
		g_raster = -1
		# Essayer d'abord les tags du dataset
		if hasattr(src, 'tags') and src.tags():
			tags = src.tags()
			if 'GTModelTypeGeoKey' in tags:
				try:
					g_model = int(tags['GTModelTypeGeoKey'])
				except (ValueError, TypeError):
					pass
			if 'GTRasterTypeGeoKey' in tags:
				try:
					g_raster = int(tags['GTRasterTypeGeoKey'])
				except (ValueError, TypeError):
					pass
		# Si pas trouvé, essayer les tags de la première bande
		if (g_model == -1 or g_raster == -1) and hasattr(src, 'tags'):
			try:
				band_tags = src.tags(1)
				if band_tags:
					if 'GTModelTypeGeoKey' in band_tags and g_model == -1:
						try:
							g_model = int(band_tags['GTModelTypeGeoKey'])
						except (ValueError, TypeError):
							pass
					if 'GTRasterTypeGeoKey' in band_tags and g_raster == -1:
						try:
							g_raster = int(band_tags['GTRasterTypeGeoKey'])
						except (ValueError, TypeError):
							pass
			except Exception:
				pass

		# Phasage : déterminer si c'est "HG" (haut-gauche) ou "CP" (centre pixel)
		phasage = "HG"  # Par défaut
		if hasattr(src, 'tags') and src.tags():
			# Vérifier les tags pour le phasage
			tags = src.tags()
			if 'AREA_OR_POINT' in tags and tags['AREA_OR_POINT'] == 'Point':
				phasage = "CP"

	return RasterInfo(
		pas_x=pas_x,
		pas_y=pas_y,
		projection=projection,
		x_0=x_0,
		x_1=x_1,
		y_0=y_0,
		y_1=y_1,
		phasage=phasage,
		nbre_col=nbre_col,
		nbre_lig=nbre_lig,
		g_model=g_model,
		g_raster=g_raster,
	)

#############################################################################################################################	
def read_as_2D_float(filename,no_data):
	
	with rasterio.open(filename, 'r') as dataset:
		data= dataset.read(1).astype(np.float32)
		data[data==no_data] = np.nan
		
	return data
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


def _same_grid(src_a, src_b):
	"""True si les deux rasters partagent taille, transform et CRS."""
	return (
		src_a.width == src_b.width
		and src_a.height == src_b.height
		and src_a.transform == src_b.transform
		and src_a.crs == src_b.crs
	)


def compute_mns_mnt_diff(chem_mns, chem_mnt, chem_out, no_data=-9999):
	"""
	Calcule la différence MNS - MNT sur la grille du MNS.

	Si le MNT n'est pas sur la même grille (résolution / transform / CRS),
	il est rééchantillonné sur la grille du MNS (bilinéaire).
	Un pixel est nodata en sortie si MNS ou MNT (resamplé) est nodata / NaN.
	"""
	with rasterio.open(chem_mns, 'r') as src_mns, rasterio.open(chem_mnt, 'r') as src_mnt:
		if src_mns.crs is None:
			raise ValueError("Le MNS n'a pas de CRS défini")
		if src_mnt.crs is None:
			raise ValueError("Le MNT n'a pas de CRS défini")

		nodata_mns = src_mns.nodata if src_mns.nodata is not None else no_data
		nodata_mnt = src_mnt.nodata if src_mnt.nodata is not None else no_data
		aligned = _same_grid(src_mns, src_mnt)

		if not aligned:
			print(
				"INFO: MNT rééchantillonné sur la grille du MNS "
				f"(MNS={src_mns.width}x{src_mns.height}, "
				f"MNT={src_mnt.width}x{src_mnt.height})"
			)

		metadata = src_mns.meta.copy()
		metadata.update({
			'dtype': 'float32',
			'count': 1,
			'nodata': no_data,
			'compress': 'lzw',
		})

		with rasterio.open(chem_out, 'w', **metadata) as dst:
			for _, window in src_mns.block_windows(1):
				mns = src_mns.read(1, window=window).astype(np.float32)

				if aligned:
					mnt = src_mnt.read(1, window=window).astype(np.float32)
					mask_mnt_invalid = (mnt == nodata_mnt) | np.isnan(mnt)
				else:
					mnt = np.full(mns.shape, no_data, dtype=np.float32)
					win_transform = rasterio.windows.transform(window, src_mns.transform)
					reproject(
						source=rasterio.band(src_mnt, 1),
						destination=mnt,
						src_transform=src_mnt.transform,
						src_crs=src_mnt.crs,
						src_nodata=nodata_mnt,
						dst_transform=win_transform,
						dst_crs=src_mns.crs,
						dst_nodata=no_data,
						resampling=Resampling.bilinear,
					)
					mask_mnt_invalid = (mnt == no_data) | np.isnan(mnt)

				mask_invalid = (
					(mns == nodata_mns) | np.isnan(mns) | mask_mnt_invalid
				)
				diff = mns - mnt
				diff[mask_invalid] = no_data
				dst.write(diff.astype(np.float32), 1, window=window)

	return chem_out
