#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Lecture / écriture rasters et utilitaires I/O."""
import os
import shutil

import numpy as np
import rasterio


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
