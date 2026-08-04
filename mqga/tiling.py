#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Découpage en tuiles et traitement parallèle."""
import os
import signal
from multiprocessing import Pool

from tqdm import tqdm

from mqga.io_raster import crop_tile
from mqga.quality_mask import create_negative_image, diff_2_mask_quality


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
