#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Point d'entrée CLI du pipeline MQGA."""
import argparse
import sys
import time

from mqga.io_raster import GetInfo, init_rep_tra
from mqga.tiling import CalculNombreDallesXY, MakeDecoupage, DoParallel
from mqga.mosaic import Make_Assemblage_FINAL
from mqga.interpolate import (
    interpolate_nodata_griddata,
    interpolate_nodata_idw,
    interpolate_nodata_idw_vectorized,
    interpolate_nodata_window,
    interpolate_nodata_with_linearnd,
    interpolate_nodata_fast,
    interpolate_nodata_hybrid,
    apply_moving_average,
)


def main(argv=None):
	if argv is None:
		argv = sys.argv[1:]
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
		
		args = parser.parse_args(argv)
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



if __name__ == "__main__":
	main()
