#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Point d'entrée CLI du pipeline MQGA."""
import argparse
import logging
import os
import sys
import time
import traceback

from mqga.io_raster import GetInfo, init_rep_tra, compute_mns_mnt_diff
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

logger = logging.getLogger("mqga")


def _build_parser():
	parser = argparse.ArgumentParser(
		description='MAKE QUALITY GREAT AGAIN ALL ZONE - Version Spatialisée & Parallélisée'
	)
	parser.add_argument("--mns", "-mns", required=True, type=str, help="Chemin vers le MNS (DSM)")
	parser.add_argument("--mnt", "-mnt", required=True, type=str, help="Chemin vers le MNT (DTM)")
	parser.add_argument("-out", required=True, type=str, help="Masque de Qualité en sortie")
	parser.add_argument("-no", type=int, default=-9999, help="Valeur de No Data")
	parser.add_argument("-per", type=float, default=0.05, help="Valeur de percentile")
	parser.add_argument("-demiwinl", type=int, default=50, help="Demie-taille en ligne de la fenêtre d'analyse")
	parser.add_argument("-demiwinc", type=int, default=50, help="Demie-taille en colonne de la fenêtre d'analyse")
	parser.add_argument("-tile", type=int, default=500, help="Tile / Taille de la tuile")
	parser.add_argument("-pad", type=int, default=50, help="Pad / Recouvrement entre tuiles")
	parser.add_argument("-RepTra", required=True, type=str, help="Répertoire de Travail")
	parser.add_argument("-cpu", required=True, type=int, help="Nombre de CPU disponibles")
	parser.add_argument(
		"-winavg", type=int, default=50,
		help="Taille de la fenêtre glissante pour la moyenne (par défaut 50x50)",
	)
	parser.add_argument(
		"-interp", type=str, default="idw",
		choices=["griddata", "idw", "idw_old", "window", "linearnd", "fast", "hybrid"],
		help="Méthode d'interpolation pour les pixels nodata (défaut: idw)",
	)
	parser.add_argument(
		"-clean", action='store_true',
		help="Supprimer le contenu du répertoire temporaire s'il existe déjà",
	)
	return parser


def _validate_args(args):
	errors = []

	if not os.path.isfile(args.mns):
		errors.append(f"Fichier MNS introuvable: {args.mns}")
	if not os.path.isfile(args.mnt):
		errors.append(f"Fichier MNT introuvable: {args.mnt}")

	if not (0.0 < args.per <= 1.0):
		errors.append(f"-per doit être dans ]0, 1], reçu: {args.per}")
	if args.tile <= 0:
		errors.append(f"-tile doit être > 0, reçu: {args.tile}")
	if args.pad < 0:
		errors.append(f"-pad doit être >= 0, reçu: {args.pad}")
	if args.pad >= args.tile:
		errors.append(f"-pad ({args.pad}) doit être strictement inférieur à -tile ({args.tile})")
	if args.cpu < 1:
		errors.append(f"-cpu doit être >= 1, reçu: {args.cpu}")
	if args.demiwinl < 1 or args.demiwinc < 1:
		errors.append("-demiwinl et -demiwinc doivent être >= 1")
	if args.winavg < 1:
		errors.append(f"-winavg doit être >= 1, reçu: {args.winavg}")

	out_dir = os.path.dirname(os.path.abspath(args.out))
	if out_dir and not os.path.isdir(out_dir):
		errors.append(f"Répertoire de sortie inexistant: {out_dir}")

	if errors:
		raise ValueError("\n".join(errors))


def _run_interpolation(interp_method, chem_out_tmp, chem_out_clean_bouchage, no_data, iNbreCPU):
	logger.info("Post-traitement: interpolation des pixels nodata (méthode: %s)...", interp_method)

	if interp_method == "griddata":
		interpolate_nodata_griddata(chem_out_tmp, chem_out_clean_bouchage, no_data)
	elif interp_method == "idw":
		interpolate_nodata_idw_vectorized(chem_out_tmp, chem_out_clean_bouchage, no_data, n_jobs=iNbreCPU)
	elif interp_method == "idw_old":
		interpolate_nodata_idw(chem_out_tmp, chem_out_clean_bouchage, no_data)
	elif interp_method == "window":
		interpolate_nodata_window(chem_out_tmp, chem_out_clean_bouchage, no_data)
	elif interp_method == "linearnd":
		interpolate_nodata_with_linearnd(chem_out_tmp, chem_out_clean_bouchage, no_data)
	elif interp_method == "fast":
		interpolate_nodata_fast(chem_out_tmp, chem_out_clean_bouchage, no_data)
	elif interp_method == "hybrid":
		interpolate_nodata_hybrid(
			chem_out_tmp, chem_out_clean_bouchage, no_data,
			connectivity=4, seuil_percent=50, poids=1, rayon=50, n=1,
		)
	else:
		logger.warning("Méthode '%s' non reconnue, utilisation de 'fast'.", interp_method)
		interpolate_nodata_fast(chem_out_tmp, chem_out_clean_bouchage, no_data)


def main(argv=None):
	if argv is None:
		argv = sys.argv[1:]

	logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

	try:
		start_time = time.time()
		args = _build_parser().parse_args(argv)
		_validate_args(args)

		chem_mns = args.mns
		chem_mnt = args.mnt
		chem_out = args.out
		no_data = args.no
		percentile = args.per
		demiwinsize_lig = args.demiwinl
		RepTra = args.RepTra
		iNbreCPU = args.cpu
		iTailleparcelle = args.tile
		iTailleRecouvrement = args.pad
		winavg_size = args.winavg
		dl = demiwinsize_lig

		init_rep_tra(RepTra, clean=args.clean)

		# Différence MNS - MNT (équivalent de l'ancien -diff)
		chem_in = os.path.join(RepTra, "diff_MNS_MNT.tif")
		logger.info("Calcul de la différence MNS - MNT (grille MNS)...")
		compute_mns_mnt_diff(chem_mns, chem_mnt, chem_in, no_data=no_data)
		logger.info("Différence écrite: %s", chem_in)

		infos = GetInfo(chem_in)
		tiles = CalculNombreDallesXY(
			infos.nbre_col, infos.nbre_lig, iTailleparcelle, iTailleRecouvrement
		)

		MakeDecoupage(
			chem_in, RepTra, tiles.nbre_dalle_x, tiles.nbre_dalle_y,
			iTailleparcelle, iTailleRecouvrement, iNbreCPU,
		)
		DoParallel(
			RepTra, tiles.nbre_dalle_x, tiles.nbre_dalle_y,
			dl, no_data, percentile, iNbreCPU,
		)

		chem_out_tmp = chem_out.replace('.tif', '_tmp.tif')
		Make_Assemblage_FINAL(chem_out_tmp, tiles.nbre_dalle_x, tiles.nbre_dalle_y, RepTra)

		chem_out_clean_bouchage = chem_out.replace('.tif', '_clean_bouchage.tif')
		_run_interpolation(args.interp, chem_out_tmp, chem_out_clean_bouchage, no_data, iNbreCPU)

		apply_moving_average(chem_out_clean_bouchage, chem_out, window_size=winavg_size, no_data=no_data)
		logger.info("Fichier final généré: %s", chem_out)

		elapsed_seconds = int(round(time.time() - start_time))
		elapsed_minutes = elapsed_seconds // 60
		elapsed_secs = elapsed_seconds % 60
		if elapsed_minutes > 0:
			logger.info(
				"Temps total d'exécution: %s minute(s) et %s seconde(s)",
				elapsed_minutes, elapsed_secs,
			)
		else:
			logger.info("Temps total d'exécution: %s seconde(s)", elapsed_seconds)

		return 0

	except SystemExit:
		raise
	except Exception as exc:
		logger.error("%s", exc)
		logger.debug(traceback.format_exc())
		# Afficher la traceback complète pour le debug opérationnel
		traceback.print_exc()
		return 1


if __name__ == "__main__":
	sys.exit(main())
