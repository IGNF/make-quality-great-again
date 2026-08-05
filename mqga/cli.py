#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Point d'entrée CLI du pipeline MQGA."""
import argparse
import os
import sys
import time
import traceback

from loguru import logger

from mqga.io_raster import GetInfo, init_rep_tra, compute_mns_mnt_diff
from mqga.tiling import CalculNombreDallesXY, MakeDecoupage, DoParallel
from mqga.mosaic import Make_Assemblage_FINAL
from mqga.interpolate import (
	interpolate_nodata_idw_vectorized,
	interpolate_nodata_hybrid,
	apply_moving_average,
)

LOG_FORMAT = "{time:YYYY-MM-DD HH:mm:ss} | {level:<7} | {message}"


def _build_parser():
	parser = argparse.ArgumentParser(
		description="""MAKE QUALITY GREAT AGAIN - Masque de qualité à partir d'un MNS et d'un MNT
Auteur: Nicolas Champion""",
		epilog="""EXEMPLE:

  python3 make_quality_great_again.py \\
    --mns /chemin/mns.tif \\
    --mnt /chemin/mnt.tif \\
    --out /chemin/masque_qualite.tif \\
    --RepTra /chemin/tmp_mqga \\
    --cpu 8 \\
    --interp hybrid \\
    --clean

Le MNT est rééchantillonné sur la grille du MNS si les résolutions diffèrent.
Un fichier de log est écrit à côté de --out (même nom, extension .log).
""",
		formatter_class=argparse.RawTextHelpFormatter,
	)
	parser.add_argument("--mns", required=True, type=str, help="Chemin vers le MNS (DSM)")
	parser.add_argument("--mnt", required=True, type=str, help="Chemin vers le MNT (DTM)")
	parser.add_argument("--out", "-out", required=True, type=str, help="Masque de qualité en sortie")
	parser.add_argument("--no", "-no", type=int, default=-9999, help="Valeur de NoData")
	parser.add_argument("--per", "-per", type=float, default=0.05, help="Percentile local (défaut: 0.05)")
	parser.add_argument(
		"--demiwinl", "-demiwinl", type=int, default=50,
		help="Demie-taille en ligne de la fenêtre d'analyse",
	)
	parser.add_argument(
		"--demiwinc", "-demiwinc", type=int, default=50,
		help="Demie-taille en colonne de la fenêtre d'analyse",
	)
	parser.add_argument("--tile", "-tile", type=int, default=500, help="Taille de la tuile")
	parser.add_argument("--pad", "-pad", type=int, default=50, help="Recouvrement entre tuiles")
	parser.add_argument("--RepTra", "-RepTra", required=True, type=str, help="Répertoire de travail")
	parser.add_argument("--cpu", "-cpu", required=True, type=int, help="Nombre de CPU disponibles")
	parser.add_argument(
		"--winavg", "-winavg", type=int, default=50,
		help="Taille de la fenêtre glissante pour la moyenne (défaut: 50)",
	)
	parser.add_argument(
		"--interp", "-interp", type=str, default="hybrid",
		choices=["hybrid", "idw"],
		help="Méthode d'interpolation des NoData: hybrid (défaut) ou idw",
	)
	parser.add_argument(
		"--clean", "-clean", action="store_true",
		help="Supprimer le contenu du répertoire temporaire s'il existe déjà",
	)
	parser.add_argument(
		"--verbose", action="store_true",
		help="Activer le niveau DEBUG sur la console",
	)
	return parser


def _setup_logging(out_path, verbose=False):
	"""Configure loguru: console + fichier à côté de --out."""
	logger.remove()
	console_level = "DEBUG" if verbose else "INFO"
	logger.add(sys.stderr, level=console_level, format=LOG_FORMAT)

	if out_path.lower().endswith(".tif"):
		log_file = out_path[:-4] + ".log"
	else:
		log_file = out_path + ".log"

	log_dir = os.path.dirname(os.path.abspath(log_file))
	if log_dir:
		os.makedirs(log_dir, exist_ok=True)

	logger.add(log_file, level="DEBUG", format=LOG_FORMAT)
	return log_file


def _validate_args(args):
	errors = []

	if not os.path.isfile(args.mns):
		errors.append(f"Fichier MNS introuvable: {args.mns}")
	if not os.path.isfile(args.mnt):
		errors.append(f"Fichier MNT introuvable: {args.mnt}")

	if not (0.0 < args.per <= 1.0):
		errors.append(f"--per doit être dans ]0, 1], reçu: {args.per}")
	if args.tile <= 0:
		errors.append(f"--tile doit être > 0, reçu: {args.tile}")
	if args.pad < 0:
		errors.append(f"--pad doit être >= 0, reçu: {args.pad}")
	if args.pad >= args.tile:
		errors.append(f"--pad ({args.pad}) doit être strictement inférieur à --tile ({args.tile})")
	if args.cpu < 1:
		errors.append(f"--cpu doit être >= 1, reçu: {args.cpu}")
	if args.demiwinl < 1 or args.demiwinc < 1:
		errors.append("--demiwinl et --demiwinc doivent être >= 1")
	if args.winavg < 1:
		errors.append(f"--winavg doit être >= 1, reçu: {args.winavg}")

	out_dir = os.path.dirname(os.path.abspath(args.out))
	if out_dir and not os.path.isdir(out_dir):
		errors.append(f"Répertoire de sortie inexistant: {out_dir}")

	if errors:
		raise ValueError("\n".join(errors))


def _run_interpolation(interp_method, chem_out_tmp, chem_out_clean_bouchage, no_data, iNbreCPU):
	logger.info("Post-traitement: interpolation des pixels nodata (méthode: {})...", interp_method)

	if interp_method == "idw":
		interpolate_nodata_idw_vectorized(
			chem_out_tmp, chem_out_clean_bouchage, no_data, n_jobs=iNbreCPU
		)
	else:
		# hybrid (défaut)
		interpolate_nodata_hybrid(
			chem_out_tmp, chem_out_clean_bouchage, no_data,
			connectivity=4, seuil_percent=50, poids=1, rayon=50, n=1,
		)


def main(argv=None):
	if argv is None:
		argv = sys.argv[1:]

	try:
		start_time = time.time()
		args = _build_parser().parse_args(argv)

		log_file = _setup_logging(args.out, verbose=args.verbose)
		_validate_args(args)

		logger.info("Programme démarré")
		logger.info("Log fichier: {}", log_file)
		logger.info("MNS: {}", args.mns)
		logger.info("MNT: {}", args.mnt)
		logger.info("OUT: {}", args.out)
		logger.debug(
			"Params: per={}, tile={}, pad={}, cpu={}, interp={}, demiwinl={}, winavg={}",
			args.per, args.tile, args.pad, args.cpu, args.interp, args.demiwinl, args.winavg,
		)

		chem_mns = args.mns
		chem_mnt = args.mnt
		chem_out = args.out
		no_data = args.no
		percentile = args.per
		RepTra = args.RepTra
		iNbreCPU = args.cpu
		iTailleparcelle = args.tile
		iTailleRecouvrement = args.pad
		winavg_size = args.winavg
		dl = args.demiwinl

		init_rep_tra(RepTra, clean=args.clean)

		# Différence MNS - MNT (grille MNS)
		chem_in = os.path.join(RepTra, "diff_MNS_MNT.tif")
		logger.info("Calcul de la différence MNS - MNT (grille MNS)...")
		compute_mns_mnt_diff(chem_mns, chem_mnt, chem_in, no_data=no_data)
		logger.info("Différence écrite: {}", chem_in)

		infos = GetInfo(chem_in)
		tiles = CalculNombreDallesXY(
			infos.nbre_col, infos.nbre_lig, iTailleparcelle, iTailleRecouvrement
		)
		logger.info(
			"Tuiles: {}×{} (image {}×{})",
			tiles.nbre_dalle_x, tiles.nbre_dalle_y,
			int(infos.nbre_col), int(infos.nbre_lig),
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
		logger.info("Fichier final généré: {}", chem_out)

		elapsed_seconds = int(round(time.time() - start_time))
		elapsed_minutes = elapsed_seconds // 60
		elapsed_secs = elapsed_seconds % 60
		if elapsed_minutes > 0:
			logger.info(
				"Temps total d'exécution: {} minute(s) et {} seconde(s)",
				elapsed_minutes, elapsed_secs,
			)
		else:
			logger.info("Temps total d'exécution: {} seconde(s)", elapsed_seconds)

		return 0

	except SystemExit:
		raise
	except Exception as exc:
		logger.error("{}", exc)
		logger.debug("{}", traceback.format_exc())
		return 1


if __name__ == "__main__":
	sys.exit(main())
