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
	DEFAULT_HOLE_ALPHA,
	DEFAULT_HOLE_LAMBDA,
	interpolate_nodata_idw_vectorized,
	interpolate_nodata_hybrid,
	apply_moving_average,
)
from mqga.decrochage import (
	DEFAULT_MORPH_RADIUS,
	DEFAULT_SEUIL_STD,
	DEFAULT_SEUIL_V,
	DEFAULT_SEUIL_Z,
	apply_mask_as_nodata,
	detect_decrochage,
)
from mqga.quality_mask import (
	DEFAULT_MAD_K,
	DEFAULT_MIN_VALID,
	DEFAULT_MIN_VALID_PCT,
	resolve_min_valid,
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
    --decrochage \\
    --clean

Le MNT est rééchantillonné sur la grille du MNS si les résolutions diffèrent.
Avec --decrochage, les zones MNS ≫ MNT suspectes sont mises en NoData
dans la carte de précision et exportées en shapefile à côté de --out.
Un fichier de log est écrit à côté de --out (même nom, extension .log).
""",
		formatter_class=argparse.RawTextHelpFormatter,
	)
	parser.add_argument("--mns", required=True, type=str, help="Chemin vers le MNS (DSM)")
	parser.add_argument("--mnt", required=True, type=str, help="Chemin vers le MNT (DTM)")
	parser.add_argument("--out", "-out", required=True, type=str, help="Masque de qualité en sortie")
	parser.add_argument("--no", "-no", type=int, default=-9999, help="Valeur de NoData")
	parser.add_argument(
		"--per", "-per", type=float, default=0.05,
		help="Percentile local si --stat percentile (défaut: 0.05)",
	)
	parser.add_argument(
		"--stat", type=str, default="mad", choices=["mad", "percentile"],
		help="Statistique locale sur la partie négative: mad (défaut) ou percentile",
	)
	parser.add_argument(
		"--mad-k", type=float, default=DEFAULT_MAD_K,
		help=f"Facteur k dans ε = |b| + k·MAD (défaut: {DEFAULT_MAD_K} ≈ LE90 gaussienne)",
	)
	parser.add_argument(
		"--bias", type=float, default=0.0,
		help="Biais systématique |b| ajouté à ε (m), défaut: 0 → ε = |b| + k·MAD",
	)
	parser.add_argument(
		"--min-valid", type=int, default=DEFAULT_MIN_VALID,
		help=(
			f"Plancher absolu d'effectif (STANAG 2215 / historique). "
			f"Seuil effectif = max(min-valid, min-valid-pct%% de la fenêtre) "
			f"[default: {DEFAULT_MIN_VALID}]"
		),
	)
	parser.add_argument(
		"--min-valid-pct", type=float, default=DEFAULT_MIN_VALID_PCT,
		help=(
			f"Taux minimal (%%) de pixels négatifs valides dans la fenêtre "
			f"(0 = désactive la contrainte relative) [default: {DEFAULT_MIN_VALID_PCT}]"
		),
	)
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
		"--hole-vcalc", type=str, default="p90",
		choices=["min", "p90"],
		help="Constante de trou en hybrid: p90=P90(ε_bord) (défaut) ou min (historique)",
	)
	parser.add_argument(
		"--hole-alpha", type=float, default=DEFAULT_HOLE_ALPHA,
		help=(
			f"V2a hybrid: pénalité α·d·Δ (m d'incertitude par m de distance au bord). "
			f"0 = rampe off (V1) [default: {DEFAULT_HOLE_ALPHA}]"
		),
	)
	parser.add_argument(
		"--hole-lambda", type=float, default=DEFAULT_HOLE_LAMBDA,
		help=(
			f"V2a hybrid: plafond ε ≤ λ·V_calc (ignoré si --hole-alpha=0) "
			f"[default: {DEFAULT_HOLE_LAMBDA}]"
		),
	)
	parser.add_argument(
		"--clean", "-clean", action="store_true",
		help="Supprimer le contenu du répertoire temporaire s'il existe déjà",
	)
	parser.add_argument(
		"--verbose", action="store_true",
		help="Activer le niveau DEBUG sur la console",
	)
	parser.add_argument(
		"--decrochage", action="store_true",
		help="Détecter les zones de décrochage MNT (MNS ≫ MNT): NoData + shapefile",
	)
	parser.add_argument(
		"--seuilZ", type=float, default=DEFAULT_SEUIL_Z,
		help=f"Seuil Z positif sur MNS-MNT pour le décrochage (m) [default: {DEFAULT_SEUIL_Z}]",
	)
	parser.add_argument(
		"--seuilV", type=float, default=DEFAULT_SEUIL_V,
		help=f"Seuil de volume (COUNT*MEAN) des zones de décrochage [default: {DEFAULT_SEUIL_V}]",
	)
	parser.add_argument(
		"--seuilSTD", type=float, default=DEFAULT_SEUIL_STD,
		help=f"Seuil d'écart-type des zones de décrochage (m) [default: {DEFAULT_SEUIL_STD}]",
	)
	parser.add_argument(
		"--morph-radius", type=int, default=DEFAULT_MORPH_RADIUS,
		help=f"Rayon (pixels) de l'opening morphologique décrochage [default: {DEFAULT_MORPH_RADIUS}]",
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
	if args.mad_k <= 0:
		errors.append(f"--mad-k doit être > 0, reçu: {args.mad_k}")
	if args.bias < 0:
		errors.append(f"--bias doit être >= 0, reçu: {args.bias}")
	if args.min_valid < 1:
		errors.append(f"--min-valid doit être >= 1, reçu: {args.min_valid}")
	if args.min_valid_pct < 0 or args.min_valid_pct > 100:
		errors.append(
			f"--min-valid-pct doit être dans [0, 100], reçu: {args.min_valid_pct}"
		)
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
	if args.seuilZ <= 0:
		errors.append(f"--seuilZ doit être > 0, reçu: {args.seuilZ}")
	if args.seuilV < 0:
		errors.append(f"--seuilV doit être >= 0, reçu: {args.seuilV}")
	if args.seuilSTD < 0:
		errors.append(f"--seuilSTD doit être >= 0, reçu: {args.seuilSTD}")
	if args.morph_radius < 0:
		errors.append(f"--morph-radius doit être >= 0, reçu: {args.morph_radius}")
	if args.hole_alpha < 0:
		errors.append(f"--hole-alpha doit être >= 0, reçu: {args.hole_alpha}")
	if args.hole_lambda <= 0:
		errors.append(f"--hole-lambda doit être > 0, reçu: {args.hole_lambda}")

	out_dir = os.path.dirname(os.path.abspath(args.out))
	if out_dir and not os.path.isdir(out_dir):
		errors.append(f"Répertoire de sortie inexistant: {out_dir}")

	if errors:
		raise ValueError("\n".join(errors))


def _run_interpolation(
	interp_method, chem_out_tmp, chem_out_clean_bouchage, no_data, iNbreCPU,
	protect_mask_path=None,
	hole_vcalc="p90",
	hole_alpha=DEFAULT_HOLE_ALPHA,
	hole_lambda=DEFAULT_HOLE_LAMBDA,
):
	logger.info("Post-traitement: interpolation des pixels nodata (méthode: {})...", interp_method)

	if interp_method == "idw":
		interpolate_nodata_idw_vectorized(
			chem_out_tmp, chem_out_clean_bouchage, no_data, n_jobs=iNbreCPU,
			protect_mask_path=protect_mask_path,
		)
	else:
		# hybrid (défaut)
		interpolate_nodata_hybrid(
			chem_out_tmp, chem_out_clean_bouchage, no_data,
			connectivity=4, seuil_percent=50, poids=1, rayon=50, n=1,
			protect_mask_path=protect_mask_path,
			vcalc_mode=hole_vcalc,
			hole_alpha=hole_alpha,
			hole_lambda=hole_lambda,
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
		n_required, win_side, win_area = resolve_min_valid(
			args.demiwinl, min_valid=args.min_valid, min_valid_pct=args.min_valid_pct,
		)
		logger.debug(
			"Params: stat={}, per={}, mad_k={}, bias={}, min_valid={}, min_valid_pct={}, "
			"tile={}, pad={}, cpu={}, interp={}, demiwinl={}, winavg={}, decrochage={}, "
			"seuilZ={}, seuilV={}, seuilSTD={}, morph_radius={}",
			args.stat, args.per, args.mad_k, args.bias, args.min_valid, args.min_valid_pct,
			args.tile, args.pad, args.cpu, args.interp, args.demiwinl, args.winavg,
			args.decrochage, args.seuilZ, args.seuilV, args.seuilSTD, args.morph_radius,
		)
		logger.info(
			"Statistique qualité: {} (mad-k={}, bias={}, min-valid={}, min-valid-pct={}%)",
			args.stat, args.mad_k, args.bias, args.min_valid, args.min_valid_pct,
		)
		logger.info(
			"Seuil effectif min-valid: {} (fenêtre {}x{}={}, abs={}, pct={}%)",
			n_required, win_side, win_side, win_area, args.min_valid, args.min_valid_pct,
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

		protect_mask_path = None
		if args.decrochage:
			chem_mask_dec = os.path.join(RepTra, "mask_decrochage.tif")
			if chem_out.lower().endswith(".tif"):
				chem_shp_dec = chem_out[:-4] + "_zones_decrochage.shp"
			else:
				chem_shp_dec = chem_out + "_zones_decrochage.shp"
			logger.info("Détection des zones de décrochage MNT (MNS ≫ MNT)...")
			dec_res = detect_decrochage(
				chem_in,
				chem_mask_dec,
				chem_shp_dec,
				no_data=no_data,
				seuil_z=args.seuilZ,
				seuil_v=args.seuilV,
				seuil_std=args.seuilSTD,
				morph_radius=args.morph_radius,
			)
			protect_mask_path = dec_res.mask_path
			if dec_res.n_zones > 0:
				apply_mask_as_nodata(chem_in, protect_mask_path, chem_out=chem_in, no_data=no_data)
				logger.info(
					"Diff masquée sur {} zone(s) de décrochage ({} pixels)",
					dec_res.n_zones, dec_res.n_pixels,
				)
			else:
				logger.info("Aucune zone de décrochage retenue — diff inchangée")

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
			stat=args.stat, mad_k=args.mad_k, bias=args.bias,
			min_valid=args.min_valid, min_valid_pct=args.min_valid_pct,
		)

		chem_out_tmp = chem_out.replace('.tif', '_tmp.tif')
		Make_Assemblage_FINAL(chem_out_tmp, tiles.nbre_dalle_x, tiles.nbre_dalle_y, RepTra)

		if protect_mask_path:
			apply_mask_as_nodata(
				chem_out_tmp, protect_mask_path, chem_out=chem_out_tmp, no_data=no_data
			)

		chem_out_clean_bouchage = chem_out.replace('.tif', '_clean_bouchage.tif')
		_run_interpolation(
			args.interp, chem_out_tmp, chem_out_clean_bouchage, no_data, iNbreCPU,
			protect_mask_path=protect_mask_path,
			hole_vcalc=args.hole_vcalc,
			hole_alpha=args.hole_alpha,
			hole_lambda=args.hole_lambda,
		)

		apply_moving_average(
			chem_out_clean_bouchage, chem_out, window_size=winavg_size, no_data=no_data,
			protect_mask_path=protect_mask_path,
		)
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
