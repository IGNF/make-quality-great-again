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
		description="""MAKE QUALITY GREAT AGAIN - Module open source d'autoqualification de MNT (GEMAUT) - Délivre une carte de précision associée au MNT 
Auteur: Nicolas Champion""",
		epilog="""EXEMPLE:

  python3 make_quality_great_again.py \\
    --mnt /chemin/mnt.tif \\
    --mns /chemin/mns.tif \\
    --out /chemin/masque_qualite.tif \\
    --RepTra /chemin/tmp_mqga \\
    --cpu 8 \\
    --interp hybrid \\
    --decrochage \\
    --clean

Par défaut la carte est calculée sur la grille du MNS.
Avec --reso, MNS et MNT sont rééchantillonnés sur une grille de travail (m).
Avec --decrochage, les zones MNS ≫ MNT suspectes sont mises en NoData
dans la carte de précision et exportées en shapefile à côté de --out.
Un fichier de log est écrit à côté de --out (même nom, extension .log).
""",
		formatter_class=argparse.RawTextHelpFormatter,
	)

	# --- Options obligatoires ---
	g_req = parser.add_argument_group("Options obligatoires")
	g_req.add_argument("--mnt", required=True, type=str, help="Chemin vers le MNT (DTM)")
	g_req.add_argument("--mns", required=True, type=str, help="Chemin vers le MNS (DSM)")
	g_req.add_argument("--out", "-out", required=True, type=str, help="Masque de qualité en sortie")
	g_req.add_argument(
		"--RepTra", "-RepTra", required=True, type=str,
		help="Répertoire de travail (fichiers temporaires)",
	)
	g_req.add_argument("--cpu", "-cpu", required=True, type=int, help="Nombre de CPU disponibles")

	# --- Options générales ---
	g_gen = parser.add_argument_group("Options générales")
	g_gen.add_argument(
		"--reso", type=float, default=None,
		help=(
			"Résolution de travail en mètres (ex. 1 ou 4).\n"
			"Si absent = Résolution du MNS en entrée."
		),
	)
	g_gen.add_argument("--no", "-no", type=int, default=-9999, help="Valeur de NoData (défaut: -9999)")
	g_gen.add_argument(
		"--clean", "-clean", action="store_true",
		help="Vider le répertoire temporaire s'il existe déjà",
	)
	g_gen.add_argument(
		"--verbose", action="store_true",
		help="Activer le niveau DEBUG sur la console",
	)

	# --- Analyse locale MNS-MNT ---
	g_stat = parser.add_argument_group(
		"Analyse MNS-MNT (partie négative) - calcul de ε_stat = statistique locale"
	)
	g_stat.add_argument(
		"--stat", type=str, default="percentile", choices=["mad", "percentile"],
		help=(
			"Statistique locale sur la partie négative de MNS-MNT:\n"
			"percentile (défaut) ou mad"
		),
	)
	g_stat.add_argument(
		"--per", "-per", type=float, default=0.10,
		help="Percentile local si --stat percentile (défaut: 0.10 ; ignoré si --stat mad)",
	)
	g_stat.add_argument(
		"--mad-k", type=float, default=DEFAULT_MAD_K,
		help=(
			f"Facteur k dans ε = |b| + k·MAD si --stat mad "
			f"(défaut: {DEFAULT_MAD_K} ≈ LE90 ; ignoré si --stat percentile)"
		),
	)
	g_stat.add_argument(
		"--bias", type=float, default=0.0,
		help="Biais |b| ajouté à ε (m) : ε ← |b| + ε_stat (défaut: 0)",
	)
	g_stat.add_argument(
		"--demiwin", "-demiwin", type=int, default=50,
		help="Demi-taille (pixels) de la fenêtre d'analyse (défaut: 50 → fenêtre ~101×101)",
	)
	g_stat.add_argument(
		"--min-valid", type=int, default=DEFAULT_MIN_VALID,
		help=(
			f"Plancher absolu d'effectif (STANAG).\n"
			f"Seuil effectif = max(min-valid, min-valid-pct%% de la fenêtre) "
			f"[défaut: {DEFAULT_MIN_VALID}]"
		),
	)
	g_stat.add_argument(
		"--min-valid-pct", type=float, default=DEFAULT_MIN_VALID_PCT,
		help=(
			f"Taux minimal (%%) de pixels négatifs (MNS-MNT) valides dans la fenêtre.\n"
			f"0 = désactive la contrainte relative [défaut: {DEFAULT_MIN_VALID_PCT}]"
		),
	)

	# --- Tuiles ---
	g_tile = parser.add_argument_group("Découpage en tuiles")
	g_tile.add_argument(
		"--tile", "-tile", type=int, default=500,
		help="Taille d'une tuile en pixels (défaut: 500)",
	)
	g_tile.add_argument(
		"--pad", "-pad", type=int, default=50,
		help="Recouvrement entre tuiles en pixels (défaut: 50 ; doit être >= --demiwin)",
	)

	# --- Interpolation ---
	g_interp = parser.add_argument_group("Interpolation des NoData (trous)")
	g_interp.add_argument(
		"--interp", "-interp", type=str, default="hybrid",
		choices=["hybrid", "idw"],
		help="Méthode d'interpolation: hybrid (défaut) ou idw",
	)
	g_interp.add_argument(
		"--hole-vcalc", type=str, default="p90",
		choices=["min", "p90"],
		help="Ancre de trou en hybrid: p90=P90(ε_bord) (défaut) ou min (historique)",
	)
	g_interp.add_argument(
		"--hole-alpha", type=float, default=DEFAULT_HOLE_ALPHA,
		help=(
			f"V2a: pénalité α·d·Δ (m d'incertitude par m de distance au bord).\n"
			f"0 = rampe off (comportement V1) [défaut: {DEFAULT_HOLE_ALPHA}]"
		),
	)
	g_interp.add_argument(
		"--hole-lambda", type=float, default=DEFAULT_HOLE_LAMBDA,
		help=(
			f"V2a: plafond ε ≤ λ·V_calc (ignoré si --hole-alpha=0) "
			f"[défaut: {DEFAULT_HOLE_LAMBDA}]"
		),
	)

	# --- Lissage ---
	g_smooth = parser.add_argument_group("Lissage final")
	g_smooth.add_argument(
		"--winavg", "-winavg", type=int, default=50,
		help="Taille (pixels) de la moyenne glissante finale (défaut: 50)",
	)

	# --- Décrochage ---
	g_dec = parser.add_argument_group("Détection des zones de décrochage MNT")
	g_dec.add_argument(
		"--decrochage", action="store_true",
		help="Activer la détection (MNS ≫ MNT): NoData + shapefile à côté de --out",
	)
	g_dec.add_argument(
		"--seuilZ", type=float, default=DEFAULT_SEUIL_Z,
		help=f"Seuil Z positif sur MNS-MNT (m) [défaut: {DEFAULT_SEUIL_Z}]",
	)
	g_dec.add_argument(
		"--seuilV", type=float, default=DEFAULT_SEUIL_V,
		help=f"Seuil de volume COUNT×MEAN [défaut: {DEFAULT_SEUIL_V}]",
	)
	g_dec.add_argument(
		"--seuilSTD", type=float, default=DEFAULT_SEUIL_STD,
		help=f"Seuil d'écart-type (m) [défaut: {DEFAULT_SEUIL_STD}]",
	)
	g_dec.add_argument(
		"--morph-radius", type=int, default=DEFAULT_MORPH_RADIUS,
		help=f"Rayon (pixels) de l'opening morphologique [défaut: {DEFAULT_MORPH_RADIUS}]",
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
	if args.reso is not None and args.reso <= 0:
		errors.append(f"--reso doit être > 0, reçu: {args.reso}")

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
	if args.demiwin < 1:
		errors.append("--demiwin doit être >= 1")
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
			args.demiwin, min_valid=args.min_valid, min_valid_pct=args.min_valid_pct,
		)
		logger.debug(
			"Params: stat={}, per={}, mad_k={}, bias={}, min_valid={}, min_valid_pct={}, "
			"tile={}, pad={}, cpu={}, interp={}, demiwin={}, winavg={}, decrochage={}, "
			"seuilZ={}, seuilV={}, seuilSTD={}, morph_radius={}",
			args.stat, args.per, args.mad_k, args.bias, args.min_valid, args.min_valid_pct,
			args.tile, args.pad, args.cpu, args.interp, args.demiwin, args.winavg,
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
		dl = args.demiwin

		init_rep_tra(RepTra, clean=args.clean)

		# Différence MNS - MNT (grille MNS, ou grille --reso)
		chem_in = os.path.join(RepTra, "diff_MNS_MNT.tif")
		if args.reso is not None:
			logger.info(
				"Calcul de la différence MNS - MNT (grille de travail {} m)...",
				args.reso,
			)
		else:
			logger.info("Calcul de la différence MNS - MNT (grille MNS)...")
		compute_mns_mnt_diff(
			chem_mns, chem_mnt, chem_in, no_data=no_data, reso=args.reso,
		)
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
