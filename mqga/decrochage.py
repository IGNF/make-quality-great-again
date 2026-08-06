#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Détection des zones de décrochage MNT (MNS ≫ MNT)."""
from __future__ import annotations

import json
import os
from dataclasses import dataclass

import numpy as np
import rasterio
from loguru import logger
from osgeo import ogr, osr
from rasterio.features import shapes
from scipy.ndimage import binary_opening, label


# Seuils par défaut (alignés sur make_quality_great_again_SERVEUR, branche POS)
DEFAULT_SEUIL_Z = 10.0
DEFAULT_SEUIL_V = 100000.0
DEFAULT_SEUIL_STD = 3.0
DEFAULT_MORPH_RADIUS = 5


@dataclass(frozen=True)
class DecrochageResult:
	"""Chemins et résumé de la détection."""
	mask_path: str
	shp_path: str
	n_zones: int
	n_pixels: int


def _box_structure(radius: int) -> np.ndarray:
	size = 2 * int(radius) + 1
	return np.ones((size, size), dtype=bool)


def _write_empty_shapefile(path_shp: str, crs) -> None:
	"""Crée un shapefile polygone vide (même schéma que les alertes)."""
	_ensure_parent(path_shp)
	driver = ogr.GetDriverByName("ESRI Shapefile")
	if os.path.exists(path_shp):
		driver.DeleteDataSource(path_shp)

	ds = driver.CreateDataSource(path_shp)
	srs = osr.SpatialReference()
	if crs is not None:
		srs.ImportFromWkt(crs.to_wkt())
	layer = ds.CreateLayer("decrochage", srs, ogr.wkbPolygon)
	layer.CreateField(ogr.FieldDefn("ID", ogr.OFTInteger))
	layer.CreateField(ogr.FieldDefn("COUNT", ogr.OFTInteger))
	for name in ("MEAN", "STD", "VOLUME"):
		field = ogr.FieldDefn(name, ogr.OFTReal)
		field.SetWidth(24)
		field.SetPrecision(8)
		layer.CreateField(field)
	ds = None


def _write_alert_shapefile(path_shp: str, labeled: np.ndarray, stats: dict, transform, crs) -> None:
	"""Écrit les polygones d'alerte avec attributs COUNT/MEAN/STD/VOLUME."""
	_ensure_parent(path_shp)
	driver = ogr.GetDriverByName("ESRI Shapefile")
	if os.path.exists(path_shp):
		driver.DeleteDataSource(path_shp)

	ds = driver.CreateDataSource(path_shp)
	srs = osr.SpatialReference()
	if crs is not None:
		srs.ImportFromWkt(crs.to_wkt())
	layer = ds.CreateLayer("decrochage", srs, ogr.wkbPolygon)
	layer.CreateField(ogr.FieldDefn("ID", ogr.OFTInteger))
	layer.CreateField(ogr.FieldDefn("COUNT", ogr.OFTInteger))
	for name in ("MEAN", "STD", "VOLUME"):
		field = ogr.FieldDefn(name, ogr.OFTReal)
		field.SetWidth(24)
		field.SetPrecision(8)
		layer.CreateField(field)

	mask = labeled > 0
	if not np.any(mask):
		ds = None
		return

	for geom_dict, value in shapes(
		labeled.astype(np.int32),
		mask=mask,
		transform=transform,
	):
		label_id = int(value)
		if label_id <= 0 or label_id not in stats:
			continue
		geom = ogr.CreateGeometryFromJson(json.dumps(geom_dict))
		if geom is None or geom.IsEmpty():
			continue
		feat = ogr.Feature(layer.GetLayerDefn())
		feat.SetGeometry(geom)
		st = stats[label_id]
		feat.SetField("ID", label_id)
		feat.SetField("COUNT", int(st["count"]))
		feat.SetField("MEAN", float(st["mean"]))
		feat.SetField("STD", float(st["std"]))
		feat.SetField("VOLUME", float(st["volume"]))
		layer.CreateFeature(feat)
		feat = None

	ds = None


def _ensure_parent(path: str) -> None:
	parent = os.path.dirname(os.path.abspath(path))
	if parent:
		os.makedirs(parent, exist_ok=True)


def detect_decrochage(
	chem_diff: str,
	chem_mask_out: str,
	chem_shp_out: str,
	no_data: float = -9999,
	seuil_z: float = DEFAULT_SEUIL_Z,
	seuil_v: float = DEFAULT_SEUIL_V,
	seuil_std: float = DEFAULT_SEUIL_STD,
	morph_radius: int = DEFAULT_MORPH_RADIUS,
) -> DecrochageResult:
	"""
	Détecte les zones où le MNT décroche (MNS ≫ MNT).

	Pipeline :
	1. seuil diff > seuil_z
	2. opening morphologique (boîte)
	3. composantes connexes
	4. filtre VOLUME > seuil_v et STD > seuil_std
	   (VOLUME = COUNT * MEAN sur la diff, comme le script SERVEUR historique)

	Returns:
		DecrochageResult (masque raster + shapefile ; shp vide si aucune zone)
	"""
	_ensure_parent(chem_mask_out)
	_ensure_parent(chem_shp_out)

	with rasterio.open(chem_diff, "r") as src:
		diff = src.read(1).astype(np.float32)
		profile = src.profile.copy()
		transform = src.transform
		crs = src.crs

	valid = (diff != no_data) & ~np.isnan(diff)
	candidates = valid & (diff > seuil_z)

	if morph_radius > 0:
		candidates = binary_opening(candidates, structure=_box_structure(morph_radius))

	labeled, n_comp = label(candidates)
	logger.info(
		"Décrochage: {} composante(s) après seuil Z={} m et morpho radius={}",
		n_comp, seuil_z, morph_radius,
	)

	alert = np.zeros(diff.shape, dtype=bool)
	alert_labeled = np.zeros(diff.shape, dtype=np.int32)
	stats: dict[int, dict] = {}
	next_id = 0

	for comp_id in range(1, n_comp + 1):
		comp_mask = labeled == comp_id
		vals = diff[comp_mask]
		count = int(vals.size)
		if count == 0:
			continue
		mean = float(np.mean(vals))
		std = float(np.std(vals))
		volume = float(count * mean)
		if volume > seuil_v and std > seuil_std:
			next_id += 1
			alert[comp_mask] = True
			alert_labeled[comp_mask] = next_id
			stats[next_id] = {
				"count": count,
				"mean": mean,
				"std": std,
				"volume": volume,
			}

	n_zones = len(stats)
	n_pixels = int(np.count_nonzero(alert))
	logger.info(
		"Décrochage: {} zone(s) d'alerte retenues ({} pixels) "
		"[seuilV={}, seuilSTD={}]",
		n_zones, n_pixels, seuil_v, seuil_std,
	)

	profile.update(dtype="uint8", count=1, nodata=0, compress="lzw")
	with rasterio.open(chem_mask_out, "w", **profile) as dst:
		dst.write(alert.astype(np.uint8), 1)

	if n_zones == 0:
		_write_empty_shapefile(chem_shp_out, crs)
	else:
		_write_alert_shapefile(chem_shp_out, alert_labeled, stats, transform, crs)

	logger.info("Masque décrochage: {}", chem_mask_out)
	logger.info("Shapefile décrochage: {}", chem_shp_out)
	return DecrochageResult(
		mask_path=chem_mask_out,
		shp_path=chem_shp_out,
		n_zones=n_zones,
		n_pixels=n_pixels,
	)


def apply_mask_as_nodata(
	chem_raster: str,
	chem_mask: str,
	chem_out: str | None = None,
	no_data: float = -9999,
) -> str:
	"""
	Met à nodata les pixels où le masque d'alerte est non nul.

	Si chem_out est None, écrase chem_raster.
	"""
	out_path = chem_out or chem_raster
	with rasterio.open(chem_raster, "r") as src_r, rasterio.open(chem_mask, "r") as src_m:
		data = src_r.read(1).astype(np.float32)
		mask = src_m.read(1)
		profile = src_r.profile.copy()
		if mask.shape != data.shape:
			raise ValueError(
				f"Masque décrochage {mask.shape} incompatible avec raster {data.shape}"
			)
		data[mask != 0] = no_data
		profile.update(dtype="float32", count=1, nodata=no_data)

	_ensure_parent(out_path)
	with rasterio.open(out_path, "w", **profile) as dst:
		dst.write(data, 1)
	return out_path
