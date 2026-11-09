"""Does past clustering improve deforestation risk allocation? A simplified VT0007-style test.

VCS Tool VT0007 (Unplanned Deforestation Allocation, v1.0, 2024), used with VM0048,
builds a *benchmark* risk map from the distance to the forest edge: 30 vulnerability
classes in a geometric series up to the negligible risk threshold (NRT, the distance
containing 99.5 % of historical deforestation), intersected with administrative
divisions into "modelling regions". The relative frequency of deforestation in each
region during a fitting period becomes a density, which is applied to the forest of
the next period and scaled to the expected quantity. An alternative model is
accepted only if its median absolute error (MedAE) is lower than the benchmark's in
both tests.

This script reproduces that logic in a simplified, illustrative way on the 90 m
Hansen aggregates and tests one alternative: the same modelling regions further
split by the Getis-Ord Gi* hot/cold spot class of the fitting period.

  test 1  fit 2001-2008 (forest 2000)  ->  predict 2009-2016 (forest 2008)
  test 2  fit 2009-2016 (forest 2008)  ->  predict 2017-2025 (forest 2016)

Simplifications (this is NOT the official tool, see Clark Labs UDef-ARP):
90 m cells instead of 30 m, a cell is forest if >= 50 % of its 30 m pixels are
forest, distances on a locally metric grid, departments as administrative
divisions, quantity scaled to the observed loss (VT0007 "testing" mode) and
MedAE computed over the 5 km hexagons.

Usage:  python scripts/05_vt0007_benchmark.py --region amazon_arc
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from pyproj import Geod
from rasterio.features import rasterize
from scipy.ndimage import distance_transform_edt

N_CLASSES = 30
NRT_Q = 99.5
REGIONS = {"amazon_arc": ["Meta", "Guaviare", "Caquetá", "Putumayo"], "choco": ["Chocó"]}
ZONES = {"cold spot": 0, "not significant": 1, "hot spot": 2}
TESTS = [
    {"name": "test 1", "fit": "p1", "predict": "p2", "fit_label": "2001-2008", "pred_label": "2009-2016"},
    {"name": "test 2", "fit": "p2", "predict": "p3", "fit_label": "2009-2016", "pred_label": "2017-2025"},
]


def pixel_area_ha(transform, height: int) -> np.ndarray:
    geod = Geod(ellps="WGS84")
    res = abs(transform.a) / 3
    out = np.empty(height, dtype="float32")
    for r in range(height):
        c = transform.f + (r + 0.5) * transform.e
        out[r] = abs(geod.polygon_area_perimeter([0, res, res, 0], [c + res / 2, c + res / 2, c - res / 2, c - res / 2])[0]) / 1e4
    return out


def edge_distance(forest_cells: np.ndarray, transform) -> np.ndarray:
    """Distance (m) from each forest cell to the nearest non-forest cell."""
    lat = transform.f + transform.e * forest_cells.shape[0] / 2
    dy = abs(transform.e) * 110_574.0
    dx = abs(transform.a) * 111_320.0 * np.cos(np.radians(lat))
    return distance_transform_edt(forest_cells, sampling=(dy, dx)).astype("float32")


def weighted_percentile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    sel = weights > 0
    if not sel.any():
        return 0.0
    v, w = values[sel], weights[sel]
    order = np.argsort(v)
    cum = np.cumsum(w[order])
    return float(v[order][np.searchsorted(cum, q / 100 * cum[-1])])


def classify(dist: np.ndarray, has_forest: np.ndarray, dmin: float, nrt: float) -> np.ndarray:
    """30 vulnerability classes: 30 = at the forest edge (distance <= one cell), 29..1 = geometric
    distance bins between one cell and the NRT, 0 = beyond the NRT or no forest."""
    upper = dmin * (nrt / dmin) ** (np.arange(1, N_CLASSES) / (N_CLASSES - 1))    # 29 upper limits, last = NRT
    idx = np.searchsorted(upper, dist, side="left")                                # 0..28 inside the NRT
    cls = np.where(dist <= dmin, N_CLASSES, (N_CLASSES - 1) - idx).astype("int16")
    cls[(dist > nrt) | ~has_forest] = 0
    return cls.astype("uint8")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--region", required=True, choices=sorted(REGIONS))
    p.add_argument("--raw", default="data/raw")
    p.add_argument("--interim", default="data/interim")
    p.add_argument("--processed", default="data/processed")
    p.add_argument("--results", default="results")
    args = p.parse_args()
    res_dir = Path(args.results)
    res_dir.mkdir(parents=True, exist_ok=True)

    with rasterio.open(Path(args.interim) / f"gfc_counts_90m_{args.region}.tif") as src:
        c = src.read()                                   # uint8 counts, kept small on purpose
        transform, shape, crs = src.transform, (src.height, src.width), src.crs
    pix = pixel_area_ha(transform, shape[0])
    forest_counts = {"p1": c[0], "p2": c[0] - c[1], "p3": c[0] - c[1] - c[2]}   # forest at the start of each period
    loss_counts = {"p1": c[1], "p2": c[2], "p3": c[3]}

    adm = gpd.read_file(Path(args.raw) / "geoBoundaries-COL-ADM1.geojson").to_crs("EPSG:4326")
    adm = adm[adm["shapeName"].isin(REGIONS[args.region])].reset_index(drop=True)
    dept = rasterize([(g, i + 1) for i, g in enumerate(adm.geometry)], out_shape=shape, transform=transform,
                     fill=0, dtype="uint8")
    n_dept = len(adm) + 1

    hexes = gpd.read_file(Path(args.processed) / f"hex5km_{args.region}_stats.gpkg")
    hex_wgs = hexes.to_crs("EPSG:4326")
    hex_idx = rasterize([(g, i + 1) for i, g in enumerate(hex_wgs.geometry)], out_shape=shape,
                        transform=transform, fill=0, dtype="int32")
    hex_idx[dept == 0] = 0
    nh = len(hexes) + 1
    dmin = float(abs(transform.e) * 110_574.0)
    chunks = [(r, min(r + 512, shape[0])) for r in range(0, shape[0], 512)]

    def region_ids(cls, zone, dep, model):
        base = cls.astype("int32") * n_dept + dep
        return base if model == "benchmark" else base * 3 + zone

    rows, hex_out = [], hexes[["hex_id", "department", "geometry"]].copy()
    for t in TESTS:
        f, q = t["fit"], t["predict"]

        # fitting: distance from the forest at the start of the fitting period
        d_fit = edge_distance(forest_counts[f] >= 5, transform)
        inside = (dept > 0) & (loss_counts[f] > 0)
        nrt = weighted_percentile(d_fit[inside], loss_counts[f][inside].astype("float64"), NRT_Q) or dmin
        cls_fit = classify(d_fit, forest_counts[f] > 0, dmin, nrt)
        del d_fit, inside
        d_pred = edge_distance(forest_counts[q] >= 5, transform)
        cls_pred = classify(d_pred, forest_counts[q] > 0, dmin, nrt)
        del d_pred

        # hot/cold spot zone of each cell, from the Gi* of the fitting period (no information from the future)
        zone_of_hex = np.concatenate([[ZONES["not significant"]],
                                      hexes[f"hotspot_{f}"].fillna("not significant").map(ZONES).to_numpy()]).astype("uint8")
        zone = zone_of_hex[hex_idx]

        obs_hex = np.zeros(nh)
        forest_hex = np.zeros(nh)
        for a, b in chunks:
            area = pix[a:b, None].astype("float64")
            obs_hex += np.bincount(hex_idx[a:b].ravel(), weights=(loss_counts[q][a:b] * area).ravel(), minlength=nh)
            forest_hex += np.bincount(hex_idx[a:b].ravel(), weights=(forest_counts[q][a:b] * area).ravel(), minlength=nh)
        obs_total = obs_hex[1:].sum()

        results = {}
        for model in ("benchmark", "alternative"):
            n = (N_CLASSES + 1) * n_dept * (1 if model == "benchmark" else 3)
            loss_r, forest_r = np.zeros(n), np.zeros(n)
            for a, b in chunks:                                   # fit densities chunk by chunk
                m = (dept[a:b] > 0) & (cls_fit[a:b] > 0)
                reg = region_ids(cls_fit[a:b], zone[a:b], dept[a:b], model)[m]
                area = np.broadcast_to(pix[a:b, None], m.shape)[m].astype("float64")
                loss_r += np.bincount(reg, weights=loss_counts[f][a:b][m] * area, minlength=n)
                forest_r += np.bincount(reg, weights=forest_counts[f][a:b][m] * area, minlength=n)
            with np.errstate(invalid="ignore", divide="ignore"):
                dens = np.where(forest_r > 0, loss_r / forest_r, 0.0)
            pred_hex = np.zeros(nh)
            for a, b in chunks:                                   # apply them to the next period
                m = (dept[a:b] > 0) & (cls_pred[a:b] > 0)
                reg = region_ids(cls_pred[a:b], zone[a:b], dept[a:b], model)[m]
                area = np.broadcast_to(pix[a:b, None], m.shape)[m].astype("float64")
                pred_hex += np.bincount(hex_idx[a:b][m], weights=dens[reg] * forest_counts[q][a:b][m] * area, minlength=nh)
            scale = obs_total / pred_hex[1:].sum() if pred_hex[1:].sum() > 0 else 0.0
            pred_hex *= scale
            valid = forest_hex[1:] > 0
            err = np.abs(pred_hex[1:] - obs_hex[1:])[valid]
            results[model] = {"medae": float(np.median(err)), "mae": float(np.mean(err)), "scale": float(scale)}
            hex_out[f"pred_{model}_{q}_ha"] = np.round(pred_hex[1:], 2)
        hex_out[f"obs_{q}_ha"] = np.round(obs_hex[1:], 2)

        if t["name"] == "test 2":
            with rasterio.open(res_dir / f"risk_class_benchmark_2016_{args.region}.tif", "w", driver="GTiff",
                               width=shape[1], height=shape[0], count=1, dtype="uint8", crs=crs,
                               transform=transform, compress="deflate", nodata=0) as dst:
                dst.write(np.where(dept > 0, cls_pred, 0).astype("uint8"), 1)
        del cls_fit, cls_pred, zone

        rows.append({
            "region": args.region, "test": t["name"], "fit": t["fit_label"], "predict": t["pred_label"],
            "nrt_m": round(nrt), "observed_loss_ha": round(float(obs_total)),
            "medae_benchmark_ha": round(results["benchmark"]["medae"], 3),
            "medae_alternative_ha": round(results["alternative"]["medae"], 3),
            "mae_benchmark_ha": round(results["benchmark"]["mae"], 2),
            "mae_alternative_ha": round(results["alternative"]["mae"], 2),
            "alternative_better": results["alternative"]["medae"] < results["benchmark"]["medae"],
        })
        print(json.dumps(rows[-1]))

    hex_out.to_file(Path(args.processed) / f"hex5km_{args.region}_vt0007.gpkg", driver="GPKG")
    out = res_dir / "vt0007_benchmark_vs_hotspot_alternative.csv"
    df = pd.DataFrame(rows)
    if out.exists():
        old = pd.read_csv(out)
        df = pd.concat([old[old["region"] != args.region], df], ignore_index=True)
    df.to_csv(out, index=False)


if __name__ == "__main__":
    main()
