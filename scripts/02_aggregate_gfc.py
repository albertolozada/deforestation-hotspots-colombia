"""Aggregate Hansen GFC (30 m) to a 90 m count raster per study region.

For every 3 x 3 block of Hansen pixels the output stores:

  band 1  forest_2000   pixels with treecover2000 >= 30 %   (0-9)
  band 2  loss_p1       of those, loss in 2001-2008          (0-9)
  band 3  loss_p2       of those, loss in 2009-2016          (0-9)
  band 4  loss_p3       of those, loss in 2017-2025          (0-9)

The 30 % canopy threshold follows the national forest definition used by IDEAM
for Colombia. Counts (not areas) are stored: pixel areas are computed later,
geodesically, from the latitude of each row.

The script is resumable: it processes horizontal strips until a time budget is
used and records progress, so it can run in several short calls.

Usage:
  python scripts/02_aggregate_gfc.py --region amazon_arc --raw data/raw --out data/interim [--budget 150]
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.transform import from_origin
from rasterio.windows import from_bounds

RES = 0.00025              # Hansen pixel size (degrees)
F = 3                      # aggregation factor: 3 x 3 pixels, about 90 m
ORES = RES * F
TILES = ["10N_080W", "00N_080W", "10N_070W"]
PERIODS = [(1, 8), (9, 16), (17, 25)]
REGIONS = {
    "amazon_arc": ["Meta", "Guaviare", "Caquetá", "Putumayo"],
    "choco": ["Chocó"],
}
STRIP_OUT_ROWS = 300       # output rows per strip (900 Hansen rows)


def region_bounds(raw: Path, region: str) -> tuple[float, float, float, float]:
    adm = gpd.read_file(raw / "geoBoundaries-COL-ADM1.geojson")
    names = REGIONS[region]
    sel = adm[adm["shapeName"].isin(names)]
    missing = set(names) - set(sel["shapeName"])
    if missing:
        raise SystemExit(f"Departments not found in boundaries: {missing}. Available: {sorted(adm.shapeName)}")
    minx, miny, maxx, maxy = sel.to_crs("EPSG:4326").total_bounds
    pad = 0.02
    # snap to the 90 m grid anchored on the Hansen tile origin (-80, 10)
    snap = lambda v, o, f: o + f(((v - o) / ORES)) * ORES  # noqa: E731
    return (snap(minx - pad, -80, math.floor), snap(miny - pad, 10, math.floor),
            snap(maxx + pad, -80, math.ceil), snap(maxy + pad, 10, math.ceil))


def read_mosaic(raw: Path, layer: str, bounds: tuple[float, float, float, float], shape: tuple[int, int]) -> np.ndarray:
    """Read a Hansen layer over `bounds` from whichever tiles intersect it."""
    out = np.zeros(shape, dtype="uint8")
    west, south, east, north = bounds
    for tile in TILES:
        path = raw / f"Hansen_GFC-2025-v1.13_{layer}_{tile}.tif"
        if not path.exists():
            continue
        with rasterio.open(path) as src:
            tw, ts, te, tn = src.bounds
            iw, i_s, ie, i_n = max(west, tw), max(south, ts), min(east, te), min(north, tn)
            if iw >= ie or i_s >= i_n:
                continue
            win = from_bounds(iw, i_s, ie, i_n, src.transform).round_offsets().round_lengths()
            data = src.read(1, window=win)
            r0 = int(round((north - i_n) / RES))
            c0 = int(round((iw - west) / RES))
            out[r0:r0 + data.shape[0], c0:c0 + data.shape[1]] = data
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--region", required=True, choices=sorted(REGIONS))
    p.add_argument("--raw", default="data/raw")
    p.add_argument("--out", default="data/interim")
    p.add_argument("--budget", type=float, default=150.0, help="seconds to work before stopping")
    args = p.parse_args()
    raw, out = Path(args.raw), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    dst_path = out / f"gfc_counts_90m_{args.region}.tif"
    prog_path = out / f"gfc_counts_90m_{args.region}.progress.json"

    west, south, east, north = region_bounds(raw, args.region)
    width = int(round((east - west) / ORES))
    height = int(round((north - south) / ORES))
    transform = from_origin(west, north, ORES, ORES)

    if not dst_path.exists():
        profile = dict(driver="GTiff", width=width, height=height, count=4, dtype="uint8", crs="EPSG:4326",
                       transform=transform, compress="deflate", predictor=2, tiled=True, blockxsize=512,
                       blockysize=512, nodata=None, BIGTIFF="IF_SAFER")
        with rasterio.open(dst_path, "w", **profile) as dst:
            dst.descriptions = ("forest_2000", "loss_2001_2008", "loss_2009_2016", "loss_2017_2025")
        prog_path.write_text(json.dumps({"next_row": 0, "height": height, "width": width}))

    prog = json.loads(prog_path.read_text())
    start = time.time()
    with rasterio.open(dst_path, "r+") as dst:
        while prog["next_row"] < height and time.time() - start < args.budget:
            r = prog["next_row"]
            n = min(STRIP_OUT_ROWS, height - r)
            s_north = north - r * ORES
            s_south = s_north - n * ORES
            shape = (n * F, width * F)
            tc = read_mosaic(raw, "treecover2000", (west, s_south, east, s_north), shape)
            ly = read_mosaic(raw, "lossyear", (west, s_south, east, s_north), shape)
            forest = tc >= 30
            bands = [forest] + [forest & (ly >= a) & (ly <= b) for a, b in PERIODS]
            agg = np.stack([b.reshape(n, F, width, F).sum(axis=(1, 3)).astype("uint8") for b in bands])
            dst.write(agg, window=rasterio.windows.Window(0, r, width, n))
            prog["next_row"] = r + n
            prog_path.write_text(json.dumps(prog))
    done = prog["next_row"] >= height
    print(json.dumps({"region": args.region, "rows_done": prog["next_row"], "rows_total": height,
                      "done": done, "seconds": round(time.time() - start, 1)}))


if __name__ == "__main__":
    main()
