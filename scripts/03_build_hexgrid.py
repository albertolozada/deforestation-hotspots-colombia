"""Aggregate the 90 m Hansen counts to a 5 km hexagon grid in an equal-area CRS.

Output: data/processed/hex5km_<region>.gpkg with, per hexagon,

  forest_2000_ha            forest (treecover2000 >= 30 %) in 2000
  loss_p1_ha, loss_p2_ha, loss_p3_ha   loss in 2001-2008, 2009-2016, 2017-2025
  forest_start_p1/p2/p3_ha  forest at the start of each period
  rate_p1/p2/p3             annual loss rate, % of the forest at the start of the period
  dens_p1/p2/p3             annual loss density, ha per km2 per year

Areas use the geodesic area of each 30 m Hansen pixel (WGS84 ellipsoid), so
nothing is measured in degrees. Hexagons with less than half of their area
inside the region are dropped; rates are left empty where the forest at the
start of a period is below MIN_FOREST_HA (a few hectares give unstable rates).

Usage:  python scripts/03_build_hexgrid.py --region amazon_arc
"""

from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from pyproj import Geod
from rasterio.features import rasterize
from shapely.geometry import Polygon

LAEA_COL = "+proj=laea +lat_0=4 +lon_0=-73.5 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
HEX_D = 5000.0                       # flat-to-flat distance (m): 21.65 km2 per hexagon
MIN_INSIDE = 0.5                     # minimum share of the hexagon inside the region
MIN_FOREST_HA = 200.0                # minimum forest at the start of a period to compute a rate
YEARS = {"p1": 8, "p2": 8, "p3": 9}  # 2001-2008, 2009-2016, 2017-2025
REGIONS = {
    "amazon_arc": ["Meta", "Guaviare", "Caquetá", "Putumayo"],
    "choco": ["Chocó"],
}


def hex_grid(bounds: tuple[float, float, float, float], d: float = HEX_D) -> gpd.GeoDataFrame:
    """Flat-topped hexagons with flat-to-flat distance d covering bounds (projected metres)."""
    s = d / np.sqrt(3)                              # side length
    minx, miny, maxx, maxy = bounds
    angles = np.radians(np.arange(0, 360, 60))
    polys, ids = [], []
    i = 0
    x = minx
    while x <= maxx + 1.5 * s:
        y = miny - (d / 2 if i % 2 else 0)
        j = 0
        while y <= maxy + d:
            polys.append(Polygon(np.column_stack([x + s * np.cos(angles), y + s * np.sin(angles)])))
            ids.append(f"{i:04d}_{j:04d}")
            y += d
            j += 1
        x += 1.5 * s
        i += 1
    return gpd.GeoDataFrame({"hex_id": ids}, geometry=polys, crs=LAEA_COL)


def pixel_area_rows(transform, height: int) -> np.ndarray:
    """Geodesic area (ha) of one 30 m Hansen pixel for each 90 m output row."""
    geod = Geod(ellps="WGS84")
    res = abs(transform.a) / 3
    areas = np.empty(height)
    for r in range(height):
        top = transform.f + r * transform.e + 1.5 * transform.e + res / 2   # centre latitude of the row
        lats = [top, top, top - res, top - res]
        lons = [0, res, res, 0]
        areas[r] = abs(geod.polygon_area_perimeter(lons, lats)[0]) / 10_000
    return areas


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--region", required=True, choices=sorted(REGIONS))
    p.add_argument("--raw", default="data/raw")
    p.add_argument("--interim", default="data/interim")
    p.add_argument("--out", default="data/processed")
    args = p.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    adm = gpd.read_file(Path(args.raw) / "geoBoundaries-COL-ADM1.geojson")
    depts = adm[adm["shapeName"].isin(REGIONS[args.region])][["shapeName", "geometry"]].to_crs(LAEA_COL)
    region = depts.union_all()

    hexes = hex_grid(region.bounds)
    hexes = hexes[hexes.intersects(region)].copy()
    inside = hexes.geometry.intersection(region).area / hexes.geometry.area
    hexes = hexes[inside >= MIN_INSIDE].reset_index(drop=True)
    hexes["share_inside"] = inside[inside >= MIN_INSIDE].round(3).values
    # department holding the largest part of each hexagon
    ov = gpd.overlay(hexes[["hex_id", "geometry"]], depts, how="intersection", keep_geom_type=True)
    ov["a"] = ov.area
    hexes["department"] = hexes["hex_id"].map(ov.sort_values("a").groupby("hex_id")["shapeName"].last())

    with rasterio.open(Path(args.interim) / f"gfc_counts_90m_{args.region}.tif") as src:
        counts = src.read().astype("float64")
        transform, shape = src.transform, (src.height, src.width)
    # hexagon index raster on the 90 m grid, restricted to the region
    hex_wgs = hexes.to_crs("EPSG:4326")
    region_wgs = gpd.GeoSeries([region], crs=LAEA_COL).to_crs("EPSG:4326").iloc[0]
    idx = rasterize([(g, i + 1) for i, g in enumerate(hex_wgs.geometry)], out_shape=shape, transform=transform,
                    fill=0, dtype="int32")
    in_region = rasterize([(region_wgs, 1)], out_shape=shape, transform=transform, fill=0, dtype="uint8")
    idx[in_region == 0] = 0

    pix_ha = pixel_area_rows(transform, shape[0])[:, None]
    n = len(hexes) + 1
    names = ["forest_2000_ha", "loss_p1_ha", "loss_p2_ha", "loss_p3_ha"]
    flat_idx = idx.ravel()
    for b, name in enumerate(names):
        hexes[name] = np.bincount(flat_idx, weights=(counts[b] * pix_ha).ravel(), minlength=n)[1:].round(2)
    hexes["area_ha"] = (hexes.geometry.intersection(region).area / 10_000).round(1)

    hexes["forest_start_p1_ha"] = hexes["forest_2000_ha"]
    hexes["forest_start_p2_ha"] = hexes["forest_start_p1_ha"] - hexes["loss_p1_ha"]
    hexes["forest_start_p3_ha"] = hexes["forest_start_p2_ha"] - hexes["loss_p2_ha"]
    for k, years in YEARS.items():
        start = hexes[f"forest_start_{k}_ha"]
        rate = hexes[f"loss_{k}_ha"] / start / years * 100
        hexes[f"rate_{k}"] = rate.where(start >= MIN_FOREST_HA).round(4)
        hexes[f"dens_{k}"] = (hexes[f"loss_{k}_ha"] / (hexes["area_ha"] / 100) / years).round(4)

    path = out / f"hex5km_{args.region}.gpkg"
    hexes.to_file(path, driver="GPKG")
    summary = {
        "region": args.region,
        "hexagons": len(hexes),
        "forest_2000_ha": round(hexes["forest_2000_ha"].sum()),
        **{f"loss_{k}_ha": round(hexes[f"loss_{k}_ha"].sum()) for k in YEARS},
        **{f"hexes_with_rate_{k}": int(hexes[f"rate_{k}"].notna().sum()) for k in YEARS},
    }
    print(pd.Series(summary).to_string())


if __name__ == "__main__":
    main()
