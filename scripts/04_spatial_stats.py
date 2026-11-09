"""Spatial autocorrelation of deforestation with gis-mcp (PySAL), confirmed with esda.

For each region and period:

1. The agent step: the gis-mcp server is launched over stdio, exactly as an AI
   host would launch it, and its PySAL tools are called through an MCP client:
   `morans_i`, `moran_local` and `getis_ord_g_local`. The data stay in the
   equal-area CRS (target_crs = LAEA Colombia) and the distance band is given in
   metres (5.5 km = the first ring of 5 km hexagons).
2. The script step: the same statistics are recomputed directly with esda to
   confirm the tool results, and then the classification used in the maps is
   produced with choices the tool does not expose: Getis-Ord Gi* (the hexagon
   itself included), binary first-ring weights, 9,999 permutations and a False
   Discovery Rate correction for multiple testing.

Usage:  python scripts/04_spatial_stats.py --gis-mcp /path/to/gis-mcp
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anyio
import esda
import geopandas as gpd
import libpysal
import numpy as np
import pandas as pd
from mcp import Client, StdioServerParameters

LAEA_COL = "+proj=laea +lat_0=4 +lon_0=-73.5 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
BAND_M = 5500.0
PERIODS = {"p1": "2001-2008", "p2": "2009-2016", "p3": "2017-2025"}
REGIONS = ["amazon_arc", "choco"]
PERMS = 9999
ALPHA = 0.05
SEED = 20261109


async def run_gis_mcp(server: StdioServerParameters, path: str, var: str) -> dict:
    """Call the three gis-mcp PySAL tools on one layer and return their results."""
    args = {"shapefile_path": path, "dependent_var": var, "target_crs": LAEA_COL, "distance_threshold": BAND_M}
    out = {}
    async with Client(server) as client:
        for tool in ("morans_i", "moran_local", "getis_ord_g_local"):
            res = await client.call_tool(tool, args)
            payload = json.loads(res.content[0].text)
            if payload.get("status") != "success":
                raise RuntimeError(f"gis-mcp {tool} failed: {payload.get('message')}")
            out[tool] = payload
    return out


def lisa_labels(lisa: esda.Moran_Local, sig: np.ndarray) -> np.ndarray:
    labels = np.array(["not significant"] * len(lisa.q), dtype=object)
    names = {1: "high-high", 2: "low-high", 3: "low-low", 4: "high-low"}
    for q, name in names.items():
        labels[(lisa.q == q) & sig] = name
    return labels


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--gis-mcp", required=True, help="path to the gis-mcp executable")
    p.add_argument("--processed", default="data/processed")
    p.add_argument("--results", default="results")
    p.add_argument("--regions", nargs="+", default=REGIONS, choices=REGIONS)
    args = p.parse_args()
    proc, res_dir = Path(args.processed), Path(args.results)
    res_dir.mkdir(parents=True, exist_ok=True)
    server = StdioServerParameters(command=args.gis_mcp, args=[])

    summary = []
    for region in args.regions:
        hexes = gpd.read_file(proc / f"hex5km_{region}.gpkg")
        for k, label in PERIODS.items():
            var = f"rate_{k}"
            sub = hexes[hexes[var].notna()].reset_index(drop=True)
            # hexagons with no neighbour within the band (isolated by hexagons without a rate) are left out,
            # so that gis-mcp and esda work on exactly the same set
            w0 = libpysal.weights.DistanceBand.from_dataframe(sub, threshold=BAND_M, binary=True, silence_warnings=True)
            isolated = len(w0.islands)
            if isolated:
                sub = sub.drop(index=w0.islands).reset_index(drop=True)
            layer = proc / f"stats_input_{region}_{k}.gpkg"
            sub[["hex_id", var, "geometry"]].to_file(layer, driver="GPKG")

            # 1. agent step: gis-mcp tools through MCP
            tool = anyio.run(run_gis_mcp, server, str(layer.resolve()), var)
            t_moran = tool["morans_i"]["result"]

            # 2. script step: confirm with esda using the same weights as gis-mcp
            y = sub[var].to_numpy(dtype=float)
            w_tool = libpysal.weights.DistanceBand.from_dataframe(sub, threshold=BAND_M, binary=False, silence_warnings=True)
            w_tool.transform = "r"
            np.random.seed(SEED)
            moran_check = esda.Moran(y, w_tool, permutations=PERMS)

            # classification used in the maps: binary first-ring weights, Gi*, FDR
            w = libpysal.weights.DistanceBand.from_dataframe(sub, threshold=BAND_M, binary=True, silence_warnings=True)
            w_b = libpysal.weights.W(w.neighbors, silence_warnings=True)   # binary copy for Gi*
            w.transform = "r"
            lisa = esda.Moran_Local(y, w, permutations=PERMS, seed=SEED, n_jobs=2)
            gstar = esda.G_Local(y, w_b, transform="B", star=True, permutations=PERMS, seed=SEED, n_jobs=2)
            fdr_l = esda.fdr(lisa.p_sim, ALPHA)
            fdr_g = esda.fdr(gstar.p_sim, ALPHA)
            sig_l = lisa.p_sim <= fdr_l
            sig_g = gstar.p_sim <= fdr_g
            hot = np.where(sig_g & (gstar.Zs > 0), "hot spot", np.where(sig_g & (gstar.Zs < 0), "cold spot", "not significant"))

            cols = pd.DataFrame({
                "hex_id": sub["hex_id"],
                f"lisa_{k}": lisa_labels(lisa, sig_l),
                f"lisa_p_{k}": lisa.p_sim,
                f"gistar_z_{k}": gstar.Zs,
                f"gistar_p_{k}": gstar.p_sim,
                f"hotspot_{k}": hot,
                f"tool_local_moran_{k}": tool["moran_local"]["result"]["Is"],
                f"tool_gi_z_{k}": tool["getis_ord_g_local"]["result"]["z_scores"],
            })
            hexes = hexes.drop(columns=[c for c in cols.columns if c != "hex_id" and c in hexes], errors="ignore")
            hexes = hexes.merge(cols, on="hex_id", how="left")

            loss = sub[f"loss_{k}_ha"].to_numpy()
            summary.append({
                "region": region,
                "period": label,
                "hexagons": len(sub),
                "isolated_hexagons_excluded": isolated,
                "morans_i_gis_mcp": round(t_moran["I"], 4),
                "morans_i_esda": round(float(moran_check.I), 4),
                "p_value_gis_mcp": t_moran["p_value"],
                "z_gis_mcp": round(t_moran["z_score"], 2),
                "local_moran_max_abs_diff": float(np.max(np.abs(np.array(tool["moran_local"]["result"]["Is"]) -
                                                     esda.Moran_Local(y, w_tool, permutations=0).Is))),
                "hot_spots": int((hot == "hot spot").sum()),
                "cold_spots": int((hot == "cold spot").sum()),
                "high_high": int((lisa_labels(lisa, sig_l) == "high-high").sum()),
                "share_of_loss_in_hot_spots_pct": round(float(loss[hot == "hot spot"].sum() / loss.sum() * 100), 1),
                "share_of_area_in_hot_spots_pct": round(float((hot == "hot spot").mean() * 100), 1),
                "fdr_p_threshold_gistar": float(fdr_g),
            })
            print(json.dumps(summary[-1]))
        hexes.to_file(proc / f"hex5km_{region}_stats.gpkg", driver="GPKG")

    out = res_dir / "spatial_autocorrelation_summary.csv"
    df = pd.DataFrame(summary)
    if out.exists():
        old = pd.read_csv(out)
        df = pd.concat([old[~old["region"].isin(args.regions)], df], ignore_index=True)
    df.to_csv(out, index=False)


if __name__ == "__main__":
    main()
