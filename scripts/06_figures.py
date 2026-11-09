"""Figures for the README: banner, hot spot migration, LISA clusters and the VT0007 test."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from matplotlib import font_manager
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

LAEA_COL = "+proj=laea +lat_0=4 +lon_0=-73.5 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
REGIONS = {"amazon_arc": ["Meta", "Guaviare", "Caquetá", "Putumayo"], "choco": ["Chocó"]}
TITLES = {"amazon_arc": "Amazon deforestation arc", "choco": "Chocó"}
PERIODS = {"p1": "2001-2008", "p2": "2009-2016", "p3": "2017-2025"}

INK, MUTED, SURFACE, LINE = "#0b0b0b", "#52514e", "#fcfcfb", "#d6d4cd"
HOT, COLD, NS, NOF = "#e34948", "#2a78d6", "#e7e5df", "#ffffff"
LISA_COLORS = {"high-high": "#c0262d", "low-low": "#1f5fae", "high-low": "#f2a3a1", "low-high": "#9cc3ef", "not significant": NS}
DARK_BG, DARK_PANEL, DARK_TEXT, DARK_MUTED = "#0b1419", "#111d24", "#f2f5f3", "#9fb0ae"
HOT_D, COLD_D, NS_D, NOF_D = "#ff5a4e", "#3d9be9", "#2a3c45", "#121e24"

for f in font_manager.findSystemFonts():
    if "Inter" in f:
        font_manager.fontManager.addfont(f)
plt.rcParams["font.family"] = "Inter"


def load(proc: Path, raw: Path):
    adm = gpd.read_file(raw / "geoBoundaries-COL-ADM1.geojson").to_crs(LAEA_COL)
    data = {r: gpd.read_file(proc / f"hex5km_{r}_stats.gpkg") for r in REGIONS}
    vt = {r: gpd.read_file(proc / f"hex5km_{r}_vt0007.gpkg") for r in REGIONS}
    return adm, data, vt


def hot_colors(series, dark=False):
    """Hexagons without a rate (little or no forest at the start of the period) get their own colour."""
    h, c, n, nf = (HOT_D, COLD_D, NS_D, NOF_D) if dark else (HOT, COLD, NS, NOF)
    return series.fillna("no forest").map({"hot spot": h, "cold spot": c, "not significant": n, "no forest": nf})


def banner(adm, data, summary, out: Path):
    fig = plt.figure(figsize=(16, 8), dpi=120)
    fig.patch.set_facecolor(DARK_BG)
    axes = {"choco": fig.add_axes([0.0, 0.1, 0.2, 0.88]), "amazon_arc": fig.add_axes([0.19, 0.1, 0.4, 0.88])}
    for r, ax in axes.items():
        g = data[r]
        g.plot(ax=ax, color=hot_colors(g["hotspot_p3"], dark=True), linewidth=0)
        adm[adm.shapeName.isin(REGIONS[r])].boundary.plot(ax=ax, color="#ffffff", linewidth=0.8, alpha=0.7)
        ax.set_axis_off()
        ax.set_aspect("equal")
        s = summary[(summary.region == r) & (summary.period == "2017-2025")].iloc[0]
        box = ax.get_position()
        fig.text(box.x0 + box.width / 2, 0.035, f"{TITLES[r]}  \u00b7  Moran's I {s.morans_i_gis_mcp:.2f}", ha="center",
                 va="bottom", color=DARK_TEXT, fontsize=12.5, fontweight="bold")

    # locator
    loc = fig.add_axes([0.215, 0.72, 0.1, 0.25])
    col = adm.dissolve()
    col.plot(ax=loc, color="#2b3b43", edgecolor="#4a5b63", linewidth=0.4)
    for r in REGIONS:
        adm[adm.shapeName.isin(REGIONS[r])].dissolve().plot(ax=loc, color="#ffa94d")
    loc.set_axis_off()
    loc.set_aspect("equal")

    tp = fig.add_axes([0.6, 0.0, 0.4, 1.0])
    tp.set_facecolor(DARK_PANEL)
    tp.set_xticks([]), tp.set_yticks([])
    for sp in tp.spines.values():
        sp.set_visible(False)
    x0 = 0.08
    tp.text(x0, 0.875, "Where does\ndeforestation cluster?", fontsize=31, fontweight="bold", color=DARK_TEXT,
            transform=tp.transAxes, va="center", linespacing=1.15)
    tp.text(x0, 0.735, "Colombia, Hansen GFC 2001-2025, 5 km hexagons.\nMoran's I and Getis-Ord Gi* run by an AI agent through\ngis-mcp, confirmed with PySAL.",
            fontsize=13, color=DARK_MUTED, transform=tp.transAxes, va="center", linespacing=1.45)

    a = summary[(summary.region == "amazon_arc") & (summary.period == "2017-2025")].iloc[0]
    c = summary[(summary.region == "choco") & (summary.period == "2017-2025")].iloc[0]
    stats = [
        (f"{a.share_of_loss_in_hot_spots_pct:.0f} %", f"of 2017-2025 tree cover loss in the\nAmazon arc falls in hot spots that\ncover {a.share_of_area_in_hot_spots_pct:.0f} % of its area", HOT_D),
        (f"{c.share_of_loss_in_hot_spots_pct:.0f} %", f"in Chocó, in hot spots that\ncover only {c.share_of_area_in_hot_spots_pct:.0f} % of its area", "#ffa94d"),
    ]
    y = 0.57
    for big, small, color in stats:
        tp.text(x0, y, big, fontsize=40, fontweight="bold", color=color, transform=tp.transAxes, va="center")
        tp.text(x0 + 0.36, y, small, fontsize=12.5, color=DARK_TEXT, transform=tp.transAxes, va="center", linespacing=1.4)
        y -= 0.15
    tp.text(x0, 0.27, "VM0048 link", fontsize=12, fontweight="bold", color="#7ee081", transform=tp.transAxes)
    tp.text(x0, 0.2, "Adding the hot spots of the previous period to a\nVT0007-style benchmark risk map roughly halves\nthe allocation error in the Amazon arc.",
            fontsize=12.5, color=DARK_TEXT, transform=tp.transAxes, va="center", linespacing=1.45)
    handles = [Patch(color=HOT_D, label="hot spot"), Patch(color=COLD_D, label="cold spot"),
               Patch(color=NS_D, label="not significant"), Patch(facecolor=NOF_D, edgecolor="#3a4b53", label="little forest")]
    tp.legend(handles=handles, loc="lower left", bbox_to_anchor=(x0 - 0.01, 0.03), ncol=4, frameon=False,
              labelcolor=DARK_TEXT, fontsize=11, title="Getis-Ord Gi*, 2017-2025, FDR 5 %", title_fontsize=10.5)
    tp.get_legend().get_title().set_color(DARK_MUTED)
    fig.savefig(out, facecolor=DARK_BG)
    plt.close(fig)


def migration(adm, data, out: Path):
    fig, axes = plt.subplots(2, 3, figsize=(14, 9.5), dpi=130, gridspec_kw={"height_ratios": [1.0, 0.85]})
    fig.patch.set_facecolor(SURFACE)
    for i, r in enumerate(["choco", "amazon_arc"]):
        g = data[r]
        for j, (k, label) in enumerate(PERIODS.items()):
            ax = axes[i, j]
            g.plot(ax=ax, color=hot_colors(g[f"hotspot_{k}"]), linewidth=0)
            adm[adm.shapeName.isin(REGIONS[r])].boundary.plot(ax=ax, color=MUTED, linewidth=0.6)
            ax.set_axis_off()
            ax.set_aspect("equal")
            ax.set_title(f"{TITLES[r]}, {label}", fontsize=12, color=INK, loc="left")
    handles = [Patch(color=HOT, label="hot spot"), Patch(color=COLD, label="cold spot"), Patch(color=NS, label="not significant"),
               Patch(facecolor=NOF, edgecolor=LINE, label="little or no forest (no rate)")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=11)
    fig.subplots_adjust(hspace=0.02, wspace=0.02, top=0.9, bottom=0.06)
    fig.suptitle("Getis-Ord Gi* hot spots of annual tree cover loss rate, by period", x=0.02, ha="left", fontsize=15,
                 fontweight="bold", color=INK)
    fig.text(0.02, 0.945, "5 km hexagons, first-ring binary weights, 9,999 permutations, False Discovery Rate 5 %.",
             fontsize=10.5, color=MUTED)
    fig.savefig(out, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)


def lisa(adm, data, out: Path):
    fig, axes = plt.subplots(1, 2, figsize=(14, 7.5), dpi=130, gridspec_kw={"width_ratios": [1, 2]})
    fig.patch.set_facecolor(SURFACE)
    for ax, r in zip(axes, ["choco", "amazon_arc"]):
        g = data[r]
        g.plot(ax=ax, color=g["lisa_p3"].map(LISA_COLORS).fillna(NOF), edgecolor=g["lisa_p3"].isna().map({True: LINE, False: "none"}), linewidth=0.2)
        adm[adm.shapeName.isin(REGIONS[r])].boundary.plot(ax=ax, color=MUTED, linewidth=0.6)
        ax.set_axis_off()
        ax.set_aspect("equal")
        ax.set_title(TITLES[r], fontsize=12, color=INK, loc="left")
    handles = [Patch(color=v, label=k) for k, v in LISA_COLORS.items()] + [Patch(facecolor=NOF, edgecolor=LINE, label="little or no forest")]
    fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False, fontsize=11)
    fig.suptitle("Local Moran's I clusters, annual tree cover loss rate 2017-2025", x=0.02, ha="left", fontsize=15,
                 fontweight="bold", color=INK)
    fig.savefig(out, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)


def vt0007(adm, vt, results: Path, out: Path):
    df = pd.read_csv(results / "vt0007_benchmark_vs_hotspot_alternative.csv")
    fig = plt.figure(figsize=(15, 7.8), dpi=130)
    fig.patch.set_facecolor(SURFACE)

    with rasterio.open(results / "risk_class_benchmark_2016_amazon_arc.tif") as src:
        cls = src.read(1, out_shape=(src.height // 4, src.width // 4), resampling=rasterio.enums.Resampling.mode)
        b = src.bounds
    ax1 = fig.add_axes([0.0, 0.12, 0.32, 0.76])
    cmap = ListedColormap(plt.get_cmap("magma_r")(np.linspace(0.05, 0.95, 30)))
    im = ax1.imshow(np.ma.masked_equal(cls, 0), extent=(b.left, b.right, b.bottom, b.top), cmap=cmap, vmin=1, vmax=30,
                    interpolation="nearest")
    cax1 = fig.add_axes([0.03, 0.07, 0.25, 0.016])
    cb = fig.colorbar(im, cax=cax1, orientation="horizontal", ticks=[1, 15, 30])
    cb.ax.set_xticklabels(["1 (near the NRT)", "15", "30 (at the edge)"], fontsize=9, color=MUTED)
    cb.outline.set_visible(False)
    fig.text(0.03, 0.015, "White: non-forest or beyond the negligible risk threshold (NRT).", fontsize=9, color=MUTED)
    adm_wgs = adm.to_crs("EPSG:4326")
    adm_wgs[adm_wgs.shapeName.isin(REGIONS["amazon_arc"])].boundary.plot(ax=ax1, color=MUTED, linewidth=0.6)
    ax1.set_axis_off()
    ax1.set_title("Benchmark vulnerability, forest 2016\n(distance to forest edge, 30 classes)", fontsize=11.5,
                  color=INK, loc="left")

    g = vt["amazon_arc"]
    ax2 = fig.add_axes([0.34, 0.12, 0.32, 0.76])
    err_b = (g["pred_benchmark_p3_ha"] - g["obs_p3_ha"])
    err_a = (g["pred_alternative_p3_ha"] - g["obs_p3_ha"])
    gain = (err_b.abs() - err_a.abs())
    lim = np.nanpercentile(np.abs(gain), 98)
    g.assign(gain=gain).plot(ax=ax2, column="gain", cmap="RdBu", vmin=-lim, vmax=lim, linewidth=0)
    cax2 = fig.add_axes([0.37, 0.07, 0.25, 0.016])
    sm = plt.cm.ScalarMappable(cmap="RdBu", norm=plt.Normalize(-lim, lim))
    cb2 = fig.colorbar(sm, cax=cax2, orientation="horizontal")
    cb2.outline.set_visible(False)
    cb2.ax.tick_params(labelsize=9, colors=MUTED)
    fig.text(0.37, 0.015, "Error reduction per hexagon (ha); red: the benchmark was closer.", fontsize=9, color=MUTED)
    adm[adm.shapeName.isin(REGIONS["amazon_arc"])].boundary.plot(ax=ax2, color=MUTED, linewidth=0.6)
    ax2.set_axis_off()
    ax2.set_aspect("equal")
    ax2.set_title("Error reduction with hot spots, 2017-2025\n(blue: alternative closer to the observed loss)",
                  fontsize=11.5, color=INK, loc="left")

    ax3 = fig.add_axes([0.79, 0.2, 0.2, 0.62])
    labels, bench, alt = [], [], []
    for _, row in df.iterrows():
        labels.append(f"{'Amazon arc' if row.region == 'amazon_arc' else 'Chocó'}\n{row.test}: {row.predict}")
        bench.append(row.medae_benchmark_ha)
        alt.append(row.medae_alternative_ha)
    ypos = np.arange(len(labels))
    ax3.barh(ypos + 0.19, bench, height=0.36, color="#9a9893", label="benchmark (distance to edge)")
    ax3.barh(ypos - 0.19, alt, height=0.36, color="#2a78d6", label="benchmark + Gi* hot spots")
    for y, v in zip(ypos + 0.19, bench):
        ax3.text(v + 0.5, y, f"{v:.1f}", va="center", fontsize=9.5, color=INK)
    for y, v in zip(ypos - 0.19, alt):
        ax3.text(v + 0.5, y, f"{v:.1f}", va="center", fontsize=9.5, color=INK)
    ax3.set_yticks(ypos, labels, fontsize=9.5)
    ax3.invert_yaxis()
    ax3.set_xlabel("median absolute error per 5 km hexagon (ha)", fontsize=10, color=MUTED)
    for sp in ("top", "right"):
        ax3.spines[sp].set_visible(False)
    ax3.spines["left"].set_color(LINE)
    ax3.spines["bottom"].set_color(LINE)
    ax3.tick_params(colors=MUTED)
    ax3.set_facecolor(SURFACE)
    ax3.legend(loc="upper center", bbox_to_anchor=(0.45, -0.12), frameon=False, fontsize=9.5, ncol=1)
    ax3.set_title("Allocation error (MedAE)", fontsize=11.5, color=INK, loc="left")

    fig.suptitle("A VT0007-style benchmark risk map, with and without the hot spots of the previous period",
                 x=0.01, ha="left", fontsize=15, fontweight="bold", color=INK)
    fig.savefig(out, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    proc, raw, results = Path("data/processed"), Path("data/raw"), Path("results")
    out = Path("outputs/maps")
    out.mkdir(parents=True, exist_ok=True)
    adm, data, vt = load(proc, raw)
    summary = pd.read_csv(results / "spatial_autocorrelation_summary.csv")
    banner(adm, data, summary, out / "hero.png")
    migration(adm, data, out / "hotspot_migration.png")
    lisa(adm, data, out / "lisa_2017_2025.png")
    vt0007(adm, vt, results, out / "vt0007_comparison.png")
    print("figures written to", out)


if __name__ == "__main__":
    main()
