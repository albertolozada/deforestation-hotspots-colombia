"""Static location maps for project READMEs.

Draws an area of interest (AOI) over a web basemap with scale bar,
north arrow and basemap attribution, and saves it as a PNG ready to
embed in a GitHub README.

Mapas de localizacion estaticos para los README de los proyectos.
"""

from __future__ import annotations

from pathlib import Path

import contextily as ctx
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar

# Basemaps that do not need an API key.
# Fondos que no necesitan API key.
BASEMAPS = {
    "satellite": ctx.providers.Esri.WorldImagery,
    # CARTO tiles now need an API key, so the light basemap is Esri.
    "light": ctx.providers.Esri.WorldGrayCanvas,
    "osm": ctx.providers.OpenStreetMap.Mapnik,
}

# Outline colour that reads well on each basemap.
_EDGE_COLOURS = {"satellite": "#ffd60a", "light": "#c1121f", "osm": "#c1121f"}


def _to_projected(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Return the layer in a metric CRS (local UTM if it is geographic)."""
    if gdf.crs is None:
        raise ValueError("The layer has no CRS. Set it with gdf.set_crs(...) first.")
    if gdf.crs.is_geographic:
        return gdf.to_crs(gdf.estimate_utm_crs())
    return gdf


def _add_north_arrow(ax, x: float = 0.94, y: float = 0.86, size: float = 0.07) -> None:
    """Draw a simple north arrow (arrow pointing up, N on top) in axes coordinates."""
    ax.annotate(
        "",
        xy=(x, y),
        xytext=(x, y - size),
        xycoords="axes fraction",
        arrowprops=dict(facecolor="black", edgecolor="white", width=4, headwidth=12),
    )
    ax.text(
        x, y + 0.035, "N",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=12,
        fontweight="bold",
        bbox=dict(boxstyle="circle,pad=0.2", fc="white", ec="none", alpha=0.85),
    )


def plot_aoi_map(
    gdf: gpd.GeoDataFrame,
    out_png: str | Path | None = None,
    basemap: str = "satellite",
    title: str | None = None,
    margin: float = 0.15,
    figsize: tuple[float, float] = (8, 8),
    dpi: int = 200,
    zoom: int | str = "auto",
    cache_dir: str | Path | None = None,
    edgecolor: str | None = None,
    linewidth: float = 2.0,
):
    """Plot an area of interest over a web basemap.

    Parameters
    ----------
    gdf : GeoDataFrame
        AOI layer (polygons, lines or points) with a defined CRS.
        If it is geographic (e.g. EPSG:4326) it is reprojected to its
        local UTM zone so the scale bar is in metres.
    out_png : path, optional
        Output PNG. If None the figure is only returned.
    basemap : {"satellite", "light", "osm"} or xyzservices provider
        Background tiles. None of the presets needs an API key.
    title : str, optional
        Map title.
    margin : float
        Extra space around the AOI as a fraction of its larger side.
    zoom : int or "auto"
        Tile zoom level. "auto" lets contextily choose.
    cache_dir : path, optional
        Folder where tiles are cached, so the map can be redrawn
        without downloading them again.

    Returns
    -------
    (fig, ax)
    """
    if cache_dir is not None:
        Path(cache_dir).mkdir(parents=True, exist_ok=True)
        ctx.set_cache_dir(str(cache_dir))

    data = _to_projected(gdf)
    key = basemap if isinstance(basemap, str) else None
    source = BASEMAPS[basemap] if key else basemap
    colour = edgecolor or _EDGE_COLOURS.get(key, "#c1121f")

    fig, ax = plt.subplots(figsize=figsize)

    geom_types = set(data.geom_type)
    if geom_types <= {"Point", "MultiPoint"}:
        data.plot(ax=ax, color=colour, edgecolor="white", markersize=40, zorder=3)
    else:
        data.plot(ax=ax, facecolor="none", edgecolor=colour, linewidth=linewidth, zorder=3)

    # Extent with margin, kept square so the map does not look stretched.
    xmin, ymin, xmax, ymax = data.total_bounds
    side = max(xmax - xmin, ymax - ymin) * (1 + 2 * margin)
    side = max(side, 1000.0)  # at least 1 km for single points
    cx, cy = (xmin + xmax) / 2, (ymin + ymax) / 2
    ax.set_xlim(cx - side / 2, cx + side / 2)
    ax.set_ylim(cy - side / 2, cy + side / 2)

    # Tiles are warped to the layer CRS, so the data stays in UTM.
    ctx.add_basemap(ax, source=source, crs=data.crs, zoom=zoom, attribution_size=6)

    # Scale bar at the top so it never covers the basemap attribution (bottom).
    ax.add_artist(ScaleBar(1, units="m", location="upper left", box_alpha=0.8))
    _add_north_arrow(ax)

    ax.set_axis_off()
    if title:
        ax.set_title(title, fontsize=13, pad=10)
    epsg = data.crs.to_epsg()
    ax.text(1.0, -0.01, f"CRS: EPSG:{epsg}" if epsg else data.crs.name,
            transform=ax.transAxes, ha="right", va="top", fontsize=7, color="0.3")

    if out_png is not None:
        out_png = Path(out_png)
        out_png.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_png, dpi=dpi, bbox_inches="tight")

    return fig, ax
