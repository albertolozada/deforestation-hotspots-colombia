"""Unit tests of the helpers behind the hexagons and the VT0007-style test."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest
from rasterio.transform import from_origin

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


hexmod = _load("03_build_hexgrid")
vt = _load("05_vt0007_benchmark")


def test_hexagons_have_the_expected_area_and_tile_the_plane():
    g = hexmod.hex_grid((0, 0, 30_000, 30_000))
    areas = g.geometry.area
    assert areas.iloc[0] == pytest.approx(np.sqrt(3) / 2 * 5000 ** 2, rel=1e-9)   # 21.65 km2
    union = g.geometry.union_all()
    assert union.area == pytest.approx(areas.sum(), rel=1e-9)                     # no overlaps
    assert union.buffer(0.01).contains(union.envelope.buffer(-6000))              # no gaps inside


def test_geodesic_pixel_area_at_the_equator():
    # a 0.00025 degree Hansen pixel at the equator is about 27.83 m x 27.64 m = 0.0769 ha
    tr = from_origin(-75, 0.00075, 0.00075, 0.00075)
    assert hexmod.pixel_area_rows(tr, 1)[0] == pytest.approx(0.0769, rel=2e-3)
    assert vt.pixel_area_ha(tr, 1)[0] == pytest.approx(0.0769, rel=2e-3)


def test_weighted_percentile():
    v = np.array([1.0, 2.0, 3.0, 4.0])
    w = np.array([1.0, 1.0, 1.0, 97.0])
    assert vt.weighted_percentile(v, w, 50) == 4.0
    assert vt.weighted_percentile(v, np.zeros(4), 99.5) == 0.0


def test_classes_are_geometric_and_bounded():
    d = np.array([0.0, 90.0, 500.0, 2000.0, 3999.0, 5000.0])
    has_forest = np.array([True, True, True, True, True, True])
    cls = vt.classify(d, has_forest, dmin=90.0, nrt=4000.0)
    assert cls[0] == 30 and cls[1] == 30          # at the edge
    assert cls[-1] == 0                           # beyond the NRT
    assert np.all(np.diff(cls[:-1].astype(int)) <= 0)   # risk decreases with distance
    assert cls[4] == 1                            # last class just inside the NRT
    assert set(vt.classify(np.linspace(91, 3999, 5000), np.ones(5000, bool), 90.0, 4000.0)) == set(range(1, 30))


def test_edge_distance_is_metric():
    forest = np.ones((50, 50), bool)
    forest[:, 0] = False
    tr = from_origin(-75, 0.0375, 0.00075, 0.00075)
    d = vt.edge_distance(forest, tr)
    assert d[25, 10] == pytest.approx(10 * 0.00075 * 111_320, rel=0.01)
