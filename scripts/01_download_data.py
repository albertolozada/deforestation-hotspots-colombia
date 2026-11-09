"""Download the public input data.

* Hansen / UMD Global Forest Change v1.13 (2000-2025), CC BY 4.0:
  `treecover2000` and `lossyear` tiles 10N_080W, 00N_080W and 10N_070W, which
  cover the Colombian Amazon deforestation arc (Meta, Guaviare, Caqueta,
  Putumayo) and Choco.
* geoBoundaries (gbOpen) administrative boundaries of Colombia, level 1
  (departments), derived from OpenStreetMap, ODbL 1.0.

If your network blocks these hosts, download the files with a browser into
data/raw (see data/README.md); the other scripts only need them to be there.

Usage:  python scripts/01_download_data.py [--out data/raw]
"""

from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

GFC = "https://storage.googleapis.com/earthenginepartners-hansen/GFC-2025-v1.13/Hansen_GFC-2025-v1.13_{layer}_{tile}.tif"
TILES = ["10N_080W", "00N_080W", "10N_070W"]
LAYERS = ["treecover2000", "lossyear"]
# geoBoundaries gbOpen COL ADM1, pinned to release commit 9469f09 (build of Dec 2023)
GEOBOUNDARIES = "https://github.com/wmgeolab/geoBoundaries/raw/9469f09/releaseData/gbOpen/COL/ADM1/geoBoundaries-COL-ADM1.geojson"


def fetch(url: str, dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        print(f"already there: {dest.name}")
        return
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"downloading {url}")
    with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as f:
        while chunk := r.read(1 << 20):
            f.write(chunk)
    tmp.rename(dest)
    print(f"  -> {dest} ({dest.stat().st_size / 1e6:.1f} MB)")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="data/raw")
    p.add_argument("--only", choices=["boundaries", "gfc"], default=None)
    args = p.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    if args.only in (None, "boundaries"):
        fetch(GEOBOUNDARIES, out / "geoBoundaries-COL-ADM1.geojson")

    if args.only in (None, "gfc"):
        for tile in TILES:
            for layer in LAYERS:
                fetch(GFC.format(layer=layer, tile=tile), out / f"Hansen_GFC-2025-v1.13_{layer}_{tile}.tif")


if __name__ == "__main__":
    sys.exit(main())
