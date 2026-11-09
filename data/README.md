# Data

## Raw data (not versioned)

`scripts/01_download_data.py` downloads them into `data/raw/`. If your network blocks these hosts, download them with a browser into the same folder.

| File | Source | Licence |
|---|---|---|
| `Hansen_GFC-2025-v1.13_treecover2000_{10N_080W, 00N_080W, 10N_070W}.tif` | [Global Forest Change v1.13](https://storage.googleapis.com/earthenginepartners-hansen/GFC-2025-v1.13/download.html), University of Maryland | CC BY 4.0 |
| `Hansen_GFC-2025-v1.13_lossyear_{10N_080W, 00N_080W, 10N_070W}.tif` | same | CC BY 4.0 |
| `geoBoundaries-COL-ADM1.geojson` | [geoBoundaries](https://www.geoboundaries.org) gbOpen, release 9469f09, derived from OpenStreetMap | ODbL 1.0 |

Tile URLs follow `https://storage.googleapis.com/earthenginepartners-hansen/GFC-2025-v1.13/Hansen_GFC-2025-v1.13_{layer}_{tile}.tif`.

## Intermediate data (not versioned)

`data/interim/gfc_counts_90m_{region}.tif`: four bands of 90 m counts (0-9 Hansen pixels): forest 2000 and loss in 2001-2008, 2009-2016 and 2017-2025. Produced by `02_aggregate_gfc.py`.

## Processed data (versioned)

| File | Content |
|---|---|
| `processed/hex5km_{region}_stats.gpkg` | 5 km hexagons (LAEA Colombia): forest and loss (ha), annual rates and densities per period, LISA class, Gi* z and class, and the gis-mcp outputs |
| `processed/hex5km_{region}_vt0007.gpkg` | observed and predicted loss per hexagon for the benchmark and the hot spot alternative |

Both derive from Hansen GFC (CC BY 4.0) and geoBoundaries (ODbL 1.0: derived databases keep the ODbL).
