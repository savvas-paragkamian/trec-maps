# trec-maps

Map scripts for the TREC expedition (EMBL, *Traversing European Coastlines*): land-sea
**profile maps** of the sampling transects, built from a DEM on land and EMODnet bathymetry
at sea, as 2D profiles (matplotlib) and 3D path-traced blocks (forge3d). Python only.

The pilot is the Athens supersite at Anavyssos: the soil transect 041-03 (HCMR) with the
Tara stations TE_176 (1.2 km offshore, 9 m) and TE_177 (131 m).

| output | what |
|---|---|
| `output/<transect>/profile.png`, `.pdf` | elevation / depth along the transect path, TREC zones, soil positions, stations with their sampling depths, a zoom of the land transect |
| `output/<transect>/block_<quality>.png`, `.pdf` | oblique path-traced block of land and sea floor with the path and stations drawn on it (vector overlay in the PDF) |

## Setup

A conda environment from free community channels only (`conda-forge`, `bioconda`,
`nodefaults`), the convention of the TREC repositories, so it works with Miniforge or with a
Miniconda whose `.condarc` avoids the Anaconda `defaults` channel. The environment lives in
conda's envs directory, never inside this Dropbox folder.

```bash
CONDA_CHANNEL_PRIORITY=strict conda env create -f environment.yml   # or: make env
conda activate trec-maps
make check          # forge3d 1.44.0 gpu True
```

forge3d renders on the local GPU (Metal on Apple silicon); the 2D profile does not need it.

## Data

Put the TREC data-portal exports in `data/` (they are not committed, see `data/README.md`):
`TREC_site_data_export_all-sites.xlsx`, `TREC_sample_data_export.xlsx` (soil and sediment
samples) and `TREC_sample_data_athens_marine_metagenome.xlsx` (marine samples). Everything
else is downloaded by the two get-scripts:

- **Land:** Copernicus DEM GLO-30, 1° tiles from the AWS open-data bucket (`scripts/00_get_dem.py`).
- **Sea floor:** EMODnet Bathymetry DTM (2024 release, 1/16 arc minute, about 115 m) through
  the WCS service, requested on its native grid for the transect's box (`scripts/01_get_bathymetry.py`).

## Run

```bash
make pilot                 # dem, bathy, grid, profile, render (preview) for T=anavyssos
make final                 # 3840 x 2160 block, 64 frames
make render Q=draft        # 2560 x 1440
python scripts/30_render_3d.py --quality preview --tag _try --set camera.azimuth_deg=200 exaggeration=4
python scripts/30_render_3d.py --quality final --recompose    # reuse the render, redo the page
```

Steps, each a thin script over `trecmaps/`:

| step | script | writes |
|---|---|---|
| 0 | `00_get_dem.py` | `data/raw/dem/*.tif` (the only DEM network code; `--dry-run`, `.part` then rename, skips existing files) |
| 1 | `01_get_bathymetry.py` | `data/raw/bathy/emodnet_mean_<transect>.tif` + `.json` (re-fetched when the grid box grows) |
| 10 | `10_build_grid.py` | `data/derived/<transect>/grid.npz`, `z.tif`, `grid.json`: DEM and bathymetry on one metric grid, merged with a land mask |
| 20 | `20_profile_2d.py` | `profile.csv`, `land_samples.csv`, `stations.csv`, `output/<transect>/profile.png/.pdf` |
| 30 | `30_render_3d.py` | `output/<transect>/block_<quality>.png/.pdf` and a render cache for `--recompose` |

## Configuration

- `config/params.yml`: paths, sources, grid (CRS, cell size, margins, land/sea merge),
  profile (step, zones, figure), render (camera, sun, exaggeration, quality presets, colours).
- `config/transects.yml`: one entry per transect, by TREC site IDs only (`land_site`,
  `marine_sites` in path order, optional `waypoints` to keep a marine leg off the beach).
  Coordinates come from the xlsx exports at run time.

To add a transect, add an entry and run `make pilot T=<name>`. `20_profile_2d.py` warns when
a marine leg crosses land; add a waypoint and rerun from `01_get_bathymetry.py`.

## How the grid is built

The path (inland end, shoreline, waypoints, stations) is boxed with margins in EPSG:2100
(GGRS87, metres; configurable) at 10 m cells. GLO-30 is resampled bilinearly, EMODnet
cubically after its land cells are set to nodata. Land is DEM > 0.5 m (speckle removed); sea
cells take the EMODnet floor, capped at −0.5 m and ramped to that cap over six cells from the
shoreline, so the 115 m bathymetry meets the 30 m coastline without a step. GLO-30 heights are
above EGM2008, EMODnet depths below LAT; the Saronic tidal range is about 0.1 m, so the offset
is ignored.

## Sources and citations

- Copernicus DEM GLO-30, ESA (produced by Airbus/DLR), https://spacedata.copernicus.eu/collections/copernicus-digital-elevation-model
- EMODnet Bathymetry Consortium (2024). *EMODnet Digital Bathymetry (DTM 2024)*. https://doi.org/10.12770/cf51df64-56f9-4a99-b1aa-36b8d7b743a1 (CC BY)
- Sample and site metadata: TREC data portal (EMBL). Publication of TREC data needs the consortium's consent.
- Rendering: forge3d (https://github.com/milos-agathon/forge3d).
