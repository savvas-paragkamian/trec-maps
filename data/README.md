# data/

Nothing in this folder is committed except this file (TREC data are shared only with
the consortium's consent; the downloads are rebuilt by the scripts).

- `TREC_sample_data_export-all-MetaBGT.xlsx` (soil and sediment samples of every LSI site) and
  `TREC_site_data_export_all-sites_marine_matagenome.xlsx` (seawater samples: Tara stations and
  LSI shoreline water): exports from the TREC data portal. Title in row 1, header in row 2;
  the site table (positions, elevation, bathymetry, station names, related sites) is derived
  from their `Site ...` columns. Older exports (`TREC_sample_data_export.xlsx`,
  `TREC_sample_data_athens_marine_metagenome.xlsx`) are no longer read.
- `raw/dem/`: Copernicus GLO-30 windows per transect and tile, `scripts/00_get_dem.py`
  (the four full Athens tiles from the first runs can be deleted).
- `raw/bathy/`: EMODnet Bathymetry WCS extracts, `scripts/01_get_bathymetry.py`.
- `raw/landcover/`: ESA WorldCover windows, `scripts/02_get_landcover.py`.
- `derived/_cache/`: the exports as pickles (rebuilt when an xlsx is newer).
- `derived/<transect>/`: merged land-sea grid (`z.tif`, `grid.npz`, `grid.json`),
  profile tables (`profile.csv`, `land_samples.csv`, `stations.csv`) and render caches.
