# Sentinel-2 + Random Forest Land Cover Classification Pipeline

Pipeline used in this thesis to classify land use and land cover from Sentinel-2
imagery, using the Brazil Data Cube (BDC), the `sits` R package for the data
cube and Random Forest classification, an OBIA (GRASS `i.segment`)
segmentation step for landscape metrics, a PRODES-based forest and non-forest masks, and a
stratified accuracy assessment following Olofsson et al. (2014).


## Pipeline steps

Run in this order. Each script's header comment lists its required packages
and the paths you need to edit (marked with placeholders like
`path/to/your/...`).

| # | Script | What it does |
|---|--------|---------------|
| 1 | `01_create_sentinel2_datacube_bdc.R` | Builds the Sentinel-2 data cube from the BDC for the study area/period, downloads it locally, and runs a spectral mixture model (MLME) to get forest/soil/water fractions. |
| 2 | `02_create_sentinel2_mosaics.py` | Mosaics the regularized tiles into one multi-band GeoTIFF per date. |
| 3 | `03_grass_segmentation.py` | Segments a mosaic date into objects (GRASS `i.segment`, region growing), tiled and parallelized. |
| 4 | `04_landscape_metrics.py` | Computes AREA/PERIM/SHAPE from the segmentation polygons and rasterizes them, replicated across every date of the data cube. |
| 5 | `05_merge_metrics_into_datacube.R` | Regularizes the landscape metrics cube and merges it into the main 16-day data cube. **Depends on objects from script 1** (`ROI`, `output_dir`, `s2_cube_local`) in the same R session. |
| 6 | `06_fast_load_train_and_classify.R` | Loads a saved cube (`.RData`), trains a Random Forest on labeled samples, classifies the cube, smooths the probabilities (Bayesian), and labels the final map. |
| 7 | `07_mosaic_classification_tiles.py` | Mosaics the classified tiles into a single raster. |
| 8 | `08_apply_prodes_mask.py` | Reclassifies secondary-vegetation pixels into FLO/NFLO using PRODES forest/non-forest masks. |
| 9 | `09_generate_stratified_sampling_points.R` | Generates stratified validation points (Olofsson/Cochran sample-size equations), eroding class borders to avoid mixed pixels. |
| 10 | `10_extract_classification_at_points.R` | Extracts the classified value at each validation point and translates class codes to names. |
| 11 | `11_olofsson_accuracy_assessment.R` | Computes the stratified accuracy assessment (overall/user's/producer's accuracy with SE and 95% CI), confusion matrices, error composition and patch-size charts. |

## Requirements

**R**: `sits`, `sf`, `terra`, `raster`, `dplyr`, `tidyr`, `ggplot2`, `stringr`,
`caret`, `readr`, `tmap`, `tibble`, `tools`, `magrittr`

**Python**: `rasterio`, `geopandas`, `numpy`, `pandas`, `rioxarray`, `xarray`,
`shapely`, `psutil`

**Other**: GRASS GIS, available in a conda environment named `grass_obia`
(used by script 3 via `conda run -n grass_obia`)

## Notes

- Every path in these scripts is a placeholder (e.g. `path/to/your/...`) —
  replace them with your own directories before running.
- Class dictionaries, example tile lists and statistical parameters (e.g.
  the Olofsson target standard error) are kept as working examples; adjust
  them to your own study area and classes.
- Scripts 1 and 5 must run in the same R session (script 5 reuses objects
  created by script 1).
