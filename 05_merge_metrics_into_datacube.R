# Required packages:
# - sits: builds/regularizes the metrics cube and merges it with the
#   spectral+MLME cube
# - magrittr: pipe operator support (%>%), used by sits internally
#
# This script assumes it runs in the same R session as script 01
# (01_create_sentinel2_datacube_bdc.R), right after it: it reuses ROI,
# output_dir and s2_cube_local from that script without redefining them here

library(sits)
library(magrittr)



# 1. PARAMETERS

period <- "P16D"  # ISO 8601 duration: regularize the metrics cube to a
                   # 16-day step, matching the main data cube's time step
output_dir_metrics <- "path/to/your/temp_dir/"
dir_metrics <- "path/to/your/metrics_dir"



if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------------------------------------------------------
# 2. LOAD THE TWO CUBES
# ---------------------------------------------------------

# Loads the rasterized metrics (AREA/PERIM/SHAPE, from step 4) as a raw local
# cube. parse_info tells sits how to read satellite/sensor/tile/band/date out
# of each file name (e.g. SENTINEL-2_MSI_<tile>_<band>_<date>.tif)

cube_metrics <- sits_cube(
  source = "BDC",
  collection = "SENTINEL-2-16D",
  data_dir = dir_metrics,
  parse_info = c("satellite", "sensor", "tile", "band", "date")
)

# Regularizes the metrics cube to the same 16-day period/resolution/ROI as
# the main data cube, overwriting dir_metrics with the regularized version

cube_metrics <- sits_regularize(
  cube = cube_metrics,
  period = period,
  roi  = ROI,
  res = 10,
  multicores = 28,
  memsize = 180,
  progress = TRUE,
  output_dir = dir_metrics
)

# ---------------------------------------------------------
# 3. MERGE THE METRICS INTO THE 16-DAY CUBE
# ---------------------------------------------------------
# s2_cube_local (from step 1) is already the 16-day target cube, so
# sits_merge is used instead of a second regularize call. sits_merge takes
# the metrics (static in time, one value per pixel) and replicates them
# across every date of the 16-day cube automatically

message("Integrando métricas estáticas na linha do tempo de 16 dias...")
reg_cube_metrics <- sits_merge(s2_cube_local, cube_metrics)



# Keeps only the bands actually needed downstream (step 6 training/
# classification): spectral bands, MLME fractions, spectral indices and the
# landscape metrics

reg_cube_metrics <- sits_select(
  data  = reg_cube_metrics,
  bands = c("B02","B03","B04","B05","B06","B07","B08","B11","B12","B8A",
            "AREA","PERIM","SHAPE","NDVI","NBR","EVI","FOREST","SOIL","WATER")
) #Set your bands of interest
