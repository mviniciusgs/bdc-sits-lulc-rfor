# Required packages:
# - sits: builds and processes the satellite image data cube
# - sf: reads the vector file used as the study area boundary (ROI)
# - tibble: builds the endmember table used by the mixture model
# - magrittr: pipe operator support (%>%), used by sits internally

library(magrittr)
library(sits)
library(tibble)
library(sf)

# Output directory where the regularized (downloaded) data cube will be stored

dir.create("path/to/your/output_dir", recursive = TRUE)
output_dir <- "path/to/your/output_dir"


# Prints the Sentinel-2 collections available from the Brazil Data Cube (BDC)
# provider, so you can confirm "SENTINEL-2-16D" (used below) is valid/active

sits_list_collections(source = "BDC")

# Study area boundary (polygon). Used to filter which images/tiles are
# downloaded, instead of manually listing tiles

ROI <- sf::st_read("path/to/your_roi.shp")

# Sentinel-2 MGRS tile identifiers covering the study area.
# Kept here for reference/documentation only: the cube below is built from
# ROI, not from this list, so TILES has no effect on the query

TILES <- c("18LXR", "18LYR", "18LZR", "18MXS", "18MXT", "18MYS", "18MYT", "18MZS", "18MZT", "19LBL", "19MBM", "19MBN", "19MCN")

# Bands/indices requested from the collection: 3 spectral indices (NDVI, EVI,
# NBR) plus the 10 Sentinel-2 surface reflectance bands used in this study

BANDAS <- c("NDVI","EVI", "NBR", "B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12")



# Queries the BDC STAC catalog and builds the (remote) cube metadata for the
# study period/area, without downloading imagery yet.
# To reuse this script for another area, replace the ROI file path above
# with your own study area boundary file

s2_cube <- sits_cube(
  source     = "BDC",
  collection = "SENTINEL-2-16D",
  bands      =  BANDAS,
  roi = ROI,
  start_date = as.Date("2024-01-01"),
  end_date   = as.Date("2024-12-31"),
  multicore = 8L
)

# Extends the default download timeout (seconds) to avoid failures on large,
# slow transfers from the BDC

options(timeout = 10000)




# Downloads the imagery listed in s2_cube, cropped to the ROI, and writes it
# to output_dir. n_tries retries each failed file download up to 10 times
reg_cube <- sits_cube_copy(
  cube       = s2_cube,
  output_dir = output_dir,
  roi        = ROI,
  n_tries = 10L,
  multicore  = 8L
)



#### MLME
# Spectral Mixture Model: decomposes each pixel into fractions of forest,
# soil and water (shadows) based on the reference spectra (endmembers) below

# Endmember reference spectra, one row per land cover type, with the
# reflectance value expected for each Sentinel-2 band (surface reflectance,
# scaled by 10000, i.e. sits/BDC convention).

em <- tibble::tribble(
  ~type, ~B02, ~B03, ~B04, ~B05, ~B06, ~B07,~B08, ~B8A, ~B11, ~B12, # CHANGE BANDS

  "forest",  164, 480,  194, 808, 2582, 3149, 3536, 3466, 1497, 598,
  "soil",    839, 1116, 1278, 1614, 2142, 2361, 2474, 2754, 3198, 2315,
  "water",   110, 193,  119, 108, 26, 35, 10, 10, 26, 24 )

# Applies the mixture model to reg_cube, producing one output band per
# endmember (forest/soil/water fraction) for every date in the cube.
# memsize is the RAM budget (GB) sits is allowed to use during processing

s2_cube_local <- sits_mixture_model(
  data = reg_cube,
  endmembers = em,
  multicores = 20,
  memsize = 150,
  output_dir = output_dir)
