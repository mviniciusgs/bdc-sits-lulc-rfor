# Required packages:
# - sits: trains the model and runs classification/smoothing/labeling
# - readr: reads the training samples CSV
# - dplyr: sample table wrangling (%>%, mutate, rename)
# - raster, tmap: loaded for interactive inspection/plotting of results;
#   not called anywhere in this script
#
# "Fast load": instead of rebuilding the data cube from steps 1-5 every time,
# this script loads a previously saved cube object (reg_cube_metrics) from an
# .RData file, so training/classification can be repeated quickly

library(sits)
library(raster)
library(readr)
library(tmap)
library(dplyr)

# Cube metadata saved from a previous run of steps 1-5 (must contain a
# reg_cube_metrics object)

load("path/to/your/datacube_metadata.RData")


# Folder where the classification outputs will be written

output_dir_class <- "path/to/your/output_dir"

# --- RUN SETTINGS ---
# version: label used to name every output file from this run, so different
# parameter sets don't overwrite each other
# START_DATE/END_DATE: time range assigned to every training sample below

version <-  "ST2A-DB_FULL-16d-rfor-2024-V3-RG"
START_DATE <- "2024-01-01"
END_DATE  <- "2024-12-31"

# Training samples table: one row per point, with at least a "class" column
# and coordinates/geometry recognized by sits_get_data

amostras_2024 <- read_csv("path/to/your/training_samples.csv")

# Adapts the samples table to what sits expects: a start/end date per sample
# (here, the same for all samples) and a "label" column (sits_get_data reads
# the class from "label", not "class")

amostras_2024_sits <- amostras_2024  %>%
  mutate(
    start_date = as.Date(START_DATE),
    end_date   = as.Date(END_DATE)
  ) %>%
  rename(
    label = class
  )


# Extracts, for each sample point, its time series from every band of
# reg_cube_metrics (loaded above from the .RData file)

all_samples_2024 <- sits_get_data(
  cube       = reg_cube_metrics,
  samples    = amostras_2024_sits,
  multicores = 28,
  memsize    = 180,
  progress   = TRUE
)


# Trains the classifier
# set.seed fixes the Random Forest's internal randomness, so results are
# reproducible across runs

set.seed(146)
model2_rfor <- sits_train(
  all_samples_2024,
  ml_method = sits_rfor(num_trees = 200)
)

# Classifies the full cube with the trained model, producing one class
# probability per pixel/date (probs_cube)

probs_cube <- sits_classify(
  data = reg_cube_metrics, # Troque aqui o tile
  ml_model = model2_rfor,
  output_dir = output_dir_class,
  version = version,
  multicores = 10,
  memsize = 180,
  progress = FALSE,
  verbose = TRUE
)

# Bayesian smoothing of the class probabilities: reduces salt-and-pepper
# noise by borrowing information from each pixel's spatial neighborhood
# before the final labeling step

probs_bayes <- sits_smooth(
  cube = probs_cube,
  window_size = 13, # Janela
  neigh_fraction = 0.5,
  smoothness = 20,
  multicores = 28,
  memsize = 180,
  output_dir = output_dir_class, # Local de saída
  version = paste(version, "w13-nf05-s10-V1", sep = "-") # Descrição do smooth
)

# Converts the smoothed probabilities into the final classified map, picking
# the most likely class per pixel/date

class_cube <- sits_label_classification(
  cube = probs_bayes,
  multicores = 28,
  memsize = 180,
  output_dir = output_dir_class, # Local de saída
  version = paste(version, "w13-nf05-s10-V1", sep = "-")
)
