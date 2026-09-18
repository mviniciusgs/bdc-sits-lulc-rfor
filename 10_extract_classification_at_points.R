# Required packages:
# - sf: reads/writes the validation points as a spatial GeoPackage
# - terra: reprojects the classified rasters and extracts values at points
# - dplyr: table wrangling (rename, mutate/across) on the extracted data
# - tools: strips file extensions to build clean raster layer names
#
# Extracts the classified value at each validation point (generated in step
# 9), for every classification raster found in raster_dir, translating the
# numeric class codes to their class names. Output feeds step 10's accuracy
# assessment (Olofsson script)

library(sf)
library(terra)
library(dplyr)
library(tools)

# ---- Settings ----
# pontos_file: validation points (step 9 output, after any manual review)
# raster_dir: folder with the masked classification raster(s) (step 8 output)
# output_file: points with the extracted class values, ready for validation

pontos_file <- "path/to/your/validation_points.gpkg"
raster_dir  <- "path/to/your/classification_dir/"
output_file <- "path/to/your/validation_points_extracted.gpkg"

# Temporary folder for the reprojected rasters, removed at the end

temp_dir <- file.path(tempdir(), "raster_reproj_temp")
if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)

# Class dictionary: maps each raster code to its class name, including the
# nodata codes used across the pipeline (0 and 255)

dicionario_classes <- c(
  "1"  = "AGPE", "2"  = "AGUA", "3"  = "AL",   "4"  = "DMC",  "5"  = "DMF",
  "6"  = "PA",   "7"  = "PH",   "8"  = "PSI",  "9"  = "SE",   "10" = "URB",
  "11" = "VSA",  "12" = "VSI",  "13" = "FLO",  "14" = "NFLO", "0"  = "NoData",
  "255"= "NoData"
)

# ---- 1. Load the points ----

message("Carregando pontos de validação...")
pontos <- sf::st_read(pontos_file, quiet = TRUE)

if ("fid" %in% names(pontos)) {
  message("Coluna 'fid' detectada nos pontos originais. Renomeando para evitar conflito com o GeoPackage...")
  pontos <- pontos %>% dplyr::rename(fid_original = fid)
}

pontos_crs <- st_crs(pontos)
cat("CRS dos pontos:", pontos_crs$input, "\n")

pontos_vect <- terra::vect(pontos)

# ---- 2. List the source rasters ----
# Only files starting with "MASK" are picked up, matching step 8's output
# naming convention

raster_paths <- list.files(
  path = raster_dir,
  pattern = "^MASK.*\\.tif$",
  full.names = TRUE
)

if (length(raster_paths) == 0) {
  stop("Nenhum arquivo raster encontrado. Verifique o caminho da pasta.")
}

message("Encontrados ", length(raster_paths), " mapeamentos.")

# ---- 3. Reproject each raster to the points' CRS ----
# terra::extract requires the raster and the points to share a CRS; rasters
# already in the right CRS are just copied, to keep the stacking step below
# uniform regardless of what was reprojected

message("Reprojetando rasters para o CRS dos pontos (", pontos_crs$input, ")...")

raster_paths_reproj <- character(length(raster_paths))

for (i in seq_along(raster_paths)) {
  r <- rast(raster_paths[i])
  nm <- basename(raster_paths[i])

  r_crs <- crs(r)
  if (r_crs == "" || is.na(r_crs)) {
    stop(sprintf("Raster '%s' não tem CRS definido -- não é possível reprojetar com segurança.", nm))
  }

  out_temp <- file.path(temp_dir, nm)

  if (st_crs(r_crs) == pontos_crs) {
    cat(sprintf("[%s] CRS já compatível -- copiando sem reprojetar.\n", nm))
    file.copy(raster_paths[i], out_temp, overwrite = TRUE)
  } else {
    cat(sprintf("[%s] Reprojetando (method = 'near', preserva classes categóricas)...\n", nm))
    r_proj <- project(r, y = crs(pontos_crs$wkt), method = "near")
    writeRaster(r_proj, out_temp, overwrite = TRUE)
  }

  raster_paths_reproj[i] <- out_temp
}

# ---- 4. Stack the reprojected rasters ----
# Each raster becomes one layer/column, named after its file (without
# extension), so multiple classification versions can be extracted at once

message("Empilhando rasters reprojetados...")
rasters <- terra::rast(raster_paths_reproj)
nomes_limpos <- tools::file_path_sans_ext(basename(raster_paths_reproj))
names(rasters) <- nomes_limpos

# ---- 5. Extract the values ----

message("Extraindo valores dos rasters para os pontos...")
extracao <- terra::extract(rasters, pontos_vect, ID = FALSE)

na_por_coluna <- sapply(extracao, function(x) sum(is.na(x)) / length(x) * 100)
cat("\n% de NA por coluna extraída:\n")
print(round(na_por_coluna, 2))

# ---- 6. Join the data and translate the codes ----
# Replaces each raster's numeric class code with its class name, using
# dicionario_classes, for every extracted column at once

message("Traduzindo códigos numéricos para nomes das classes...")
pontos_final <- dplyr::bind_cols(pontos, extracao) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(nomes_limpos),
      ~ dicionario_classes[as.character(.)]
    )
  )

# ---- 7. Save the result ----

message("Salvando novo GeoPackage...")
if (!dir.exists(dirname(output_file))) {
  dir.create(dirname(output_file), recursive = TRUE)
}

sf::st_write(pontos_final, output_file,
             append = FALSE,
             delete_dsn = TRUE,
             quiet = TRUE)

message("GeoPackage salvo com sucesso!")

# ---- 8. Remove the temporary reprojected rasters ----

message("Removendo rasters temporários reprojetados...")
unlink(temp_dir, recursive = TRUE, force = TRUE)

if (!dir.exists(temp_dir)) {
  message("Pasta temporária removida com sucesso: ", temp_dir)
} else {
  warning("Falha ao remover a pasta temporária: ", temp_dir)
}

# ---- Final check: preview the first extracted/translated columns ----

num_colunas_print <- min(6, length(nomes_limpos))
print(head(pontos_final %>% sf::st_drop_geometry() %>% dplyr::select(dplyr::all_of(nomes_limpos[1:num_colunas_print]))))
