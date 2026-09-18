# Required packages:
# - geopandas: reads the segmentation polygons and computes their geometry
# - rasterio: reads the reference grid and rasterizes (burns) each metric
# - pandas / numpy: array handling and Inf/NaN checks
# - re / datetime: detects and validates dates from file names
#
# math and rasterio.enums.MergeAlg are imported but not used in this script

import os
import re
import geopandas as gpd
import pandas as pd
import numpy as np
import rasterio
from rasterio import features
from rasterio.enums import MergeAlg
import math
from datetime import datetime

# 1. PATH SETTINGS
# path_gpkg: segmentation output from the previous step (step 3)
# dir_referencia: folder with the regularized per-date/per-tile bands from
# the data cube step (step 1) — used only to read each tile's raster grid
# (extent/resolution) and to detect which dates the metrics must be
# duplicated for
# dir_final: folder where the rasterized metrics (one file per tile, metric
# and date) will be written

path_gpkg = os.path.expanduser("path/to/your/segmentation.gpkg")
dir_referencia = os.path.expanduser("path/to/your/datacube_dir/")
dir_final = os.path.expanduser("path/to/your/output_dir")
os.makedirs(dir_final, exist_ok=True)

# Tile identifiers covering the study area, matching the reference file names

tiles = ["001014", "002015", "002014", "002013", "003013", "003015", "003014"]

# Metrics computed from the segmentation polygons and rasterized below

metricas_alvo = ["AREA", "PERIM", "SHAPE"]

# 2. AUTOMATIC DATE DETECTION FROM THE REFERENCE FOLDER
# Matches a date written as YYYY-MM-DD or YYYYMMDD inside a file/folder name

PADRAO_DATA = re.compile(r"(\d{4}-\d{2}-\d{2}|\d{8})")

def normaliza_data(match_str):
    # Converts a matched date string to 'YYYY-MM-DD' and validates that it
    # is a real calendar date

    fmt_entrada = "%Y-%m-%d" if "-" in match_str else "%Y%m%d"
    try:
        return datetime.strptime(match_str, fmt_entrada).strftime("%Y-%m-%d")
    except ValueError:
        return None

def extrai_datas(diretorio, tiles=None):
    # Walks diretorio recursively and returns every distinct date found in
    # file/folder names, optionally restricted to the given tile ids. These
    # are the dates the segmentation metrics will be replicated for, since
    # the segmentation itself is static (from a single mosaic date) but the
    # data cube needs one metric value per time step

    datas_encontradas = set()
    for raiz, dirs, arquivos in os.walk(diretorio):
        for nome in dirs + arquivos:
            if tiles and not any(t in nome for t in tiles):
                continue
            for m in PADRAO_DATA.findall(nome):
                d = normaliza_data(m)
                if d:
                    datas_encontradas.add(d)
    return sorted(datas_encontradas)

datas_16d = extrai_datas(dir_referencia, tiles=tiles)

if not datas_16d:
    raise RuntimeError(
        f"Nenhuma data detectada em {dir_referencia}. "
        "Verifique se o padrao de nomenclatura dos arquivos bate com o regex."
    )

print(f"{len(datas_16d)} datas detectadas em dir_referencia:")
for d in datas_16d:
    print(" ", d)

# 3. LOAD THE VECTOR AND COMPUTE THE METRICS

print("Carregando vetor e calculando métricas espaciais...")
gdf = gpd.read_file(path_gpkg)

area_m2 = gdf.geometry.area

# Drops geometries with no real area (LineString, MultiLineString, Point, or
# degenerate polygons) — otherwise they cause a division by zero -> Infinity
# in SHAPE below, which corrupts the whole raster of the tile it falls into

n_antes = len(gdf)
gdf = gdf[area_m2 > 0].copy()
n_removidos = n_antes - len(gdf)
if n_removidos > 0:
    print(f"Aviso: {n_removidos} geometrias sem área removidas de {n_antes} ({100*n_removidos/n_antes:.3f}%)")

area_m2 = gdf.geometry.area
perim_metros = gdf.geometry.length

# Area in hectares

gdf['AREA'] = area_m2 / 10000

# Shape index, using the standard formula (perimeter in meters, area in m²):
# 1.0 for a circle-like compact patch, increasing with elongation/irregularity

gdf['SHAPE'] = perim_metros / (4 * np.sqrt(area_m2))

# Perimeter saved in kilometers so it fits in int16/float32 without
# overflowing downstream in sits

gdf['PERIM'] = perim_metros / 1000

# Confirms no Inf/NaN remains in any of the three metrics after cleaning

for col in ['AREA', 'PERIM', 'SHAPE']:
    n_inf = np.isinf(gdf[col]).sum()
    n_nan = gdf[col].isna().sum()
    if n_inf > 0 or n_nan > 0:
        print(f"Aviso: coluna {col} ainda tem {n_inf} Inf e {n_nan} NaN após limpeza")

print(f"Vetor carregado com {len(gdf)} polígonos válidos.")
print(f"   AREA  -> min: {gdf['AREA'].min():.2f}  | max: {gdf['AREA'].max():.2f}")
print(f"   PERIM -> min: {gdf['PERIM'].min():.2f}  | max: {gdf['PERIM'].max():.2f}")
print(f"   SHAPE -> min: {gdf['SHAPE'].min():.2f}  | max: {gdf['SHAPE'].max():.2f}")

# 4. RASTERIZE AND DUPLICATE PER DATE
# sat/sen/data_ref_orig identify one reference file per tile (band B02),
# used only to copy its grid (transform/shape) for the output rasters.
# data_ref_orig must match a date that actually exists in dir_referencia —
# here it is the same mosaic date the segmentation (step 3) was built from

sat, sen, data_ref_orig = "SENTINEL-2", "MSI", "2023-08-29"
for tile_id in tiles:
    ref_name = f"{sat}_{sen}_{tile_id}_B02_{data_ref_orig}.tif"
    ref_path = os.path.join(dir_referencia, ref_name)
    if not os.path.exists(ref_path):
        print(f"Aviso: Tile {tile_id} não encontrado em {ref_path}")
        continue
    print(f"\n>> Processando Tile: {tile_id}")
    with rasterio.open(ref_path) as src:
        meta = src.meta.copy()
        out_transform = src.transform
        out_shape = src.shape
        meta.update(dtype=rasterio.float32, count=1, compress='lzw', tiled=True)
        for met_nome in metricas_alvo:
            print(f"   Rasterizando métrica: {met_nome}...")

            # Burns each polygon's metric value into the tile's raster grid

            shapes = ((geom, value) for geom, value in zip(gdf.geometry, gdf[met_nome]))
            burned = features.rasterize(
                shapes=shapes,
                out_shape=out_shape,
                transform=out_transform,
                fill=0,
                all_touched=False,
                dtype=rasterio.float32
            )

            # Writes the same rasterized metric once per detected date, so
            # the metric can be merged into every time step of the data cube

            for data_atual in datas_16d:
                nome_final = f"{sat}_{sen}_{tile_id}_{met_nome}_{data_atual}.tif"
                path_save = os.path.join(dir_final, nome_final)
                with rasterio.open(path_save, 'w', **meta) as dst:
                    dst.write(burned, 1)
            del burned
print("\n" + "-"*50)
print("PROCESSO CONCLUÍDO!")
print(f"Arquivos gerados em: {dir_final}")
