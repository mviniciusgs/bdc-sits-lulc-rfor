# Required packages:
# - rasterio: reads the input mosaic and writes the per-tile crops
# - geopandas / shapely: reads, clips and merges the segmented polygons
# - pandas: concatenates the per-tile results into a single GeoDataFrame
# - numpy: normalizes pixel values before segmentation
# GRASS GIS is also required, installed in a conda environment named
# "grass_obia" (see the "conda run -n grass_obia" calls below). This script
# drives GRASS's i.segment (region growing) through the command line instead
# of a Python GRASS binding

import os
import subprocess
import shutil
import rasterio
from rasterio.windows import Window
import geopandas as gpd
import pandas as pd
import numpy as np
import warnings
from shapely.geometry import box
from shapely.errors import ShapelyDeprecationWarning
from concurrent.futures import ProcessPoolExecutor, as_completed

# Silences noisy internal Geopandas/Shapely warnings unrelated to the result

warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# =========================================================================
# 1. SEGMENTATION PARAMETERS
# =========================================================================
# The input mosaic is split into square tiles and each tile is segmented
# independently (in parallel), then the resulting polygons are merged back
# into a single vector file

tile_size = 3000       # tile size in pixels (before adding the buffer)
buffer_size = 500      # extra margin (pixels) read around each tile, so
                        # segments are not cut artificially at tile edges
vMemory = 4000         # RAM limit per GRASS worker, in MB
max_workers = 28       # number of tiles segmented in parallel
vBands = [1, 2, 3]     # mosaic band indices (1-based) used for segmentation (I used ["B11", "B8A", "B02"])

vScale = 0.7           # i.segment similarity threshold: higher merges more,
                        # producing larger/fewer segments
vMinSize = 20          # minimum segment size, in pixels
# =========================================================================

# 2. Paths
# Replace path_img with your mosaic file and dir_out with where the
# segmentation output (and temporary GRASS data) should be written

path_img = os.path.expanduser("path/to/your/mosaic.tif")
dir_out = os.path.expanduser("path/to/your/output_dir/")
os.makedirs(dir_out, exist_ok=True)

str_bandas = "_".join(map(str, vBands))
nome_arquivo = f"AGPE_Ab1_Standard_vScale_{vScale}_Bandas_{str_bandas}.gpkg"
path_out_final = os.path.join(dir_out, nome_arquivo)

dir_temp = os.path.join(dir_out, "temp_tiles")
os.makedirs(dir_temp, exist_ok=True)
grass_db = os.path.join(dir_temp, "grassdata")

def process_tile_with_overlap(row_off, col_off, tile_size, buffer_size, vBands, path_img, dir_temp, grass_db):
    # Runs the full segmentation of a single tile: crop -> GRASS location ->
    # i.segment -> vectorize -> clip back to the tile's real (unbuffered)
    # extent. Executed once per tile, in a separate worker process

    idx = f"{row_off}_{col_off}"

    path_crop = os.path.join(dir_temp, f"crop_{idx}.tif")
    path_gpkg = os.path.join(dir_temp, f"vec_{idx}.gpkg")
    loc_temp = os.path.join(grass_db, f"loc_{idx}")

    try:
        with rasterio.open(path_img) as src:
            w, h = src.width, src.height

            # 2.1 Exact tile bounds (NO buffer), used later to clip the
            # segmented output back to this tile's real extent

            width_unbuffered = min(tile_size, w - col_off)
            height_unbuffered = min(tile_size, h - row_off)
            window_unbuffered = Window(col_off, row_off, width_unbuffered, height_unbuffered)
            transform_unbuffered = src.window_transform(window_unbuffered)

            tl_x, tl_y = transform_unbuffered * (0, 0)
            br_x, br_y = transform_unbuffered * (width_unbuffered, height_unbuffered)
            unbuffered_box = box(min(tl_x, br_x), min(tl_y, br_y), max(tl_x, br_x), max(tl_y, br_y))

            # 2.2 Expanded bounds (WITH buffer): this is the window actually
            # read and segmented, so segments near the tile border are formed
            # using real neighboring pixels instead of an artificial edge

            row_start = max(0, row_off - buffer_size)
            row_end = min(h, row_off + height_unbuffered + buffer_size)
            col_start = max(0, col_off - buffer_size)
            col_end = min(w, col_off + width_unbuffered + buffer_size)

            buf_width = col_end - col_start
            buf_height = row_end - row_start
            window_buffered = Window(col_start, row_start, buf_width, buf_height)
            transform_buffered = src.window_transform(window_buffered)

            # Reads the raster window

            crop_data = src.read(vBands, window=window_buffered)
            nodata_val = src.nodata if src.nodata is not None else -9999

            # Skips tiles that are entirely nodata or blank

            if (crop_data == nodata_val).all() or crop_data.max() <= 0:
                return None

            # Global normalization to 0.0-1.0, assuming input reflectance
            # scaled by 10000 (sits/BDC convention)

            crop_data = crop_data.astype('float32')
            crop_data[crop_data < 0] = 0.0
            crop_data = np.clip(crop_data / 10000.0, 0.0, 1.0)

            kwargs = src.meta.copy()
            kwargs.update({
                'height': buf_height,
                'width': buf_width,
                'transform': transform_buffered,
                'count': len(vBands),
                'dtype': 'float32',
                'nodata': 0.0
            })

            with rasterio.open(path_crop, 'w', **kwargs) as dst:
                dst.write(crop_data)

        if os.path.exists(loc_temp):
            shutil.rmtree(loc_temp)

        # Creates a throwaway GRASS location/mapset for this tile, matching
        # the cropped raster's projection

        subprocess.run([
            "conda", "run", "-n", "grass_obia",
            "grass", "-c", path_crop, loc_temp, "-e"
        ], capture_output=True, text=True, check=True)

        # Segmentation with Queen contiguity (-d) and normalized band weights
        # (-w): imports the crop, groups its bands, runs i.segment, then
        # converts the resulting segment raster to vector polygons

        comandos_grass = f"""
        r.import input="{path_crop}" output=matriz_{idx} --overwrite
        i.group group=grupo_{idx} input=$(g.list type=raster pattern="matriz_{idx}.*" separator=",") --overwrite
        i.segment -w -d group=grupo_{idx} output=seg_{idx} threshold={vScale} minsize={vMinSize} memory={vMemory} --overwrite
        r.to.vect input=seg_{idx} output=vetor_{idx} type=area --overwrite
        v.out.ogr input=vetor_{idx} output="{path_gpkg}" format=GPKG --overwrite
        """

        mapset_path = os.path.join(loc_temp, "PERMANENT")
        subprocess.run([
            "conda", "run", "-n", "grass_obia",
            "grass", mapset_path, "--exec", "sh", "-c", comandos_grass
        ], capture_output=True, text=True, check=True)

        # Clips the segmented polygons back to the tile's unbuffered extent,
        # discarding the buffer area so neighboring tiles don't overlap

        if os.path.exists(path_gpkg):
            gdf = gpd.read_file(path_gpkg)
            if not gdf.empty:
                bbox_gdf = gpd.GeoDataFrame(geometry=[unbuffered_box], crs=gdf.crs)
                gdf_filt = gpd.clip(gdf, bbox_gdf)
                return gdf_filt if not gdf_filt.empty else None

        return None

    except subprocess.CalledProcessError as e:
        print(f"\n[ERRO GRASS - Bloco {idx}] Falha no motor C.\nLog:\n{e.stderr}\n")
        return None
    except Exception as e:
        print(f"\n[ERRO PYTHON - Bloco {idx}]: {str(e)}\n")
        return None
    finally:
        # Always cleans up this tile's temporary GRASS location and files,
        # regardless of success or failure

        if os.path.exists(loc_temp): shutil.rmtree(loc_temp)
        if os.path.exists(path_crop): os.remove(path_crop)
        if os.path.exists(path_gpkg): os.remove(path_gpkg)


if __name__ == '__main__':
    print("==================================================================")
    print(" INICIANDO OBIA STANDARD: NORMALIZAÇÃO GLOBAL + BUFFER DE RECORTES")
    print("==================================================================")

    # Builds the list of tile origins (row/col offsets) covering the whole
    # mosaic, based on tile_size

    tiles_to_process = []
    with rasterio.open(path_img) as src:
        w, h = src.width, src.height
        for row_off in range(0, h, tile_size):
            for col_off in range(0, w, tile_size):
                tiles_to_process.append((row_off, col_off))

    lista_geometrias = []
    total_tiles = len(tiles_to_process)
    tiles_concluidos = 0

    # Segments every tile in parallel (max_workers processes), collecting
    # each tile's polygons as its worker finishes

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(process_tile_with_overlap, r, c, tile_size, buffer_size, vBands, path_img, dir_temp, grass_db): (r, c)
            for r, c in tiles_to_process
        }

        for future in as_completed(futures):
            r, c = futures[future]
            tiles_concluidos += 1
            try:
                gdf_result = future.result()
                if gdf_result is not None:
                    lista_geometrias.append(gdf_result)
                    print(f"[{tiles_concluidos}/{total_tiles}] Segmentação concluída: Bloco Linha {r}, Coluna {c}")
                else:
                    print(f"[{tiles_concluidos}/{total_tiles}] Ignorado (NoData/Fundo): Bloco Linha {r}, Coluna {c}")
            except Exception as exc:
                print(f"[{tiles_concluidos}/{total_tiles}] Interrupção crítica no bloco {r}_{c}: {exc}")

    print("\nIniciando junção vetorial e exportação...")

    # Merges every tile's polygons into the single final GeoPackage

    if lista_geometrias:
        gdf_final = pd.concat(lista_geometrias, ignore_index=True)
        gdf_final.to_file(path_out_final, driver="GPKG")
        print(f"-> Operação concluída. Malha vetorial salva: {nome_arquivo}")
    else:
        print("-> Operação abortada. Nenhuma matriz válida gerou geometrias.")

    # Removes the temporary tile crops and GRASS locations/mapsets

    if os.path.exists(dir_temp):
        shutil.rmtree(dir_temp)
    print("\nProcessamento finalizado.")
