# Required packages:
# - os, re, datetime: filesystem access, filename pattern matching and date parsing
# - psutil: reports RAM usage during processing, for memory monitoring
# - xarray: stacks the per-band mosaics into a single multi-band array
# - rioxarray / rasterio: read, merge, reproject and write georeferenced rasters

import os
import re
import psutil
from datetime import datetime

import xarray as xr
import rioxarray
import rasterio
from rioxarray.merge import merge_arrays

# =============================================================================
# 1. INITIAL SETTINGS
# =============================================================================

# Folder with the regularized Sentinel-2 tiles produced by the create data cube step
# (one file per tile/band/date), and the folder where the daily mosaics will
# be written

diretorio_origem = os.path.expanduser("path/to/your/datacube_dir")
output_dir = os.path.expanduser("path/to/your/output_dir")
os.makedirs(output_dir, exist_ok=True)

# BDC tile identifiers covering the study area. Only files whose name
# contains one of these ids are considered when scanning diretorio_origem

tiles = ["001014", "002015", "002014", "002013", "003013", "003015", "003014"] #These sequence tiles were used in my case, but you need to put your tiles here.

# Matches a date written as YYYY-MM-DD or YYYYMMDD inside a file/folder name.
# Adjust the pattern if the actual naming convention is different

PADRAO_DATA = re.compile(r"(\d{4}-\d{2}-\d{2}|\d{8})")

def normaliza_data(match_str):
    """Converts a matched date string to 'YYYY-MM-DD' and validates that it
    is a real calendar date."""

    if "-" in match_str:
        fmt_entrada = "%Y-%m-%d"
    else:
        fmt_entrada = "%Y%m%d"
    try:
        return datetime.strptime(match_str, fmt_entrada).strftime("%Y-%m-%d")
    except ValueError:
        return None  # discards invalid dates (e.g. 2024-02-30)

def extrai_datas(diretorio, tiles=None):
    """Walks diretorio recursively and returns every distinct date found in
    file/folder names, optionally restricted to the given tile ids."""

    datas_encontradas = set()
    for raiz, dirs, arquivos in os.walk(diretorio):
        nomes = dirs + arquivos
        for nome in nomes:
            # Restricts the scan to the tiles of interest
            if tiles and not any(t in nome for t in tiles):
                continue
            for m in PADRAO_DATA.findall(nome):
                data_norm = normaliza_data(m)
                if data_norm:
                    datas_encontradas.add(data_norm)
    return sorted(datas_encontradas)

datas = extrai_datas(diretorio_origem, tiles=tiles)

print(f"{len(datas)} datas detectadas:")
for d in datas:
    print(" ", d)

# Bands used to build the mosaic (B02 is the reference resolution, 10 m)

bandas_alvo = ["B11", "B8A", "B02"] # used these bands for my visualization; you can select other bands from your data cube.

# Regex fragment combining all tile ids into a single alternation, reused by
# listar_arquivos below

pattern_tiles = "|".join(tiles)


# =============================================================================
# FILE LISTING FUNCTION (recursive, equivalent to R's list.files)
# =============================================================================

def listar_arquivos(diretorio, tile_pattern, banda, data):
    """Returns every file under diretorio whose name matches one of the
    tiles, the given band and the given date."""

    arquivos_encontrados = []
    regex_pattern = re.compile(rf"({tile_pattern}).*_{banda}_.*{data}")

    for root, _, files in os.walk(diretorio):
        for file in files:
            if regex_pattern.search(file):
                arquivos_encontrados.append(os.path.join(root, file))
    return arquivos_encontrados


def log_memoria(rotulo):
    """Logs the process's resident RAM. Called before/after each iteration
    to check memory usage empirically instead of assuming that closing the
    datasets was enough."""

    processo = psutil.Process()
    rss_gb = processo.memory_info().rss / 1e9
    print(f"   [RAM] {rotulo}: {rss_gb:.2f} GB")


# =============================================================================
# 2. LOOP OVER DATES
# =============================================================================
# For each date: merge each band's tiles into one mosaic, resample every
# band to match B02's grid, stack the bands and write a single multi-band
# GeoTIFF

for d in datas:
    print(f"\n--- Processando data: {d} ---")
    log_memoria("início da iteração")

    # Lists the files available for each band on this specific date

    f_list = {}
    falta_banda = False

    for b in bandas_alvo:
        files = listar_arquivos(diretorio_origem, pattern_tiles, b, d)
        f_list[b] = files

        if len(files) == 0:
            print(f"[Aviso] Pulando data {d} pois faltam arquivos para a banda: {b}")
            falta_banda = True
            break

    if falta_banda:
        continue

    # 3. Build the individual mosaics and standardize resolution

    mosaicos_res = {}

    # Reference band (B02, 10 m). merge_arrays mosaics all tiles of this
    # band into a single array

    lista_b02_rasters = [rioxarray.open_rasterio(f, masked=True) for f in f_list["B02"]]
    m_b02 = merge_arrays(lista_b02_rasters)
    mosaicos_res["B02"] = m_b02

    # Closes the individual tile datasets right after the merge — otherwise
    # each one stays open in memory until the garbage collector eventually
    # runs, with no guaranteed timing

    for r in lista_b02_rasters:
        r.close()
    del lista_b02_rasters

    # Processes every other band the same way

    for b in bandas_alvo:
        if b == "B02":
            continue

        lista_banda_rasters = [rioxarray.open_rasterio(f, masked=True) for f in f_list[b]]
        m_banda = merge_arrays(lista_banda_rasters)

        for r in lista_banda_rasters:
            r.close()
        del lista_banda_rasters

        # Resamples to B02's resolution/extent (nearest neighbor preserves
        # discrete/categorical values)

        mosaicos_res[b] = m_banda.rio.reproject_match(
            m_b02, resampling=rasterio.enums.Resampling.nearest
        )

        # m_banda (before reproject_match) is no longer needed

        m_banda.close()

    log_memoria("após merges e reprojeções")

    # 4. Stack the bands and write the output raster

    lista_para_stack = [mosaicos_res[b] for b in bandas_alvo]
    mosaico_final = xr.concat(lista_para_stack, dim="band")
    mosaico_final = mosaico_final.assign_coords(band=bandas_alvo)

    data_formatada = d.replace("-", "_")
    out_name = f"MOSAICO_MLME_{data_formatada}.tif"
    out_path = os.path.join(output_dir, out_name)

    # windowed=True writes the raster in blocks instead of holding the full
    # array in RAM during the write — relevant together with bigtiff="YES",
    # since the output file can be large

    mosaico_final.rio.to_raster(
        out_path,
        tiled=True,
        compress="LZW",
        bigtiff="YES",
        windowed=True,
    )

    # Closes datasets to free RAM before moving to the next date

    m_b02.close()
    for b in mosaicos_res:
        mosaicos_res[b].close()
    mosaico_final.close()

    log_memoria("fim da iteração")
    print(f"Salvo com sucesso: {out_name}")

print("\nProcessamento concluído para todas as datas.")
