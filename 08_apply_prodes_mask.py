# Required packages:
# - rioxarray / xarray: reads, reprojects and writes georeferenced rasters
# - numpy: builds the reclassification matrix
# - os / glob: path handling and file discovery

import os
import glob
import xarray as xr
import rioxarray
import numpy as np

# 1. PATHS
# diretorio_mosaicos: folder with the classification mosaics (step 7)
# diretorio_mascaras: folder with the two PRODES reference masks
# diretorio_saida: where the masked ("MASK_...") mosaics will be written

diretorio_mosaicos = os.path.expanduser("path/to/your/mosaic_dir")
diretorio_mascaras = os.path.expanduser("path/to/your/prodes_masks_dir")
diretorio_saida    = os.path.expanduser("path/to/your/output_dir")

# Ensures the output folder exists

os.makedirs(diretorio_saida, exist_ok=True)

# 2. Load the original PRODES masks
# mask_flo_path: PRODES forest mask (pixels = 1 mean "forest")
# mask_nf_path: PRODES non-forest mask (pixel code 3 = deforested/non-forest,
# per the PRODES legend)

mask_flo_path = os.path.join(diretorio_mascaras, "your_forest_mask.tif")
mask_nf_path  = os.path.join(diretorio_mascaras, "your_non_forest_mask.tif")

print("Carregando máscaras PRODES...")
mask_floresta_raw = rioxarray.open_rasterio(mask_flo_path, masked=False)
mask_nao_flo_raw  = rioxarray.open_rasterio(mask_nf_path, masked=False)

# 3. List every classification mosaic to process
arquivos_mosaico = glob.glob(os.path.join(diretorio_mosaicos, "MOSAICO_*.tif"))

if not arquivos_mosaico:
    print(f"Nenhum mosaico encontrado em: {diretorio_mosaicos}")
    exit()

print(f"Total de mosaicos encontrados para processar: {len(arquivos_mosaico)}")

# ==============================================================================
# OPTIMIZATION: crop the PRODES masks to the combined extent of all mosaics
# up front, once, instead of reprojecting the full PRODES scene (which can
# cover an entire state) inside the loop for every mosaic
# ==============================================================================
print("Calculando limites geográficos para corte otimizado...")
min_x, min_y, max_x, max_y = float("inf"), float("inf"), float("-inf"), float("-inf")
for cam in arquivos_mosaico:
    with rioxarray.open_rasterio(cam) as r:
        b = r.rio.bounds()
        min_x, min_y = min(min_x, b[0]), min(min_y, b[1])
        max_x, max_y = max(max_x, b[2]), max(max_y, b[3])

# Crops the raw masks down to the study area's bounding box

crs_alvo = rioxarray.open_rasterio(arquivos_mosaico[0]).rio.crs
mask_floresta_crop = mask_floresta_raw.rio.clip_box(min_x, min_y, max_x, max_y, crs=crs_alvo)
mask_nao_flo_crop  = mask_nao_flo_raw.rio.clip_box(min_x, min_y, max_x, max_y, crs=crs_alvo)
# ==============================================================================

# 4. Processing loop
for caminho_mosaico in sorted(arquivos_mosaico):
    nome_base = os.path.basename(caminho_mosaico)

    # Skips files that are already a masking output, in case this script is
    # re-run over the same output folder

    if nome_base.startswith("MASK_"):
        continue

    print(f"\n--- Processando com Regras Customizadas (NDWI): {nome_base} ---")

    mosaico = rioxarray.open_rasterio(caminho_mosaico, masked=False)
    mosaico.rio.write_nodata(0, inplace=True)

    # Reprojects/resamples the (already cropped) masks onto this mosaic's
    # exact grid. resampling=0 is nearest neighbor, required for categorical
    # data

    print("Alinhando máscaras PRODES...")
    m_flo_res = mask_floresta_crop.rio.reproject_match(mosaico, resampling=0) # 0 = nearest
    m_nf_res  = mask_nao_flo_crop.rio.reproject_match(mosaico, resampling=0)

    # --- RECLASSIFICATION RULES ---
    # Classes 11 and 12 are the secondary vegetation classes (VSA/VSI) in
    # this study's class dictionary. Any such pixel that falls inside PRODES
    # forest becomes 13 (FLO); any pixel inside PRODES non-forest (mask code
    # 3) becomes 14 (NFLO), regardless of its original class. Everything
    # else keeps its original classified value

    print("Aplicando matriz de decisão NumPy...")
    is_secondary_veg = (mosaico.values == 11) | (mosaico.values == 12)

    mosaico_final_values = np.where(
        m_nf_res.values == 3, 14,
        np.where(
            (is_secondary_veg) & (m_flo_res.values == 1), 13,
            mosaico.values
        )
    )

    # Builds a new DataArray with the reclassified values, keeping the
    # original raster's metadata (CRS, transform, etc.)

    mosaico_final = mosaico.copy(data=mosaico_final_values)
    mosaico_final.rio.write_nodata(0, inplace=True)
    mosaico_final = mosaico_final.astype(np.uint8)

    # 5. Save to the output folder

    nome_saida = f"MASK_{nome_base}"
    caminho_final = os.path.join(diretorio_saida, nome_saida)

    print(f"Gravando arquivo compactado: {nome_saida}")
    mosaico_final.rio.to_raster(
        caminho_final,
        dtype="uint8",
        compress="lzw",
        tiled=True,
        windowed=True # Evita estouro de RAM ao escrever no disco
    )

    # Frees this iteration's rasters before moving to the next mosaic

    mosaico.close()
    m_flo_res.close()
    m_nf_res.close()
    print(f"Sucesso!")

# Closes the base mask datasets

mask_floresta_raw.close()
mask_nao_flo_raw.close()

print("\nProcesso concluído com sucesso para a pasta NDWI!")
