# Required packages:
# - rasterio: reads the classified tiles and merges them into one mosaic
# - pathlib / contextlib: path handling and safe closing of open rasters

import contextlib
from pathlib import Path
import rasterio
from rasterio.merge import merge

# 1. MAIN PATHS
# diretorio_raiz: base folder for this run's outputs (step 6)
# diretorio_saida: where the merged mosaic(s) will be written
# pastas_alvo: one entry per folder to mosaic; each entry produces its own
# "MOSAICO_<folder name>.tif". Add more paths here to mosaic several
# classification runs/versions in the same execution

diretorio_raiz = Path("path/to/your/classification_dir")
diretorio_saida = diretorio_raiz / "MOSAIC"

pastas_alvo = [
    diretorio_raiz / "CLASS",
]

diretorio_saida.mkdir(parents=True, exist_ok=True)

# 2. MAIN LOOP - one mosaic per folder in pastas_alvo
for pasta in pastas_alvo:
    if not pasta.exists():
        print(f"Aviso: Pasta não encontrada -> {pasta}")
        continue

    print(f"--- Processando: {pasta.name} ---")

    # Selects the classified tiles to merge: any .tif with "class" in the name

    arquivos_class = list(pasta.glob("*class*.tif"))

    if not arquivos_class:
        print(f"  Aviso: Nenhum arquivo 'class*.tif' encontrado em {pasta}")
        continue

    print(f"  Encontrados {len(arquivos_class)} arquivos. Criando mosaico...")

    nome_saida = f"MOSAICO_{pasta.name}.tif"
    caminho_final_arquivo = diretorio_saida / nome_saida

    # 3. Mosaic operation with safe memory handling (ExitStack)
    # ExitStack keeps every tile open only for as long as merge() needs them,
    # then closes all of them together as soon as the block exits
    try:
        with contextlib.ExitStack() as stack:
            rasters_abertos = [stack.enter_context(rasterio.open(f)) for f in arquivos_class]

            mosaico_final, out_transform = merge(rasters_abertos)

            out_meta = rasters_abertos[0].meta.copy()

        out_meta.update({
            "driver": "GTiff",
            "height": mosaico_final.shape[1],
            "width": mosaico_final.shape[2],
            "transform": out_transform
        })

        with rasterio.open(caminho_final_arquivo, "w", **out_meta) as dest:
            dest.write(mosaico_final)

        print(f"  Salvo com sucesso em: {caminho_final_arquivo}\n")

    except Exception as e:
        print(f"  Erro ao processar o mosaico {nome_saida}: {e}\n")

print("Processamento de todos os mosaicos concluído!")
