# Required packages:
# - terra: reads the classified raster and draws the sample points
# - sf: builds/writes the output points as a spatial GeoPackage
# - dplyr: table wrangling for the per-class weights/allocation
#
# Generates stratified sampling points (Olofsson et al., 2014) for accuracy
# validation of the masked classification (step 8 output), EXCLUDING the FLO
# (13) and NFLO (14) mask classes. The remaining classes are treated as 100%
# of the area of interest for the weight/sample-size calculations below

library(terra)
library(sf)
library(dplyr)

# ============================================================
# CONFIG
# ============================================================

RASTER_PATH <- "path/to/your/masked_classification.tif"
OUTPUT_GPKG <- "path/to/your/sampling_points.gpkg"

NODATA_VAL   <- 255
CLASSES_EXCLUIR <- c(13, 14)  # FLO e NFLO

# Class dictionary with the expected user's accuracy (Ui) per class, used
# only to size the sample (Cochran/Olofsson equation below) — Ui here is a
# planning assumption, not a measured value

dicionario <- data.frame(
  id_raster = c(1,2,3,4,5,6,7,8,9,10,11,12),
  Classe    = c("AGPE","AGUA","AL","DMC","DMF","PA","PH","PSI","SE","URB","VSA","VSI"),
  Ui        = c(0.80, 0.90, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.90, 0.90)
)

SE_ALVO       <- 0.01   # target standard error for the overall accuracy estimate
MIN_RARA      <- 75     # minimum number of points enforced for rare classes
LIMIAR_RARA   <- 0.05   # a class is "rare" when its area weight Wi < 5%
EXCEDENTE_PCT <- 0.50   # extra points drawn on top of the statistical minimum (50%)
BORDA_PX      <- 5      # erosion buffer, in pixels, to avoid points near nodata/class borders

SEED <- 42
# ============================================================

set.seed(SEED)

dir.create(dirname(OUTPUT_GPKG), recursive = TRUE, showWarnings = FALSE)

cat("Lendo raster...\n")
r <- rast(RASTER_PATH)
r[r == NODATA_VAL] <- NA
r[r %in% CLASSES_EXCLUIR] <- NA   # exclui FLO e NFLO do calculo e do sorteio

# ------------------------------------------------------------
# Wi (area weight) recomputed only over the remaining classes (= 100%)
# ------------------------------------------------------------

freq_df <- as.data.frame(freq(r))
freq_df <- freq_df[!is.na(freq_df$value), ]

total_px <- sum(freq_df$count)

df_wi <- freq_df %>%
  rename(id_raster = value, contagem = count) %>%
  left_join(dicionario, by = "id_raster") %>%
  filter(!is.na(Classe)) %>%
  mutate(Wi = contagem / total_px)

cat("\nPesos (Wi) recalculados apos exclusao de FLO/NFLO:\n")
print(df_wi %>% select(Classe, contagem, Wi, Ui))

# ------------------------------------------------------------
# Total sample size (Cochran Eq. 5.25 / Olofsson Eq. 13)
# ------------------------------------------------------------

df_wi <- df_wi %>% mutate(Si = sqrt(Ui * (1 - Ui)))
soma_WiSi <- sum(df_wi$Wi * df_wi$Si)
n_total <- ceiling((soma_WiSi / SE_ALVO)^2)

cat(sprintf("\nn_total recomendado (Eq. 13, SE_alvo=%.3f): %d\n", SE_ALVO, n_total))

# ------------------------------------------------------------
# Allocation: MIN_RARA floor for rare classes, remaining points split
# proportionally to Wi among the common classes
# ------------------------------------------------------------

df_wi <- df_wi %>%
  mutate(tipo = ifelse(Wi < LIMIAR_RARA, "rara", "comum"))

n_raras_total <- MIN_RARA * sum(df_wi$tipo == "rara")
n_restante    <- max(n_total - n_raras_total, 0)
soma_Wi_comuns <- sum(df_wi$Wi[df_wi$tipo == "comum"])

df_wi <- df_wi %>%
  mutate(
    n_min = ifelse(tipo == "rara", MIN_RARA,
                   round(n_restante * Wi / soma_Wi_comuns)),
    # Global floor: no class, not even a "common" one, ends up with fewer
    # points than the floor set for rare classes (MIN_RARA)
    n_min = pmax(n_min, MIN_RARA),
    n_max = ceiling(n_min * (1 + EXCEDENTE_PCT))
  )

cat("\n=== Alocacao final (n minimo = alvo estatistico Olofsson | n maximo = com excedente de 50%) ===\n")
print(df_wi %>% select(Classe, Wi, Ui, tipo, n_min, n_max))
cat(sprintf("\nSoma n_min: %d | Soma n_max: %d\n", sum(df_wi$n_min), sum(df_wi$n_max)))

# ------------------------------------------------------------
# Draws n_max points per class (the allocation including the surplus),
# sampled from inside that class's eroded mask
# ------------------------------------------------------------
# ------------------------------------------------------------
# Mask erosion: a pixel is only eligible if its entire (2*BORDA_PX+1)
# neighborhood belongs to the same class. This automatically excludes
# pixels near nodata (incomplete neighborhood = NA) AND pixels near any
# other class (including FLO/NFLO, already set to NA in the raster) —
# it keeps sample points away from class boundaries, where classification
# error and mixed pixels are most likely
# ------------------------------------------------------------

erodir_mascara <- function(r, id_classe, buffer_px) {
  m <- as.numeric(r == id_classe)  # 1 / 0 / NA
  janela <- matrix(1, nrow = 2 * buffer_px + 1, ncol = 2 * buffer_px + 1)
  eroded <- focal(m, w = janela, fun = min, na.rm = FALSE)
  eroded[eroded == 0] <- NA
  eroded
}

sortear_pontos_classe <- function(r, id_classe, n, nome_classe, buffer_px) {
  m <- erodir_mascara(r, id_classe, buffer_px)
  celulas_validas <- which(!is.na(values(m, mat = FALSE)))

  cat(sprintf("  %s: %d celulas disponiveis apos erosao de %d px\n",
              nome_classe, length(celulas_validas), buffer_px))

  if (length(celulas_validas) < n) {
    warning(sprintf("%s: apenas %d celulas disponiveis, mas %d solicitados.",
                    nome_classe, length(celulas_validas), n))
    n <- length(celulas_validas)
  }

  celulas <- sample(celulas_validas, size = n, replace = FALSE)
  xy <- xyFromCell(m, celulas)
  pts <- st_as_sf(as.data.frame(xy), coords = c("x", "y"), crs = crs(r))
  pts$Classe <- nome_classe
  pts
}

cat("\nSorteando pontos (n_max, com excedente) por classe, com erosao de borda...\n")
lista_pontos <- lapply(seq_len(nrow(df_wi)), function(i) {
  sortear_pontos_classe(r, df_wi$id_raster[i], df_wi$n_max[i], df_wi$Classe[i], BORDA_PX)
})

pontos_finais <- do.call(rbind, lista_pontos)

cat(sprintf("\nTotal de pontos gerados: %d\n", nrow(pontos_finais)))

# ------------------------------------------------------------
# Compares the points actually generated per class against what was
# planned (n_min = statistical target, n_max = target plus surplus)
# ------------------------------------------------------------

contagem_final <- as.data.frame(table(pontos_finais$Classe))
names(contagem_final) <- c("Classe", "n_gerado")

resumo_contagem <- df_wi %>%
  select(Classe, n_min, n_max) %>%
  left_join(contagem_final, by = "Classe")

cat("\n=== Contagem final de pontos gerados por classe ===\n")
print(resumo_contagem)

st_write(pontos_finais, OUTPUT_GPKG, delete_dsn = TRUE)
cat(sprintf("Salvo em: %s\n", OUTPUT_GPKG))
