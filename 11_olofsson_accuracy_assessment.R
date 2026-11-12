# Required packages:
# - sf: reads the validation points (step 10 output)
# - dplyr / tidyr: table wrangling, joins and reshaping
# - ggplot2: every chart produced by this script
# - stringr: string helpers used by dplyr pipelines below
# - caret: builds the confusion matrix (confusionMatrix)
# - terra: reads the classification raster(s) to compute class areas/weights
#
# Computes the Olofsson et al. (2014) stratified accuracy assessment (overall,
# user's and producer's accuracy, with standard errors and 95% confidence
# intervals) for one classification raster/year, plus confusion matrices,
# error composition and patch-size charts

# ==============================================================================
# 1. REQUIRED PACKAGES AND DIRECTORIES
# ==============================================================================
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(caret)
library(terra)


# ------------------------------------------------------------------------------
# PATHS — set these manually
# ------------------------------------------------------------------------------
# pasta_validacao: folder where every chart/table produced below is written
# caminho_gpkg: validation points with extracted class values (step 10 output)
# raster_dir / nome_raster: folder and file name (without .tif) of the
# classification raster being evaluated

pasta_validacao <- "path/to/your/validation_output_dir/" # <<< AJUSTAR A SAIDA
caminho_gpkg    <- "path/to/your/validation_points_extracted.gpkg" # <<< AJUSTAR OS PONTOS
raster_dir      <- "path/to/your/classification_dir/" # <<< AJUSTAR
nome_raster      <- "your_classification_raster"  # <<< AJUSTAR (sem .tif)

if (!dir.exists(pasta_validacao)) {
  dir.create(pasta_validacao, recursive = TRUE)
}


# ------------------------------------------------------------------------------
# YEAR IDENTIFIER (set manually, not taken from the file name)
# ------------------------------------------------------------------------------

Ano <- 2024   # <<< AJUSTAR: valor entre 2018 e 2024

if (!Ano %in% 2018:2024) {
  stop(sprintf("Ano = %s fora do intervalo esperado (2018 a 2024).", Ano))
}

# Color per year — dark tones chosen so the white bar labels stay readable.
# 2024 is kept fixed as requested; the other years use darker/more muted
# tones than a default ggplot palette would give

cores_anos <- c(
  "2018" = "#1B3A4B",  # azul petroleo escuro
  "2019" = "#6B1E3C",  # vinho escuro
  "2020" = "#2D4A1E",  # verde floresta escuro
  "2021" = "#3D2A5C",  # violeta escuro
  "2022" = "#0F4C4C",  # teal escuro
  "2023" = "#4A3B1E",  # bronze/oliva escuro
  "2024" = "#a64901"   # laranja queimado (fixo, nao alterado)
)

cor_ano    <- cores_anos[[as.character(Ano)]]
rotulo_ano <- sprintf("Ano %d", Ano)


# ==============================================================================
# 2. EXPERIMENT AND CLASS DEFINITION (MASK CLASSES EXCLUDED)
# ==============================================================================
# A single "experiment", identified by Ano rather than by an algorithm/method
# name — this script evaluates one classification at a time

lista_experimentos <- c(nome_raster)

# Classes 13 (FLO) and 14 (NFLO) are removed: the assessment treats these
# remaining 12 classes as 100% of the area of interest, same convention used
# when generating the sample points (step 9)

dicionario_classes <- data.frame(
  id_raster = c(1,2,3,4,5,6,7,8,9,10,11,12),
  Classe = c(
    "AGPE","AGUA","AL","DMC","DMF","PA","PH",
    "PSI","SE","URB","VSA","VSI"
  )
)

classes_validas <- dicionario_classes$Classe

class_colors <- c(
  "AGPE" = "#FF7F00",
  "AGUA" = "#0000CD",
  "AL"   = "#000080",
  "DMC"  = "#FF0000",
  "DMF"  = "#800000",
  "PA"   = "#F5DEB3",
  "PH"   = "#FFFF00",
  "PSI"  = "#00FFFF",
  "SE"   = "#800080",
  "URB"  = "#EE82EE",
  "VSA"  = "#228B22",
  "VSI"  = "#00FF00",
  "FLO"  = "#004600",
  "NFLO" = "#808080",
  "TOTAL"= "white"
)

# ---------------------------------------------------------------------------
# Experiment label — now based on Ano instead of the raster file name
# ---------------------------------------------------------------------------

mapa_nomes_experimentos <- c(rotulo_ano)
names(mapa_nomes_experimentos) <- nome_raster

rotular_algoritmo <- function(x) {
  unname(mapa_nomes_experimentos[x])
}

niveis_exp <- c(rotulo_ano)

cores_mapeamentos <- setNames(cor_ano, rotulo_ano)

# ==============================================================================
# 3. AUTOMATIC EXTRACTION OF THE WEIGHTS (Wi)
# ==============================================================================
# Wi is each class's share of the total mapped area, read directly from the
# classification raster's pixel counts — this is the "map area" weight
# required by the Olofsson stratified estimator (used in the function below)

pesos_por_experimento <- list()

cat("\nExtraindo áreas e calculando pesos (Wi)...\n")

for (exp_name in lista_experimentos) {

  caminho_raster <- file.path(raster_dir, paste0(exp_name, ".tif"))

  if (!file.exists(caminho_raster)) {
    stop(sprintf("Raster não encontrado: %s", caminho_raster))
  }

  r <- rast(caminho_raster)

  freq_df <- as.data.frame(freq(r))

  df_peso_exp <- freq_df %>%
    rename(id_raster = value, contagem = count) %>%
    left_join(dicionario_classes, by = "id_raster") %>%
    filter(!is.na(Classe)) %>%
    mutate(
      Area_ha = contagem * (res(r)[1] * res(r)[2]) / 10000,
      Peso = contagem / sum(contagem)
    ) %>%
    select(Classe, Peso, Area_ha)

  pesos_por_experimento[[exp_name]] <- df_peso_exp
}

cat("\nPesos dinâmicos calculados com sucesso.\n")

# ==============================================================================
# 4. READING THE VALIDATION DATA
# ==============================================================================

cat("\nCarregando pontos de validação...\n")

if (!file.exists(caminho_gpkg)) {
  stop(sprintf("Arquivo GPKG não encontrado: %s", caminho_gpkg))
}

dados_df <- st_read(caminho_gpkg, quiet = TRUE) %>%
  st_drop_geometry()

# ==============================================================================
# 5. STRATIFIED OLOFSSON FUNCTION (WITH VARIANCES - Eqs. 5, 6 AND 7)
# ==============================================================================
# Implements Olofsson et al. (2014): builds the area-weighted error matrix
# (p_ij) from the sample's Map x Reference cross-tabulation and the map-area
# weights (Wi), then derives overall/user's/producer's accuracy and their
# standard errors and 95% confidence intervals (SE * 1.96)

formatar_4_casas <- function(x) {
  round(x * 100, 4)
}

calcular_olofsson_estratificado <- function(df, col_pred, col_ref, df_pesos) {

  df[[col_pred]] <- factor(df[[col_pred]], levels = classes_validas)
  df[[col_ref]]  <- factor(df[[col_ref]],  levels = classes_validas)

  df_valido <- df %>%
    filter(!is.na(.data[[col_pred]]), !is.na(.data[[col_ref]]))

  # n_i: number of validation samples that fall in each mapped class,
  # needed to turn area weights into per-sample weights below

  n_i_df <- df_valido %>%
    group_by(Map = .data[[col_pred]]) %>%
    summarise(n_i = n(), .groups = "drop")

  # w_u: per-sample weight (Wi / n_i), the building block of the
  # area-weighted error matrix p_ij

  df_pesado <- df_valido %>%
    mutate(Map = .data[[col_pred]]) %>%
    left_join(n_i_df, by = "Map") %>%
    left_join(df_pesos %>% select(Classe, Peso), by = c("Map" = "Classe")) %>%
    mutate(
      w_u = ifelse(is.na(Peso) | n_i == 0, 0, Peso / n_i)
    )

  # p_ij: area-weighted proportion of the map in cell (Map=i, Reference=j).
  # complete() fills in combinations absent from the sample with 0, so the
  # matrix always has every class in both dimensions

  p_ij_df <- df_pesado %>%
    group_by(Map = .data[[col_pred]], Ref = .data[[col_ref]]) %>%
    summarise(p_ij = sum(w_u), .groups = "drop") %>%
    complete(Map = classes_validas, Ref = classes_validas, fill = list(p_ij = 0))

  p_ij_mat <- xtabs(p_ij ~ Map + Ref, data = p_ij_df)

  soma_total <- sum(p_ij_mat)
  cat("\n", col_pred, " | Soma matriz = ", round(soma_total, 8), "\n")

  # W_i: row totals of p_ij (equivalent to the map-area weights)
  # U_i: user's accuracy per class (diagonal / row total)
  # P_j: producer's accuracy per class (diagonal / column total)
  # OA: overall accuracy (sum of the diagonal)

  W_i <- rowSums(p_ij_mat)

  U_i <- diag(p_ij_mat) / W_i
  U_i[is.nan(U_i)]      <- 0
  U_i[is.infinite(U_i)] <- 0

  P_j_den <- colSums(p_ij_mat)

  P_j <- diag(p_ij_mat) / P_j_den
  P_j[is.nan(P_j)]      <- 0
  P_j[is.infinite(P_j)] <- 0

  OA <- sum(diag(p_ij_mat))

  # n_ij: raw sample counts (not area-weighted), needed by the variance
  # formulas below

  n_ij_df <- df_valido %>%
    group_by(Map = .data[[col_pred]], Ref = .data[[col_ref]]) %>%
    summarise(n_ij = n(), .groups = "drop") %>%
    complete(Map = classes_validas, Ref = classes_validas, fill = list(n_ij = 0))

  n_ij_mat <- xtabs(n_ij ~ Map + Ref, data = n_ij_df)

  n_i_vec <- setNames(rep(0, length(classes_validas)), classes_validas)
  n_i_vec[n_i_df$Map] <- n_i_df$n_i
  n_i_vec <- n_i_vec[classes_validas]

  # Variance of the overall accuracy (Olofsson Eq. 5): a weighted sum of
  # each class's binomial variance, undefined (set to 0) for classes with
  # 0 or 1 sample

  termo_O <- (W_i^2) * U_i * (1 - U_i) / pmax(n_i_vec - 1, 1)
  termo_O[n_i_vec <= 1] <- 0
  V_O  <- sum(termo_O, na.rm = TRUE)
  SE_O <- sqrt(V_O)

  # Variance of the user's accuracy per class (Olofsson Eq. 6): a simple
  # binomial variance of U_i within that class's own sample

  V_Ui <- U_i * (1 - U_i) / pmax(n_i_vec - 1, 1)
  V_Ui[n_i_vec <= 1] <- 0
  V_Ui[is.nan(V_Ui) | is.infinite(V_Ui)] <- 0
  SE_Ui <- sqrt(V_Ui)

  # Variance of the producer's accuracy per class (Olofsson Eq. 7): unlike
  # U_i and OA, this depends on every other class's confusion with class j,
  # hence the inner loop over "i != j" below

  SE_Pj <- setNames(rep(0, length(classes_validas)), classes_validas)

  for (j in classes_validas) {
    nj     <- n_i_vec[j]
    Wj     <- W_i[j]
    Uj     <- U_i[j]
    Pj     <- P_j[j]
    pj_hat <- P_j_den[j]

    if (is.na(pj_hat) || pj_hat == 0 || nj <= 1) next

    termo1 <- (Wj^2) * ((1 - Pj)^2) * Uj * (1 - Uj) / (nj - 1)

    termo2 <- 0
    for (i in classes_validas) {
      if (i == j) next
      ni <- n_i_vec[i]
      if (ni <= 1) next
      frac   <- n_ij_mat[i, j] / ni
      termo2 <- termo2 + (W_i[i]^2) * frac * (1 - frac) / (ni - 1)
    }

    V_Pj     <- (1 / pj_hat^2) * (termo1 + (Pj^2) * termo2)
    SE_Pj[j] <- sqrt(V_Pj)
  }

  # Assembles the one-row result: overall accuracy plus, for every class,
  # user's and producer's accuracy, each with SE and 95% CI (SE * 1.96)

  resultado <- data.frame(
    Experiment  = col_pred,
    Global      = formatar_4_casas(OA),
    Global_SE   = formatar_4_casas(SE_O),
    Global_IC95 = formatar_4_casas(1.96 * SE_O)
  )

  for (classe in classes_validas) {
    resultado[[paste0("Usuario_", classe)]]          <- formatar_4_casas(U_i[classe])
    resultado[[paste0("Usuario_", classe, "_SE")]]    <- formatar_4_casas(SE_Ui[classe])
    resultado[[paste0("Usuario_", classe, "_IC95")]]  <- formatar_4_casas(1.96 * SE_Ui[classe])

    resultado[[paste0("Produtor_", classe)]]          <- formatar_4_casas(P_j[classe])
    resultado[[paste0("Produtor_", classe, "_SE")]]    <- formatar_4_casas(SE_Pj[classe])
    resultado[[paste0("Produtor_", classe, "_IC95")]]  <- formatar_4_casas(1.96 * SE_Pj[classe])
  }

  return(resultado)
}

# ==============================================================================
# 6. EXECUTION
# ==============================================================================
# Runs the function above for every entry in lista_experimentos (a single
# raster/year in this script), comparing its predicted classes against the
# "Revisor" (reviewer) column in the validation points

resultados_list <- lapply(lista_experimentos, function(exp) {
  df_pesos_especifico <- pesos_por_experimento[[exp]]
  if (is.null(df_pesos_especifico)) {
    stop(sprintf("Pesos não encontrados para %s", exp))
  }
  calcular_olofsson_estratificado(
    df       = dados_df,
    col_pred = exp,
    col_ref  = "Revisor",
    df_pesos = df_pesos_especifico
  )
})

df_limpo <- do.call(rbind, resultados_list)

# ==============================================================================
# 7. NAME STANDARDIZATION
# ==============================================================================

df_limpo <- df_limpo %>%
  mutate(Experiment_Group = rotular_algoritmo(Experiment))

df_limpo$Experiment_Group <- factor(df_limpo$Experiment_Group, levels = niveis_exp)

# ==============================================================================
# 8. ACCURACY BAR CHARTS (AGPE CLASS) — WITH 95% CI
# ==============================================================================
# AGPE is this study's main class of interest, so it gets its own dedicated
# charts (global, producer's and user's accuracy) in addition to the
# all-classes charts further below

gerar_plot_barra <- function(dados_plot, coluna_y, coluna_ic, titulo_grafico, sufixo_arquivo) {

  dados_plot <- dados_plot %>%
    mutate(
      ymin_ic = pmax(0,   !!sym(coluna_y) - !!sym(coluna_ic)),
      ymax_ic = pmin(100, !!sym(coluna_y) + !!sym(coluna_ic)),
      y_label_dentro = pmax(!!sym(coluna_y) * 0.92, 4)
    )

  p <- ggplot(dados_plot, aes(x = Experiment_Group, y = !!sym(coluna_y), fill = Experiment_Group)) +
    geom_bar(stat = "identity", color = "black", width = 0.4) +
    geom_errorbar(aes(ymin = ymin_ic, ymax = ymax_ic), width = 0.1, linewidth = 0.6, color = "grey20") +
    geom_text(aes(y = ymax_ic, label = paste0("±", round(!!sym(coluna_ic), 2))),
              vjust = -0.8, size = 2.8, fontface = "italic", color = "grey30") +
    geom_text(aes(y = y_label_dentro, label = paste0(!!sym(coluna_y), "%")),
              vjust = 1, size = 3.8, fontface = "bold", color = "white") +
    scale_fill_manual(values = cores_mapeamentos) +
    ylim(0, 110) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 10, color = "black"),
      axis.text.y = element_text(face = "bold", size = 10, color = "black"),
      axis.title = element_text(face = "bold", size = 11),
      legend.position = "none",
      panel.grid.major.x = element_blank()
    ) +
    labs(
      title = titulo_grafico,
      y = "Acurácia (%)",
      x = "",
      caption = "Bigode = Intervalo de Confiança de 95% (±1.96 x EP), conforme Olofsson et al. (2014)"
    )

  ggsave(filename = paste0(pasta_validacao, sufixo_arquivo), plot = p, device = "png", width = 6, height = 6, dpi = 600)
  return(p)
}

print("Gerando gráficos de acurácia (AGPE)...")
p1 <- gerar_plot_barra(df_limpo, "Global", "Global_IC95", paste("Acurácia Global - Classificação Ano", Ano), "Comparativo_Acuracia_Global.png")
p2 <- gerar_plot_barra(df_limpo, "Produtor_AGPE", "Produtor_AGPE_IC95", paste("Acurácia do Produtor AGPE - Ano", Ano), "Comparativo_Acuracia_Produtor_AGPE.png")
p3 <- gerar_plot_barra(df_limpo, "Usuario_AGPE", "Usuario_AGPE_IC95", paste("Acurácia do Usuário AGPE - Ano", Ano), "Comparativo_Acuracia_Usuario_AGPE.png")

# Combined AGPE chart (producer's + user's accuracy side by side)

gerar_plot_duplo_agpe <- function(dados_agpe, sufixo_arquivo) {
  df_prod <- dados_agpe %>% transmute(Experiment_Group, Tipo_Acuracia = "Acurácia do Produtor", Valor_Acuracia = Produtor_AGPE, IC95 = Produtor_AGPE_IC95)
  df_user <- dados_agpe %>% transmute(Experiment_Group, Tipo_Acuracia = "Acurácia do Usuário", Valor_Acuracia = Usuario_AGPE, IC95 = Usuario_AGPE_IC95)

  df_long_agpe <- bind_rows(df_prod, df_user) %>%
    mutate(
      ymin_ic = pmax(0,   Valor_Acuracia - IC95),
      ymax_ic = pmin(100, Valor_Acuracia + IC95),
      y_label_dentro = pmax(Valor_Acuracia * 0.92, 4)
    )

  p_duplo <- ggplot(df_long_agpe, aes(x = Experiment_Group, y = Valor_Acuracia, fill = Experiment_Group)) +
    geom_bar(stat = "identity", color = "black", width = 0.4) +
    geom_errorbar(aes(ymin = ymin_ic, ymax = ymax_ic), width = 0.1, linewidth = 0.6, color = "grey20") +
    geom_text(aes(y = ymax_ic, label = paste0("±", round(IC95, 2))), vjust = -0.8, size = 2.8, fontface = "italic", color = "grey30") +
    geom_text(aes(y = y_label_dentro, label = paste0(Valor_Acuracia, "%")), vjust = 1, size = 3.8, fontface = "bold", color = "white") +
    scale_fill_manual(values = cores_mapeamentos) +
    ylim(0, 110) +
    facet_wrap(~ Tipo_Acuracia, ncol = 1, scales = "free_y") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 10, color = "black"),
      axis.text.y = element_text(face = "bold", size = 10, color = "black"),
      axis.title = element_text(face = "bold", size = 11),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey90", color = "grey70", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 11, color = "black"),
      panel.spacing = unit(2, "lines")
    ) +
    labs(title = paste("Desempenho da Classe AGPE - Ano", Ano), y = "Acurácia (%)", x = "", caption = "Bigode = Intervalo de Confiança de 95% (±1.96 x EP)")

  png(filename = paste0(pasta_validacao, sufixo_arquivo), width = 6, height = 8, units = "in", res = 600)
  print(p_duplo)
  dev.off()
  return(p_duplo)
}

p_agpe_combinado <- gerar_plot_duplo_agpe(df_limpo, "Comparativo_Acuracias_Combinadas_AGPE.png")

# ==============================================================================
# 8B. PRODUCER'S AND USER'S ACCURACY FOR EVERY CLASS
# ==============================================================================
print("Extraindo acurácia de Produtor e Usuário para todas as classes...")

extrair_valores_classe <- function(df_row, prefixo, sufixo = "") {
  sapply(classes_validas, function(cl) {
    col <- paste0(prefixo, cl, sufixo)
    if (col %in% names(df_row)) as.numeric(df_row[[col]][1]) else NA_real_
  })
}

df_todas_classes <- bind_rows(
  data.frame(
    Classe         = classes_validas,
    Valor_Acuracia = extrair_valores_classe(df_limpo, "Produtor_"),
    IC95           = extrair_valores_classe(df_limpo, "Produtor_", "_IC95"),
    Tipo_Acuracia  = "Acurácia do Produtor"
  ),
  data.frame(
    Classe         = classes_validas,
    Valor_Acuracia = extrair_valores_classe(df_limpo, "Usuario_"),
    IC95           = extrair_valores_classe(df_limpo, "Usuario_", "_IC95"),
    Tipo_Acuracia  = "Acurácia do Usuário"
  )
)

df_todas_classes$Classe <- factor(df_todas_classes$Classe, levels = classes_validas)

gerar_plot_barra_classes <- function(dados_plot, titulo_grafico, sufixo_arquivo) {

  # Classes with light fill colors get black bar labels instead of white,
  # so the accuracy value stays readable against the bar

  classes_claras <- c("PH", "PA", "PSI", "VSI")

  dados_plot <- dados_plot %>%
    mutate(
      ymin_ic = pmax(0,   Valor_Acuracia - IC95),
      ymax_ic = pmin(100, Valor_Acuracia + IC95),
      y_label_dentro = pmax(Valor_Acuracia * 0.92, 4),
      cor_label = ifelse(Classe %in% classes_claras, "black", "white")
    )

  p <- ggplot(dados_plot, aes(x = Classe, y = Valor_Acuracia, fill = Classe)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    geom_errorbar(aes(ymin = ymin_ic, ymax = ymax_ic), width = 0.2, linewidth = 0.5, color = "grey20") +
    geom_text(aes(y = ymax_ic, label = paste0("±", round(IC95, 1))),
              vjust = -0.6, size = 2.4, fontface = "italic", color = "grey30") +
    geom_text(aes(y = y_label_dentro, label = paste0(round(Valor_Acuracia, 1), "%"), color = cor_label),
              vjust = 1, size = 3, fontface = "bold") +
    scale_fill_manual(values = class_colors, name = "Classe") +
    scale_color_identity() +
    ylim(0, 115) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 9, color = "black"),
      axis.text.y = element_text(face = "bold", size = 10, color = "black"),
      axis.title = element_text(face = "bold", size = 11),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      title = titulo_grafico,
      y = "Acurácia (%)",
      x = "Classe",
      caption = "Bigode = Intervalo de Confiança de 95% (±1.96 x EP), conforme Olofsson et al. (2014)"
    )

  ggsave(filename = paste0(pasta_validacao, sufixo_arquivo), plot = p, device = "png", width = 11, height = 6, dpi = 600)
  return(p)
}

p_prod_todas <- gerar_plot_barra_classes(
  df_todas_classes %>% filter(Tipo_Acuracia == "Acurácia do Produtor"),
  paste("Acurácia do Produtor por Classe - Ano", Ano),
  "Comparativo_Acuracia_Produtor_TodasClasses.png"
)

p_user_todas <- gerar_plot_barra_classes(
  df_todas_classes %>% filter(Tipo_Acuracia == "Acurácia do Usuário"),
  paste("Acurácia do Usuário por Classe - Ano", Ano),
  "Comparativo_Acuracia_Usuario_TodasClasses.png"
)

# Combined chart (producer's + user's accuracy, every class)

gerar_plot_duplo_todas_classes <- function(dados_todas, sufixo_arquivo) {

  classes_claras <- c("PH", "PA", "PSI", "VSI")

  dados_todas <- dados_todas %>%
    mutate(
      ymin_ic = pmax(0,   Valor_Acuracia - IC95),
      ymax_ic = pmin(100, Valor_Acuracia + IC95),
      y_label_dentro = pmax(Valor_Acuracia * 0.92, 4),
      cor_label = ifelse(Classe %in% classes_claras, "black", "white")
    )

  p_duplo <- ggplot(dados_todas, aes(x = Classe, y = Valor_Acuracia, fill = Classe)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    geom_errorbar(aes(ymin = ymin_ic, ymax = ymax_ic), width = 0.2, linewidth = 0.5, color = "grey20") +
    geom_text(aes(y = ymax_ic, label = paste0("±", round(IC95, 1))), vjust = -0.6, size = 2.2, fontface = "italic", color = "grey30") +
    geom_text(aes(y = y_label_dentro, label = paste0(round(Valor_Acuracia, 1), "%"), color = cor_label), vjust = 1, size = 2.8, fontface = "bold") +
    scale_fill_manual(values = class_colors, name = "Classe") +
    scale_color_identity() +
    ylim(0, 115) +
    facet_wrap(~ Tipo_Acuracia, ncol = 1) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 9, color = "black"),
      axis.text.y = element_text(face = "bold", size = 10, color = "black"),
      axis.title = element_text(face = "bold", size = 11),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey90", color = "grey70", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 11, color = "black"),
      panel.spacing = unit(2, "lines")
    ) +
    labs(
      title = paste("Desempenho por Classe - Ano", Ano),
      y = "Acurácia (%)",
      x = "Classe",
      caption = "Bigode = Intervalo de Confiança de 95% (±1.96 x EP)"
    )

  png(filename = paste0(pasta_validacao, sufixo_arquivo), width = 11, height = 9, units = "in", res = 600)
  print(p_duplo)
  dev.off()
  return(p_duplo)
}

p_todas_classes_combinado <- gerar_plot_duplo_todas_classes(df_todas_classes, "Comparativo_Acuracias_Combinadas_TodasClasses.png")

print("Gráficos de acurácia por classe (Produtor e Usuário) gerados com sucesso.")

# ==============================================================================
# 9. EXPORT THE FULL RESULTS TABLE (accuracies + SE + 95% CI)
# ==============================================================================
write.csv2(
  df_limpo,
  file = paste0(pasta_validacao, "Tabela_Resultados_Olofsson_com_Variancia.csv"),
  row.names = FALSE
)

print("Processamento concluído: acurácias, erros padrão e IC95% calculados e exportados.")

# ==============================================================================
# 9B/9C. PAIRED STATISTICAL TEST (BOOTSTRAP) — REMOVED IN THIS VERSION
# ==============================================================================
# With a single Ano (equivalent to the previous version's "single
# experiment"), combn(lista_experimentos, 2) returns no pairs, so a paired
# comparison would not be statistically meaningful here. Left removed, as in
# the original reference script (individual Region Growing run)
# ==============================================================================

# ==============================================================================
# 10. CONFUSION MATRIX (MOSAIC - RAW / SAMPLE FREQUENCY)
# ==============================================================================
print("Gerando matriz de confusão...")
matrizes_longas_lista <- list()

for (exp_name in lista_experimentos) {
  pred_vec <- factor(dados_df[[exp_name]], levels = classes_validas)
  ref_vec  <- factor(dados_df$Revisor, levels = classes_validas)

  cm <- confusionMatrix(pred_vec, ref_vec)
  df_cm <- as.data.frame(cm$table)

  # Adds TOTAL row/column and grand total, so the matrix chart below can
  # render row/column sums alongside the per-cell counts

  total_ref    <- df_cm %>% group_by(Reference) %>% summarise(Freq = sum(Freq), .groups = 'drop') %>% mutate(Prediction = "TOTAL")
  total_pred   <- df_cm %>% group_by(Prediction) %>% summarise(Freq = sum(Freq), .groups = 'drop') %>% mutate(Reference = "TOTAL")
  grande_total <- data.frame(Prediction = "TOTAL", Reference = "TOTAL", Freq = sum(df_cm$Freq))

  df_cm$Prediction <- as.character(df_cm$Prediction)
  df_cm$Reference  <- as.character(df_cm$Reference)

  df_completa_exp <- bind_rows(df_cm, total_ref, total_pred, grande_total)
  exp_label <- rotular_algoritmo(exp_name)

  df_completa_exp <- df_completa_exp %>%
    mutate(
      Experiment  = exp_name,
      Exp_Label   = exp_label,
      Eh_Diagonal = if_else(Prediction == Reference & Prediction != "TOTAL", as.character(Prediction), NA_character_),
      Eh_Total    = if_else(Prediction == "TOTAL" | Reference == "TOTAL", TRUE, FALSE)
    )

  matrizes_longas_lista[[exp_name]] <- df_completa_exp
}

df_mosaico_completo <- bind_rows(matrizes_longas_lista)

classes_com_total_X <- c(classes_validas, "TOTAL")
classes_com_total_Y <- c(classes_validas, "TOTAL")

num_cols_mosaico <- length(lista_experimentos)

p_matrizes <- ggplot(df_mosaico_completo, aes(x = Reference, y = Prediction)) +
  geom_tile(fill = "white", color = "grey90", linewidth = 0.3) +
  geom_tile(data = filter(df_mosaico_completo, Eh_Total == TRUE), fill = "grey95", color = "grey80", linewidth = 0.4) +
  geom_tile(data = filter(df_mosaico_completo, !is.na(Eh_Diagonal)), aes(fill = Eh_Diagonal), color = "grey90", linewidth = 0.3) +
  geom_text(aes(
    label = Freq,
    color = case_when(
      !is.na(Eh_Diagonal) & Eh_Diagonal %in% c("AL", "DMF", "SE", "VSA") ~ "white",
      !is.na(Eh_Diagonal) ~ "black",
      Eh_Total == TRUE ~ "black",
      Freq == 0 ~ "grey75",
      TRUE ~ "black"
    ),
    fontface = if_else(Freq > 0 | Eh_Total == TRUE, "bold", "plain")
  ), size = 3.5) +
  scale_fill_manual(values = class_colors, guide = "none", na.value = "transparent") +
  scale_color_identity() +
  scale_y_discrete(limits = rev(classes_com_total_Y)) +
  scale_x_discrete(limits = classes_com_total_X, position = "top") +
  facet_wrap(~ Exp_Label, ncol = num_cols_mosaico) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 15)),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 9, color = "black"),
    axis.text.y = element_text(face = "bold", size = 9, color = "black"),
    axis.title = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "navy", color = "black"),
    strip.text = element_text(face = "bold", size = 11, color = "white"),
    panel.grid = element_blank()
  ) +
  labs(
    title = paste("Matriz de Confusão (Contagem Amostral Bruta) - Ano", Ano),
    x = "Referência (Revisor)",
    y = "Predição (Modelo)"
  )

ggsave(filename = paste0(pasta_validacao, "Mosaico_Matrizes_Confusao.png"), plot = p_matrizes, device = "png", width = 8, height = 7, dpi = 600)

print("Matriz de confusão gerada com sucesso.")

# ==============================================================================
# 11. CONFUSION MATRIX (INDIVIDUAL PLOT)
# ==============================================================================
# Same matrix as section 10, redrawn one experiment at a time (useful when
# lista_experimentos has more than one entry)

print("Gerando gráfico individual da matriz de confusão...")

for (exp_name in names(matrizes_longas_lista)) {

  df_exp    <- matrizes_longas_lista[[exp_name]]
  exp_label <- unique(df_exp$Exp_Label)

  p_individual <- ggplot(df_exp, aes(x = Reference, y = Prediction)) +
    geom_tile(fill = "white", color = "grey90", linewidth = 0.3) +
    geom_tile(data = filter(df_exp, Eh_Total == TRUE), fill = "grey95", color = "grey80", linewidth = 0.4) +
    geom_tile(data = filter(df_exp, !is.na(Eh_Diagonal)), aes(fill = Eh_Diagonal), color = "grey90", linewidth = 0.3) +
    geom_text(aes(
      label = Freq,
      color = case_when(
        !is.na(Eh_Diagonal) & Eh_Diagonal %in% c("AL", "AGUA", "DMF", "SE", "VSA") ~ "white",
        !is.na(Eh_Diagonal) ~ "black",
        Eh_Total == TRUE ~ "black",
        Freq == 0 ~ "grey75",
        TRUE ~ "black"
      ),
      fontface = if_else(Freq > 0 | Eh_Total == TRUE, "bold", "plain")
    ), size = 4) +
    scale_fill_manual(values = class_colors, guide = "none", na.value = "transparent") +
    scale_color_identity() +
    scale_y_discrete(limits = rev(classes_com_total_Y)) +
    scale_x_discrete(limits = classes_com_total_X, position = "top") +
    theme_bw() +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 20)),
      axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0, face = "bold", size = 10, color = "black"),
      axis.text.y = element_text(face = "bold", size = 10, color = "black"),
      axis.title.x = element_text(face = "bold", size = 12, margin = margin(b = 15)),
      axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 15)),
      panel.grid = element_blank(),
      plot.margin = margin(t = 20, r = 20, b = 10, l = 10)
    ) +
    labs(
      title = paste("Matriz de Confusão -", exp_label),
      x = "Referência (Revisor)",
      y = "Predição (Modelo)"
    )

  nome_arquivo <- paste0(pasta_validacao, "Matriz_Confusao_Individual_", exp_label, ".png")
  ggsave(filename = nome_arquivo, plot = p_individual, device = "png", width = 8, height = 7, dpi = 600)
}

print("Matriz individual salva com sucesso.")

# ==============================================================================
# 12. ERROR COMPOSITION (STACKED BARS) - OMISSION AND COMMISSION (AGPE)
# ==============================================================================
# Breaks down AGPE's omission error (producer's error: AGPE points the model
# mapped as something else) and commission error (user's error: other
# points the model mapped as AGPE) by which other class they were confused
# with

print("Calculando proporções e gerando gráfico de composição de erros...")

df_erros_lista <- list()

for (exp_name in lista_experimentos) {

  linha_metrics <- df_limpo %>% filter(Experiment == exp_name)
  if (nrow(linha_metrics) == 0) next

  erro_produtor_ajustado <- 100 - as.numeric(linha_metrics$Produtor_AGPE)
  erro_usuario_ajustado  <- 100 - as.numeric(linha_metrics$Usuario_AGPE)

  df_exp <- matrizes_longas_lista[[exp_name]]

  df_prod <- df_exp %>%
    filter(Reference == "AGPE", Prediction != "AGPE", Prediction != "TOTAL") %>%
    mutate(Tipo_Erro = "Producer Error", Classe_Confusa = as.character(Prediction), Total_Pts = sum(Freq), Erro_Ajustado = erro_produtor_ajustado, Experiment = exp_name)

  df_user <- df_exp %>%
    filter(Prediction == "AGPE", Reference != "AGPE", Reference != "TOTAL") %>%
    mutate(Tipo_Erro = "User Error", Classe_Confusa = as.character(Reference), Total_Pts = sum(Freq), Erro_Ajustado = erro_usuario_ajustado, Experiment = exp_name)

  df_erros_lista[[paste0(exp_name, "_PROD")]] <- df_prod
  df_erros_lista[[paste0(exp_name, "_USER")]] <- df_user
}

df_erros_combinado <- bind_rows(df_erros_lista)

# Keeps, per algorithm/error type, only the label text for the 3 classes
# that contribute the most error (Rank <= 3), to avoid cluttering the chart

df_erros_combinado <- df_erros_combinado %>%
  mutate(Algoritmo = rotular_algoritmo(Experiment)) %>%
  filter(Freq > 0) %>%
  group_by(Algoritmo, Tipo_Erro) %>%
  mutate(
    Perc       = (Freq / sum(Freq)) * 100,
    Rank       = rank(-Perc, ties.method = "first"),
    Label_Text = ifelse(Rank <= 3, sprintf("%s\n%.1f%%", Classe_Confusa, Perc), "")
  ) %>%
  ungroup()

df_erros_combinado$Tipo_Erro      <- factor(df_erros_combinado$Tipo_Erro, levels = c("Producer Error", "User Error"))
df_erros_combinado$Algoritmo      <- factor(df_erros_combinado$Algoritmo, levels = niveis_exp)
df_erros_combinado$Classe_Confusa <- factor(df_erros_combinado$Classe_Confusa, levels = classes_validas)

p_barras_erros <- ggplot(df_erros_combinado, aes(x = Algoritmo, y = Perc, fill = Classe_Confusa)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.4, width = 0.5) +
  geom_text(
    aes(label = Label_Text, color = ifelse(Classe_Confusa %in% c("PH", "PA"), "black", "white"), group = Classe_Confusa),
    position = position_stack(vjust = 0.5), size = 3.2, fontface = "bold", lineheight = 0.8
  ) +
  scale_fill_manual(values = class_colors) +
  scale_color_identity() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = function(x) paste0(x, "%")) +
  facet_wrap(~ Tipo_Erro, ncol = 2) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 15)),
    axis.text.x = element_text(face = "bold", size = 11, color = "black"),
    axis.text.y = element_text(face = "bold", size = 10, color = "black"),
    axis.title = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "grey85", color = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 12, color = "black", margin = margin(t = 8, b = 8)),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(15, 15, 15, 15)
  ) +
  labs(
    title = paste("Composição dos Erros de Omissão e Comissão - Classe AGPE - Ano", Ano),
    x = "",
    y = "Proporção do Erro (%)"
  )

arquivo_saida_barras <- paste0(pasta_validacao, "Composicao_Erros_AGPE_Barras.png")
ggsave(filename = arquivo_saida_barras, plot = p_barras_erros, device = "png", width = 9, height = 7, dpi = 600)

print(paste("Gráfico de composição de erros gerado com sucesso:", arquivo_saida_barras))

# ==============================================================================
# 13. TOTAL MAPPED AREA CHARTS (INCLUDING FLO AND NFLO)
# ==============================================================================
# Unlike sections 3-12 (which excluded FLO/NFLO from the accuracy
# assessment), this section reports mapped area for all 14 raster classes,
# since area itself is a valid summary even for the mask classes

print("Extraindo área de TODAS as 14 classes (incluindo máscaras) para plotagem...")

df_areas_lista <- list()

dicionario_todas_classes <- data.frame(
  id_raster = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
  Classe = c("AGPE", "AGUA", "AL", "DMC", "DMF", "PA", "PH", "PSI", "SE", "URB", "VSA", "VSI", "FLO", "NFLO")
)

for (exp_name in lista_experimentos) {
  caminho_raster <- file.path(raster_dir, paste0(exp_name, ".tif"))

  if (file.exists(caminho_raster)) {
    r <- rast(caminho_raster)
    freq_df <- freq(r)

    df_temp <- freq_df %>%
      rename(id_raster = value, contagem = count) %>%
      left_join(dicionario_todas_classes, by = "id_raster") %>%
      filter(!is.na(Classe)) %>%
      mutate(Area_ha = contagem * (res(r)[1] * res(r)[2]) / 10000, Experiment_Group = rotular_algoritmo(exp_name))

    df_areas_lista[[exp_name]] <- df_temp
  }
}

df_areas <- bind_rows(df_areas_lista)

df_areas$Experiment_Group <- factor(df_areas$Experiment_Group, levels = niveis_exp)
df_areas$Classe <- factor(df_areas$Classe, levels = dicionario_todas_classes$Classe)

formatar_numero_br <- function(x) {
  formatC(x, format = "f", big.mark = ".", decimal.mark = ",", digits = 0)
}

p_areas_tudo <- ggplot(df_areas, aes(x = Classe, y = Area_ha, fill = Experiment_Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +
  geom_text(aes(label = formatar_numero_br(Area_ha)), vjust = -0.2, hjust = 0, angle = 45, size = 2.5, fontface = "bold") +
  scale_fill_manual(values = cores_mapeamentos) +
  scale_y_continuous(labels = formatar_numero_br, expand = expansion(mult = c(0, 0.15))) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 15)),
    axis.text.x = element_text(face = "bold", size = 10, color = "black"),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  labs(title = paste("Área Total Mapeada por Classe - Ano", Ano), x = "Classe", y = "Área (Hectares)")

arquivo_saida_tudo <- paste0(pasta_validacao, "Comparativo_Areas_TUDO.png")
ggsave(filename = arquivo_saida_tudo, plot = p_areas_tudo, device = "png", width = 12, height = 7, dpi = 600)

# Same chart, zoomed in on the three classes most relevant to this study
# (PH, AGPE, PA)

df_areas_filtrado <- df_areas %>% filter(Classe %in% c("PH", "AGPE", "PA"))
df_areas_filtrado$Classe <- factor(df_areas_filtrado$Classe, levels = c("PH", "AGPE", "PA"))

p_areas_filtrado <- ggplot(df_areas_filtrado, aes(x = Classe, y = Area_ha, fill = Experiment_Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +
  geom_text(aes(label = formatar_numero_br(Area_ha)), vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = cores_mapeamentos) +
  scale_y_continuous(labels = formatar_numero_br, expand = expansion(mult = c(0, 0.15))) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 15)),
    axis.text.x = element_text(face = "bold", size = 12, color = "black"),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  labs(title = paste("Área Total Mapeada (Classes Selecionadas) - Ano", Ano), x = "Classe", y = "Área (Hectares)")

arquivo_saida_filtrado <- paste0(pasta_validacao, "Comparativo_Areas_Filtrado.png")
ggsave(filename = arquivo_saida_filtrado, plot = p_areas_filtrado, device = "png", width = 8, height = 7, dpi = 600)

print("Gráficos de área gerados e salvos com sucesso.")

# ==============================================================================
# 14. PATCH-SIZE ANALYSIS (PA, PH, AGPE)
# ==============================================================================
# For each target class, isolates it as a binary raster, delineates
# connected patches (terra::patches) and bins them by area, to show whether
# a class is dominated by large, contiguous patches or by many small,
# fragmented ones

print("Extraindo conectividade espacial e calculando tamanho dos fragmentos...")

terraOptions(memfrac = 0.8, tempdir = tempdir())

classes_alvo <- c("AGPE", "PA", "PH")
df_fragmentos_list <- list()

for (exp_name in lista_experimentos) {
  caminho_raster <- file.path(raster_dir, paste0(exp_name, ".tif"))

  if (file.exists(caminho_raster)) {
    r <- rast(caminho_raster)
    fator_ha <- (res(r)[1] * res(r)[2]) / 10000

    for (classe in classes_alvo) {
      print(sprintf("Processando fragmentos: %s | Classe %s", rotulo_ano, classe))

      id_classe <- dicionario_todas_classes$id_raster[dicionario_todas_classes$Classe == classe]
      matriz_reclass <- matrix(c(id_classe, 1), ncol = 2, byrow = TRUE)
      r_bin <- classify(r, matriz_reclass, others = NA)

      r_patches    <- patches(r_bin, directions = 8)
      freq_patches <- freq(r_patches)

      if (nrow(freq_patches) > 0) {
        freq_patches$Area_ha <- freq_patches$count * fator_ha

        resumo_tamanho <- freq_patches %>%
          mutate(Classe_Tamanho = cut(Area_ha, breaks = c(0, 2, 5, 10, 50, Inf), labels = c("0 a 2 ha", "2 a 5 ha", "5 a 10 ha", "10 a 50 ha", "> 50 ha"), include.lowest = TRUE)) %>%
          group_by(Classe_Tamanho) %>%
          summarise(Numero_Fragmentos = n(), .groups = "drop") %>%
          mutate(Experiment_Group = rotular_algoritmo(exp_name), Classe_Uso = classe)

        df_fragmentos_list[[paste(exp_name, classe)]] <- resumo_tamanho
      }

      rm(r_bin, r_patches, freq_patches)
      gc(verbose = FALSE)
    }

    rm(r)
    gc(verbose = FALSE)
  }
}

terraOptions(memfrac = 0.6)

df_fragmentos <- do.call(rbind, df_fragmentos_list)

df_fragmentos_pct <- df_fragmentos %>%
  group_by(Experiment_Group, Classe_Uso) %>%
  mutate(Total_Fragmentos = sum(Numero_Fragmentos), Pct_Fragmentos = (Numero_Fragmentos / Total_Fragmentos) * 100) %>%
  ungroup() %>%
  mutate(
    Classe_Tamanho   = factor(Classe_Tamanho, levels = rev(c("0 a 2 ha", "2 a 5 ha", "5 a 10 ha", "10 a 50 ha", "> 50 ha"))),
    Experiment_Group = factor(Experiment_Group, levels = niveis_exp)
  )

cores_tamanho <- c(
  "0 a 2 ha"    = "#FEE391",
  "2 a 5 ha"    = "#FEC44F",
  "5 a 10 ha"   = "#FE9929",
  "10 a 50 ha"  = "#D95F0E",
  "> 50 ha"     = "#993404"
)

p_tamanho <- ggplot(df_fragmentos_pct, aes(x = Experiment_Group, y = Pct_Fragmentos, fill = Classe_Tamanho)) +
  geom_bar(stat = "identity", color = "white", width = 0.5) +
  geom_text(aes(label = if_else(Pct_Fragmentos > 3, paste0(round(Pct_Fragmentos, 1), "%"), "")), position = position_stack(vjust = 0.5), size = 3.5, fontface = "bold", color = "black") +
  scale_fill_manual(values = cores_tamanho) +
  facet_wrap(~ Classe_Uso, ncol = 3) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 10)),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10, color = "black"),
    strip.background = element_rect(fill = "navy", color = "black"),
    strip.text = element_text(face = "bold", size = 12, color = "white"),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  ) +
  labs(
    title = paste("Distribuição do Número de Fragmentos por Classe de Tamanho - Ano", Ano),
    x = "",
    y = "Proporção de Fragmentos (%)",
    fill = "Tamanho do\nFragmento"
  )

caminho_salvamento_tamanho <- paste0(pasta_validacao, "Comparativo_Tamanho_Fragmentos_AGPE_PA_PH.png")
ggsave(filename = caminho_salvamento_tamanho, plot = p_tamanho, device = "png", width = 12, height = 8, dpi = 600)

print("Gráfico de distribuição de tamanhos gerado e salvo com sucesso.")

print(sprintf("Script completo executado com sucesso para %s.", rotulo_ano))
