# ============================================================
# v2_12_sensibilidad_umbrales.R — Sensibilidad al umbral de la regla de saltos
# ------------------------------------------------------------
# Compara las conclusiones centrales entre la serie sin ajustar y las series ajustadas
# con umbrales de 200, 300 (principal) y 500 trabajadores (v2_00_datos.R, SALTO_MIN).
# Requiere correr v2_02, v2_03, v2_04, v2_08 y v2_10 con AJUSTE_SALTOS=1 y SALTO_MIN=200/500.
# ============================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readxl); library(mclust) })
V <- c("Sin ajustar" = "Resultados_v2", "Umbral 200" = "Resultados_v2_ajustado_200",
       "Umbral 300" = "Resultados_v2_ajustado", "Umbral 500" = "Resultados_v2_ajustado_500")
CFG <- list("Sin ajustar" = c(A = "0", U = "300"), "Umbral 200" = c(A = "1", U = "200"),
            "Umbral 300" = c(A = "1", U = "300"), "Umbral 500" = c(A = "1", U = "500"))
URB <- c("Centro regional", "Centro subregional")
tip_p <- readRDS(file.path(V[[1]], "Tipologia_funcional/tipologia_funcional.rds"))$ind

fila <- function(v) {
  dir <- V[[v]]
  # datos de la versión (saltos marcados y participación primaria)
  e <- new.env()
  Sys.setenv(AJUSTE_SALTOS = CFG[[v]][["A"]], SALTO_MIN = CFG[[v]][["U"]])
  suppressMessages(suppressWarnings(source("v2_00_datos.R", local = e)))
  prim <- e$ancha %>% mutate(p = 100 * Primario / total) %>% group_by(anio) %>% summarise(n40 = sum(p > 40))
  tip <- readRDS(file.path(dir, "Tipologia_funcional/tipologia_funcional.rds"))$ind
  pas <- readRDS(file.path(dir, "ShiftShare/shiftshare_v2.rds"))$pasos
  area <- pas %>% group_by(periodo) %>% summarise(IS = sum(IS), RS = sum(RS)) %>%
    left_join(pas %>% filter(anio == as.integer(substr(periodo, 1, 4))) %>% group_by(periodo) %>% summarise(E0 = sum(E)), by = "periodo") %>%
    mutate(IS = 100 * IS / E0, RS = 100 * RS / E0)
  d1 <- read_excel(file.path(dir, "Capa3/capa3.xlsx"), "d1_CR_aleatorizacion")
  d2 <- read_excel(file.path(dir, "Capa3/capa3.xlsx"), "d2_conjunta_aleatorizacion")
  d4 <- read_excel(file.path(dir, "Capa3/capa3.xlsx"), "d4_eta2_holm")
  rel <- read_excel(file.path(dir, "Relacion/relacion.xlsx"), "modelos")
  nav <- read_excel(file.path(dir, "Cambio_estructural/cambio_estructural.xlsx"), "por_periodo")
  comun <- inner_join(tip_p %>% select(id, a = cluster), tip %>% select(id, b = cluster), by = "id")
  tibble(
    version = v,
    saltos_marcados = if (CFG[[v]][["A"]] == "1") nrow(e$saltos) else 0L,
    ARI_vs_sin_ajustar = adjustedRandIndex(comun$a, comun$b),
    urbano_2010 = sum(tip$anio == 2010 & tip$tipo %in% URB),
    urbano_2024 = sum(tip$anio == 2024 & tip$tipo %in% URB),
    urbano_2024_comunas = paste(sort(tip$comuna[tip$anio == 2024 & tip$tipo %in% URB]), collapse = ", "),
    primario_40_2006 = prim$n40[prim$anio == 2006], primario_40_2024 = prim$n40[prim$anio == 2024],
    tipo_primario_2006 = sum(tip$anio == 2006 & tip$tipo == "Primario especializado"),
    tipo_primario_2024 = sum(tip$anio == 2024 & tip$tipo == "Primario especializado"),
    IS_area = paste(sprintf("%.1f", area$IS), collapse = " / "),
    RS_area = paste(sprintf("%.1f", area$RS), collapse = " / "),
    IS_negativo_todos = all(area$IS < 0),
    p_conjunta_min = min(d2$p_aleatorizacion),
    eta2_p_holm_min = min(d4$p_holm),
    CR_vs_intermedia = paste(sprintf("%.1f", d1$estimado), collapse = " / "),
    CR_p_colas_iguales = paste(sprintf("%.3f", d1$p_colas_iguales), collapse = " / "),
    relacion_efecto_area = rel$efecto_1de[1], relacion_p_area = rel$p[1],
    nav_prop_local = paste(sprintf("%.2f", range(nav$prop_local_med)), collapse = "-")
  )
}
res <- bind_rows(lapply(names(V), fila))
Sys.unsetenv(c("AJUSTE_SALTOS", "SALTO_MIN"))
print(as.data.frame(t(res %>% mutate(across(where(is.double), ~ round(.x, 3))))), right = FALSE)
dir.create("Resultados_v2/Sensibilidad_saltos", showWarnings = FALSE, recursive = TRUE)
writexl::write_xlsx(res, "Resultados_v2/Sensibilidad_saltos/sensibilidad_umbrales.xlsx")
