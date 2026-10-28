# ============================================================
# v2_06_sensibilidad_saltos.R — Versión principal vs. serie ajustada por saltos de registro
# ------------------------------------------------------------
# Requiere correr v2_02-v2_04 dos veces: sin variable (Resultados_v2/) y con
# AJUSTE_SALTOS=1 (Resultados_v2_ajustado/). Ver regla de detección en v2_00_datos.R.
# ============================================================
Sys.unsetenv("AJUSTE_SALTOS")
source("v2_00_datos.R")
suppressPackageStartupMessages(library(mclust))
OUT <- "Resultados_v2/Sensibilidad_saltos"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
V <- c(principal = "Resultados_v2", ajustada = "Resultados_v2_ajustado")
leer <- function(v, f) readRDS(file.path(V[[v]], f))
URBANO <- c("Centro regional", "Centro subregional")

tp <- leer("principal", "Tipologia_funcional/tipologia_funcional.rds")$ind
ta <- leer("ajustada", "Tipologia_funcional/tipologia_funcional.rds")$ind
cruce <- inner_join(tp %>% select(id, comuna, anio, tipo_principal = tipo, cl_p = cluster),
                    ta %>% select(id, tipo_ajustada = tipo, cl_a = cluster), by = "id")
ari <- adjustedRandIndex(cruce$cl_p, cruce$cl_a)

nivel_urbano <- cruce %>% group_by(anio) %>% summarise(
  principal = paste(sort(comuna[tipo_principal %in% URBANO]), collapse = ", "),
  ajustada  = paste(sort(comuna[tipo_ajustada %in% URBANO]), collapse = ", "))

area <- lapply(names(V), function(v) {
  s <- leer(v, "ShiftShare/shiftshare_v2.rds")$pasos
  s %>% group_by(periodo) %>% summarise(across(c(dE, NS, IS, RS), sum)) %>%
    left_join(s %>% filter(anio == as.integer(substr(periodo, 1, 4))) %>% group_by(periodo) %>%
                summarise(E0 = sum(E)), by = "periodo") %>%
    transmute(version = v, periodo, CO = 100 * dE / E0, NS = 100 * NS / E0, IS = 100 * IS / E0, RS = 100 * RS / E0)
}) %>% bind_rows()

hoja <- function(h) bind_rows(lapply(names(V), function(v)
  readxl::read_excel(file.path(V[[v]], "Capa3/capa3.xlsx"), h) %>% mutate(version = v, .before = 1)))
cr_vs_ref <- hoja("M2_diferencias") %>% filter(tipo == "Centro regional")
wald      <- hoja("M2_wald")
logit     <- hoja("b_logit")
tamano    <- hoja("c_robustez_tamano")
ri_cr     <- hoja("d1_CR_aleatorizacion")
ri_conj   <- hoja("d2_conjunta_aleatorizacion")
ri_trans  <- hoja("d3_transiciones_permutacion")
eta <- bind_rows(lapply(names(V), function(v) leer(v, "ShiftShare/shiftshare_v2.rds")$cap2_var %>% mutate(version = v, .before = 1)))

cat("Saltos marcados:", nrow(saltos), "| ARI entre tipologías:", round(ari, 2), "\n\nNivel urbano:\n"); print(as.data.frame(nivel_urbano))
cat("\nShift-share del área (%):\n"); print(as.data.frame(area %>% mutate(across(where(is.double), ~ round(.x, 1)))))
cat("\nCentro regional vs intermedia — aleatorización:\n"); print(as.data.frame(ri_cr %>% mutate(across(where(is.double), ~ round(.x, 3)))))
cat("\nCentro regional vs intermedia (M2):\n"); print(as.data.frame(cr_vs_ref %>% select(version, periodo, estimado, ic_inf, ic_sup, p) %>% mutate(across(where(is.double), ~ round(.x, 3)))))

writexl::write_xlsx(list(saltos_marcados = saltos, cruce_tipologias = cruce, nivel_urbano = nivel_urbano,
                         shiftshare_area = area, centro_regional_M2 = cr_vs_ref, wald_M2 = wald,
                         eta2 = eta, logit = logit, tamano = tamano,
                         CR_aleatorizacion = ri_cr, conjunta_aleatorizacion = ri_conj, transiciones_permutacion = ri_trans,
                         resumen = tibble(indicador = "ARI tipologías (principal vs ajustada)", valor = ari)),
                    file.path(OUT, "sensibilidad_saltos.xlsx"))
