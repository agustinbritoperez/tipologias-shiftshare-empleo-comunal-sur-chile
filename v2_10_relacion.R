# ============================================================
# v2_10_relacion.R — Diversificación relacionada: densidad de relación entre rubros
# ------------------------------------------------------------
# Hidalgo et al. (2007); Neffke et al. (2011); Boschma (2017).
#  - Especialización: RCA_cj = (E_cj/E_c)/(E_j/E) >= 1, con todas las comunas del país
#    (>= 100 trabajadores) al inicio de cada periodo.
#  - Proximidad entre rubros j,k: phi_jk = sum_c M_cj M_ck / max(sum_c M_cj, sum_c M_ck).
#  - Densidad de relación de la comuna c con el rubro j: sum_k phi_jk M_ck / sum_k phi_jk (k != j).
# Prueba (modelo de probabilidad lineal): entrada = pasar de RCA < 1 a RCA >= 1 en el periodo,
# en función de la densidad inicial, con efectos fijos comuna-periodo y rubro-periodo;
# errores agrupados por comuna. Se estima para el área de estudio y, por potencia, para el país.
# Rubros: 15 (A-N por separado; O-S en un bloque por los quiebres de registro).
# ============================================================
source("v2_00_datos.R")
suppressPackageStartupMessages({ library(sandwich); library(lmtest) })
OUT <- file.path(RES, "Relacion")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CORTES <- c(2006, 2010, 2015, 2020, 2024)

emp15 <- base %>%
  filter(anio %in% CORTES, !grepl("^Sin informaci", rubro, ignore.case = TRUE),
         !grepl("^Sin informaci", comuna, ignore.case = TRUE)) %>%
  mutate(letra = substr(rubro, 1, 1)) %>% filter(letra %in% names(GRUPOS)) %>%
  mutate(r15 = ifelse(letra %in% c("O", "P", "Q", "R", "S"), "OS", letra), area = region %in% REGIONES) %>%
  group_by(anio, comuna, area, r15) %>% summarise(E = sum(trabajadores), .groups = "drop") %>%
  complete(nesting(anio, comuna, area), r15, fill = list(E = 0)) %>%
  group_by(anio, comuna) %>% filter(sum(E) >= MIN_TRAB) %>% ungroup()

rca <- emp15 %>% group_by(anio, comuna) %>% mutate(Ec = sum(E)) %>%
  group_by(anio, r15) %>% mutate(Ej = sum(E)) %>% ungroup() %>%
  group_by(anio) %>% mutate(Et = sum(E)) %>% ungroup() %>%
  mutate(RCA = (E / Ec) / (Ej / Et), M = as.numeric(RCA >= 1))

densidad <- function(a) {
  d <- rca %>% filter(anio == a)
  Mx <- xtabs(M ~ comuna + r15, d); Mx <- matrix(Mx, nrow(Mx), dimnames = dimnames(Mx))
  ub <- colSums(Mx); co <- crossprod(Mx)
  phi <- co / outer(ub, ub, pmax); diag(phi) <- 0
  dens <- (Mx %*% phi) / matrix(colSums(phi), nrow(Mx), ncol(phi), byrow = TRUE)
  list(phi = phi, dens = as_tibble(as.data.frame.table(dens, responseName = "densidad")) %>%
         mutate(comuna = as.character(comuna), r15 = as.character(r15), anio = a))
}
dens_all <- lapply(CORTES[-length(CORTES)], densidad)
dens <- bind_rows(lapply(dens_all, `[[`, "dens"))

panel <- bind_rows(lapply(seq_len(length(CORTES) - 1), function(i) {
  a0 <- CORTES[i]; a1 <- CORTES[i + 1]
  rca %>% filter(anio == a0) %>% select(comuna, area, r15, M0 = M, RCA0 = RCA) %>%
    inner_join(rca %>% filter(anio == a1) %>% select(comuna, r15, M1 = M), by = c("comuna", "r15")) %>%
    left_join(dens %>% filter(anio == a0) %>% select(comuna, r15, densidad), by = c("comuna", "r15")) %>%
    mutate(periodo = paste0(a0, "–", a1))
})) %>% mutate(densidad_z = (densidad - mean(densidad)) / sd(densidad))

estimar <- function(d, etiqueta) {
  d <- d %>% filter(M0 == 0) %>% mutate(entrada = M1, cp = paste(comuna, periodo), rp = paste(r15, periodo))
  m <- lm(entrada ~ densidad_z + factor(cp) + factor(rp), data = d)
  ct <- coeftest(m, vcov. = vcovCL(m, cluster = ~ comuna, type = "HC1"))["densidad_z", ]
  tibble(muestra = etiqueta, celdas_candidatas = nrow(d), entradas = sum(d$entrada),
         tasa_entrada = mean(d$entrada), efecto_1de = ct[1], ee = ct[2], p = ct[4],
         efecto_relativo = ct[1] / mean(d$entrada))
}
res <- bind_rows(estimar(panel %>% filter(area), "Área de estudio (42 comunas)"),
                 estimar(panel, "Chile (todas las comunas)"),
                 estimar(panel %>% filter(!area), "Chile sin el área de estudio"))

# Descriptivo para el área: ¿en qué rubros entran las comunas que dejan la especialización primaria?
entradas_area <- panel %>% filter(area, M0 == 0, M1 == 1) %>% count(r15, name = "entradas") %>% arrange(desc(entradas))
phi_2006 <- dens_all[[1]]$phi
prox_A <- sort(phi_2006["A", ], decreasing = TRUE)

cat("Entrada en nuevos rubros (RCA >= 1) según densidad de relación (1 d.e.):\n")
print(as.data.frame(res %>% mutate(across(where(is.double), ~ round(.x, 3)))))
cat("\nRubros más próximos al primario (A) en 2006:", paste(names(prox_A)[1:5], round(prox_A[1:5], 2), collapse = " | "), "\n")
cat("Entradas en el área por rubro:", paste(entradas_area$r15, entradas_area$entradas, collapse = " | "), "\n")
writexl::write_xlsx(list(modelos = res, panel = panel, entradas_area = entradas_area,
                         proximidad_2006 = as_tibble(phi_2006, rownames = "rubro")),
                    file.path(OUT, "relacion.xlsx"))
