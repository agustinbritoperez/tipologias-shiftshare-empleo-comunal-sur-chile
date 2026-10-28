# ============================================================
# v2_07_sensibilidad_k.R — ¿Dependen las conclusiones del número de grupos?
# ------------------------------------------------------------
# Usa las matrices de consenso ya calculadas en v2_02 (k = 2..8) y evalúa, para
# k = 3..7, las conclusiones que dependen de la tipología:
#   (1) persistencia del nivel urbano, (2) retroceso del grupo primario,
#   (3) movilidad entre tipos, (4) la tipología no explica el RS comunal (η², permutación).
# Sin etiquetas: el grupo urbano es el de mayor empleo mediano y el primario el de
# mayor participación primaria media.
# ============================================================
Sys.unsetenv("AJUSTE_SALTOS")
source("v2_00_datos.R")
suppressPackageStartupMessages(library(mclust))
OUT <- file.path(RES, "Sensibilidad_k")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
set.seed(2026)

tf  <- readRDS(file.path(RES, "Tipologia_funcional/tipologia_funcional.rds"))
ss  <- readRDS(file.path(RES, "ShiftShare/shiftshare_v2.rds"))
ind <- tf$ind %>% bind_cols(as_tibble(P[tf$ind$id, , drop = FALSE]))
cl5 <- cutree(hclust(as.dist(1 - tf$cons[["5"]]), "average"), 5)

eta2 <- function(y, g) { m <- tapply(y, g, mean); n <- tapply(y, g, length)
  sum(n * (m - mean(y))^2) / sum((y - mean(y))^2) }
inicio <- c("2006–2010" = 2006, "2010–2015" = 2010, "2015–2020" = 2015, "2020–2024" = 2020)

res <- bind_rows(lapply(3:7, function(k) {
  C  <- tf$cons[[as.character(k)]]
  cl <- cutree(hclust(as.dist(1 - C), "average"), k)
  d  <- ind %>% mutate(g = cl)
  perfil <- d %>% group_by(g) %>% summarise(emp = median(total), prim = mean(Primario), n = n(),
                                            cons = { m <- C[cl == first(g), cl == first(g)]; mean(m[upper.tri(m)]) })
  g_urb  <- perfil$g[which.max(perfil$emp)]
  g_prim <- perfil$g[which.max(perfil$prim)]
  cambios <- d %>% arrange(comuna, anio) %>% group_by(comuna) %>%
    summarise(c = sum(g[-1] != g[-n()]), .groups = "drop")
  # η² del RS comunal según el grupo al inicio de cada periodo
  eta <- lapply(names(inicio), function(per) {
    x <- ss$cap2 %>% filter(periodo == per) %>% select(comuna, RS_pct) %>%
      inner_join(d %>% filter(anio == inicio[[per]]) %>% select(comuna, g), by = "comuna")
    obs <- eta2(x$RS_pct, x$g); perm <- replicate(999, eta2(x$RS_pct, sample(x$g)))
    c(obs, (1 + sum(perm >= obs)) / 1000)
  })
  tibble(k = k, ARI_vs_k5 = adjustedRandIndex(cl, cl5),
         consenso_min = min(perfil$cons), tamano_min = min(perfil$n),
         urbano_2006 = paste(sort(d$comuna[d$g == g_urb & d$anio == 2006]), collapse = ", "),
         urbano_2024 = paste(sort(d$comuna[d$g == g_urb & d$anio == 2024]), collapse = ", "),
         consenso_urbano = perfil$cons[perfil$g == g_urb],
         primario_pct_medio = round(100 * perfil$prim[perfil$g == g_prim], 1),
         primario_n_2006 = sum(d$g == g_prim & d$anio == 2006),
         primario_n_2024 = sum(d$g == g_prim & d$anio == 2024),
         comunas_sin_cambio = sum(cambios$c == 0),
         eta2 = paste(sprintf("%.2f", sapply(eta, `[`, 1)), collapse = " / "),
         eta2_p = paste(sprintf("%.2f", sapply(eta, `[`, 2)), collapse = " / "),
         eta2_min_p = min(sapply(eta, `[`, 2)))
}))
print(as.data.frame(res %>% mutate(across(where(is.double), ~ round(.x, 2)))), right = FALSE)
writexl::write_xlsx(res, file.path(OUT, "sensibilidad_k.xlsx"))
