# ============================================================
# v2_08_cambio_estructural.R — Cambio estructural continuo (índice NAV)
# ------------------------------------------------------------
# Sustituye el cambio discreto de tipología (pocos eventos) por una medida continua:
# cuánto se desplaza la composición sectorial de cada comuna en cada periodo, con el
# índice de norma de valores absolutos NAV = ½ Σ|s1 - s0| (proporción del empleo
# reasignada entre rubros). Se prefiere a la distancia de Aitchison porque esta
# sobrerrepresenta los rubros diminutos (log de participaciones cercanas a cero).
#   d_total : NAV entre la composición inicial y la final
#   d_nac   : desplazamiento esperado si cada rubro creciera al ritmo nacional
#             (exposición sectorial: composición contrafactual)
#   d_loc   : distancia entre la composición contrafactual y la observada (dinámica local)
# Se usan los mismos pasos anuales que el shift-share (sin los quiebres de registro), de
# modo que la composición final es E0 + suma de ΔE de esos pasos.
# Pruebas: (a) peso relativo de d_nac y d_loc; (b) ¿qué tipos cambian más y por qué?;
# (c) ¿el RS/IS previo anticipa la magnitud del cambio posterior? (Spearman, permutación).
# ============================================================
source("v2_00_datos.R")
suppressPackageStartupMessages(library(purrr))
OUT <- file.path(RES, "Cambio_estructural")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
set.seed(2026)

ss <- readRDS(file.path(RES, "ShiftShare/shiftshare_v2.rds"))
pasos <- ss$pasos

nav <- function(a, b) 0.5 * sum(abs(a / sum(a) - b / sum(b)))
factor_nac <- pasos %>% distinct(periodo, anio, rubro8, g_j) %>%
  group_by(periodo, rubro8) %>% summarise(F_j = prod(1 + g_j), .groups = "drop")

celdas <- pasos %>% group_by(periodo, comuna, tipo, rubro8) %>%
  summarise(E0 = E[anio == min(anio)], E1 = pmax(E0 + sum(dE), 0), .groups = "drop") %>%
  left_join(factor_nac, by = c("periodo", "rubro8")) %>%
  mutate(E_cf = E0 * F_j)

dist <- celdas %>% group_by(periodo, comuna, tipo) %>%
  filter(sum(E0) >= MIN_TRAB, sum(E1) >= MIN_TRAB) %>%
  summarise(d_total = nav(E0, E1), d_nac = nav(E0, E_cf), d_loc = nav(E_cf, E1), .groups = "drop") %>%
  mutate(prop_local = d_loc / (d_nac + d_loc)) %>%
  left_join(ss$cap2 %>% select(periodo, comuna, RS_pct, IS_pct), by = c("periodo", "comuna"))

# (a) Peso de la exposición sectorial frente a la dinámica local
resumen_a <- dist %>% group_by(periodo) %>%
  summarise(n = n(), d_total_med = median(d_total), d_nac_med = median(d_nac),
            d_loc_med = median(d_loc), prop_local_med = median(prop_local), .groups = "drop")
# (b) Por tipología (al inicio del periodo)
resumen_b <- dist %>% group_by(tipo) %>%
  summarise(n = n(), d_total_med = median(d_total), d_nac_med = median(d_nac),
            d_loc_med = median(d_loc), prop_local_med = median(prop_local), .groups = "drop")
kw_tipo <- kruskal.test(d_total ~ factor(tipo), data = dist)

# (c) ¿El desempeño previo anticipa la magnitud del cambio estructural siguiente?
siguiente <- c("2006–2010" = "2010–2015", "2010–2015" = "2015–2020", "2015–2020" = "2020–2024")
rez <- dist %>% select(periodo, comuna, RS_prev = RS_pct, IS_prev = IS_pct) %>%
  mutate(periodo = siguiente[periodo]) %>% filter(!is.na(periodo)) %>%
  inner_join(dist %>% select(periodo, comuna, d_total, d_loc), by = c("periodo", "comuna"))
perm_spearman <- function(x, y, estrato, B = 9999) {
  obs  <- cor(x, y, method = "spearman")
  perm <- replicate(B, cor(x, ave(y, estrato, FUN = sample), method = "spearman"))
  c(rho = obs, p = (1 + sum(abs(perm) >= abs(obs))) / (1 + B))
}
prueba_c <- bind_rows(lapply(c("RS_prev", "IS_prev"), function(v) bind_rows(lapply(c("d_total", "d_loc"), function(y) {
  r <- perm_spearman(rez[[v]], rez[[y]], rez$periodo)
  tibble(predictor = v, resultado = y, n = nrow(rez), rho = r[["rho"]], p_permutacion = r[["p"]])
})))) %>% mutate(p_holm = p.adjust(p_permutacion, "holm"))

# Validación: ¿la distancia continua coincide con los cambios de tipo observados?
tip <- readRDS(file.path(RES, "Tipologia_funcional/tipologia_funcional.rds"))$ind
fin <- c("2006–2010" = 2010, "2010–2015" = 2015, "2015–2020" = 2020, "2020–2024" = 2024)
valid <- dist %>% mutate(anio_fin = fin[periodo]) %>%
  left_join(tip %>% select(comuna, anio_fin = anio, tipo_fin = tipo), by = c("comuna", "anio_fin")) %>%
  filter(!is.na(tipo_fin)) %>% mutate(cambia = as.character(tipo) != as.character(tipo_fin))
val_w <- wilcox.test(d_total ~ cambia, data = valid)

cat("(a) Índice NAV por periodo (medianas):\n"); print(as.data.frame(resumen_a %>% mutate(across(where(is.double), ~ round(.x, 3)))))
cat("\n(b) Por tipología al inicio:\n"); print(as.data.frame(resumen_b %>% mutate(across(where(is.double), ~ round(.x, 3)))))
cat("Kruskal-Wallis d_total ~ tipología: p =", signif(kw_tipo$p.value, 3), "\n")
cat("\n(c) RS/IS previos y cambio estructural siguiente (Spearman, permutación dentro de periodo):\n"); print(as.data.frame(prueba_c %>% mutate(across(where(is.double), ~ round(.x, 3)))))
cat("\nValidación: mediana d_total en comunas que cambian de tipo vs no:",
    round(median(valid$d_total[valid$cambia]), 2), "vs", round(median(valid$d_total[!valid$cambia]), 2),
    "| Wilcoxon p =", signif(val_w$p.value, 3), "\n")

writexl::write_xlsx(list(distancias = dist, por_periodo = resumen_a, por_tipologia = resumen_b,
                         rezago = prueba_c, kruskal = tibble(p = kw_tipo$p.value),
                         validacion = tibble(med_cambia = median(valid$d_total[valid$cambia]),
                                             med_no = median(valid$d_total[!valid$cambia]), p = val_w$p.value)),
                    file.path(OUT, "cambio_estructural.xlsx"))
