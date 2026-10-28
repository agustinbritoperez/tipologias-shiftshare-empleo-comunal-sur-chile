# ============================================================
# v2_09_nacional.R — ¿Es nacional el menor desempeño competitivo de las comunas grandes?
# ------------------------------------------------------------
# En el área de estudio solo hay 3 comunas de más de 100 000 habitantes, lo que impide
# una inferencia potente. Se replica el shift-share dinámico (mismos 8 rubros y mismos
# pasos neutralizados) para todas las comunas de Chile excepto la Región Metropolitana
# (sede de la mayoría de las casas matrices) y se estima, dentro de cada región
# (efectos fijos regionales), la relación entre el RS comunal y el tamaño poblacional.
# ============================================================
source("v2_00_datos.R")
suppressPackageStartupMessages({ library(sandwich); library(lmtest); library(readxl) })
OUT <- file.path(RES, "Nacional")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
set.seed(2026)
PERIODOS <- list(c(2006, 2010), c(2010, 2015), c(2015, 2020), c(2020, 2024))
PASO_NEUTRALIZADO <- c(2014, 2017)          # igual que v2_03
GRANDE <- 100000

norm <- function(x) gsub("[^a-z]", "", iconv(tolower(x), "UTF-8", "ASCII//TRANSLIT"))
EQUIV <- c(paihuano = "paiguano", lacalera = "calera", marchigue = "marchihue", trehuaco = "treguaco",
           aisen = "aysen", coihaique = "coyhaique", titil = "tiltil")
pob_nac <- read_excel("estimaciones-y-proyecciones-2002-2035-comunas.xlsx", sheet = "Est. y Proy. de Pob. Comunal") %>%
  select(comuna = `Nombre Comuna`, starts_with("Poblacion ")) %>%
  pivot_longer(-comuna, names_to = "anio", names_prefix = "Poblacion ", values_to = "pob") %>%
  mutate(anio = as.integer(anio), clave = norm(comuna)) %>%
  group_by(clave, anio) %>% summarise(pob = sum(pob), .groups = "drop")

emp <- base %>%
  filter(anio >= 2006, !grepl("^Sin informaci", rubro, ignore.case = TRUE),
         !grepl("^Sin informaci", region, ignore.case = TRUE), !grepl("^Sin informaci", comuna, ignore.case = TRUE)) %>%
  mutate(letra = substr(rubro, 1, 1)) %>% filter(letra %in% names(GRUPOS)) %>%
  mutate(rubro8 = GRUPOS[letra]) %>%
  group_by(region, comuna, anio, rubro8) %>% summarise(E = sum(trabajadores), .groups = "drop") %>%
  complete(nesting(region, comuna), anio, rubro8, fill = list(E = 0))

# Referencia: Chile completo
ref  <- emp %>% group_by(anio, rubro8) %>% summarise(E = sum(E), .groups = "drop") %>%
  arrange(rubro8, anio) %>% group_by(rubro8) %>% mutate(g_j = lead(E) / E - 1) %>% ungroup()
gnac <- ref %>% group_by(anio) %>% summarise(E = sum(E)) %>% arrange(anio) %>% mutate(g = lead(E) / E - 1)

loc <- emp %>% filter(!grepl("Metropolitana", region)) %>%
  arrange(region, comuna, rubro8, anio) %>% group_by(region, comuna, rubro8) %>%
  mutate(dE = lead(E) - E) %>% ungroup() %>%
  left_join(ref %>% select(anio, rubro8, g_j), by = c("anio", "rubro8")) %>%
  left_join(gnac %>% select(anio, g), by = "anio")

rs <- bind_rows(lapply(PERIODOS, function(p) {
  loc %>% filter(anio >= p[1], anio < p[2], !anio %in% PASO_NEUTRALIZADO) %>%
    mutate(RS = dE - E * g_j) %>%
    group_by(region, comuna) %>%
    summarise(RS = sum(RS), .groups = "drop") %>%
    left_join(emp %>% filter(anio == p[1]) %>% group_by(comuna) %>% summarise(E0 = sum(E)), by = "comuna") %>%
    mutate(periodo = paste0(p[1], "–", p[2]), anio0 = p[1])
})) %>%
  mutate(clave = recode(norm(comuna), !!!EQUIV), RS_pct = 100 * RS / E0,
         area_estudio = region %in% REGIONES) %>%
  left_join(pob_nac, by = c("clave", "anio0" = "anio")) %>%
  filter(E0 >= MIN_TRAB)
stopifnot(!anyNA(rs$pob))
rs <- rs %>% mutate(grande = as.numeric(pob >= GRANDE), log_pob = log(pob))

# Modelos por periodo, con efectos fijos regionales (comparación dentro de cada región)
ajuste <- function(d, f, ponderar = FALSE) {
  d$w_ <- if (ponderar) d$E0 else 1
  m  <- lm(f, data = d, weights = w_)
  ct <- coeftest(m, vcov. = vcovHC(m, type = "HC1"))
  ct
}
res <- bind_rows(lapply(unique(rs$periodo), function(per) {
  d <- filter(rs, periodo == per)
  a <- ajuste(d, RS_pct ~ log_pob + factor(region))
  b <- ajuste(d, RS_pct ~ log_pob + factor(region), ponderar = TRUE)
  c <- ajuste(d, RS_pct ~ grande + factor(region))
  e <- ajuste(d, RS_pct ~ grande + factor(region), ponderar = TRUE)
  # permutación de la condición "grande" dentro de cada región (robusta a pocos casos)
  obs  <- coef(lm(RS_pct ~ grande + factor(region), d))["grande"]
  perm <- replicate(4999, coef(lm(RS_pct ~ ave(grande, region, FUN = sample) + factor(region), d))[2])
  tibble(periodo = per, n_comunas = nrow(d), n_grandes = sum(d$grande),
         b_logpob = a["log_pob", 1], p_logpob = a["log_pob", 4],
         b_logpob_pond = b["log_pob", 1], p_logpob_pond = b["log_pob", 4],
         b_grande = c["grande", 1], p_grande = c["grande", 4],
         b_grande_pond = e["grande", 1], p_grande_pond = e["grande", 4],
         p_grande_permutacion = (1 + sum(abs(perm) >= abs(obs))) / 5000,
         rho_spearman = cor(d$RS_pct, d$pob, method = "spearman"))
})) %>% mutate(p_holm_logpob = p.adjust(p_logpob, "holm"), p_holm_grande = p.adjust(p_grande_permutacion, "holm"))

cat("Comunas por periodo (sin RM, >= 100 trabajadores):", paste(table(rs$periodo), collapse = " / "), "\n")
print(as.data.frame(res %>% mutate(across(where(is.double), ~ round(.x, 3)))))
writexl::write_xlsx(list(modelos = res, datos = rs %>% select(-clave)), file.path(OUT, "nacional_tamano.xlsx"))
