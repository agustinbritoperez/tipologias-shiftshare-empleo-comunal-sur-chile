# ============================================================
# v2_03_shiftshare.R — Shift-share dinámico por tipología (capas 1 y 2)
# ------------------------------------------------------------
# Capa 1: shift-share DINÁMICO (Barff & Knight, 1988): la descomposición se
#   calcula año a año y se suma dentro de cada periodo, evitando depender de
#   los años extremos (p.ej. 2020) y de la estructura sectorial inicial.
#   Se agrega por tipología FIJA al inicio de cada periodo (v2_02), para que
#   el crecimiento de un tipo no incluya comunas que entran o salen de él.
#   NS = E_ij,t * g_t ; IS = E_ij,t * (g_j,t - g_t) ; RS = ΔE_ij,t - NS - IS
#   (RS como residuo: la identidad CO = NS + IS + RS se cumple exactamente y
#   el empleo en rubros que parten en cero no se pierde).
#   Referencia: Chile SIN el área de estudio ni registros sin comuna.
# Capa 2: RS anidado. RS_ij = RS_tipo_ij + RS_comuna_ij, con
#   RS_tipo_ij = E_ij,t * (g_Tj,t - g_j,t) (ventaja compartida por el tipo T) y
#   RS_comuna_ij = E_ij,t * (g_ij,t - g_Tj,t) (desempeño de la comuna frente a sus pares).
#   Se reporta qué parte de la varianza del RS comunal explica la tipología
#   (η²) y su significancia por permutación.
# ============================================================
source("v2_00_datos.R")

OUT      <- file.path(RES, "ShiftShare")
PERIODOS <- list(c(2006, 2010), c(2010, 2015), c(2015, 2020), c(2020, 2024))
N_PERM   <- 9999
# Pasos anuales neutralizados por discontinuidades nacionales de registro:
#   2014->2015: A cae ~14% y N sube ~47% a escala nacional (reclasificación/externalización);
#   2017->2018: O y P suben ~32% a escala nacional (cambio administrativo del registro).
# Los periodos 2010-2015 y 2015-2020 suman los otros 4 pasos; las versiones completas
# se reportan como robustez.
PASO_NEUTRALIZADO <- c(2014, 2017)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

TIPOS  <- c("Centro regional", "Centro subregional", "Intermedia agro-comercial",
            "Primario especializado", "Rural dependiente de servicios públicos",
            "Periferia cordillerana-turística")
PALETA <- setNames(c("#7a0177", "#dd3497", "#fe9929", "#4d9221", "#4575b4", "#74add1"), TIPOS)

# ---------------- Datos anuales ----------------
anual <- base %>%
  filter(anio >= 2006, !grepl("^Sin informaci", rubro, ignore.case = TRUE)) %>%
  mutate(letra = substr(rubro, 1, 1)) %>%
  filter(letra %in% names(GRUPOS)) %>%
  mutate(rubro8 = GRUPOS[letra], area = region %in% REGIONES)

# Referencia nacional: resto de Chile (sin área de estudio ni "Sin Información")
ref <- anual %>%
  filter(!area, !grepl("^Sin informaci", region, ignore.case = TRUE)) %>%
  group_by(anio, rubro8) %>% summarise(E = sum(trabajadores), .groups = "drop") %>%
  arrange(rubro8, anio) %>% group_by(rubro8) %>%
  mutate(g_j = lead(E) / E - 1) %>% ungroup()
g_nac <- ref %>% group_by(anio) %>% summarise(E = sum(E), .groups = "drop") %>%
  arrange(anio) %>% mutate(g = lead(E) / E - 1) %>% select(anio, g)

loc <- anual %>% filter(area) %>%
  group_by(comuna, anio, rubro8) %>% summarise(E = sum(trabajadores), .groups = "drop") %>%
  complete(comuna, anio, rubro8, fill = list(E = 0)) %>%
  arrange(comuna, rubro8, anio) %>% group_by(comuna, rubro8) %>%
  mutate(dE = lead(E) - E) %>% ungroup()

# ---------------- Tipología al inicio de cada periodo ----------------
tip <- readRDS(file.path(RES, "Tipologia_funcional/tipologia_funcional.rds"))$ind %>%
  select(comuna, anio, tipo) %>% mutate(tipo = as.character(tipo))
tipo_inicio <- bind_rows(lapply(PERIODOS, function(p) {
  # comuna-año excluidas en v2_02 (<100 trabajadores): se usa el siguiente corte disponible
  tip %>% filter(anio >= p[1]) %>% arrange(comuna, anio) %>% group_by(comuna) %>%
    slice(1) %>% ungroup() %>%
    transmute(comuna, periodo = paste0(p[1], "–", p[2]), tipo, tipo_de = anio)
}))
imputadas <- tipo_inicio %>% filter(tipo_de != as.integer(substr(periodo, 1, 4)))
cat("Tipologías tomadas de un corte posterior (comuna-año excluidas):\n"); print(imputadas)

# ---------------- Descomposición anual ----------------
pasos <- bind_rows(lapply(PERIODOS, function(p) {
  loc %>% filter(anio >= p[1], anio < p[2]) %>% mutate(periodo = paste0(p[1], "–", p[2]))
})) %>%
  left_join(tipo_inicio %>% select(comuna, periodo, tipo), by = c("comuna", "periodo")) %>%
  left_join(ref %>% select(anio, rubro8, g_j), by = c("anio", "rubro8")) %>%
  left_join(g_nac, by = "anio")
stopifnot(!anyNA(pasos$tipo), !anyNA(pasos$g_j), !anyNA(pasos$dE))

# Crecimiento sectorial de cada tipología (miembros fijos en el periodo)
g_tipo <- pasos %>% group_by(periodo, tipo, anio, rubro8) %>%
  summarise(E_T = sum(E), dE_T = sum(dE), .groups = "drop") %>%
  mutate(g_Tj = ifelse(E_T > 0, dE_T / E_T, NA_real_))

pasos <- pasos %>%
  left_join(g_tipo %>% select(periodo, tipo, anio, rubro8, g_Tj), by = c("periodo", "tipo", "anio", "rubro8")) %>%
  mutate(g_Tj = coalesce(g_Tj, g_j),        # tipo sin empleo en el rubro: sin efecto tipo
         NS = E * g, IS = E * (g_j - g), RS = dE - NS - IS,
         RS_tipo = E * (g_Tj - g_j), RS_comuna = RS - RS_tipo)

pasos_completo <- pasos
pasos <- pasos_completo %>% filter(!anio %in% PASO_NEUTRALIZADO)

# Controles de consistencia
stopifnot(isTRUE(all.equal(sum(pasos$dE), sum(pasos$NS + pasos$IS + pasos$RS))))
chk <- pasos %>% group_by(periodo, tipo, anio, rubro8) %>%
  summarise(s = sum(RS_comuna), E_T = sum(E), .groups = "drop") %>% filter(E_T > 0)
stopifnot(max(abs(chk$s)) < 1e-6)            # RS_comuna suma 0 dentro de cada tipo

# ---------------- Capa 1: por tipología ----------------
E0 <- function(df) df %>% filter(anio == as.integer(substr(periodo, 1, 4)))
cap1 <- pasos %>% group_by(periodo, tipo) %>%
  summarise(across(c(dE, NS, IS, RS), sum), .groups = "drop") %>%
  left_join(E0(pasos) %>% group_by(periodo, tipo) %>% summarise(E0 = sum(E), n_comunas = n_distinct(comuna), .groups = "drop"),
            by = c("periodo", "tipo")) %>%
  rename(CO = dE) %>%
  mutate(across(c(CO, NS, IS, RS), ~ 100 * .x / E0, .names = "{.col}_pct"),
         tipo = factor(tipo, TIPOS)) %>%
  arrange(periodo, tipo)
cat("\nCapa 1 — Shift-share dinámico por tipología (% del empleo inicial):\n")
print(as.data.frame(cap1 %>% select(periodo, tipo, n_comunas, E0, ends_with("_pct")) %>%
                      mutate(across(ends_with("_pct"), ~ round(.x, 1)))))

area_tot <- pasos %>% group_by(periodo) %>% summarise(across(c(dE, NS, IS, RS), sum)) %>%
  left_join(E0(pasos) %>% group_by(periodo) %>% summarise(E0 = sum(E)), by = "periodo") %>%
  mutate(across(c(dE, NS, IS, RS), ~ round(100 * .x / E0, 1)))
cat("\nÁrea de estudio completa (%):\n"); print(area_tot)

# RS por rubro dentro de cada tipología (qué sectores explican la ventaja/desventaja)
cap1_rubro <- pasos %>% group_by(periodo, tipo, rubro8) %>%
  summarise(across(c(dE, NS, IS, RS), sum), .groups = "drop") %>%
  left_join(cap1 %>% select(periodo, tipo, E0) %>% mutate(tipo = as.character(tipo)), by = c("periodo", "tipo")) %>%
  mutate(RS_pct = 100 * RS / E0, IS_pct = 100 * IS / E0, tipo = factor(tipo, TIPOS))

# ---------------- Capa 2: RS anidado por comuna ----------------
cap2 <- pasos %>% group_by(periodo, comuna, tipo) %>%
  summarise(across(c(dE, NS, IS, RS, RS_tipo, RS_comuna), sum), .groups = "drop") %>%
  left_join(E0(pasos) %>% group_by(periodo, comuna) %>% summarise(E0 = sum(E), .groups = "drop"),
            by = c("periodo", "comuna")) %>%
  mutate(across(c(dE, NS, IS, RS, RS_tipo, RS_comuna), ~ 100 * .x / E0, .names = "{.col}_pct"),
         tipo = factor(tipo, TIPOS))
# Mismo criterio que v2_00 (MIN_TRAB): con empleo inicial muy bajo la tasa carece de sentido
# (p.ej. Hualaihué 2010: 14 trabajadores -> RS ~ 1.400%). Se excluyen solo de la capa 2.
excl_cap2 <- cap2 %>% filter(E0 < MIN_TRAB) %>% select(periodo, comuna, E0)
cat("\nExcluidas de la capa 2 (empleo inicial <", MIN_TRAB, "):\n"); print(excl_cap2)
cap2 <- cap2 %>% filter(E0 >= MIN_TRAB)

# η² (ponderado opcionalmente por el empleo inicial, para que no dominen comunas pequeñas y volátiles)
eta2 <- function(y, g, w = rep(1, length(y))) {
  mu <- weighted.mean(y, w); m <- tapply(y * w, g, sum) / tapply(w, g, sum)
  sum(tapply(w, g, sum) * (m - mu)^2) / sum(w * (y - mu)^2)
}
set.seed(2026)
cap2_var <- cap2 %>% group_by(periodo) %>% group_modify(function(d, ...) {
  obs  <- eta2(d$RS_pct, d$tipo)
  perm <- replicate(N_PERM, eta2(d$RS_pct, sample(d$tipo)))
  obs_w  <- eta2(d$RS_pct, d$tipo, d$E0)
  perm_w <- replicate(N_PERM, eta2(d$RS_pct, sample(d$tipo), d$E0))
  tibble(eta2_tipologia = obs, p_permutacion = (1 + sum(perm >= obs)) / (N_PERM + 1),
         eta2_ponderado = obs_w, p_perm_ponderado = (1 + sum(perm_w >= obs_w)) / (N_PERM + 1),
         eta2_esperado_H0 = (nlevels(droplevels(d$tipo)) - 1) / (nrow(d) - 1),
         n_comunas = nrow(d))
}) %>% ungroup()
cat("\nCapa 2 — Varianza del RS comunal (% del empleo inicial) explicada por la tipología:\n")
print(cap2_var %>% mutate(across(where(is.double), ~ round(.x, 3))))

# ---------------- Robustez: quiebre de registro 2017->2018 ----------------
# El empleo nacional en O y P salta ~32% en 2018 por un cambio administrativo.
# Se compara 2015-2020 con y sin ese paso (la versión principal lo excluye).
rob_2018 <- pasos_completo %>% filter(periodo == "2015–2020") %>%
  mutate(version = ifelse(anio == 2017, "excluido", "incluido")) %>%
  { bind_rows(mutate(., escenario = "Con paso 2017-18"),
              filter(., anio != 2017) %>% mutate(escenario = "Sin paso 2017-18")) } %>%
  group_by(escenario, tipo) %>% summarise(across(c(dE, NS, IS, RS), sum), .groups = "drop") %>%
  left_join(cap1 %>% filter(periodo == "2015–2020") %>% transmute(tipo = as.character(tipo), E0), by = "tipo") %>%
  mutate(across(c(dE, NS, IS, RS), ~ round(100 * .x / E0, 1)), tipo = factor(tipo, TIPOS)) %>%
  arrange(tipo, escenario)
cat("\nRobustez 2015–2020 sin el paso 2017→2018 (% del empleo inicial):\n"); print(as.data.frame(rob_2018))

# Discontinuidad nacional 2014->2015 entre A (primario) y N (servicios de apoyo):
# A cae de ~430 000 a ~372 000 y N sube de ~320 000 a ~470 000 en un año (reclasificación
# o externalización). Se compara 2010-2015 con y sin ese paso (la versión principal lo excluye).
rob_2015 <- pasos_completo %>% filter(periodo == "2010–2015") %>%
  { bind_rows(mutate(., escenario = "Con paso 2014-15"),
              filter(., anio != 2014) %>% mutate(escenario = "Sin paso 2014-15")) } %>%
  group_by(escenario, tipo) %>% summarise(across(c(dE, NS, IS, RS), sum), .groups = "drop") %>%
  left_join(cap1 %>% filter(periodo == "2010–2015") %>% transmute(tipo = as.character(tipo), E0), by = "tipo") %>%
  mutate(across(c(dE, NS, IS, RS), ~ round(100 * .x / E0, 1)), tipo = factor(tipo, TIPOS)) %>%
  arrange(tipo, escenario)
cat("\nRobustez 2010–2015 sin el paso 2014→2015 (% del empleo inicial):\n"); print(as.data.frame(rob_2015))

# ---------------- Figuras ----------------
p1 <- cap1 %>% select(periodo, tipo, NS_pct, IS_pct, RS_pct, CO_pct) %>%
  pivot_longer(c(NS_pct, IS_pct, RS_pct), names_to = "efecto", values_to = "valor") %>%
  mutate(efecto = factor(efecto, c("NS_pct", "IS_pct", "RS_pct"),
                         c("Nacional (NS)", "Industrial (IS)", "Regional (RS)"))) %>%
  ggplot(aes(x = periodo, y = valor, fill = efecto)) +
  geom_col(position = "stack", width = .7) +
  geom_point(aes(y = CO_pct, shape = "Crecimiento observado"), size = 2.3, colour = "black", fill = "white") +
  geom_hline(yintercept = 0, linewidth = .3) +
  facet_wrap(~ tipo, ncol = 3) +
  scale_fill_manual(values = c("#bdbdbd", "#e6550d", "#3182bd")) +
  scale_shape_manual(values = 21) +
  labs(x = NULL, y = "% del empleo inicial del tipo", fill = NULL, shape = NULL,
       title = "Shift-share dinámico por tipología funcional",
       subtitle = "Tipología fija al inicio de cada periodo; referencia: resto de Chile; 2010–2015 sin el paso 2014→2015 y 2015–2020 sin el paso 2017→2018 (quiebres de registro)") +
  theme_minimal() + theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUT, "01_shiftshare_tipologia.png"), p1, width = 12, height = 8, dpi = 300, bg = "white")

lim <- max(abs(cap1_rubro$RS_pct))
p2 <- ggplot(cap1_rubro, aes(periodo, rubro8, fill = RS_pct)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.1f", RS_pct)), size = 2.6) +
  facet_wrap(~ tipo, ncol = 3) +
  scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac", limits = c(-lim, lim),
                       name = "RS (% empleo\ninicial del tipo)") +
  labs(x = NULL, y = NULL, title = "Efecto regional por rubro y tipología") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUT, "02_RS_rubro_tipologia.png"), p2, width = 13, height = 8, dpi = 300, bg = "white")

p3 <- ggplot(cap2, aes(x = tipo, y = RS_pct, colour = tipo)) +
  geom_hline(yintercept = 0, linewidth = .3) +
  geom_boxplot(outlier.shape = NA, alpha = .2) +
  geom_jitter(width = .15, size = 1.6) +
  geom_text(data = cap2_var, aes(x = 0.6, y = Inf, label = sprintf("η² = %.2f (p = %.3f) | ponderado: %.2f (p = %.3f)", eta2_tipologia, p_permutacion, eta2_ponderado, p_perm_ponderado)),
            inherit.aes = FALSE, hjust = 0, vjust = 1.3, size = 3) +
  facet_wrap(~ periodo, ncol = 2) +
  scale_colour_manual(values = PALETA, guide = "none") +
  labs(x = NULL, y = "RS comunal (% del empleo inicial)",
       title = "¿La ventaja competitiva local es propia del tipo o de cada comuna?",
       subtitle = "η²: proporción de la varianza del RS comunal explicada por la tipología (p por permutación)") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(OUT, "03_RS_comunal_por_tipo.png"), p3, width = 12, height = 9, dpi = 300, bg = "white")

write_xlsx(list(
  capa1_tipologia = cap1, capa1_rubro = cap1_rubro, area_total = area_tot,
  capa2_comuna = cap2, capa2_varianza = cap2_var, robustez_2018 = rob_2018, robustez_2015 = rob_2015, tipologias_imputadas = imputadas, excluidas_capa2 = excl_cap2
), file.path(OUT, "shiftshare_v2.xlsx"))
saveRDS(list(pasos = pasos, cap1 = cap1, cap1_rubro = cap1_rubro, cap2 = cap2, cap2_var = cap2_var),
        file.path(OUT, "shiftshare_v2.rds"))

# ---------------- Material suplementario: comuna x rubro x periodo ----------------
# Detalle para consulta (no para el cuerpo del artículo): valores absolutos y en %
# del empleo inicial de la celda, con marca de base pequeña.
BASE_PEQUENA <- 10
suplem <- pasos %>%
  group_by(periodo, comuna, tipo, rubro8) %>%
  summarise(E0 = E[anio == min(anio)], CO = sum(dE), NS = sum(NS), IS = sum(IS), RS = sum(RS),
            RS_tipo = sum(RS_tipo), RS_comuna = sum(RS_comuna), .groups = "drop") %>%
  mutate(across(c(CO, NS, IS, RS, RS_tipo, RS_comuna), ~ ifelse(E0 > 0, 100 * .x / E0, NA_real_), .names = "{.col}_pct"),
         base_pequena = E0 < BASE_PEQUENA,
         tipo = factor(tipo, TIPOS)) %>%
  arrange(periodo, comuna, rubro8)
notas <- tibble(nota = c(
  "Shift-share dinámico: suma de descomposiciones anuales dentro de cada periodo (Barff & Knight, 1988).",
  "Referencia: resto de Chile (excluye Los Ríos, Los Lagos y registros sin comuna).",
  "2015–2020 excluye el paso 2017→2018 (quiebre de registro en O y P); por eso CO no equivale a E2020 − E2015.",
  "RS = CO − NS − IS; RS = RS_tipo (ventaja común de la tipología) + RS_comuna (desempeño frente a sus pares).",
  "Tipología: fija al inicio de cada periodo (v2_02_tipologia_funcional.R).",
  paste0("base_pequena = TRUE si el empleo inicial de la celda es < ", BASE_PEQUENA, ": las tasas (%) son poco confiables."),
  "Columnas _pct: % del empleo inicial de la celda (NA si E0 = 0)."))
write_xlsx(list(notas = notas, comuna_rubro_periodo = suplem),
           file.path(OUT, "S1_shiftshare_comuna_rubro_periodo.xlsx"))

# ---------------- Mapa del RS comunal por periodo ----------------
geo <- st_read("comunasloslagoslosrios.gpkg", quiet = TRUE) |> st_zm() |> st_make_valid()
TOPE <- 50   # la escala se satura en ±50% para que pocos casos extremos no dominen el color
mapa_rs <- geo %>% select(Comuna) %>%
  inner_join(pasos %>% distinct(periodo, comuna) %>%
               left_join(cap2 %>% select(periodo, comuna, RS_pct), by = c("periodo", "comuna")) %>%
               rename(Comuna = comuna),
             by = "Comuna", relationship = "many-to-many")
p4 <- ggplot(mapa_rs) +
  geom_sf(aes(fill = RS_pct), colour = "grey35", linewidth = .1) +
  scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac", midpoint = 0,
                       limits = c(-TOPE, TOPE), oob = scales::squish, na.value = "grey80",
                       name = paste0("RS (% del empleo\ninicial; saturado en ±", TOPE, ")")) +
  facet_wrap(~ periodo, nrow = 1) +
  labs(title = "Efecto regional (RS) comunal por periodo",
       subtitle = "Gris: empleo inicial < 100 trabajadores (sin tasa); 2010–2015 y 2015–2020 sin los pasos con quiebre de registro") +
  theme_void() + theme(legend.position = "right")
ggsave(file.path(OUT, "04_mapa_RS_comunal.png"), p4, width = 14, height = 6, dpi = 300, bg = "white")
