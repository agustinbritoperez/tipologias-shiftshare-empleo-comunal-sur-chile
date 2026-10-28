# ============================================================
# v2_11_perfiles_tipologias.R — Perfiles sectoriales de las tipologías
# ------------------------------------------------------------
# (1) Perfil conjunto (todos los años): media y mediana de la participación de cada
#     gran sector, e indicadores de escala y función urbana, por tipología.
# (2) Observación representativa: comuna-año con mayor certeza de consenso dentro de
#     su tipo, excluyendo comuna-año con un salto de registro marcado (v2_00) en los
#     cuatro años previos al corte.
# (3) Material suplementario: perfiles por tipología y año, y estabilidad temporal
#     (correlación de cada perfil anual con el perfil conjunto de su tipo).
# ============================================================
Sys.unsetenv("AJUSTE_SALTOS")
source("v2_00_datos.R")
OUT <- file.path(RES, "Perfiles_tipologias")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SECTORES <- unique(GRUPOS)

ind <- readRDS(file.path(RES, "Tipologia_funcional/tipologia_funcional.rds"))$ind %>%
  bind_cols(as_tibble(100 * P[readRDS(file.path(RES, "Tipologia_funcional/tipologia_funcional.rds"))$ind$id, ]))

# (1) Perfil conjunto
larga <- ind %>% select(tipo, anio, comuna, all_of(SECTORES)) %>%
  pivot_longer(all_of(SECTORES), names_to = "sector", values_to = "pct")
perfil_sector <- larga %>% group_by(tipo, sector) %>%
  summarise(media = mean(pct), mediana = median(pct), .groups = "drop")
perfil_otros <- ind %>% group_by(tipo) %>%
  summarise(n = n(), empleo_mediana = median(total), pob_mediana = median(pob),
            shannon = mean(shannon), emp_pc = mean(emp_pc), certeza = mean(certeza), .groups = "drop")

# (2) Representativa: mayor certeza, sin saltos marcados en (anio-4, anio]
con_salto <- saltos %>% mutate(anio_salto = as.integer(substr(paso, 7, 10))) %>% distinct(comuna, anio_salto)
afectada <- ind %>% select(id, comuna, anio) %>%
  left_join(con_salto, by = "comuna", relationship = "many-to-many") %>%
  group_by(id) %>% summarise(afectada = any(!is.na(anio_salto) & anio_salto > anio - 4 & anio_salto <= anio))
repr <- ind %>% left_join(afectada, by = "id") %>% filter(!afectada) %>%
  group_by(tipo) %>% slice_max(certeza, n = 1, with_ties = FALSE) %>% ungroup() %>%
  transmute(tipo, representativa = paste(comuna, anio), certeza_repr = certeza)

# (3) Por año y estabilidad temporal
perfil_anual <- larga %>% group_by(tipo, anio, sector) %>%
  summarise(n = n(), media = mean(pct), mediana = median(pct), .groups = "drop")
estabilidad <- perfil_anual %>% left_join(perfil_sector %>% select(tipo, sector, media_conjunta = media),
                                          by = c("tipo", "sector")) %>%
  group_by(tipo, anio) %>% filter(first(n) >= 2) %>%
  summarise(r = cor(media, media_conjunta), .groups = "drop")

cat("Perfil sectorial conjunto (media / mediana, %):\n")
print(as.data.frame(perfil_sector %>% mutate(txt = sprintf("%.1f (%.1f)", media, mediana)) %>%
                      select(tipo, sector, txt) %>% pivot_wider(names_from = tipo, values_from = txt)), right = FALSE)
cat("\nOtros indicadores:\n"); print(as.data.frame(perfil_otros %>% mutate(across(where(is.double), ~ round(.x, 2)))))
cat("\nObservación representativa:\n"); print(as.data.frame(repr %>% mutate(certeza_repr = round(certeza_repr, 2))))
cat("\nEstabilidad temporal (correlación perfil anual vs conjunto; tipos con >= 2 comunas en el año):\n")
print(as.data.frame(estabilidad %>% group_by(tipo) %>% summarise(r_min = round(min(r), 3), r_media = round(mean(r), 3))))

writexl::write_xlsx(list(perfil_conjunto = perfil_sector, indicadores = perfil_otros, representativas = repr),
                    file.path(OUT, "perfiles_tipologias.xlsx"))
writexl::write_xlsx(list(
  notas = tibble(nota = c("Material suplementario S2: participación del empleo por gran sector (%), media y mediana por tipología y año.",
                          "Tipologías del agrupamiento conjunto 2006-2024 (v2_02_tipologia_funcional.R).",
                          "estabilidad: correlación entre el perfil medio anual y el perfil medio conjunto de la tipología (años con >= 2 comunas).")),
  perfil_por_anio = perfil_anual, estabilidad = estabilidad),
  file.path(OUT, "S2_perfiles_por_anio.xlsx"))
