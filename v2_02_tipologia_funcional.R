# ============================================================
# v2_02_tipologia_funcional.R — Tipología funcional comunal (robusta)
# ------------------------------------------------------------
# La tipología por composición sectorial (v2_01) no distingue jerarquía
# urbana: las proporciones no miden tamaño ni centralidad. Aquí cada
# comuna-año se describe con tres bloques, equilibrados mediante Análisis
# Factorial Múltiple (MFA; Escofier & Pagès, 1994):
#   1. Composición: CLR de 8 grandes rubros (v2_00_datos.R)
#   2. Escala relativa: log(peso en el empleo del área), log(peso en la población)
#      (relativos al año -> no capturan el crecimiento general del periodo)
#   3. Función urbana: diversidad de Shannon (15 rubros) y empleo formal por habitante
# Robustez: clustering de consenso (Monti et al., 2003) con 4 algoritmos
# (Ward, PAM, k-medias, GMM) x 100 submuestras del 80% de las comunas (con todos sus años). k se elige por
# mínimo PAC (Șenbabaoğlu et al., 2014). No hay anclas ni etiquetas a priori.
# ============================================================
source("v2_00_datos.R")
suppressPackageStartupMessages({ library(FactoMineR); library(mclust) })

# Variante de sensibilidad (TIPO_SIN_DOMICILIO=1): sin variables basadas en el empleo
# total (participación en el empleo y empleo por habitante), sesgadas por el domicilio
# tributario; la escala se mide solo con población.
SIN_DOMICILIO <- Sys.getenv("TIPO_SIN_DOMICILIO", "0") == "1"
OUT   <- file.path(RES, if (SIN_DOMICILIO) "Tipologia_sin_domicilio" else "Tipologia_funcional")
B     <- 100
FRAC  <- 0.8
K_MAX <- 8
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---------------- Indicadores ----------------
# Diversidad: Shannon sobre 15 rubros (A-N por separado; O-S en un bloque por los
# quiebres de registro descritos en v2_00_datos.R)
shannon <- base %>%
  filter(region %in% REGIONES, anio %in% ANIOS, !grepl("^Sin informaci", rubro, ignore.case = TRUE)) %>%
  mutate(letra = substr(rubro, 1, 1)) %>%
  filter(letra %in% names(GRUPOS)) %>%
  mutate(r15 = ifelse(letra %in% c("O", "P", "Q", "R", "S"), "OS", letra)) %>%
  group_by(anio, comuna, r15) %>% summarise(t = sum(trabajadores), .groups = "drop") %>%
  group_by(anio, comuna) %>% mutate(s = t / sum(t)) %>%
  summarise(shannon = -sum(s[s > 0] * log(s[s > 0])), .groups = "drop")

ind <- ancha %>%
  select(id, comuna, region, anio, total) %>%
  left_join(pob, by = c("comuna", "anio")) %>%
  left_join(shannon, by = c("comuna", "anio")) %>%
  group_by(anio) %>%
  mutate(tam_rel = log(total / sum(total)), pob_rel = log(pob / sum(pob))) %>%
  ungroup() %>%
  mutate(emp_pc = total / pob)
stopifnot(!anyNA(ind), identical(ind$id, rownames(CLR)))

VARS_ESC <- if (SIN_DOMICILIO) "pob_rel" else c("tam_rel", "pob_rel")
VARS_FUN <- if (SIN_DOMICILIO) "shannon" else c("shannon", "emp_pc")
df_mfa <- data.frame(CLR, ind[, c(VARS_ESC, VARS_FUN)], row.names = ind$id, check.names = FALSE)
mfa <- MFA(df_mfa, group = c(ncol(CLR), length(VARS_ESC), length(VARS_FUN)), type = c("c", "s", "s"),
           name.group = c("Composicion", "Escala", "Funcion_urbana"),
           ncp = ncol(df_mfa), graph = FALSE)
ncp <- which(mfa$eig[, 3] >= 90)[1]                       # componentes hasta 90% de inercia
Y   <- mfa$ind$coord[, 1:ncp]
cat("MFA: componentes retenidos =", ncp, "(", round(mfa$eig[ncp, 3], 1), "% inercia )\n")

# ---------------- Consenso ----------------
agrupar <- function(Ys, k, metodo) {
  switch(metodo,
    ward   = cutree(hclust(dist(Ys), "ward.D2"), k),
    pam    = cluster::pam(Ys, k, cluster.only = TRUE),
    kmeans = kmeans(Ys, k, nstart = 25)$cluster,
    gmm    = tryCatch(Mclust(Ys, G = k, modelNames = "VVI", verbose = FALSE)$classification,
                      error = function(e) NULL))
}
METODOS <- c("ward", "pam", "kmeans", "gmm")

COMUNA_OBS <- ind$comuna; COMUNAS_U <- unique(COMUNA_OBS)
consenso <- function(k) {
  n <- nrow(Y); juntos <- muestreados <- matrix(0, n, n)
  for (b in seq_len(B)) {
    # Submuestreo por COMUNA (con todos sus años): las 5 observaciones de una comuna
    # no son independientes, por lo que la comuna es la unidad de remuestreo.
    sel <- sample(COMUNAS_U, floor(FRAC * length(COMUNAS_U)))
    idx <- which(COMUNA_OBS %in% sel)
    for (m in METODOS) {
      cl <- agrupar(Y[idx, ], k, m)
      if (is.null(cl)) next
      muestreados[idx, idx] <- muestreados[idx, idx] + 1
      juntos[idx, idx]      <- juntos[idx, idx] + outer(cl, cl, "==")
    }
  }
  C <- juntos / pmax(muestreados, 1); dimnames(C) <- list(rownames(Y), rownames(Y)); C
}

set.seed(2026)
cons <- lapply(2:K_MAX, consenso); names(cons) <- 2:K_MAX
pac  <- sapply(cons, function(C) { v <- C[upper.tri(C)]; mean(v > 0.1 & v < 0.9) })
diag_k <- tibble(k = 2:K_MAX, PAC = pac)
print(diag_k %>% mutate(PAC = round(PAC, 3)))
# El PAC decrece monótonamente (no hay k "natural"); se fija k = 5 por
# interpretabilidad y porque todos los clústeres superan 0,57 de consenso interno.
K <- as.integer(Sys.getenv("K_CLUST", "5"))
cat("k elegido:", K, "\n")

C  <- cons[[as.character(K)]]
hc <- hclust(as.dist(1 - C), "average")
cl <- cutree(hc, K)
certeza <- sapply(seq_along(cl), function(i) mean(C[i, cl == cl[i] & seq_along(cl) != i]))

# Concordancia de cada algoritmo (datos completos) con el consenso
ari <- sapply(METODOS, function(m) { a <- agrupar(Y, K, m); if (is.null(a)) NA else adjustedRandIndex(a, cl) })
cat("ARI de cada algoritmo vs consenso:", paste(names(ari), round(ari, 2), collapse = " | "), "\n")

ind <- ind %>% mutate(cluster = cl, certeza = certeza)

# ---------------- Etiquetado a posteriori (reglas sobre los perfiles) ----------------
# Reglas explícitas para que el nombre no dependa de la numeración del árbol.
car <- ind %>% bind_cols(as_tibble(P)) %>% group_by(cluster) %>%
  summarise(emp = median(total), prim = mean(Primario), aloj = mean(`Alojamiento y comidas`),
            pub = mean(`Serv. públicos, sociales y personales`))
# Asignación secuencial: cada regla elige entre los clústeres aún sin nombre (orden fijo).
stopifnot(K == 5)
et <- setNames(character(K), car$cluster)
reglas <- list("Centro urbano" = "emp", "Primario especializado" = "prim",
               "Periferia cordillerana-turística" = "aloj", "Rural dependiente de servicios públicos" = "pub")
for (nombre in names(reglas)) {
  libres <- car %>% filter(et[as.character(cluster)] == "")
  et[as.character(libres$cluster[which.max(libres[[reglas[[nombre]]]])])] <- nombre
}
et[et == ""] <- "Intermedia agro-comercial"
cat("\nPerfiles usados para etiquetar:\n"); print(car %>% mutate(tipo = et[as.character(cluster)]))

# Subdivisión del nivel urbano por criterio externo explícito (no por clúster):
# el consenso no separa a las mayores ciudades del resto de centros (ver informe),
# por lo que la jerarquía se define por umbral poblacional.
UMBRAL_POB <- 100000
TIPOS <- c("Centro regional", "Centro subregional", "Intermedia agro-comercial",
           "Primario especializado", "Rural dependiente de servicios públicos",
           "Periferia cordillerana-turística")
ind <- ind %>% mutate(
  tipo = et[as.character(cluster)],
  tipo = case_when(tipo == "Centro urbano" & pob >= UMBRAL_POB ~ "Centro regional",
                   tipo == "Centro urbano" ~ "Centro subregional",
                   TRUE ~ tipo),
  tipo = factor(tipo, levels = TIPOS))
stopifnot(!anyNA(ind$tipo))
cat("\nCentros regionales por año:\n")
print(ind %>% filter(tipo == "Centro regional") %>% group_by(anio) %>% summarise(comunas = paste(sort(comuna), collapse = ", ")))

# ---------------- Resultados ----------------
perfil <- ind %>%
  bind_cols(as_tibble(P)) %>%
  group_by(tipo) %>%
  summarise(n = n(), empleo_mediana = median(total), pob_mediana = median(pob),
            emp_pc = mean(emp_pc), shannon = mean(shannon), certeza_media = mean(certeza),
            across(all_of(colnames(P)), ~ 100 * mean(.x)))
cat("\nPerfil por clúster:\n"); print(as.data.frame(perfil %>% mutate(across(where(is.double), ~ round(.x, 2)))))

cat("\nComunas por tipo (2024):\n")
ind %>% filter(anio == 2024) %>% group_by(tipo) %>%
  summarise(comunas = paste(sort(comuna), collapse = ", ")) %>%
  { for (i in seq_len(nrow(.))) cat(as.character(.$tipo[i]), ":", .$comunas[i], "\n") }
cat("\nComunas por tipo y año:\n"); print(table(ind$anio, ind$tipo))
cat("\nComuna-año con certeza < 0,6:", sum(ind$certeza < 0.6), "de", nrow(ind), "\n")

# Transiciones
wide_cl <- ind %>% select(comuna, anio, tipo) %>% mutate(tipo = as.character(tipo)) %>%
  pivot_wider(names_from = anio, values_from = tipo, names_prefix = "a")
trans <- bind_rows(lapply(seq_len(length(ANIOS) - 1), function(i) {
  a0 <- paste0("a", ANIOS[i]); a1 <- paste0("a", ANIOS[i + 1])
  wide_cl %>% filter(!is.na(.data[[a0]]), !is.na(.data[[a1]])) %>%
    dplyr::count(desde = .data[[a0]], hacia = .data[[a1]]) %>% mutate(intervalo = paste0(ANIOS[i], "–", ANIOS[i + 1]))
}))
cols_a <- paste0("a", ANIOS)
wide_cl$cambios <- apply(wide_cl[, cols_a], 1, function(z) { z <- z[!is.na(z)]; sum(z[-1] != z[-length(z)]) })
cat("\nN° de cambios por comuna:\n"); print(table(wide_cl$cambios))

# Figuras
PALETA <- setNames(c("#7a0177", "#dd3497", "#fe9929", "#4d9221", "#4575b4", "#74add1"), TIPOS)
p_pac <- ggplot(diag_k, aes(k, PAC)) + geom_line() + geom_point() +
  geom_point(data = filter(diag_k, k == K), colour = "firebrick", size = 3) +
  labs(title = "Selección de k por consenso", subtitle = "PAC: proporción de pares ambiguos (menor = más estable)") +
  theme_minimal()
ggsave(file.path(OUT, "00_PAC.png"), p_pac, width = 7, height = 4.5, dpi = 300, bg = "white")

orden <- hc$order
Cm <- as_tibble(C[orden, orden], rownames = "i") %>% pivot_longer(-i, names_to = "j") %>%
  mutate(i = factor(i, rownames(C)[orden]), j = factor(j, rownames(C)[orden]))
p_cons <- ggplot(Cm, aes(i, j, fill = value)) + geom_raster() +
  scale_fill_gradient(low = "white", high = "#08306b", name = "Co-asignación") +
  labs(x = NULL, y = NULL, title = paste0("Matriz de consenso (k = ", K, ")")) +
  theme_minimal() + theme(axis.text = element_blank(), panel.grid = element_blank())
ggsave(file.path(OUT, "01_matriz_consenso.png"), p_cons, width = 7, height = 6, dpi = 300, bg = "white")

p_mfa <- as_tibble(mfa$ind$coord[, 1:2], rownames = "id") %>%
  left_join(ind %>% select(id, comuna, anio, tipo), by = "id") %>%
  ggplot(aes(Dim.1, Dim.2, colour = tipo)) + geom_point(alpha = .7) +
  ggrepel::geom_text_repel(data = ~ filter(.x, anio == 2024), aes(label = comuna), size = 2.5, max.overlaps = 30, show.legend = FALSE) +
  scale_colour_manual(values = PALETA) + labs(x = paste0("Dim 1 (", round(mfa$eig[1, 2], 1), "%)"), y = paste0("Dim 2 (", round(mfa$eig[2, 2], 1), "%)"),
       colour = "Tipología", title = "Comuna-año en el espacio MFA (etiquetas: 2024)") + theme_minimal()
ggsave(file.path(OUT, "02_mfa_individuos.png"), p_mfa, width = 10, height = 8, dpi = 300, bg = "white")

alluv <- ind %>% select(comuna, anio, tipo)
p_alluv <- ggplot(alluv, aes(x = factor(anio), stratum = tipo, alluvium = comuna, fill = tipo)) +
  geom_flow(alpha = .5) + geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(count)), size = 3) + scale_fill_manual(values = PALETA) +
  labs(x = NULL, y = "N° de comunas", fill = "Tipología", title = "Transiciones entre tipologías funcionales") + theme_minimal()
ggsave(file.path(OUT, "03_aluvial.png"), p_alluv, width = 9, height = 6, dpi = 300, bg = "white")

geo <- st_read("comunasloslagoslosrios.gpkg", quiet = TRUE) |> st_zm() |> st_make_valid()
mapa <- geo %>% select(Comuna) %>%
  inner_join(ind %>% select(Comuna = comuna, anio, tipo, certeza), by = "Comuna", relationship = "many-to-many")
p_mapa <- ggplot(mapa) + geom_sf(aes(fill = tipo, alpha = certeza), colour = "grey30", linewidth = .1) +
  scale_alpha(range = c(.35, 1), guide = "none") + scale_fill_manual(values = PALETA) +
  facet_wrap(~ anio, nrow = 1) + labs(fill = "Tipología", title = "Tipología funcional por año",
                                      subtitle = "Transparencia = menor certeza de asignación") +
  theme_void() + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 2))
ggsave(file.path(OUT, "04_mapas.png"), p_mapa, width = 14, height = 6, dpi = 300, bg = "white")

write_xlsx(list(
  asignacion = ind, perfil = perfil, diagnostico_k = diag_k,
  ari_algoritmos = tibble(metodo = names(ari), ARI = ari),
  cambios = wide_cl, transiciones = trans,
  mfa_inercia = as_tibble(mfa$eig, rownames = "dim"),
  mfa_grupos = as_tibble(mfa$group$contrib, rownames = "grupo")
), file.path(OUT, "tipologia_funcional.xlsx"))
saveRDS(list(ind = ind, mfa = mfa, C = C, cons = cons, diag_k = diag_k), file.path(OUT, "tipologia_funcional.rds"))
