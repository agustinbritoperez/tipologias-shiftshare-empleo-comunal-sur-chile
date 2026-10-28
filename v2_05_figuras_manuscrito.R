# ============================================================
# v2_05_figuras_manuscrito.R — Figuras y material suplementario
# ------------------------------------------------------------
# Formato: título y fuente FUERA de la imagen (van en el pie del
# manuscrito); mapas con coordenadas en grados-minutos-segundos, escala gráfica
# y leyenda explícita, sin rosa de los vientos; sin datos de autoría.
# Requiere haber corrido v2_02, v2_03 y v2_04.
# ============================================================
source("v2_00_datos.R")
suppressPackageStartupMessages({ library(ggspatial); library(patchwork) })

OUT <- "manuscrito/figuras"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CRS_UTM <- 32718   # WGS 84 / UTM 18S: escala gráfica correcta en metros

TIPOS  <- c("Centro regional", "Centro subregional", "Intermedia agro-comercial",
            "Primario especializado", "Rural dependiente de servicios públicos",
            "Periferia cordillerana-turística")
PALETA <- setNames(c("#7a0177", "#dd3497", "#fe9929", "#4d9221", "#4575b4", "#74add1"), TIPOS)
tema   <- theme_minimal(base_size = 10) + theme(legend.position = "bottom")

# Coordenadas en grados, minutos y segundos
dms <- function(x, pos, neg) {
  s <- ifelse(x < 0, neg, pos); a <- abs(x)
  g <- floor(a); m <- floor((a - g) * 60); sec <- round(((a - g) * 60 - m) * 60)
  sprintf("%d°%02d'%02d\" %s", g, m, sec, s)
}
ejes_dms <- list(
  scale_x_continuous(breaks = c(-74, -73, -72), labels = function(x) dms(x, "E", "O")),
  scale_y_continuous(breaks = seq(-44, -40, 1), labels = function(y) dms(y, "N", "S")),
  annotation_scale(location = "bl", width_hint = .3, text_cex = .6, style = "ticks"),
  theme(axis.text = element_text(size = 6), axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_line(colour = "grey85", linewidth = .2), axis.title = element_blank())
)
geo <- st_read("comunasloslagoslosrios.gpkg", layer = "comunas", quiet = TRUE) |>
  st_zm() |> st_make_valid() |> st_transform(CRS_UTM)
guardar <- function(p, nombre, w, h) ggsave(file.path(OUT, nombre), p, width = w, height = h, dpi = 300, bg = "white")

# Abreviaturas de las comunas nombradas en el texto (clave en la Figura 1;
# etiquetas en los mapas de las Figuras 3 y 6)
ABREV <- c("Valdivia" = "VAL", "Osorno" = "OSO", "Puerto Montt" = "PMO", "Puerto Varas" = "PVA",
           "Castro" = "CAS", "Ancud" = "ANC", "Frutillar" = "FRU", "La Unión" = "LUN",
           "Llanquihue" = "LLA", "Chonchi" = "CHO", "Curaco de Vélez" = "CUR",
           "San Juan de la Costa" = "SJC", "Lago Ranco" = "LRA", "Máfil" = "MAF",
           "Río Bueno" = "RBU", "Maullín" = "MAU", "Futaleufú" = "FUT", "Puyehue" = "PUY")
stopifnot(all(names(ABREV) %in% geo$Comuna))
etq_abrev <- geo %>% filter(Comuna %in% names(ABREV)) %>% st_point_on_surface() %>%
  mutate(ab = ABREV[Comuna])
capa_abrev <- ggrepel::geom_text_repel(
  data = etq_abrev, aes(label = ab, geometry = geom), stat = "sf_coordinates",
  size = 1.7, fontface = "bold", colour = "grey5", bg.colour = "white", bg.r = .12,
  segment.size = .15, min.segment.length = .2, box.padding = .1, max.overlaps = Inf, seed = 1)

# Población censada 2024 (INE, cuadro D1, hoja 2) para la Figura 1
censo24 <- read.csv("censo2024_poblacion_comunal.csv", encoding = "UTF-8") %>%
  select(Comuna = comuna, pob = pob2024)

# ---- Figura 1: área de estudio ----
# Regiones diferenciadas, nombres de las 42 comunas (con abreviatura si aparece
# en las Figuras 3 y 6), población censada 2024 en círculos proporcionales y
# recuadro de localización en Chile.
regiones <- geo %>% group_by(Region) %>% summarise(.groups = "drop")
centros <- geo %>% st_point_on_surface() %>%
  left_join(censo24, by = "Comuna") %>%
  mutate(Comuna = ifelse(Comuna %in% names(ABREV), paste0(Comuna, " (", ABREV[Comuna], ")"), Comuna)) %>%
  arrange(desc(pob))                                   # círculos grandes debajo
stopifnot(!anyNA(centros$pob))
COL_REG <- c("Región de Los Ríos" = "#cfe3c4", "Región de Los Lagos" = "#cddff0")
etq_reg <- data.frame(x = c(-74.35, -74.75), y = c(-39.75, -42.95),
                      lab = c("REGIÓN DE\nLOS RÍOS", "REGIÓN DE\nLOS LAGOS")) |>
  st_as_sf(coords = c("x", "y"), crs = 4326)
f1_main <- ggplot() +
  geom_sf(data = geo, aes(fill = Region), colour = "grey55", linewidth = .15, show.legend = FALSE) +
  geom_sf(data = regiones, fill = NA, colour = "grey20", linewidth = .35) +
  geom_sf(data = centros, aes(size = pob), shape = 21, fill = alpha("#c2410c", .45), colour = "#7c2d12", stroke = .3) +
  ggrepel::geom_text_repel(data = centros, aes(label = Comuna, geometry = geom), stat = "sf_coordinates",
                           size = 2.1, colour = "grey10", segment.size = .15, segment.colour = "grey40",
                           min.segment.length = .3, box.padding = .18, point.padding = .1,
                           max.overlaps = Inf, seed = 1, bg.colour = "white", bg.r = .1) +
  geom_sf_text(data = etq_reg, aes(label = lab), size = 2.6, fontface = "bold.italic",
               colour = "grey25", lineheight = .9) +
  scale_fill_manual(values = COL_REG) +
  scale_size_area(max_size = 13, breaks = c(5000, 25000, 100000, 275000),
                  labels = function(x) format(x, big.mark = " ", scientific = FALSE),
                  name = "Población censada 2024") +
  coord_sf(crs = CRS_UTM, datum = 4326, xlim = c(-76.6, -71.45), ylim = c(-44.15, -39.2),
           default_crs = 4326, expand = FALSE) +
  scale_x_continuous(breaks = -76:-72, labels = function(x) dms(x, "E", "O")) +
  scale_y_continuous(breaks = -44:-40, labels = function(y) dms(y, "N", "S")) +
  annotation_scale(location = "br", width_hint = .25, text_cex = .6, style = "ticks") +
  guides(size = guide_legend(override.aes = list(fill = alpha("#c2410c", .45)))) +
  theme_minimal(base_size = 9) +
  theme(legend.position = c(.13, .97), legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", .85), colour = NA),
        legend.title = element_text(size = 7.5), legend.text = element_text(size = 6.5),
        panel.background = element_rect(fill = "#f4f8fb", colour = NA),
        axis.text = element_text(size = 6), axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_line(colour = "grey85", linewidth = .2), axis.title = element_blank())
# Recuadro: Chile continental y países vecinos (Natural Earth)
sa <- rnaturalearth::ne_countries(scale = 50, continent = "South America", returnclass = "sf")
area <- st_union(geo)
f1_ctx <- ggplot() +
  geom_sf(data = sa, fill = "grey92", colour = "grey60", linewidth = .1) +
  geom_sf(data = sa %>% filter(admin == "Chile"), fill = "grey75", colour = "grey35", linewidth = .15) +
  geom_sf(data = area, fill = "#c2410c", colour = "#7c2d12", linewidth = .2) +
  annotate("text", x = -64.5, y = -33, label = "ARGENTINA", size = 1.7, colour = "grey45") +
  annotate("text", x = -79, y = -38, label = "Océano\nPacífico", size = 1.7, colour = "grey45",
           fontface = "italic", lineheight = .85) +
  annotate("text", x = -76.5, y = -24, label = "CHILE", size = 2, fontface = "bold", colour = "grey20") +
  coord_sf(xlim = c(-83, -62), ylim = c(-56, -17), expand = FALSE) +
  theme_void() +
  theme(panel.background = element_rect(fill = "#f4f8fb", colour = "grey30", linewidth = .4))
f1 <- f1_main + inset_element(f1_ctx, left = .01, bottom = .02, right = .27, top = .52, align_to = "panel")
guardar(f1, "figura1_area_estudio.png", 7, 8.6)

# ---- Figura 2: tipología funcional 2006, 2015, 2024 ----
tip <- readRDS(file.path(RES, "Tipologia_funcional/tipologia_funcional.rds"))$ind
m2 <- geo %>% select(Comuna) %>%
  inner_join(tip %>% filter(anio %in% c(2006, 2015, 2024)) %>% select(Comuna = comuna, anio, tipo, certeza),
             by = "Comuna", relationship = "many-to-many")
f2 <- ggplot(m2) +
  geom_sf(data = geo, fill = "grey80", colour = NA) +              # comuna-año excluidas
  geom_sf(aes(fill = tipo), colour = "grey30", linewidth = .1) +
  capa_abrev +
  facet_wrap(~ anio, nrow = 1) +
  scale_fill_manual(values = PALETA, name = NULL, drop = FALSE) +
  coord_sf(crs = CRS_UTM, datum = 4326) + tema + ejes_dms +
  guides(fill = guide_legend(nrow = 2))
guardar(f2, "figura3_mapa_tipologias.png", 10, 6.5)

# ---- Figura 3: transiciones (aluvial) ----
f3 <- ggplot(tip %>% select(comuna, anio, tipo),
             aes(x = factor(anio), stratum = tipo, alluvium = comuna, fill = tipo)) +
  ggalluvial::geom_flow(alpha = .5) + ggalluvial::geom_stratum(linewidth = .2) +
  geom_text(stat = ggalluvial::StatStratum, aes(label = after_stat(count)), size = 2.8) +
  scale_fill_manual(values = PALETA, name = NULL) +
  labs(x = NULL, y = "Número de comunas") + tema + guides(fill = guide_legend(nrow = 2))
guardar(f3, "figura4_transiciones.png", 8, 5.5)

# ---- Figura 4: shift-share dinámico por tipología ----
ss <- readRDS(file.path(RES, "ShiftShare/shiftshare_v2.rds"))
f4 <- ss$cap1 %>%
  select(periodo, tipo, NS_pct, IS_pct, RS_pct, CO_pct) %>%
  pivot_longer(c(NS_pct, IS_pct, RS_pct), names_to = "efecto", values_to = "valor") %>%
  mutate(efecto = factor(efecto, c("NS_pct", "IS_pct", "RS_pct"),
                         c("Efecto nacional (NS)", "Efecto industrial (IS)", "Efecto regional (RS)"))) %>%
  ggplot(aes(periodo, valor, fill = efecto)) +
  geom_col(width = .7) + geom_hline(yintercept = 0, linewidth = .3) +
  geom_point(aes(y = CO_pct, shape = "Crecimiento observado"), size = 1.8, fill = "white") +
  facet_wrap(~ tipo, ncol = 3, labeller = label_wrap_gen(28)) +
  scale_fill_manual(values = c("#bdbdbd", "#e6550d", "#3182bd"), name = NULL) +
  scale_shape_manual(values = 21, name = NULL) +
  labs(x = NULL, y = "% del empleo inicial de la tipología") + tema +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
guardar(f4, "figura5_shiftshare_tipologia.png", 9, 6.5)

# ---- Figura 6: gradiente nacional de tamaño del efecto regional ----
# RS comunal centrado en la media de su región (equivale a efectos fijos regionales)
# frente a la población inicial; Chile sin la Región Metropolitana (v2_09).
nac <- readxl::read_excel(file.path(RES, "Nacional/nacional_tamano.xlsx"), "datos") %>%
  group_by(periodo, region) %>% mutate(RS_intra = RS_pct - mean(RS_pct)) %>% ungroup() %>%
  mutate(grupo = case_when(area_estudio & pob >= 100000 ~ "Centro regional del área de estudio",
                           area_estudio ~ "Otras comunas del área de estudio",
                           TRUE ~ "Resto de Chile (sin RM)"))
f5 <- ggplot(nac, aes(pob, RS_intra)) +
  geom_hline(yintercept = 0, linewidth = .3) +
  geom_point(aes(colour = grupo, size = grupo), alpha = .7) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black", linewidth = .6, se = TRUE) +
  facet_wrap(~ periodo, ncol = 2) +
  scale_x_log10(labels = scales::label_number(big.mark = " ")) +
  coord_cartesian(ylim = c(-60, 60)) +
  scale_colour_manual(values = c("#7a0177", "#fe9929", "grey65"), name = NULL) +
  scale_size_manual(values = c(2.4, 1.6, 1), name = NULL) +
  labs(x = "Población al inicio del periodo (escala logarítmica)",
       y = "RS respecto de la media regional (% del empleo inicial)") + tema
guardar(f5, "figura7_gradiente_nacional.png", 9, 6.5)

# ---- Figura 5: RS comunal por periodo ----
TOPE <- 50
m6 <- geo %>% select(Comuna) %>%
  inner_join(ss$pasos %>% distinct(periodo, comuna) %>%
               left_join(ss$cap2 %>% select(periodo, comuna, RS_pct), by = c("periodo", "comuna")) %>%
               rename(Comuna = comuna), by = "Comuna", relationship = "many-to-many")
f6 <- ggplot(m6) +
  geom_sf(aes(fill = RS_pct), colour = "grey35", linewidth = .1) +
  scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac", midpoint = 0,
                       limits = c(-TOPE, TOPE), oob = scales::squish, na.value = "grey80",
                       name = "RS (% del empleo inicial)") +
  capa_abrev +
  facet_wrap(~ periodo, nrow = 1) +
  coord_sf(crs = CRS_UTM, datum = 4326) + tema + ejes_dms +
  theme(legend.key.width = unit(1.2, "cm"))
guardar(f6, "figura6_mapa_RS.png", 12, 6)
cat("Figuras guardadas en", OUT, "\n")

# ---- Figura 2: matriz de consenso ordenada por tipología ----
tf  <- readRDS(file.path(RES, "Tipologia_funcional/tipologia_funcional.rds"))
Cm  <- tf$C
ind2 <- tf$ind
# orden: por tipología y, dentro de cada una, por el árbol de consenso
orden_arbol <- hclust(as.dist(1 - Cm), "average")$order
pos <- setNames(seq_along(orden_arbol), rownames(Cm)[orden_arbol])
ord <- ind2 %>% mutate(p = pos[id]) %>% arrange(tipo, p) %>% pull(id)
dm <- as_tibble(Cm[ord, ord], rownames = "i") %>% pivot_longer(-i, names_to = "j", values_to = "v") %>%
  mutate(i = factor(i, ord), j = factor(j, rev(ord)))
barra <- ind2 %>% select(id, tipo) %>% mutate(i = factor(id, ord))
f_cons <- ggplot() +
  geom_raster(data = dm, aes(i, j, fill = v)) +
  scale_fill_gradient(low = "white", high = "#08306b", limits = c(0, 1), name = "Frecuencia de\nco-asignación") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = barra, aes(x = i, y = -4, fill = tipo), height = 8) +
  scale_fill_manual(values = PALETA, name = NULL) +
  coord_cartesian(clip = "off") +
  labs(x = "Observaciones comuna-año (ordenadas por tipología)", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text = element_blank(), panel.grid = element_blank(), legend.position = "right") +
  guides(fill = guide_legend(ncol = 1))
guardar(f_cons, "figura2_matriz_consenso.png", 9, 7)

# ---- Material suplementario S3: validación de la tipología ----
SUP <- "manuscrito/suplementario"; dir.create(SUP, showWarnings = FALSE)
diag_k <- tf$diag_k
ari <- readxl::read_excel(file.path(RES, "Tipologia_funcional/tipologia_funcional.xlsx"), "ari_algoritmos") %>%
  mutate(metodo = recode(metodo, ward = "Ward", pam = "PAM (medoides)", kmeans = "k-medias", gmm = "Mezclas gaussianas"))
s3a <- ggplot(diag_k, aes(k, PAC)) + geom_line() + geom_point() +
  geom_point(data = filter(diag_k, k == 5), colour = "firebrick", size = 3) +
  labs(x = "Número de grupos (k)", y = "PAC (menor = más estable)", title = "A. Proporción de pares ambiguos") + tema
s3b <- ggplot(ari, aes(reorder(metodo, ARI), ARI)) + geom_col(fill = "grey45", width = .6) +
  geom_text(aes(label = sprintf("%.2f", ARI)), hjust = -.2, size = 3) + coord_flip(ylim = c(0, 1)) +
  labs(x = NULL, y = "Índice de Rand ajustado respecto del consenso", title = "B. Concordancia de cada algoritmo") + tema
cert <- geo %>% select(Comuna) %>%
  inner_join(ind2 %>% group_by(comuna) %>% summarise(certeza = mean(certeza)) %>% rename(Comuna = comuna), by = "Comuna")
s3c <- ggplot(cert) + geom_sf(aes(fill = certeza), colour = "grey35", linewidth = .1) +
  scale_fill_distiller(palette = "Greens", direction = 1, limits = c(0, 1), name = "Certeza media") +
  coord_sf(crs = CRS_UTM, datum = 4326) + tema + ejes_dms + labs(title = "C. Certeza de asignación (media 2006-2024)")
s3 <- (s3a / s3b) | s3c
ggsave(file.path(SUP, "S3_validacion_tipologia.png"), s3, width = 12, height = 8, dpi = 300, bg = "white")
writexl::write_xlsx(list(PAC = diag_k, ARI_algoritmos = ari,
                         certeza_comunal = st_drop_geometry(cert)),
                    file.path(SUP, "S3_validacion_tipologia.xlsx"))
# Reunir el material suplementario
file.copy(file.path(RES, "ShiftShare/S1_shiftshare_comuna_rubro_periodo.xlsx"), SUP, overwrite = TRUE)
file.copy(file.path(RES, "Perfiles_tipologias/S2_perfiles_por_anio.xlsx"), SUP, overwrite = TRUE)
cat("Figura 2 y material suplementario (S1-S3) guardados\n")
