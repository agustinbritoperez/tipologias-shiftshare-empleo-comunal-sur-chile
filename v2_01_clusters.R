# ============================================================
# v2_01_clusters.R — Tipologías comunales solo por composición (sensibilidad)
# ------------------------------------------------------------
# Decisiones del procedimiento:
#  1. Clúster AGRUPADO (pooled): todas las comuna-año en un solo análisis,
#     de modo que las tipologías son fijas y las transiciones son reales
#     (no artefactos de re-etiquetar clústeres independientes por año).
#  2. Datos composicionales: reemplazo multiplicativo de ceros (0,5
#     trabajador) + transformación CLR; Ward sobre distancia de Aitchison
#     usando TODAS las dimensiones (el PCA solo se usa para visualizar).
#  3. 19 rubros CIIU -> 8 grandes rubros (reduce ceros estructurales).
#  4. k elegido con Silhouette, Calinski-Harabasz y estabilidad bootstrap
#     (Jaccard, fpc::clusterboot); etiquetas asignadas a posteriori.
# ============================================================
source("v2_00_datos.R")
K         <- as.integer(Sys.getenv("K_CLUST", "4"))  # ver diagnostico_k.xlsx; pendiente de confirmar
OUT       <- file.path(RES, "Clusters")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

D   <- dist(CLR)                                             # distancia de Aitchison
hc  <- hclust(D, method = "ward.D2")

# ---------------- Diagnóstico de k ----------------
boot_jaccard <- function(k) {
  cb <- clusterboot(CLR, B = 200, clustermethod = hclustCBI, method = "ward.D2",
                    k = k, count = FALSE, showplots = FALSE)
  cb$bootmean
}
diag_k <- bind_rows(lapply(2:8, function(k) {
  cl <- cutree(hc, k)
  jac <- boot_jaccard(k)
  tibble(k = k,
         silhouette = mean(silhouette(cl, D)[, "sil_width"]),
         calinski_harabasz = calinhara(CLR, cl),
         jaccard_min = min(jac), jaccard_medio = mean(jac),
         tamano_min = min(table(cl)))
}))

cat("\nComuna-año excluidas (<", MIN_TRAB, "trabajadores):\n"); print(excluidas)
cat("\nObservaciones en el análisis:", nrow(CLR), "\n\nDiagnóstico de k:\n")
print(diag_k %>% mutate(across(where(is.double), ~ round(.x, 3))))
write_xlsx(list(diagnostico_k = diag_k, excluidas = excluidas), file.path(OUT, "diagnostico_k.xlsx"))

p_diag <- diag_k %>%
  select(k, silhouette, jaccard_medio, jaccard_min) %>%
  pivot_longer(-k) %>%
  ggplot(aes(k, value, colour = name)) + geom_line() + geom_point() +
  geom_hline(yintercept = 0.75, linetype = "dotted") +
  labs(title = "Selección del número de tipologías",
       subtitle = "Línea punteada: Jaccard 0,75 (umbral usual de clúster estable)",
       y = NULL, colour = NULL) + theme_minimal()
ggsave(file.path(OUT, "00_diagnostico_k.png"), p_diag, width = 8, height = 5, dpi = 300, bg = "white")


# ---------------- Solución final ----------------
ancha$cluster <- cutree(hc, K)
jac_final <- boot_jaccard(K)
cat("\nJaccard bootstrap por clúster (k =", K, "):", round(jac_final, 2), "\n")
sil_obs <- silhouette(ancha$cluster, D)
ancha$silhouette <- sil_obs[, "sil_width"]

# Perfil de cada clúster: media de proporciones y razón respecto del área (tipo cociente de localización)
prop_df <- as_tibble(P) %>% mutate(cluster = ancha$cluster)
area    <- colSums(X) / sum(X)
perfil  <- prop_df %>% group_by(cluster) %>% summarise(across(everything(), mean), n = n())
perfil_long <- perfil %>% select(-n) %>% pivot_longer(-cluster, names_to = "rubro", values_to = "prop") %>%
  mutate(razon_area = prop / area[rubro])
cat("\nPerfil medio (% del empleo) por clúster:\n")
print(perfil %>% mutate(across(-c(cluster, n), ~ round(100 * .x, 1))) %>% t())

p_perfil <- ggplot(perfil_long, aes(rubro, factor(cluster), fill = log2(razon_area))) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.0f%%", 100 * prop)), size = 3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       name = "log2(prop. / área)") +
  labs(x = NULL, y = "Clúster", title = "Perfil sectorial de cada clúster",
       subtitle = "Texto: % medio del empleo; color: sobre/sub-representación respecto del área de estudio") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(OUT, "01_perfil_clusters.png"), p_perfil, width = 10, height = 4, dpi = 300, bg = "white")

# Biplot PCA (solo visual)
pca <- prcomp(CLR, scale. = FALSE)
var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
ind <- as_tibble(pca$x[, 1:2]) %>% mutate(cluster = factor(ancha$cluster), anio = ancha$anio)
cargas <- as_tibble(pca$rotation[, 1:2], rownames = "rubro") %>% mutate(across(c(PC1, PC2), ~ .x * 4))
p_pca <- ggplot(ind, aes(PC1, PC2)) +
  geom_point(aes(colour = cluster, shape = factor(anio)), alpha = .8) +
  geom_segment(data = cargas, aes(0, 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(.15, "cm")), colour = "grey40") +
  geom_text(data = cargas, aes(PC1, PC2, label = rubro), size = 3, vjust = -.4) +
  labs(x = paste0("CP1 (", var_exp[1], "%)"), y = paste0("CP2 (", var_exp[2], "%)"),
       colour = "Clúster", shape = "Año", title = "Comuna-año en el espacio CLR (PCA)") + theme_minimal()
ggsave(file.path(OUT, "02_pca_biplot.png"), p_pca, width = 9, height = 7, dpi = 300, bg = "white")

# Composición por año
cat("\nComunas por clúster y año:\n"); print(table(ancha$anio, ancha$cluster))

# Transiciones
wide_cl <- ancha %>% select(comuna, anio, cluster) %>%
  pivot_wider(names_from = anio, values_from = cluster, names_prefix = "a")
trans <- bind_rows(lapply(seq_len(length(ANIOS) - 1), function(i) {
  a0 <- paste0("a", ANIOS[i]); a1 <- paste0("a", ANIOS[i + 1])
  wide_cl %>% filter(!is.na(.data[[a0]]), !is.na(.data[[a1]])) %>%
    count(desde = .data[[a0]], hacia = .data[[a1]]) %>%
    mutate(intervalo = paste0(ANIOS[i], "–", ANIOS[i + 1]))
}))
cat("\nMatrices de transición:\n")
for (iv in unique(trans$intervalo)) {
  cat("\n", iv, "\n"); print(xtabs(n ~ desde + hacia, trans %>% filter(intervalo == iv)))
}
cols_a <- paste0("a", ANIOS)
wide_cl$cambios <- apply(wide_cl[, cols_a], 1, function(z) { z <- z[!is.na(z)]; sum(z[-1] != z[-length(z)]) })
cat("\nN° de cambios por comuna:\n"); print(table(wide_cl$cambios))

alluv <- ancha %>% select(comuna, anio, cluster) %>% mutate(cluster = factor(cluster))
p_alluv <- ggplot(alluv, aes(x = factor(anio), stratum = cluster, alluvium = comuna, fill = cluster)) +
  geom_flow(alpha = .5) + geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(count)), size = 3) +
  labs(x = NULL, y = "N° de comunas", fill = "Clúster", title = "Transiciones entre tipologías") + theme_minimal()
ggsave(file.path(OUT, "03_aluvial.png"), p_alluv, width = 9, height = 6, dpi = 300, bg = "white")

# Mapas
geo <- st_read("comunasloslagoslosrios.gpkg", quiet = TRUE) |> st_zm() |> st_make_valid()
faltan <- setdiff(unique(ancha$comuna), geo$Comuna)
if (length(faltan)) warning("Comunas sin geometría: ", paste(faltan, collapse = ", "))
mapa <- geo %>% select(Comuna) %>%
  inner_join(ancha %>% select(Comuna = comuna, anio, cluster), by = "Comuna", relationship = "many-to-many") %>%
  mutate(cluster = factor(cluster))
p_mapa <- ggplot(mapa) + geom_sf(aes(fill = cluster), colour = "grey30", linewidth = .1) +
  facet_wrap(~ anio, nrow = 1) + labs(fill = "Clúster", title = "Tipología comunal por año") +
  theme_void() + theme(legend.position = "bottom")
ggsave(file.path(OUT, "04_mapas.png"), p_mapa, width = 14, height = 6, dpi = 300, bg = "white")

write_xlsx(list(
  asignacion   = ancha %>% select(comuna, region, anio, total, cluster, silhouette),
  perfil       = perfil,
  transiciones = trans,
  cambios      = wide_cl,
  pca_varianza = tibble(cp = seq_along(var_exp), var_exp = var_exp)
), file.path(OUT, "clusters_v2.xlsx"))
saveRDS(list(ancha = ancha, CLR = CLR, hc = hc, diag_k = diag_k), file.path(OUT, "clusters_v2.rds"))
