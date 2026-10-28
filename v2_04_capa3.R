# ============================================================
# v2_04_capa3.R — Capa 3: inferencia sobre el shift-share
# ------------------------------------------------------------
# (a) Regresión shift-share (Patterson, 1991; Knudsen, 2000; Blien & Wolf, 2002)
#     Unidad: comuna x rubro x periodo. Variable dependiente: RS dinámico de la
#     celda en % de su empleo inicial (crecimiento local menos crecimiento
#     nacional del rubro). Ponderada por empleo inicial; errores agrupados por comuna.
#       M1: RS ~ tipología x periodo                  (RS medio por tipo, con IC)
#       M2: M1 + efectos fijos rubro x periodo        (neto de shocks sectoriales
#           del área, p.ej. la crisis salmonera, que afectan a todos los tipos)
# (b) ¿El RS anticipa el cambio de tipología?
#     RS del periodo p -> transición en el intervalo siguiente (p.ej. RS 2006-10
#     -> cambio 2010-15). El RS contemporáneo se reporta solo como contraste: su
#     asociación es en parte mecánica (el mismo cambio de empleo define ambos).
# ============================================================
source("v2_00_datos.R")
suppressPackageStartupMessages({ library(sandwich); library(lmtest) })

OUT <- file.path(RES, "Capa3")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ss  <- readRDS(file.path(RES, "ShiftShare/shiftshare_v2.rds"))
TIPOS <- levels(ss$cap1$tipo)
REF   <- "Intermedia agro-comercial"   # tipo más numeroso: categoría de referencia en M2

# ---------------- (a) Regresión shift-share ----------------
celdas <- ss$pasos %>%
  group_by(periodo, comuna, tipo, rubro8) %>%
  summarise(RS = sum(RS), E0 = E[anio == min(anio)], .groups = "drop") %>%
  mutate(tipo = factor(tipo, TIPOS), periodo = factor(periodo),
         rs_pct = 100 * RS / E0)
sin_base <- celdas %>% filter(E0 == 0)
cat("Celdas con empleo inicial 0 (excluidas de la regresión):", nrow(sin_base),
    "| su RS equivale al", round(100 * sum(sin_base$RS) / sum(celdas$E0), 2), "% del empleo inicial total\n")
celdas <- celdas %>% filter(E0 > 0)

ajustar <- function(formula) {
  m <- lm(formula, data = celdas, weights = E0)
  list(modelo = m, vc = vcovCL(m, cluster = ~ comuna, type = "HC1"))
}
tabla_coef <- function(aj, patron) {
  ct <- coeftest(aj$modelo, vcov. = aj$vc)
  ci <- coefci(aj$modelo, vcov. = aj$vc)
  keep <- grepl(patron, rownames(ct))
  tibble(termino = rownames(ct)[keep], estimado = ct[keep, 1], ee = ct[keep, 2],
         p = ct[keep, 4], ic_inf = ci[keep, 1], ic_sup = ci[keep, 2])
}
wald <- function(aj, idx) {   # H0: todos los coeficientes idx = 0
  b <- coef(aj$modelo)[idx]; V <- aj$vc[idx, idx]
  W <- as.numeric(t(b) %*% solve(V, b)); df <- length(b)
  # ponytail: chi-cuadrado asintótico con 42 clústeres; bootstrap wild si un revisor lo pide
  c(W = W, df = df, p = pchisq(W, df, lower.tail = FALSE))
}

# M1: RS medio ponderado por tipo y periodo (reproduce la capa 1, ahora con IC)
m1 <- ajustar(rs_pct ~ 0 + periodo:tipo)
coef_m1 <- tabla_coef(m1, "periodo") %>%
  mutate(periodo = sub("periodo(.*):tipo.*", "\\1", termino),
         tipo = factor(sub(".*:tipo", "", termino), TIPOS))

# M2: diferencias entre tipos netas de shocks rubro x periodo del área
# Dummies explícitas (tipo != REF) x periodo: evita la parametrización ambigua de
# interacciones sin efectos principales en lm().
OTROS <- setdiff(TIPOS, REF)
for (i in seq_along(OTROS)) celdas[[paste0("D", i)]] <- as.numeric(celdas$tipo == OTROS[i])
f_m2 <- as.formula(paste("rs_pct ~ 0 + periodo:rubro8 +",
                         paste0("periodo:D", seq_along(OTROS), collapse = " + ")))
m2 <- ajustar(f_m2)
na_d <- names(which(is.na(coef(m2$modelo)[grepl(":D", names(coef(m2$modelo)))])))
if (length(na_d)) message("Términos no estimables (tipo sin comunas en el periodo): ", paste(na_d, collapse = ", "))
coef_m2 <- tabla_coef(m2, ":D") %>%
  mutate(periodo = sub("periodo(.*):D.*", "\\1", termino),
         tipo = factor(OTROS[as.integer(sub(".*:D", "", termino))], TIPOS))
wald_m2 <- bind_rows(lapply(levels(celdas$periodo), function(per) {
  idx <- grep(paste0("^periodo", per, ":D"), names(coef(m2$modelo)))
  idx <- idx[!is.na(coef(m2$modelo)[idx])]
  tibble::as_tibble_row(wald(m2, idx)) %>% mutate(periodo = per, .before = 1)
}))
cat("\nM2 — ¿difiere el RS entre tipologías (neto de shocks sectoriales)? Wald por periodo:\n")
print(wald_m2 %>% mutate(across(where(is.double), ~ round(.x, 3))))
cat("\nM2 — Diferencia de RS vs", REF, "(p.p.), con IC 95% agrupado por comuna:\n")
print(as.data.frame(coef_m2 %>% select(periodo, tipo, estimado, ic_inf, ic_sup, p) %>%
                      mutate(across(where(is.double), ~ round(.x, 2)))))

# ---------------- (b) ¿El RS anticipa el cambio de tipología? ----------------
# Jerarquía funcional para definir "ascenso" (de menor a mayor diversificación/centralidad)
RANGO <- c("Primario especializado" = 1, "Rural dependiente de servicios públicos" = 1,
           "Periferia cordillerana-turística" = 1, "Intermedia agro-comercial" = 2,
           "Centro subregional" = 3, "Centro regional" = 4)
tip <- readRDS(file.path(RES, "Tipologia_funcional/tipologia_funcional.rds"))$ind %>%
  transmute(comuna, anio, tipo = as.character(tipo))
cortes <- sort(unique(tip$anio))
trans <- bind_rows(lapply(seq_len(length(cortes) - 1), function(i) {
  a <- tip %>% filter(anio == cortes[i]) %>% select(comuna, t0 = tipo)
  b <- tip %>% filter(anio == cortes[i + 1]) %>% select(comuna, t1 = tipo)
  inner_join(a, b, by = "comuna") %>%
    mutate(intervalo = paste0(cortes[i], "–", cortes[i + 1]),
           cambia = as.integer(t0 != t1),
           asciende = as.integer(RANGO[t1] > RANGO[t0]),
           rango0 = RANGO[t0])
}))

rs_com <- ss$cap2 %>% transmute(periodo, comuna, RS_pct, IS_pct, E0)
siguiente <- setNames(c("2010–2015", "2015–2020", "2020–2024", NA), c("2006–2010", "2010–2015", "2015–2020", "2020–2024"))
datos_b <- bind_rows(
  rs_com %>% mutate(intervalo = siguiente[periodo], especificacion = "Rezagado (RS p -> cambio p+1)"),
  rs_com %>% mutate(intervalo = periodo, especificacion = "Contemporáneo (mecánico)")
) %>% filter(!is.na(intervalo)) %>%
  inner_join(trans, by = c("comuna", "intervalo"))

modelos_b <- datos_b %>% group_by(especificacion) %>% group_modify(function(d, ...) {
  bind_rows(lapply(c("cambia", "asciende"), function(y) {
    d2 <- if (y == "asciende") filter(d, rango0 < 4) else d
    f  <- as.formula(paste(y, "~ I(RS_pct / 10) + I(IS_pct / 10) + factor(rango0) + factor(intervalo)"))
    m  <- glm(f, family = binomial, data = d2)
    ct <- coeftest(m, vcov. = vcovCL(m, cluster = ~ comuna, type = "HC0"))
    tibble(resultado = y, n = nrow(d2), eventos = sum(d2[[y]]),
           OR_RS_10pp = exp(ct["I(RS_pct/10)", 1]), p_RS = ct["I(RS_pct/10)", 4],
           OR_IS_10pp = exp(ct["I(IS_pct/10)", 1]), p_IS = ct["I(IS_pct/10)", 4])
  }))
}) %>% ungroup()
cat("\n(b) Logit: odds ratio por +10 p.p. de RS/IS comunal (errores agrupados por comuna):\n")
print(as.data.frame(modelos_b %>% mutate(across(where(is.double), ~ round(.x, 3)))))

# ---------------- Figuras ----------------
PALETA <- setNames(c("#7a0177", "#dd3497", "#fe9929", "#4d9221", "#4575b4", "#74add1"), TIPOS)
p_m1 <- ggplot(coef_m1, aes(x = estimado, y = tipo, colour = tipo)) +
  geom_vline(xintercept = 0, linewidth = .3) +
  geom_errorbar(aes(xmin = ic_inf, xmax = ic_sup), width = .25, orientation = "y") + geom_point(size = 2) +
  facet_wrap(~ periodo, ncol = 2) + scale_y_discrete(limits = rev) +
  scale_colour_manual(values = PALETA, guide = "none") +
  labs(x = "RS (% del empleo inicial), IC 95% agrupado por comuna", y = NULL,
       title = "Efecto regional por tipología (regresión shift-share ponderada)") + theme_minimal()
ggsave(file.path(OUT, "01_RS_tipologia_IC.png"), p_m1, width = 11, height = 7, dpi = 300, bg = "white")

p_m2 <- ggplot(coef_m2, aes(x = estimado, y = tipo, colour = tipo)) +
  geom_vline(xintercept = 0, linewidth = .3) +
  geom_errorbar(aes(xmin = ic_inf, xmax = ic_sup), width = .25, orientation = "y") + geom_point(size = 2) +
  geom_text(data = wald_m2, aes(x = Inf, y = -Inf, label = sprintf("Wald p = %.3f", p)),
            inherit.aes = FALSE, hjust = 1.1, vjust = -.6, size = 3) +
  facet_wrap(~ periodo, ncol = 2) + scale_y_discrete(limits = rev) +
  scale_colour_manual(values = PALETA, guide = "none") +
  labs(x = paste0("Diferencia de RS vs ", REF, " (p.p.), neta de efectos rubro × periodo"), y = NULL,
       title = "¿Difiere la competitividad local entre tipologías?") + theme_minimal()
ggsave(file.path(OUT, "02_RS_diferencias_M2.png"), p_m2, width = 11, height = 7, dpi = 300, bg = "white")


# ---------------- (c) Robustez: ¿tipología o tamaño? ----------------
# Los centros regionales son las 3 comunas más pobladas: se evalúa si el RS comunal
# se asocia al tamaño (log de la población al inicio del periodo) y si la diferencia
# de los centros regionales se sostiene al controlar por tamaño.
tam <- ss$cap2 %>% mutate(anio = as.integer(substr(periodo, 1, 4))) %>%
  left_join(pob, by = c("comuna", "anio")) %>% mutate(CR = as.numeric(tipo == "Centro regional"))
fila <- function(m, term) {
  ct <- coeftest(m, vcov. = vcovHC(m, type = "HC1"))
  if (!term %in% rownames(ct)) return(c(NA, NA))
  ct[term, c(1, 4)]
}
robustez_tamano <- bind_rows(lapply(unique(tam$periodo), function(per) {
  d <- filter(tam, periodo == per)
  m_s <- lm(RS_pct ~ log(pob), d); m_w <- lm(RS_pct ~ log(pob), d, weights = E0)
  m_cr <- lm(RS_pct ~ CR, d, weights = E0); m_crp <- lm(RS_pct ~ CR + log(pob), d, weights = E0)
  sp <- cor.test(d$RS_pct, d$pob, method = "spearman", exact = FALSE)
  tibble(periodo = per, n = nrow(d),
         b_logpob = fila(m_s, "log(pob)")[1], p_logpob = fila(m_s, "log(pob)")[2],
         b_logpob_pond = fila(m_w, "log(pob)")[1], p_logpob_pond = fila(m_w, "log(pob)")[2],
         rho_spearman = unname(sp$estimate), p_spearman = sp$p.value, r2 = summary(m_s)$r.squared,
         b_CR = fila(m_cr, "CR")[1], p_CR = fila(m_cr, "CR")[2],
         b_CR_ctrl = fila(m_crp, "CR")[1], p_CR_ctrl = fila(m_crp, "CR")[2],
         b_logpob_ctrl = fila(m_crp, "log(pob)")[1], p_logpob_ctrl = fila(m_crp, "log(pob)")[2])
}))
cat("\n(c) RS comunal ~ log(población inicial) y centros regionales controlando por tamaño:\n")
print(as.data.frame(robustez_tamano %>% mutate(across(where(is.double), ~ round(.x, 3)))))

# ---------------- (d) Inferencia robusta con pocos grupos ----------------
# Los errores agrupados son poco confiables cuando el grupo de interés tiene muy pocas
# comunas (3 centros regionales). Se usa inferencia por aleatorización a nivel de comuna
# (MacKinnon & Webb, 2020). Se reemplaza además el logit de (b) —pocos eventos por
# parámetro— por pruebas de permutación, y se ajustan los p por Holm.
set.seed(2026)
N_PERM <- 9999
per_split <- split(celdas, celdas$periodo)

# (d1) Centro regional vs REF: enumeración exacta de todos los tríos posibles dentro de CR ∪ REF
coef_cr <- function(d, cr) {
  g <- ifelse(d$comuna %in% cr, "CR",
              ifelse(d$tipo %in% c("Centro regional", REF), "REF", as.character(d$tipo)))
  g <- factor(g, levels = c("REF", setdiff(unique(g), "REF")))
  unname(coef(lm(rs_pct ~ 0 + rubro8 + g, data = d, weights = E0))["gCR"])
}
ri_cr <- bind_rows(lapply(names(per_split), function(per) {
  d <- per_split[[per]]
  cr_obs <- unique(d$comuna[d$tipo == "Centro regional"])
  pool   <- unique(d$comuna[d$tipo %in% c("Centro regional", REF)])
  obs    <- coef_cr(d, cr_obs)
  combos <- combn(pool, length(cr_obs), simplify = FALSE)
  dist   <- vapply(combos, function(cr) coef_cr(d, cr), numeric(1))
  # Distribución asimétrica (tríos de comunas pequeñas generan valores extremos positivos):
  # se usa el p de colas iguales en lugar del basado en |coeficiente|.
  tibble(periodo = per, estimado = obs, n_combinaciones = length(combos),
         percentil_obs = mean(dist <= obs + 1e-9),
         p_colas_iguales = min(1, 2 * min(mean(dist <= obs + 1e-9), mean(dist >= obs - 1e-9))),
         p_valor_absoluto = mean(abs(dist) >= abs(obs) - 1e-9))
}))

# (d2) ¿Difieren las tipologías? F ponderado de tipología (sobre efectos fijos de rubro),
# con etiquetas de tipología permutadas entre comunas (tamaños de grupo fijos)
f_tipo <- function(d, tipo_com) {
  d$g <- factor(tipo_com[d$comuna])
  if (nlevels(d$g) < 2) return(NA_real_)
  anova(lm(rs_pct ~ 0 + rubro8, d, weights = E0), lm(rs_pct ~ 0 + rubro8 + g, d, weights = E0))$F[2]
}
ri_conjunta <- bind_rows(lapply(names(per_split), function(per) {
  d <- per_split[[per]]
  tc <- d %>% distinct(comuna, tipo) %>% { setNames(as.character(.$tipo), .$comuna) }
  obs  <- f_tipo(d, tc)
  perm <- replicate(1999, f_tipo(d, setNames(sample(tc), names(tc))))
  tibble(periodo = per, F_obs = obs, p_aleatorizacion = (1 + sum(perm >= obs)) / (1 + length(perm)))
})) %>% mutate(p_holm = p.adjust(p_aleatorizacion, "holm"))

# (d3) ¿El RS/IS previo difiere entre comunas que cambian (o ascienden) y las que no?
# Diferencia de medias; el resultado se permuta dentro de cada intervalo.
perm_dif <- function(x, y, estrato) {
  obs <- mean(x[y == 1]) - mean(x[y == 0])
  perm <- replicate(N_PERM, { yp <- ave(y, estrato, FUN = sample); mean(x[yp == 1]) - mean(x[yp == 0]) })
  c(diferencia = obs, p = (1 + sum(abs(perm) >= abs(obs))) / (1 + N_PERM))
}
ri_transiciones <- datos_b %>% group_by(especificacion) %>% group_modify(function(d, ...) {
  bind_rows(lapply(c("cambia", "asciende"), function(y) {
    d2 <- if (y == "asciende") filter(d, rango0 < 4) else d
    bind_rows(lapply(c("RS_pct", "IS_pct"), function(v) {
      r <- perm_dif(d2[[v]], d2[[y]], d2$intervalo)
      tibble(resultado = y, variable = v, n = nrow(d2), eventos = sum(d2[[y]]),
             diferencia_pp = r[["diferencia"]], p_permutacion = r[["p"]])
    }))
  }))
}) %>% ungroup() %>% group_by(especificacion) %>% mutate(p_holm = p.adjust(p_permutacion, "holm")) %>% ungroup()

# (d4) Holm sobre los η² de la capa 2 (4 periodos)
eta_holm <- ss$cap2_var %>% mutate(p_holm = p.adjust(p_permutacion, "holm"),
                                   p_holm_ponderado = p.adjust(p_perm_ponderado, "holm"))

cat("\n(d1) Centro regional vs", REF, "— inferencia por aleatorización:\n"); print(as.data.frame(ri_cr %>% mutate(across(where(is.double), ~ round(.x, 3)))))
cat("\n(d2) Diferencias conjuntas entre tipologías — permutación:\n"); print(as.data.frame(ri_conjunta %>% mutate(across(where(is.double), ~ round(.x, 3)))))
cat("\n(d3) Transiciones — permutación:\n"); print(as.data.frame(ri_transiciones %>% mutate(across(where(is.double), ~ round(.x, 3)))))
cat("\n(d4) η² con Holm:\n"); print(as.data.frame(eta_holm %>% select(periodo, eta2_tipologia, p_permutacion, p_holm, eta2_ponderado, p_perm_ponderado, p_holm_ponderado) %>% mutate(across(where(is.double), ~ round(.x, 3)))))

write_xlsx(list(M1_RS_por_tipo = coef_m1, M2_diferencias = coef_m2, M2_wald = wald_m2,
                b_logit = modelos_b, b_datos = datos_b, transiciones = trans,
                c_robustez_tamano = robustez_tamano,
                d1_CR_aleatorizacion = ri_cr, d2_conjunta_aleatorizacion = ri_conjunta,
                d3_transiciones_permutacion = ri_transiciones, d4_eta2_holm = eta_holm),
           file.path(OUT, "capa3.xlsx"))
