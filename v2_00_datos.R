# ============================================================
# v2_00_datos.R — Preparación común de datos (SII + INE) para el re-análisis
# Genera: base, ancha (comuna-año x 8 grandes rubros), X, P, CLR, excluidas, pob
# ============================================================
Sys.setlocale("LC_ALL", "en_US.UTF-8")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(cluster)
  library(fpc); library(sf); library(tmap); library(ggalluvial); library(writexl)
})

ANIOS     <- c(2006, 2010, 2015, 2020, 2024)
REGIONES  <- c("Región de Los Ríos", "Región de Los Lagos")
MIN_TRAB  <- 100   # comuna-año con menos trabajadores se excluye (dato no confiable)
ZERO_REPL <- 0.5   # trabajadores imputados a celdas cero antes de CLR
set.seed(2026)

# B, D y E se fusionan con C: por sí solos tienen demasiados ceros y, tras CLR,
# dominaban la distancia (análisis de sensibilidad: Jaccard mín. 0,56 -> 0,66 con k=3).
# O, P y Q se fusionan: el SII registra al personal municipal de educación y salud
# bajo O o bajo P/Q según el RUT declarante (p.ej. Osorno tiene O = 0 en todos los
# años con 47.727 trabajadores; Valdivia pasa de O = 274 a 3.763 entre 2015 y 2020).
# Separarlos generaba clústeres definidos por la forma de registro, no por la estructura.
# R y S se suman al mismo bloque (agregado CIIU rev.3 "servicios comunales, sociales
# y personales"): S muestra quiebres de registro (área: 10.480 en 2005 -> 3.537 en 2008;
# Queilén 156 -> 0 entre 2010 y 2015) que generaban un clúster espurio 2006-2015.
GRUPOS <- c(A = "Primario", B = "Industria y utilities", C = "Industria y utilities",
            D = "Industria y utilities", E = "Industria y utilities", F = "Construcción",
            G = "Comercio", H = "Transporte", I = "Alojamiento y comidas",
            J = "Serv. empresariales", K = "Serv. empresariales", L = "Serv. empresariales",
            M = "Serv. empresariales", N = "Serv. empresariales",
            O = "Serv. públicos, sociales y personales", P = "Serv. públicos, sociales y personales",
            Q = "Serv. públicos, sociales y personales", R = "Serv. públicos, sociales y personales",
            S = "Serv. públicos, sociales y personales")

# ---------------- Datos ----------------
lin  <- iconv(readLines("sii_trabajadores_comuna_ciiu_2005_2024.csv", encoding = "bytes"), from = "MACROMAN", to = "UTF-8")
base <- read.csv(text = lin, sep = ";", check.names = FALSE)
names(base)[1] <- "anio"
stopifnot(is.numeric(base$trabajadores), !anyNA(base$trabajadores))

# ---------------- Saltos de registro (sensibilidad) ----------------
# El SII asigna el empleo al domicilio tributario: cambios de domicilio o de registro
# producen saltos abruptos en una comuna-rubro (p.ej. Enseñanza en Llanquihue: 178 -> 2 096
# entre 2020 y 2021). Regla explícita y simétrica, por comuna x sección CIIU x paso anual:
#   |ΔE| >= 300, |ΔE| >= 5 % del empleo comunal del año inicial y el rubro al menos se
#   duplica o se reduce a la mitad (|log((E1+1)/(E0+1))| >= log 2).
# Además, el paso 2017->2018 de O, P y Q se trata como quiebre de registro (salto nacional).
# Con AJUSTE_SALTOS=1 los saltos marcados se restan de forma acumulada del nivel posterior
# (se elimina el salto, se conserva la dinámica restante) y los resultados se escriben en
# Resultados_v2_ajustado/. Sin la variable, solo se marcan (versión principal).
SALTO_MIN <- as.numeric(Sys.getenv("SALTO_MIN", "300")); SALTO_PCT <- 0.05   # umbral parametrizable (sensibilidad: 200 y 500)
AJUSTE <- Sys.getenv("AJUSTE_SALTOS", "0") == "1"
RES    <- if (!AJUSTE) "Resultados_v2" else if (SALTO_MIN == 300) "Resultados_v2_ajustado" else paste0("Resultados_v2_ajustado_", SALTO_MIN)
LETRAS <- c(LETTERS[1:19])   # A-S (T, U y "Sin información" no se usan)

celdas_area <- base %>%
  filter(region %in% c("Región de Los Ríos", "Región de Los Lagos"),
         !grepl("^Sin informaci", rubro, ignore.case = TRUE)) %>%
  mutate(letra = substr(rubro, 1, 1)) %>% filter(letra %in% LETRAS) %>%
  group_by(comuna, region, letra, anio) %>% summarise(E = sum(trabajadores), .groups = "drop") %>%
  complete(nesting(comuna, region), letra, anio = min(base$anio):max(base$anio), fill = list(E = 0)) %>%
  group_by(comuna, anio) %>% mutate(tot = sum(E)) %>%
  group_by(comuna, letra) %>% arrange(anio, .by_group = TRUE) %>%
  mutate(dE = lead(E) - E,
         salto = !is.na(dE) & abs(dE) >= SALTO_MIN & abs(dE) >= SALTO_PCT * tot &
                 abs(log((lead(E) + 1) / (E + 1))) >= log(2),
         quiebre_2018 = anio == 2017 & letra %in% c("O", "P", "Q"),
         quitar = ifelse(salto | quiebre_2018, dE, 0),
         E_ajust = pmax(E - lag(cumsum(coalesce(quitar, 0)), default = 0), 0)) %>%
  ungroup()
saltos <- celdas_area %>% filter(salto) %>%
  transmute(comuna, rubro = letra, paso = paste0(anio, "->", anio + 1), E_inicial = E, dE,
            pct_empleo_comunal = round(100 * dE / tot, 1))

if (AJUSTE) {
  message("AJUSTE_SALTOS=1: se restan ", nrow(saltos), " saltos marcados y el quiebre 2017-2018 de O, P y Q")
  base <- bind_rows(
    base %>% filter(!(region %in% c("Región de Los Ríos", "Región de Los Lagos") &
                        substr(rubro, 1, 1) %in% LETRAS & !grepl("^Sin informaci", rubro, ignore.case = TRUE))),
    celdas_area %>% transmute(anio, comuna, region, rubro = letra, trabajadores = E_ajust)
  )
}

larga <- base %>%
  filter(region %in% REGIONES, anio %in% ANIOS, !grepl("^Sin informaci", rubro, ignore.case = TRUE)) %>%
  mutate(letra = substr(rubro, 1, 1)) %>%
  filter(letra %in% names(GRUPOS)) %>%                      # excluye T y U (marginales)
  mutate(grupo = GRUPOS[letra]) %>%
  group_by(anio, comuna, region, grupo) %>%
  summarise(t = sum(trabajadores), .groups = "drop") %>%
  complete(nesting(anio, comuna, region), grupo, fill = list(t = 0))

ancha <- larga %>% pivot_wider(names_from = grupo, values_from = t)
ancha$total <- rowSums(ancha[, unique(GRUPOS)])
excluidas <- ancha %>% filter(total < MIN_TRAB) %>% select(anio, comuna, total)
ancha <- ancha %>% filter(total >= MIN_TRAB) %>% mutate(id = paste(comuna, anio, sep = "_"))

X <- as.matrix(ancha[, unique(GRUPOS)]); rownames(X) <- ancha$id
P <- X / rowSums(X)                                          # proporciones (para describir)

# Reemplazo multiplicativo de ceros + CLR
Z <- t(apply(X, 1, function(v) {
  p <- v / sum(v); d <- ZERO_REPL / sum(v); nz <- p == 0
  p[nz] <- d; p[!nz] <- p[!nz] * (1 - d * sum(nz)); p
}))
CLR <- log(Z) - rowMeans(log(Z))

# Población comunal (INE, estimaciones y proyecciones 2002-2035)
pob <- readxl::read_excel("estimaciones-y-proyecciones-2002-2035-comunas.xlsx",
                          sheet = "Est. y Proy. de Pob. Comunal") %>%
  select(comuna = `Nombre Comuna`, starts_with("Poblacion ")) %>%
  pivot_longer(-comuna, names_to = "anio", names_prefix = "Poblacion ", values_to = "pob") %>%
  mutate(anio = as.integer(anio)) %>%
  filter(anio %in% ANIOS, comuna %in% ancha$comuna) %>%
  group_by(comuna, anio) %>% summarise(pob = sum(pob), .groups = "drop")
stopifnot(nrow(pob) == length(unique(ancha$comuna)) * length(ANIOS))
