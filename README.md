# Tipologías funcionales y shift-share dinámico del empleo formal comunal

Código en R y datos para reproducir los resultados, cuadros, figuras y pruebas de sensibilidad de un análisis del empleo formal en las 42 comunas de las regiones de Los Ríos y Los Lagos (Chile), 2006-2024.

## Datos

| Archivo | Contenido | Fuente |
|---|---|---|
| `sii_trabajadores_comuna_ciiu_2005_2024.csv` | Trabajadores dependientes ponderados por meses trabajados, por comuna, sección CIIU y año (2005-2024), todas las comunas del país | Servicio de Impuestos Internos, Estadísticas de empresas por comuna y actividad económica (https://www.sii.cl/sobre_el_sii/estadisticas_de_empresas.html) |
| `estimaciones-y-proyecciones-2002-2035-comunas.xlsx` | Población comunal anual | Instituto Nacional de Estadísticas, Estimaciones y proyecciones de la población de Chile 2002-2035 a nivel comunal |
| `censo2024_poblacion_comunal.csv` | Población censada 2024 de las 42 comunas (solo para la Figura 1) | Instituto Nacional de Estadísticas, Censo 2024, cuadro D1, hoja 2 (https://censo2024.ine.gob.cl/estadisticas/) |
| `comunasloslagoslosrios.gpkg` | Límites de las 42 comunas (capa `comunas`) | Biblioteca del Congreso Nacional de Chile, División comunal (https://www.bcn.cl/siit/mapas_vectoriales) |

## Scripts

| Script | Resultado |
|---|---|
| `v2_00_datos.R` | Lectura de datos, agregación en 8 grandes sectores y detección de saltos de registro (lo llaman los demás scripts) |
| `v2_01_clusters.R` | Tipología solo por composición (sensibilidad) |
| `v2_02_tipologia_funcional.R` | Tipología funcional: análisis factorial múltiple y agrupamiento de consenso |
| `v2_03_shiftshare.R` | Shift-share dinámico por tipología |
| `v2_04_capa3.R` | Inferencia: aleatorización, permutaciones, Holm y prueba de tamaño |
| `v2_05_figuras_manuscrito.R` | Figuras 1-7 y material suplementario S1-S3 |
| `v2_06_sensibilidad_saltos.R` | Comparación entre la serie sin ajustar y la ajustada |
| `v2_07_sensibilidad_k.R` | Sensibilidad al número de grupos (3 a 7) |
| `v2_08_cambio_estructural.R` | Índice de cambio estructural (NAV) y su descomposición |
| `v2_09_nacional.R` | Gradiente nacional del efecto regional según tamaño (Chile sin la Región Metropolitana) |
| `v2_10_relacion.R` | Densidad de relación y entrada en nuevos rubros |
| `v2_11_perfiles_tipologias.R` | Perfiles de las tipologías (Cuadro 3) y material suplementario S2 |
| `v2_12_sensibilidad_umbrales.R` | Sensibilidad al umbral de la regla de saltos (200, 300 y 500 trabajadores) |

## Ejecución

Desde esta carpeta, en este orden:

```sh
# 1. Serie principal y series ajustadas por saltos (umbral 300, 200 y 500)
for s in v2_02 v2_03 v2_04 v2_08 v2_10; do
  Rscript ${s}_*.R
  AJUSTE_SALTOS=1 Rscript ${s}_*.R
  AJUSTE_SALTOS=1 SALTO_MIN=200 Rscript ${s}_*.R
  AJUSTE_SALTOS=1 SALTO_MIN=500 Rscript ${s}_*.R
done

# 2. Análisis complementarios, figuras y sensibilidad (serie principal)
for s in v2_05 v2_06 v2_07 v2_09 v2_11 v2_12; do Rscript ${s}_*.R; done

# Variante de la tipología sin variables basadas en el empleo total
TIPO_SIN_DOMICILIO=1 Rscript v2_02_tipologia_funcional.R
```

Resultados: `Resultados_v2/` (serie principal), `Resultados_v2_ajustado/` (umbral 300), `Resultados_v2_ajustado_200/` y `Resultados_v2_ajustado_500/`. Figuras: `manuscrito/figuras/`.

## Requisitos

R 4.4.3 con los paquetes (versiones usadas): dplyr 1.2.0, tidyr 1.3.1, ggplot2 4.0.2, cluster 2.1.8, fpc 2.2.14, sf 1.0.20, tmap 4.1, ggalluvial 0.12.5, writexl 1.5.4, readxl 1.4.5, FactoMineR 2.11, mclust 6.1.2, sandwich 3.1.1, lmtest 0.9.40, purrr 1.0.4, ggspatial 1.1.10, patchwork 1.3.0, ggrepel 0.9.6, ggnewscale 0.5.2, rnaturalearth 1.0.1, scales 1.4.0, tibble.
