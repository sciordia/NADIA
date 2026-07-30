# NADIA

**Missing Value-Aware DIA Proteomics Analysis** — pipeline en R para el análisis
de expresión diferencial de proteínas, con especial atención al tratamiento de
los valores ausentes (**NA**) característicos de la adquisición independiente de
datos (**DIA**).

## Qué cubre

Del report crudo a la figura interactiva:

1. **Preprocesado** — Spectronaut/DIA-NN, TMT (Proteome Discoverer) y LFQ
   (Proteome Discoverer); todos producen el mismo objeto S3 `proteomics_data`.
2. **Procesado** — normalización (13 métodos), corrección de lote opcional
   (HarmonizR/BERT/ComBat + diagnóstico PVCA), imputación (19 métodos, incluidos
   los híbridos MAR/MNAR `combo` y `softHybrid`, y el modelo probabilístico
   `limpa`) y expresión diferencial (`limma` o `limpa`).
3. **Métricas y benchmarking** — evaluación de normalización (PCV/PMAD/PEV,
   correlación intragrupo, separación de grupos) e imputación (marco NAguideR:
   NRMSE, SOR, PSS, ACC_OI), más benchmarking con datasets *spike-in* y ranking
   OpDEA de combinaciones normalización × imputación.
4. **Visualización** — interactiva con Highcharts (boxplots, volcano, PCA,
   perfiles de clúster) y estática con ggplot2/ComplexHeatmap, además de tablas
   interactivas con reactable.

## Instalación

NADIA es un paquete de R y requiere R ≥ 4.4. Todavía no está en Bioconductor, así
que se instala desde el repositorio:

```r
# install.packages("remotes")
remotes::install_github("sciordia/NADIA")
library(NADIA)
```

Las dependencias imprescindibles (campo `Imports:`) se instalan solas. Las
**opcionales** (`Suggests:`) solo hacen falta si se usa el método que las
invoca — `mice` únicamente con `imp_method = "mice"`, `pROC` para las métricas
AUC/pAUC del benchmarking, `Mfuzz` para el Pattern Profiler. Cuando falta alguna,
la función lo indica con un mensaje explícito. Para instalarlas todas de golpe:

```r
source("install_dependencies.R")         # desde un clon del repositorio
install_nadia_deps(dry_run = TRUE)       # solo informa de lo que falta
install_nadia_deps(optional = FALSE)     # solo lo imprescindible
```

No se usa `renv`.

## Un primer ejemplo

```r
library(NADIA)

# Dataset de ejemplo ya preprocesado: 3 condiciones x 4 réplicas, 2.000 proteínas
data(nadia_dia)

res <- process_proteomics(nadia_dia,
                          norm_method = "cycloess",
                          imp_method  = "combo",
                          de_method   = "limma")
head(res$DEPs_results)

# O partiendo del report crudo
prep <- preprocess_spectronaut(
  system.file("extdata", "nadia_dia_report.tsv.gz", package = "NADIA"),
  condition_order = c("A", "B", "D"))
```

`process_proteomics()` no escribe nada en disco a menos que se le pase
`export_dir`.

## Estructura del repositorio

- `R/` — código del paquete: 21 archivos, 102 funciones exportadas.
- `man/`, `NAMESPACE` — generados con roxygen2; no editar a mano.
- `inst/extdata/` — reports recortados de Spectronaut, TMT y LFQ para los
  ejemplos; `inst/scripts/make_extdata.R` documenta cómo se obtuvieron.
- `data/` — el dataset de ejemplo `nadia_dia`.
- `data-raw/`, `results/` — datos completos y salidas de análisis reales. No
  forman parte del paquete (`.Rbuildignore`).
- `example_workflow.R` y los scripts numerados — recorridos de principio a fin
  sobre los datos completos.
- `CLAUDE.md` — descripción detallada de la arquitectura y de cada módulo.
- `CODE_REVIEW_*.md` — revisiones de código y análisis de impacto de sus
  correcciones.

## Licencia

MIT © 2025 Sergio Ciordia
