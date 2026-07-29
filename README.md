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

## Estructura

- `R/` — módulos independientes: cada uno puede usarse por separado con
  `source()`; las dependencias pesadas están protegidas con `requireNamespace()`.
- `example_workflow.R` — ejemplo completo de principio a fin.
- `data/`, `results/` — datos de entrada y salidas generadas.
- `CLAUDE.md` — descripción detallada de la arquitectura y de cada módulo.
- `CODE_REVIEW_*.md` — revisiones de código y análisis de impacto de sus
  correcciones.

Gestión de dependencias con `renv`.

## Licencia

MIT © 2025 Sergio Ciordia
