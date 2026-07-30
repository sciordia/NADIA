# NADIA 0.99.0

Primera versión preparada para su envío a Bioconductor. El proyecto pasa de ser
una colección de módulos que se consumían con `source()` a un paquete de R
instalable, sin cambios en los resultados que produce el pipeline.

## Nuevo

* Estructura de paquete: `DESCRIPTION`, `NAMESPACE` y páginas de ayuda generadas
  con roxygen2. El flujo de trabajo pasa a ser `library(NADIA)`.
* 102 funciones exportadas, que cubren el pipeline completo: preprocesado de
  Spectronaut/DIA-NN y de Proteome Discoverer (TMT y LFQ), normalización (13
  métodos), corrección de lote, imputación (19 métodos), expresión diferencial
  con limma o limpa, métricas y benchmarking, y visualización interactiva y
  estática.
* Dataset de ejemplo `nadia_dia` y reports recortados en `inst/extdata/`, con su
  procedencia documentada en `inst/scripts/`.

## Cambios

* Las dependencias se declaran en `Imports:` y `Suggests:`; ya no se llama a
  `library()` desde el código. Los paquetes que solo hacen falta para un método
  concreto se comprueban en el punto de uso.
* `process_proteomics()` y `pattern_profiler_analysis()` ya no escriben en disco
  a menos que se les indique una ruta de salida: `export_dir` y `output_file`
  tienen ahora `NULL` por defecto. `pattern_profiler_analysis()` devuelve
  `long_output` en su lista de resultados.
* Las funciones que usan `set.seed()` restauran el generador de números
  aleatorios al salir, de modo que no alteran la reproducibilidad del código que
  se ejecute después. `.nm_hopkins()` acepta la semilla como argumento.

## Correcciones

* Ocho definiciones duplicadas de helpers de color y filtrado (`hex_to_rgba`,
  `darken_hex`, `normalize_hex`, `get_feature_ids` y otros) repartidas por seis
  módulos se unifican en una sola. Al vivir todas en el entorno global, la copia
  activa dependía del orden de carga.
* `get_feature_ids()` acepta ahora `mode = "target"` y `mode = "specific"` como
  sinónimos. Antes, cargar `Heatmap_tidyHeatmap.R` después de
  `PCA_Highcharts_Final.R` hacía que `pca_highchart_list(modes = "specific")`
  abortara.
* El filtrado por `mode = "any"` ya no puede devolver `FeatureID` `NA` cuando la
  columna `sig_any` contiene valores ausentes.
* `summary_list_widget()` funciona con metadatos de LFQ y TMT, que carecen de las
  columnas de recuento de identificaciones propias de Spectronaut.
