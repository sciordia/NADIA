# =============================================================================
# Boxplot Estilizado para Datos de Proteómica
# =============================================================================

# Librerías necesarias
librerias <- c("ggplot2", "RColorBrewer", "paletteer", "ggtext")
# renv::install(librerias)

library(ggplot2)
library(RColorBrewer)

# -----------------------------------------------------------------------------
# Tema personalizado para boxplots de proteómica
# -----------------------------------------------------------------------------

theme_proteomics_boxplot <- function(
    base_size = 12,
    base_family = "sans",
    strip_background_color = "#F8F9FA",
    grid_color = "#E9ECEF"
) {

  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Fondo general
      plot.background = element_rect(fill = "#FFFFFF", color = NA),
      panel.background = element_rect(fill = "#FFFFFF", color = NA),

      # Grid
      panel.grid.major.x = element_line(color = grid_color, linewidth = 0.4, linetype = "dotted"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),

      # Bordes del panel
      panel.border = element_rect(color = "#DEE2E6", fill = NA, linewidth = 0.5),

      # Ejes
      axis.line = element_blank(),
      axis.ticks = element_line(color = "#ADB5BD", linewidth = 0.3),
      axis.ticks.length = unit(3, "pt"),
      axis.text.x = element_text(color = "#495057", size = rel(0.9), margin = margin(t = 4)),
      axis.text.y = element_text(color = "#495057", size = rel(0.85), margin = margin(r = 4)),
      axis.title.x = element_text(color = "#212529", size = rel(1), face = "bold",
                                   margin = margin(t = 10)),
      axis.title.y = element_text(color = "#212529", size = rel(1), face = "bold",
                                   margin = margin(r = 10), angle = 90),

      # Título y subtítulo
      plot.title = element_text(color = "#1D3557", size = rel(1.3), face = "bold",
                                 hjust = 0, margin = margin(b = 8)),
      plot.subtitle = element_text(color = "#6C757D", size = rel(0.95),
                                    hjust = 0, margin = margin(b = 12)),
      plot.caption = element_text(color = "#ADB5BD", size = rel(0.75),
                                   hjust = 1, margin = margin(t = 10)),

      # Facetas
      strip.background = element_rect(fill = strip_background_color, color = "#DEE2E6"),
      strip.text = element_text(color = "#1D3557", size = rel(0.95), face = "bold",
                                 margin = margin(t = 6, b = 6, l = 8, r = 8)),

      # Leyenda
      legend.position = "bottom",
      legend.background = element_rect(fill = "#FFFFFF", color = NA),
      legend.key = element_rect(fill = "#FFFFFF", color = NA),
      legend.title = element_text(color = "#212529", size = rel(0.9), face = "bold"),
      legend.text = element_text(color = "#495057", size = rel(0.85)),
      legend.key.size = unit(0.9, "lines"),
      legend.spacing.x = unit(0.3, "cm"),
      legend.margin = margin(t = 10),

      # Márgenes generales
      plot.margin = margin(t = 15, r = 15, b = 10, l = 10)
    )
}


# -----------------------------------------------------------------------------
# Funciones auxiliares (fallbacks independientes de PRONE)
# -----------------------------------------------------------------------------

#' Validar y obtener assays del SummarizedExperiment
check_input_assays <- function(se, ain = NULL) {
  available <- SummarizedExperiment::assayNames(se)
  if (is.null(ain)) {
    return(available)
  }
  valid <- ain[ain %in% available]
  if (length(valid) == 0) {
    warning("Ninguno de los assays especificados existe. Disponibles: ",
            paste(available, collapse = ", "))
    return(NULL)
  }
  if (length(valid) < length(ain)) {
    missing <- setdiff(ain, available)
    warning("Assays no encontrados (ignorados): ", paste(missing, collapse = ", "))
  }
  valid
}

#' Obtener variable de color del colData
get_color_value <- function(se, color_by = NULL) {
  if (is.null(color_by)) return(NULL)
  cd_names <- names(SummarizedExperiment::colData(se))
  if (!(color_by %in% cd_names)) {
    warning("'", color_by, "' no está en colData. Columnas disponibles: ",
            paste(cd_names, collapse = ", "))
    return(NULL)
  }
  color_by
}

#' Obtener variable de etiqueta del colData
get_label_value <- function(se, label_by = NULL) {
  cd_names <- names(SummarizedExperiment::colData(se))
  if (is.null(label_by) || !(label_by %in% cd_names)) {
    # Por defecto usar nombres de columna
    return(list(TRUE, "Column"))
  }
  list(TRUE, label_by)
}


# -----------------------------------------------------------------------------
# Fallback para obtener datos del SummarizedExperiment
# -----------------------------------------------------------------------------

get_complete_dt_safe <- function(se, ain = NULL) {
  stopifnot(inherits(se, "SummarizedExperiment"))
  if (is.null(ain)) ain <- SummarizedExperiment::assayNames(se)

  cd <- as.data.frame(SummarizedExperiment::colData(se), check.names = FALSE)
  if (!("Column" %in% names(cd))) cd$Column <- rownames(cd)

  long_list <- lapply(ain, function(an) {
    if (!(an %in% SummarizedExperiment::assayNames(se))) {
      stop("Assay '", an, "' no está en 'se'.")
    }
    X <- SummarizedExperiment::assay(se, an)
    if (!is.matrix(X)) X <- as.matrix(as.data.frame(X, check.names = FALSE))
    storage.mode(X) <- "double"

    data.frame(
      Assay     = an,
      Column    = rep(colnames(X), each = nrow(X)),
      Intensity = as.numeric(X),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  melted <- do.call(rbind, long_list)
  out <- merge(melted, cd, by = "Column", all.x = TRUE, sort = FALSE)
  out <- out[is.finite(out$Intensity), , drop = FALSE]
  out
}


# -----------------------------------------------------------------------------
# Función principal: Boxplots estilizados para proteómica
# -----------------------------------------------------------------------------

#' Boxplots Estilizados para Datos de Proteómica
#'
#' @param se SummarizedExperiment con los datos
#' @param ain Vector de assays a incluir
#' @param color_by Variable para colorear (ej: "Condition")
#' @param label_by Variable para etiquetas del eje Y (ej: "Column")
#' @param facet_norm Logical, usar facetas para cada normalización
#' @param ncol Número de columnas en facetas
#' @param assay_order Orden de los assays en facetas
#' @param group_order Orden de los grupos/condiciones
#' @param palette Paleta de colores: "ggsci::palette", "brewer:Name", o vector
#' @param title Título del gráfico (opcional)
#' @param subtitle Subtítulo del gráfico (opcional)
#' @param show_points Mostrar puntos individuales (jitter)
#' @param point_alpha Transparencia de los puntos jitter (0-1)
#' @param box_alpha Transparencia del relleno de los boxplots (0-1)
#' @param outlier_shape Forma de outliers (NA para ocultar)
#' @param notch Añadir muesca de intervalo de confianza
#'
#' @return ggplot object o lista de ggplots
#'
plot_boxplots_proteomics <- function(
    se,
    ain = NULL,
    color_by = NULL,
    label_by = NULL,
    facet_norm = TRUE,
    ncol = 2,
    assay_order = NULL,
    group_order = NULL,
    palette = NULL,
    title = NULL,
    subtitle = NULL,
    show_points = FALSE,
    point_alpha = 0.4,
    box_alpha = 0.75,
    outlier_shape = 21,
    notch = FALSE
) {

  # -------------------------------------------------------------------------
  # 1) Validación de inputs
  # -------------------------------------------------------------------------
  ain <- check_input_assays(se, ain)
  if (is.null(ain)) return(NULL)

  color_by <- get_color_value(se, color_by)
  tmp <- get_label_value(se, label_by)
  show_sample_names <- tmp[[1]]
  label_by <- tmp[[2]]

  # -------------------------------------------------------------------------
  # 2) Obtener datos en formato largo
  # -------------------------------------------------------------------------
  melted_dt <- if (exists("get_complete_dt")) {
    get_complete_dt(se, ain = ain)
  } else {
    get_complete_dt_safe(se, ain)
  }

  # -------------------------------------------------------------------------
  # 3) Ordenar Assays (facetas)
  # -------------------------------------------------------------------------
  if (!is.null(assay_order)) {
    melted_dt$Assay <- factor(melted_dt$Assay, levels = assay_order)
    ain <- intersect(assay_order, ain)
  } else {
    melted_dt$Assay <- factor(melted_dt$Assay, levels = unique(melted_dt$Assay))
  }

  # Etiquetas más descriptivas para los assays
  assay_labels <- c(
    "log2" = "Log₂ Intensity",
    "LoessCyc" = "LOESS Cyclic Normalization",
    "raw" = "Raw Intensity",
    "vsn" = "VSN Normalized",
    "quantile" = "Quantile Normalized"
  )

  melted_dt$Assay_Label <- ifelse(
    as.character(melted_dt$Assay) %in% names(assay_labels),
    assay_labels[as.character(melted_dt$Assay)],
    as.character(melted_dt$Assay)
  )
  melted_dt$Assay_Label <- factor(melted_dt$Assay_Label,
                                   levels = unique(melted_dt$Assay_Label))

  # -------------------------------------------------------------------------
  # 4) Ordenar grupos (color_by)
  # -------------------------------------------------------------------------
  if (!is.null(color_by)) {
    if (!is.null(group_order)) {
      melted_dt[[color_by]] <- factor(melted_dt[[color_by]], levels = group_order)
    } else {
      melted_dt[[color_by]] <- factor(melted_dt[[color_by]],
                                       levels = unique(melted_dt[[color_by]]))
    }
  }

  # -------------------------------------------------------------------------
  # 5) Ordenar muestras por condición y réplica
  # -------------------------------------------------------------------------
  if ("Replicate" %in% names(melted_dt)) {
    rep_num <- suppressWarnings(as.integer(melted_dt$Replicate))
  } else {
    rep_num <- suppressWarnings(
      as.integer(sub(".*?(\\d+)\\s*$", "\\1", as.character(melted_dt$Column)))
    )
  }
  rep_num[is.na(rep_num)] <- 1L

  ord_idx <- order(
    if (!is.null(color_by)) melted_dt[[color_by]] else factor(1),
    rep_num,
    as.character(melted_dt$Column)
  )

  labels_in_order <- unique(
    melted_dt[[if (show_sample_names) label_by else "Column"]][ord_idx]
  )
  label_levels <- rev(labels_in_order)

  if (show_sample_names) {
    melted_dt$Label <- factor(melted_dt[[label_by]], levels = label_levels)
  } else {
    melted_dt$Label <- factor(melted_dt$Column, levels = label_levels)
  }

  # -------------------------------------------------------------------------
  # 6) Configuración de paleta de colores
  # -------------------------------------------------------------------------
  col_values <- NULL
  if (!is.null(color_by)) {
    lvls <- levels(melted_dt[[color_by]])
    if (is.null(lvls)) lvls <- unique(melted_dt[[color_by]])

    # Paleta por defecto moderna
    default_palette <- c(
      "#457B9D",  # Azul acero
      "#E63946",  # Rojo coral
      "#2A9D8F",  # Verde azulado
      "#E9C46A",  # Amarillo mostaza
      "#9B5DE5",  # Púrpura
      "#F4A261",  # Naranja melocotón
      "#264653",  # Azul oscuro
      "#00BBF9"   # Cian brillante
    )

    if (is.null(palette)) {
      pal <- default_palette
      if (length(pal) < length(lvls)) pal <- rep(pal, length.out = length(lvls))
      col_values <- stats::setNames(pal[seq_along(lvls)], lvls)

    } else if (is.character(palette) && length(palette) == 1 && grepl("::", palette)) {
      if (!requireNamespace("paletteer", quietly = TRUE)) {
        stop("Para usar paletteer, instala 'paletteer'.")
      }
      pal <- as.character(paletteer::paletteer_d(palette))
      if (length(pal) < length(lvls)) pal <- rep(pal, length.out = length(lvls))
      col_values <- stats::setNames(pal[seq_along(lvls)], lvls)

    } else if (is.character(palette) && length(palette) == 1 && startsWith(palette, "brewer:")) {
      nm <- sub("^brewer:", "", palette)
      maxc <- RColorBrewer::brewer.pal.info[nm, "maxcolors"]
      pal <- RColorBrewer::brewer.pal(maxc, nm)
      if (length(pal) < length(lvls)) pal <- rep(pal, length.out = length(lvls))
      col_values <- stats::setNames(pal[seq_along(lvls)], lvls)

    } else if (is.character(palette) && length(palette) > 1) {
      if (length(palette) < length(lvls)) palette <- rep(palette, length.out = length(lvls))
      col_values <- stats::setNames(palette[seq_along(lvls)], lvls)

    } else {
      palette <- unlist(palette)
      if (is.null(names(palette))) {
        stop("Si pasas un vector para 'palette', debe estar nombrado por niveles de '",
             color_by, "'.")
      }
      missing <- setdiff(lvls, names(palette))
      if (length(missing)) {
        stop("Faltan colores para niveles: ", paste(missing, collapse = ", "))
      }
      col_values <- palette[lvls]
    }
  }

  # -------------------------------------------------------------------------
  # 7) Construcción del gráfico
  # -------------------------------------------------------------------------

  build_boxplot <- function(data, use_facet = FALSE) {

    # Iniciar ggplot
    if (is.null(color_by)) {
      p <- ggplot(data, aes(x = Intensity, y = Label))
    } else {
      p <- ggplot(data, aes(x = Intensity, y = Label, fill = .data[[color_by]]))
    }

    # Añadir puntos jitter si se solicita
    if (show_points) {
      if (is.null(color_by)) {
        p <- p + geom_point(
          position = position_jitter(width = 0, height = 0.2),
          alpha = point_alpha,
          size = 0.8,
          color = "#6C757D",
          shape = 16
        )
      } else {
        p <- p + geom_point(
          aes(color = .data[[color_by]]),
          position = position_jitter(width = 0, height = 0.2),
          alpha = point_alpha,
          size = 0.8,
          shape = 16
        )
      }
    }

    # Boxplot principal
    if (is.null(color_by)) {
      p <- p + geom_boxplot(
        fill = "#457B9D",
        color = "#1D3557",
        alpha = box_alpha,
        outlier.shape = outlier_shape,
        outlier.size = 1.5,
        outlier.alpha = 0.6,
        outlier.fill = "#ADB5BD",
        outlier.color = "#495057",
        notch = notch,
        linewidth = 0.5,
        width = 0.7
      )
    } else {
      p <- p + geom_boxplot(
        color = "#495057",
        alpha = box_alpha,
        outlier.shape = outlier_shape,
        outlier.size = 1.5,
        outlier.alpha = 0.6,
        notch = notch,
        linewidth = 0.5,
        width = 0.7
      )
    }

    # Bigotes con estilo
    p <- p + stat_boxplot(
      geom = "errorbar",
      width = 0.3,
      linewidth = 0.4,
      color = "#495057"
    )

    # Escalas de color
    if (!is.null(color_by) && !is.null(col_values)) {
      p <- p +
        scale_fill_manual(name = color_by, values = col_values) +
        scale_color_manual(name = color_by, values = col_values, guide = "none")
    }

    # Etiquetas
    p <- p + labs(
      x = expression(bold(log[2]~Intensity)),
      y = NULL
    )

    # Títulos si se proporcionan
    if (!is.null(title) || !is.null(subtitle)) {
      p <- p + labs(title = title, subtitle = subtitle)
    }

    # Facetas
    if (use_facet) {
      p <- p + facet_wrap(
        ~Assay_Label,
        scales = "free_x",
        ncol = ncol
      )
    }

    # Tema
    p <- p + theme_proteomics_boxplot()

    # Ocultar etiquetas del eje Y si no se muestran nombres de muestra
    if (!show_sample_names) {
      p <- p + theme(axis.text.y = element_blank())
    }

    # Guías de leyenda
    p <- p + guides(
      fill = guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        nrow = 1,
        override.aes = list(alpha = 1)
      )
    )

    p
  }

  # -------------------------------------------------------------------------
  # 8) Generar output
  # -------------------------------------------------------------------------

  if (isTRUE(facet_norm)) {
    p <- build_boxplot(melted_dt, use_facet = TRUE)
    return(p)

  } else {
    pl <- list()
    for (method in ain) {
      dt <- melted_dt[melted_dt$Assay == method, , drop = FALSE]
      assay_title <- if (method %in% names(assay_labels)) {
        assay_labels[method]
      } else {
        method
      }

      tmp <- build_boxplot(dt, use_facet = FALSE) +
        labs(title = assay_title)

      pl[[method]] <- tmp
    }
    return(pl)
  }
}


# =============================================================================
# EJEMPLOS DE USO
# =============================================================================

# --- Ejemplo básico ---
# p <- plot_boxplots_proteomics(
#   se          = se_proc,
#   ain         = c("log2", "LoessCyc"),
#   color_by    = "Condition",
#   label_by    = "Column",
#   facet_norm  = TRUE,
#   ncol        = 2,
#   assay_order = c("log2", "LoessCyc"),
#   group_order = c("A", "B", "C", "D")
# )
# p

# --- Con puntos jitter para ver distribución ---
# p <- plot_boxplots_proteomics(
#   se          = se_proc,
#   ain         = c("log2", "LoessCyc"),
#   color_by    = "Condition",
#   label_by    = "Column",
#   show_points = TRUE,
#   point_alpha = 0.3,
#   box_alpha   = 0.7
# )
# p

# --- Con título y paleta personalizada ---
# p <- plot_boxplots_proteomics(
#   se          = se_proc,
#   ain         = c("log2", "LoessCyc"),
#   color_by    = "Condition",
#   label_by    = "Column",
#   title       = "Distribución de Intensidades por Condición",
#   subtitle    = "Comparación antes y después de normalización LOESS",
#   palette     = "ggsci::nrc_npg"
# )
# p

# --- Con muescas de intervalo de confianza ---
# p <- plot_boxplots_proteomics(
#   se          = se_proc,
#   ain         = c("log2", "LoessCyc"),
#   color_by    = "Condition",
#   notch       = TRUE,
#   outlier_shape = NA  # Ocultar outliers
# )
# p

# --- Paleta manual ---
# mis_colores <- c(
#   "A" = "#264653",
#   "B" = "#2A9D8F",
#   "C" = "#E9C46A",
#   "D" = "#E76F51"
# )
# p <- plot_boxplots_proteomics(
#   se          = se_proc,
#   ain         = c("log2", "LoessCyc"),
#   color_by    = "Condition",
#   palette     = mis_colores
# )
# p
