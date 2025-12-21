# =============================================================================
# Boxplot Estilizado para Datos de Proteómica
# =============================================================================

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
# Función principal: Boxplots estilizados para proteómica
# -----------------------------------------------------------------------------

#' Boxplots Estilizados para Datos de Proteómica
#'
#' @param data Data frame con columnas: Column, Assay, Intensity, Condition, Replicate
#' @param assays Vector de assays a incluir (NULL = todos)
#' @param color_by Columna para colorear (default: "Condition")
#' @param label_by Columna para etiquetas del eje Y (default: "Column")
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
    data,
    assays = NULL,
    color_by = "Condition",
    label_by = "Column",
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
  required_cols <- c("Column", "Assay", "Intensity")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Columnas requeridas faltantes: ", paste(missing_cols, collapse = ", "))
  }

  # Verificar columna de color
  if (!is.null(color_by) && !(color_by %in% names(data))) {
    warning("'", color_by, "' no está en el data frame. Se ignorará color_by.")
    color_by <- NULL
  }

  # Verificar columna de etiqueta
  if (!(label_by %in% names(data))) {
    label_by <- "Column"
  }

  # Copiar datos
 melted_dt <- data

  # Filtrar assays si se especifican
  if (!is.null(assays)) {
    melted_dt <- melted_dt[melted_dt$Assay %in% assays, , drop = FALSE]
  }

  # Eliminar valores no finitos
  melted_dt <- melted_dt[is.finite(melted_dt$Intensity), , drop = FALSE]

  # -------------------------------------------------------------------------
  # 2) Ordenar Assays (facetas)
  # -------------------------------------------------------------------------
  if (!is.null(assay_order)) {
    melted_dt$Assay <- factor(melted_dt$Assay, levels = assay_order)
  } else {
    melted_dt$Assay <- factor(melted_dt$Assay, levels = unique(melted_dt$Assay))
  }

  # Etiquetas más descriptivas para los assays
  assay_labels <- c(
    "log2" = "Log\u2082 Intensity",
    "LoessCyc" = "LOESS Cyclic",
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
  # 3) Ordenar grupos (color_by)
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
  # 4) Ordenar muestras por condición y réplica
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

  labels_in_order <- unique(melted_dt[[label_by]][ord_idx])
  label_levels <- rev(labels_in_order)

  melted_dt$Label <- factor(melted_dt[[label_by]], levels = label_levels)

  # -------------------------------------------------------------------------
  # 5) Configuración de paleta de colores
  # -------------------------------------------------------------------------
  col_values <- NULL
  if (!is.null(color_by)) {
    lvls <- levels(melted_dt[[color_by]])
    if (is.null(lvls)) lvls <- unique(melted_dt[[color_by]])

    # Paleta por defecto moderna
    default_palette <- c(
      "#457B9D",
      "#E63946",
      "#2A9D8F",
      "#E9C46A",
      "#9B5DE5",
      "#F4A261",
      "#264653",
      "#00BBF9"
    )

    if (is.null(palette)) {
      pal <- default_palette
      if (length(pal) < length(lvls)) pal <- rep(pal, length.out = length(lvls))
      col_values <- stats::setNames(pal[seq_along(lvls)], lvls)

    } else if (is.character(palette) && length(palette) == 1 && grepl("::", palette)) {
      if (!requireNamespace("paletteer", quietly = TRUE)) {
        stop("Para usar paletteer, instala con: install.packages('paletteer')")
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
  # 6) Construcción del gráfico
  # -------------------------------------------------------------------------

  build_boxplot <- function(dt, use_facet = FALSE) {

    # Iniciar ggplot
    if (is.null(color_by)) {
      p <- ggplot(dt, aes(x = Intensity, y = Label))
    } else {
      p <- ggplot(dt, aes(x = Intensity, y = Label, fill = .data[[color_by]]))
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

    # Guías de leyenda
    if (!is.null(color_by)) {
      p <- p + guides(
        fill = guide_legend(
          title.position = "top",
          title.hjust = 0.5,
          nrow = 1,
          override.aes = list(alpha = 1)
        )
      )
    }

    p
  }

  # -------------------------------------------------------------------------
  # 7) Generar output
  # -------------------------------------------------------------------------

  current_assays <- levels(melted_dt$Assay)

  if (isTRUE(facet_norm)) {
    p <- build_boxplot(melted_dt, use_facet = TRUE)
    return(p)

  } else {
    pl <- list()
    for (method in current_assays) {
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
#   data        = mi_dataframe,
#   assays      = c("log2", "LoessCyc"),
#   color_by    = "Condition",
#   label_by    = "Column",
#   facet_norm  = TRUE,
#   ncol        = 2,
#   group_order = c("A", "B", "C", "D")
# )
# p

# --- Con puntos jitter ---
# p <- plot_boxplots_proteomics(
#   data        = mi_dataframe,
#   assays      = c("log2", "LoessCyc"),
#   color_by    = "Condition",
#   show_points = TRUE,
#   point_alpha = 0.3
# )
# p

# --- Con título y subtítulo ---
# p <- plot_boxplots_proteomics(
#   data     = mi_dataframe,
#   color_by = "Condition",
#   title    = "Distribución de Intensidades",
#   subtitle = "Comparación de métodos de normalización"
# )
# p

# --- Con paleta personalizada ---
# p <- plot_boxplots_proteomics(
#   data     = mi_dataframe,
#   color_by = "Condition",
#   palette  = c("#264653", "#2A9D8F", "#E9C46A", "#E76F51")
# )
# p
