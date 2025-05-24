library(tidyverse)
library(pheatmap)
library(cowplot)
#library(ggpubr)

library(RColorBrewer)
library(ggsci)


jjVolcano <- function(diffData = NULL, myMarkers = NULL, order.by = c("avg_log2FC"),
                      log2FC.cutoff = 0.25, pvalue.cutoff = 0.05, adjustP.cutoff = 0.01,
                      topGeneN = 5, col.type = "updown", back.col = "grey93",
                      pSize = 0.75, aesCol = c("#0099CC", "#CC3333"), legend.position = c(
                        0.7,
                        0.9
                      ), base_size = 14, tile.col = jjAnno::useMyCol("paired",
                        n = 9
                      ), cluster.order = NULL, polar = FALSE, expand = c(
                        -1,
                        1
                      ), flip = FALSE, ...) {
  diff.marker <- diffData %>% dplyr::filter(abs(avg_log2FC) >=
    log2FC.cutoff & p_val < pvalue.cutoff)
  diff.marker <- diff.marker %>%
    dplyr::mutate(type = ifelse(avg_log2FC >=
      log2FC.cutoff, "sigUp", "sigDown")) %>%
    dplyr::mutate(type2 = ifelse(p_val_adj <
      adjustP.cutoff, paste("adjust Pvalue < ", adjustP.cutoff,
      sep = ""
    ), paste("adjust Pvalue >= ", adjustP.cutoff,
      sep = ""
    )))
  if (!is.null(cluster.order)) {
    diff.marker$cluster <- factor(diff.marker$cluster, levels = cluster.order)
  }
  back.data <- purrr::map_df(
    unique(diff.marker$cluster),
    function(x) {
      tmp <- diff.marker %>% dplyr::filter(cluster ==
        x)
      new.tmp <- data.frame(cluster = x, min = min(tmp$avg_log2FC) -
        0.2, max = max(tmp$avg_log2FC) + 0.2)
      return(new.tmp)
    }
  )
  top.marker.tmp <- diff.marker %>% dplyr::group_by(cluster)
  top.marker.max <- top.marker.tmp %>% dplyr::slice_max(
    n = topGeneN,
    order_by = get(order.by)
  )
  top.marker.min <- top.marker.tmp %>% dplyr::slice_min(
    n = topGeneN,
    order_by = get(order.by)
  )
  top.marker <- rbind(top.marker.max, top.marker.min)
  if (!is.null(myMarkers)) {
    top.marker <- diff.marker %>% dplyr::filter(gene %in%
      myMarkers)
  } else {
    top.marker <- top.marker
  }
  p1 <- ggplot2::ggplot(diff.marker, ggplot2::aes(
    x = cluster,
    y = avg_log2FC
  )) +
    ggplot2::geom_col(
      data = back.data,
      ggplot2::aes(x = cluster, y = min), fill = back.col
    ) +
    ggplot2::geom_col(data = back.data, ggplot2::aes(
      x = cluster,
      y = max
    ), fill = back.col)
  if (col.type == "updown") {
    p2 <- p1 + ggplot2::geom_jitter(ggplot2::aes(color = type),
      size = pSize
    ) + ggplot2::scale_color_manual(values = c(
      sigDown = aesCol[1],
      sigUp = aesCol[2]
    ))
  } else if (col.type == "adjustP") {
    p2 <- p1 + ggplot2::geom_jitter(ggplot2::aes(color = type2),
      size = pSize
    ) + ggplot2::scale_color_manual(values = c(
      aesCol[2],
      aesCol[1]
    ))
  }
  p3 <- p2 + ggplot2::scale_y_continuous(n.breaks = 6) + ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = legend.position, legend.title = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank()
    ) +
    ggplot2::xlab("Clusters") + ggplot2::ylab("Average log2FoldChange") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5)))
  p4 <- p3 + ggplot2::geom_tile(ggplot2::aes(
    x = cluster,
    y = 0, fill = cluster
  ), color = "black", height = log2FC.cutoff *
    2, alpha = 0.3, show.legend = F) + ggplot2::scale_fill_manual(values = tile.col) +
    ggrepel::geom_text_repel(
      data = top.marker, ggplot2::aes(
        x = cluster,
        y = avg_log2FC, label = gene
      ), max.overlaps = 50,
      ...
    )
  if (polar == TRUE) {
    p5 <- p4 + geomtextpath::geom_textpath(ggplot2::aes(
      x = cluster,
      y = 0, label = cluster
    )) + ggplot2::scale_y_continuous(
      n.breaks = 6,
      expand = ggplot2::expansion(mult = expand)
    ) + ggplot2::theme_void(base_size = base_size) +
      ggplot2::theme(
        legend.position = legend.position,
        legend.title = ggplot2::element_blank()
      ) + ggplot2::coord_polar(
        clip = "off",
        theta = "x"
      )
  } else {
    if (flip == TRUE) {
      p5 <- p4 + ggplot2::scale_y_continuous(n.breaks = 6) +
        ggplot2::geom_label(ggplot2::aes(
          x = cluster,
          y = 0, label = cluster
        )) + ggplot2::theme(
          axis.line.y = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()
        ) +
        ggplot2::coord_flip()
    } else {
      p5 <- p4 + ggplot2::scale_y_continuous(n.breaks = 6) +
        ggplot2::geom_text(ggplot2::aes(
          x = cluster,
          y = 0, label = cluster
        )) + ggplot2::theme(
          axis.line.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()
        )
    }
  }
  return(p5)
}



#####################################################
# basic
#####################################################

# get ann_col for pheatmap color bar
getAnnCol <- function(col_ann = NULL) {
  # color
  ann_col <- list()
  for (col in colnames(col_ann)) {
    if (!is.factor(col_ann[, col])) col_ann[, col] <- factor(col_ann[, col])
    num <- length(levels(col_ann[, col]))
    colors <-
      if (num <= 9) {
        if (num < 3) {
          num <- 3
        }
        RColorBrewer::brewer.pal(num, "Set1")
      } else {
        gg_color_hue(num)
      }
    ann_col[[col]] <- stats::setNames(colors, levels(col_ann[, col]))
  }
  # if num < 3, remove na
  ann_col <- lapply(ann_col, FUN = function(x) {
    x[!is.na(names(x))]
  })
  return(ann_col)
}

# modified from pheatmap:::scale_rows, not use sd
scale_rows <- function(x, scale = FALSE) {
  m <- apply(x, 1, mean, na.rm = T)
  if (scale) {
    s <- apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
  } else {
    return(x - m)
  }
}

range01_rows <- function(x) {
  range01 <- function(a) {
    (a - min(a)) / (max(a) - min(a))
  }
  m <- apply(x, 1, range01)
  return(t(m))
}

# get point size for geom_point
get_pt_size <- function(x) {
  ### Adjust point size
  if (x < 200) {
    ptsize <- 4
  } else if (x > 200 & x < 1000) {
    ptsize <- 3.5
  } else if (x > 1000 & x < 2000) {
    ptsize <- 3
  } else if (x > 2000 & x < 5000) {
    ptsize <- 2
  } else if (x > 5000 & x < 10000) {
    ptsize <- 1
  } else if (x > 10000) {
    ptsize <- 0.5
  }
  return(ptsize)
}


## extract the max value of the y axis
extract_max <- function(p) {
  ymax <- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

#####################################################
# rast DimPlot
#####################################################

# Set a default value if an object is NOT null
#
# @param lhs An object to set if it's NOT null
# @param rhs The value to provide if x is NOT null
#
# @return lhs if lhs is null, else rhs
#
# @author Hadley Wickham
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
#
`%iff%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(rhs)
  } else {
    return(lhs)
  }
}

SingleDimPlot_rast <- function(data, dims, col.by = NULL, cols = NULL, pt.size = NULL,
                               shape.by = NULL, order = NULL, label = FALSE, repel = FALSE,
                               label.size = 4, cells.highlight = NULL, cols.highlight = "#DE2D26",
                               sizes.highlight = 1, na.value = "grey50") {
  pt.size <- pt.size %||% Seurat:::AutoPointSize(data = data)
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  if (!is.data.frame(x = data)) {
    data <- as.data.frame(x = data)
  }
  if (is.character(x = dims) && !all(dims %in% colnames(x = data))) {
    stop("Cannot find dimensions to plot in data")
  } else if (is.numeric(x = dims)) {
    dims <- colnames(x = data)[dims]
  }
  if (!is.null(x = cells.highlight)) {
    highlight.info <- SetHighlight(
      cells.highlight = cells.highlight,
      cells.all = rownames(x = data), sizes.highlight = sizes.highlight %||%
        pt.size, cols.highlight = cols.highlight, col.base = cols[1] %||%
        "#C3C3C3", pt.size = pt.size
    )
    order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- "highlight"
    pt.size <- highlight.info$size
    cols <- highlight.info$color
  }
  if (!is.null(x = order) && !is.null(x = col.by)) {
    if (typeof(x = order) == "logical") {
      if (order) {
        data <- data[order(data[, col.by]), ]
      }
    } else {
      order <- rev(x = c(order, setdiff(x = unique(x = data[
        ,
        col.by
      ]), y = order)))
      data[, col.by] <- factor(x = data[, col.by], levels = order)
      new.order <- order(x = data[, col.by])
      data <- data[new.order, ]
      if (length(x = pt.size) == length(x = new.order)) {
        pt.size <- pt.size[new.order]
      }
    }
  }
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find ", col.by, " in plotting data, not coloring plot")
    col.by <- NULL
  } else {
    col.index <- match(x = col.by, table = colnames(x = data))
    if (grepl(pattern = "^\\d", x = col.by)) {
      col.by <- paste0("x", col.by)
    } else if (grepl(pattern = "-", x = col.by)) {
      col.by <- gsub(
        pattern = "-", replacement = ".",
        x = col.by
      )
    }
    colnames(x = data)[col.index] <- col.by
  }
  if (!is.null(x = shape.by) && !shape.by %in% colnames(x = data)) {
    warning("Cannot find ", shape.by, " in plotting data, not shaping plot")
  }
  plot <- ggplot(data = data) +
    ggrastr::geom_point_rast(
      mapping = aes_string(
        x = dims[1],
        y = dims[2], color = paste0("`", col.by, "`"), shape = shape.by
      ),
      size = pt.size
    ) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(color = NULL)
  if (label && !is.null(x = col.by)) {
    plot <- LabelClusters(
      plot = plot, id = col.by, repel = repel,
      size = label.size
    )
  }
  if (!is.null(x = cols)) {
    if (length(x = cols) == 1 && (is.numeric(x = cols) ||
      cols %in% rownames(x = brewer.pal.info))) {
      scale <- scale_color_brewer(palette = cols, na.value = na.value)
    } else if (length(x = cols) == 1 && (cols %in% c(
      "alphabet",
      "alphabet2", "glasbey", "polychrome", "stepped"
    ))) {
      colors <- DiscretePalette(length(unique(data[[col.by]])),
        palette = cols
      )
      scale <- scale_color_manual(values = colors, na.value = na.value)
    } else {
      scale <- scale_color_manual(values = cols, na.value = na.value)
    }
    plot <- plot + scale
  }
  plot <- plot + cowplot::theme_cowplot(
    font_size = 12,
    rel_small = 10 / 12, rel_large = 12 / 12, rel_tiny = 9 / 14,
    font_family = "Arial"
  )
  return(plot)
}

DimPlot_rast <- function(object, dims = c(1, 2), cells = NULL, cols = NULL,
                         pt.size = NULL, reduction = NULL, group.by = NULL, split.by = NULL,
                         shape.by = NULL, order = NULL, label = FALSE, label.size = 4,
                         repel = FALSE, cells.highlight = NULL, cols.highlight = "#DE2D26",
                         sizes.highlight = 1, na.value = "grey50", ncol = NULL, combine = TRUE) {
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  reduction <- reduction %||% DefaultDimReduc(object = object)
  cells <- cells %||% colnames(x = object)
  data <- Embeddings(object = object[[reduction]])[cells, dims]
  data <- as.data.frame(x = data)
  dims <- paste0(Key(object = object[[reduction]]), dims)
  object[["ident"]] <- Idents(object = object)
  orig.groups <- group.by
  group.by <- group.by %||% "ident"
  data[, group.by] <- object[[group.by]][cells, , drop = FALSE]
  for (group in group.by) {
    if (!is.factor(x = data[, group])) {
      data[, group] <- factor(x = data[, group])
    }
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  if (!is.null(x = split.by)) {
    data[, split.by] <- object[[split.by, drop = TRUE]]
  }
  plots <- lapply(X = group.by, FUN = function(x) {
    plot <- SingleDimPlot_rast(
      data = data[, c(
        dims, x, split.by,
        shape.by
      )], dims = dims, col.by = x, cols = cols,
      pt.size = pt.size, shape.by = shape.by, order = order,
      label = FALSE, cells.highlight = cells.highlight,
      cols.highlight = cols.highlight, sizes.highlight = sizes.highlight,
      na.value = na.value
    )
    if (label) {
      plot <- LabelClusters(
        plot = plot, id = x, repel = repel,
        size = label.size, split.by = split.by
      )
    }
    if (!is.null(x = split.by)) {
      plot <- plot + FacetTheme() + facet_wrap(
        facets = vars(!!sym(x = split.by)),
        ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
          length(x = unique(x = data[, split.by]))
        } else {
          ncol
        }
      )
    }
    return(plot)
  })
  if (!is.null(x = split.by)) {
    ncol <- 1
  }
  if (combine) {
    plots <- wrap_plots(plots, ncol = orig.groups %iff% ncol)
  }
  return(plots)
}

#####################################################
# StackedVlnPlot
# from: CellChat
#####################################################

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot <- function(object,
                           features,
                           idents = NULL,
                           split.by = NULL,
                           cols = NULL,
                           show.median = FALSE,
                           median.size = 1,
                           show.text.y = TRUE,
                           line.size = NULL,
                           pt.size = 0,
                           plot.margin = margin(0, 0, 0, 0, "cm"),
                           ...) {
  options(warn = -1)
  p <- Seurat::VlnPlot(object, features = features, cols = cols, pt.size = pt.size, idents = idents, split.by = split.by, ...) +
    xlab("") + ylab(features) + ggtitle("")
  if (show.median) {
    p <- p + stat_summary(fun.y = median, geom = "point", shape = 3, size = median.size)
  }
  p <- p + theme(text = element_text(size = 6)) + theme(axis.line = element_line(size = line.size)) +
    theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.line.x = element_line(colour = "black", size = line.size), axis.line.y = element_line(colour = "black", size = line.size))
  # theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
  p <- p + theme(
    plot.title = element_blank(), legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = rel(1), angle = 0),
    axis.text.y = element_text(size = rel(1)),
    plot.margin = plot.margin
  ) +
    theme(axis.text.y = element_text(size = 6))
  p <- p + theme(element_line(size = line.size))

  if (!show.text.y) {
    p <- p + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
  }
  return(p)
}

## extract the max value of the y axis
extract_max <- function(p) {
  ymax <- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

## main function
StackedVlnPlot <- function(object, features, idents = NULL, split.by = NULL,
                           color.use = NULL, colors.ggplot = FALSE, show.median = FALSE, median.size = 1,
                           angle.x = 60, vjust.x = NULL, hjust.x = NULL, show.text.y = TRUE, line.size = NULL,
                           pt.size = 0,
                           plot.margin = margin(0, 0, 0, 0, "cm"),
                           ...) {
  options(warn = -1)
  if (is.null(color.use)) {
    numCluster <- length(levels(Seurat::Idents(object)))
    if (colors.ggplot) {
      color.use <- NULL
    } else {
      color.use <- scPalette(numCluster)
    }
  }
  if (is.null(vjust.x) | is.null(hjust.x)) {
    angle <- c(0, 45, 60, 90)
    hjust <- c(0, 1, 1, 1)
    vjust <- c(0, 1, 1, 0.5)
    vjust.x <- vjust[angle == angle.x]
    hjust.x <- hjust[angle == angle.x]
  }

  plot_list <- purrr::map(features, function(x) {
    modify_vlnplot(
      object = object, features = x, idents = idents, split.by = split.by, cols = color.use, show.median = show.median, median.size = median.size, pt.size = pt.size,
      show.text.y = show.text.y, line.size = line.size, ...
    )
  })

  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
    theme(axis.text.x = element_text(), axis.ticks.x = element_line()) +
    theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x, vjust = vjust.x)) +
    theme(axis.text.x = element_text(size = 6))

  # change the y-axis tick to only max value
  ymaxs <- purrr::map_dbl(plot_list, extract_max)
  plot_list <- purrr::map2(plot_list, ymaxs, function(x, y) {
    x +
      scale_y_continuous(breaks = c(y)) +
      expand_limits(y = y)
  })

  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1) + patchwork::plot_layout(guides = "collect")
  return(p)
}


#####################################################
# custom VlnPlot: no legend, no point, custom color
#####################################################

# plot gene expression using VlnPlot: one gene, across group
customVlnPlot <- function(seurat_obj = NULL, features = NULL, groupby = "seurat_clusters",
                          color_maps = cluster_color_maps, ...) {
  p <- VlnPlot(seurat_obj, features = features, group.by = groupby, combine = FALSE, ...)
  p <- lapply(X = p, FUN = function(x) {
    # change font size
    x <- x + FontSize(
      main = 6, x.title = 6, y.title = 6,
      x.text = 6, y.text = 6
    )
    x +
      labs(y = "Exp.") +
      scale_fill_manual(values = color_maps) +
      NoLegend()
  })
  p
  # patchwork::wrap_plots(p, ncol = ncol)
}

#####################################################
# custom FeaturePlot
# no legend, no point, custom color
# rast
#####################################################

FeaturePlot_rast <- function(
    object, features, dims = c(1, 2), cells = NULL, cols = if (blend) {
      c("lightgrey", "#ff0000", "#00ff00")
    } else {
      c("lightgrey", "blue")
    }, pt.size = NULL, order = FALSE, min.cutoff = NA, max.cutoff = NA,
    reduction = NULL, split.by = NULL, shape.by = NULL, slot = "data",
    blend = FALSE, blend.threshold = 0.5, label = FALSE, label.size = 4,
    repel = FALSE, ncol = NULL, coord.fixed = FALSE, by.col = TRUE,
    sort.cell = NULL, combine = TRUE) {
  if (!is.null(x = sort.cell)) {
    warning("The sort.cell parameter is being deprecated. Please use the order ",
      "parameter instead for equivalent functionality.",
      call. = FALSE, immediate. = TRUE
    )
    if (isTRUE(x = sort.cell)) {
      order <- sort.cell
    }
  }
  no.right <- theme(
    axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank(), axis.title.y.right = element_text(
      face = "bold",
      size = 14, margin = margin(r = 7)
    )
  )
  reduction <- reduction %||% Seurat:::DefaultDimReduc(object = object)
  if (length(x = dims) != 2 || !is.numeric(x = dims)) {
    stop("'dims' must be a two-length integer vector")
  }
  if (blend && length(x = features) != 2) {
    stop("Blending feature plots only works with two features")
  }
  if (blend) {
    default.colors <- eval(expr = formals(fun = FeaturePlot)$cols)
    cols <- switch(EXPR = as.character(x = length(x = cols)),
      `0` = {
        warning("No colors provided, using default colors",
          call. = FALSE, immediate. = TRUE
        )
        default.colors
      },
      `1` = {
        warning("Only one color provided, assuming specified is double-negative and augmenting with default colors",
          call. = FALSE, immediate. = TRUE
        )
        c(cols, default.colors[2:3])
      },
      `2` = {
        warning("Only two colors provided, assuming specified are for features and agumenting with '",
          default.colors[1], "' for double-negatives",
          call. = FALSE, immediate. = TRUE
        )
        c(default.colors[1], cols)
      },
      `3` = cols,
      {
        warning("More than three colors provided, using only first three",
          call. = FALSE, immediate. = TRUE
        )
        cols[1:3]
      }
    )
  }
  if (blend && length(x = cols) != 3) {
    stop("Blending feature plots only works with three colors; first one for negative cells")
  }
  dims <- paste0(Key(object = object[[reduction]]), dims)
  cells <- cells %||% colnames(x = object)
  data <- FetchData(object = object, vars = c(
    dims, "ident",
    features
  ), cells = cells, slot = slot)
  if (ncol(x = data) < 4) {
    stop("None of the requested features were found: ", paste(features,
      collapse = ", "
    ), " in slot ", slot, call. = FALSE)
  } else if (!all(dims %in% colnames(x = data))) {
    stop("The dimensions requested were not found", call. = FALSE)
  }
  features <- colnames(x = data)[4:ncol(x = data)]
  min.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = min(data[
      ,
      feature
    ]), no = cutoff))
  }, cutoff = min.cutoff, feature = features)
  max.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = max(data[
      ,
      feature
    ]), no = cutoff))
  }, cutoff = max.cutoff, feature = features)
  check.lengths <- unique(x = vapply(X = list(
    features, min.cutoff,
    max.cutoff
  ), FUN = length, FUN.VALUE = numeric(length = 1)))
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols, ]$maxcolors, no = length(x = cols))
  data[, 4:ncol(x = data)] <- sapply(
    X = 4:ncol(x = data),
    FUN = function(index) {
      data.feature <- as.vector(x = data[, index])
      min.use <- Seurat:::SetQuantile(cutoff = min.cutoff[index -
        3], data.feature)
      max.use <- Seurat:::SetQuantile(cutoff = max.cutoff[index -
        3], data.feature)
      data.feature[data.feature < min.use] <- min.use
      data.feature[data.feature > max.use] <- max.use
      if (brewer.gran == 2) {
        return(data.feature)
      }
      data.cut <- if (all(data.feature == 0)) {
        0
      } else {
        as.numeric(x = as.factor(x = cut(
          x = as.numeric(x = data.feature),
          breaks = brewer.gran
        )))
      }
      return(data.cut)
    }
  )
  colnames(x = data)[4:ncol(x = data)] <- features
  rownames(x = data) <- cells
  data$split <- if (is.null(x = split.by)) {
    Seurat:::RandomName()
  } else {
    switch(EXPR = split.by,
      ident = Idents(object = object)[cells],
      object[[split.by, drop = TRUE]][cells]
    )
  }
  if (!is.factor(x = data$split)) {
    data$split <- factor(x = data$split)
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  plots <- vector(mode = "list", length = ifelse(test = blend,
    yes = 4, no = length(x = features) * length(x = levels(x = data$split))
  ))
  xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[
    ,
    dims[1]
  ])))
  ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[
    ,
    dims[2]
  ])))
  if (blend) {
    ncol <- 4
    color.matrix <- BlendMatrix(
      two.colors = cols[2:3], col.threshold = blend.threshold,
      negative.color = cols[1]
    )
    cols <- cols[2:3]
    colors <- list(
      color.matrix[, 1], color.matrix[1, ],
      as.vector(x = color.matrix)
    )
  }
  for (i in 1:length(x = levels(x = data$split))) {
    ident <- levels(x = data$split)[i]
    data.plot <- data[as.character(x = data$split) == ident, ,
      drop = FALSE
    ]
    if (blend) {
      features <- features[1:2]
      no.expression <- features[colMeans(x = data.plot[
        ,
        features
      ]) == 0]
      if (length(x = no.expression) != 0) {
        stop("The following features have no value: ",
          paste(no.expression, collapse = ", "),
          call. = FALSE
        )
      }
      data.plot <- cbind(
        data.plot[, c(dims, "ident")],
        BlendExpression(data = data.plot[, features[1:2]])
      )
      features <- colnames(x = data.plot)[4:ncol(x = data.plot)]
    }
    for (j in 1:length(x = features)) {
      feature <- features[j]
      if (blend) {
        cols.use <- as.numeric(x = as.character(x = data.plot[
          ,
          feature
        ])) + 1
        cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
      } else {
        cols.use <- NULL
      }
      data.single <- data.plot[, c(
        dims, "ident", feature,
        shape.by
      )]
      plot <- SingleDimPlot_rast(
        data = data.single, dims = dims,
        col.by = feature, order = order, pt.size = pt.size,
        cols = cols.use, shape.by = shape.by, label = FALSE
      ) +
        scale_x_continuous(limits = xlims) + scale_y_continuous(limits = ylims) +
        theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
      if (label) {
        plot <- LabelClusters(
          plot = plot, id = "ident",
          repel = repel, size = label.size
        )
      }
      if (length(x = levels(x = data$split)) > 1) {
        plot <- plot + theme(panel.border = element_rect(
          fill = NA,
          colour = "black"
        ))
        plot <- plot + if (i == 1) {
          labs(title = feature)
        } else {
          labs(title = NULL)
        }
        if (j == length(x = features) && !blend) {
          suppressMessages(expr = plot <- plot + scale_y_continuous(
            sec.axis = dup_axis(name = ident),
            limits = ylims
          ) + no.right)
        }
        if (j != 1) {
          plot <- plot + theme(
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(), axis.text.y = element_blank(),
            axis.title.y.left = element_blank()
          )
        }
        if (i != length(x = levels(x = data$split))) {
          plot <- plot + theme(
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(), axis.text.x = element_blank(),
            axis.title.x = element_blank()
          )
        }
      } else {
        plot <- plot + labs(title = feature)
      }
      if (!blend) {
        plot <- plot + guides(color = NULL)
        cols.grad <- cols
        if (length(x = cols) == 1) {
          plot <- plot + scale_color_brewer(palette = cols)
        } else if (length(x = cols) > 1) {
          unique.feature.exp <- unique(data.plot[, feature])
          if (length(unique.feature.exp) == 1) {
            warning(
              "All cells have the same value (",
              unique.feature.exp, ") of ", feature, "."
            )
            if (unique.feature.exp == 0) {
              cols.grad <- cols[1]
            } else {
              cols.grad <- cols
            }
          }
          plot <- suppressMessages(expr = plot + scale_color_gradientn(
            colors = cols.grad,
            guide = "colorbar"
          ))
        }
      }
      if (coord.fixed) {
        plot <- plot + coord_fixed()
      }
      plot <- plot
      plots[[(length(x = features) * (i - 1)) + j]] <- plot
    }
  }
  if (blend) {
    blend.legend <- BlendMap(color.matrix = color.matrix)
    for (ii in 1:length(x = levels(x = data$split))) {
      suppressMessages(expr = plots <- append(
        x = plots,
        values = list(blend.legend + scale_y_continuous(
          sec.axis = dup_axis(name = ifelse(test = length(x = levels(x = data$split)) >
            1, yes = levels(x = data$split)[ii], no = "")),
          expand = c(0, 0)
        ) + labs(
          x = features[1], y = features[2],
          title = if (ii == 1) {
            paste("Color threshold:", blend.threshold)
          } else {
            NULL
          }
        ) + no.right), after = 4 * ii - 1
      ))
    }
  }
  plots <- Filter(f = Negate(f = is.null), x = plots)
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = features) == 1) {
      ncol <- 1
    }
    if (length(x = features) > 6) {
      ncol <- 3
    }
    if (length(x = features) > 9) {
      ncol <- 4
    }
  }
  ncol <- ifelse(test = is.null(x = split.by) || blend, yes = ncol,
    no = length(x = features)
  )
  legend <- if (blend) {
    "none"
  } else {
    split.by %iff% "none"
  }
  if (combine) {
    if (by.col && !is.null(x = split.by) && !blend) {
      plots <- lapply(X = plots, FUN = function(x) {
        return(suppressMessages(expr = x + theme_cowplot() +
          ggtitle("") + scale_y_continuous(
            sec.axis = dup_axis(name = ""),
            limits = ylims
          ) + no.right))
      })
      nsplits <- length(x = levels(x = data$split))
      idx <- 1
      for (i in (length(x = features) * (nsplits - 1) +
        1):(length(x = features) * nsplits)) {
        plots[[i]] <- suppressMessages(expr = plots[[i]] +
          scale_y_continuous(
            sec.axis = dup_axis(name = features[[idx]]),
            limits = ylims
          ) + no.right)
        idx <- idx + 1
      }
      idx <- 1
      for (i in which(x = 1:length(x = plots) %% length(x = features) ==
        1)) {
        plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) +
          theme(plot.title = element_text(hjust = 0.5))
        idx <- idx + 1
      }
      idx <- 1
      if (length(x = features) == 1) {
        for (i in 1:length(x = plots)) {
          plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) +
            theme(plot.title = element_text(hjust = 0.5))
          idx <- idx + 1
        }
        ncol <- 1
        nrow <- nsplits
      } else {
        nrow <- split.by %iff% length(x = levels(x = data$split))
      }
      plots <- plots[c(do.call(what = rbind, args = split(
        x = 1:length(x = plots),
        f = ceiling(x = seq_along(along.with = 1:length(x = plots)) / length(x = features))
      )))]
      plots <- wrap_plots(plots, ncol = nrow, nrow = ncol)
      if (!is.null(x = legend) && legend == "none") {
        plots <- plots & NoLegend()
      }
    } else {
      plots <- wrap_plots(plots, ncol = ncol, nrow = split.by %iff%
        length(x = levels(x = data$split)))
    }
    if (!is.null(x = legend) && legend == "none") {
      plots <- plots & NoLegend()
    }
  }
  return(plots)
}

# plot gene expression using FeaturePlot
customFeaturePlot <- function(seurat_obj = NULL,
                              features = NULL,
                              no_legend = TRUE,
                              no_axis = TRUE,
                              add_border = TRUE, ...) {
  p <- FeaturePlot(seurat_obj,
    features = features,
    combine = FALSE, ...
  )
  p <- lapply(X = p, FUN = function(x) {
    # change font size
    x <- x + FontSize(
      main = 6, x.title = 6, y.title = 6,
      x.text = 6, y.text = 6
    )
    # no axis
    if (no_axis) {
      x <- x +
        NoAxes()
    }
    # no legend
    if (no_legend) {
      x <- x + NoLegend()
    }
    # add border
    if (add_border) {
      x <- x +
        theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
    }
    return(x)
  })
  p
  # patchwork::wrap_plots(p, ncol = ncol)
}

#####################################################

# plot several genes using ComplexHeatmap

#####################################################

remove_nan <- function(
    features = NULL,
    seurat_obj = NULL) {
  .helper_remove_nan <- function(
      features = NULL,
      seurat_obj = NULL) {
    # set idents
    seuratObj <- seurat_obj
    Idents(seuratObj) <- "seurat_clusters"
    # get average exp
    var_ave_byCluster <- AverageExpression(seuratObj,
      features = features,
      assays = "RNA", slot = "data", verbose = F
    )$RNA

    # log1p
    mat <- log1p(as.matrix(var_ave_byCluster))
    mat <- pheatmap:::scale_rows(mat)
    # print(mat)
    # remove NaN
    mat <- mat[complete.cases(mat), ]
    return(rownames(mat))
  }

  if (typeof(features) == "character") {
    .helper_remove_nan(features = features, seurat_obj = seurat_obj)
  } else {
    lapply(features, FUN = .helper_remove_nan, seurat_obj = seurat_obj)
  }
}

customDoHeatmap <- function(seurat_obj = NULL,
                            features = NULL,
                            z_score = FALSE,
                            cluster_columns = TRUE,
                            cluster_rows = TRUE,
                            show_row_names = TRUE, show_column_names = TRUE,
                            clustering_method_rows = "average",
                            clustering_method_cols = "average",
                            idents = "seurat_clusters",
                            assay = "RNA",
                            slot = "data") {
  # set idents
  seuratObj <- seurat_obj
  Idents(seuratObj) <- idents
  # get average exp
  var_ave_byCluster <- AverageExpression(seuratObj,
    features = unlist(features),
    assays = assay, slot = slot, verbose = F
  )$RNA

  # log1p
  mat <- log1p(as.matrix(var_ave_byCluster))

  if (z_score) {
    mat <- pheatmap:::scale_rows(mat)
    # print(mat)
    # remove NaN
    mat <- mat[complete.cases(mat), ]
    x <- ceiling(max(abs(mat)))
    col_fun <- circlize::colorRamp2(c(-x, 0, x), c("#476fa9", "#ffffff", "#ca3226"))
    legend_name <- "z-score"
  } else {
    x <- ceiling(max(abs(mat)))
    col_fun <- circlize::colorRamp2(c(0, x), c("#ffffff", "#ca3226"))
    legend_name <- "Abs. exp (log1p)"
  }

  # draw
  if (typeof(features) == "character") {
    ht <- ComplexHeatmap::Heatmap(mat,
      name = legend_name,
      col = col_fun,
      cluster_columns = cluster_columns, cluster_rows = cluster_rows,
      show_row_dend = FALSE, show_column_dend = FALSE,
      row_names_side = "right",
      row_names_gp = gpar(fontsize = 6, fontface = "italic"),
      column_names_gp = gpar(fontsize = 6),
      column_names_rot = 60,
      row_title_gp = gpar(fontsize = 6),
      heatmap_legend_param = list(
        # title_position = "topcenter",
        # direction = 'horizontal',
        labels_gp = gpar(fontsize = 6),
        title_gp = gpar(fontsize = 6),
        grid_height = unit(2, "mm"),
        grid_width = unit(2, "mm")
      )
    )
  } else {
    ht <- ComplexHeatmap::Heatmap(mat,
      name = legend_name,
      col = col_fun,
      cluster_columns = cluster_columns, cluster_rows = cluster_rows,
      show_row_dend = FALSE, show_column_dend = FALSE,
      show_row_names = show_row_names, show_column_names = show_column_names,
      row_names_side = "right",
      row_names_gp = gpar(fontsize = 6, fontface = "italic"),
      column_names_gp = gpar(fontsize = 6),
      column_names_rot = 60,
      row_title_gp = gpar(fontsize = 6),
      row_split = factor(
        rep(names(features), sapply(features, length)),
        levels = names(features)
      ),
      heatmap_legend_param = list(
        # title_position = "topcenter",
        # direction = 'horizontal',
        labels_gp = gpar(fontsize = 6),
        title_gp = gpar(fontsize = 6),
        grid_height = unit(2, "mm"),
        grid_width = unit(2, "mm")
      )
    )
  }

  # draw
  # ht_list = ht1 + ht2
  return(ht)
}

#####################################################

# gene_BoxPlot2

#####################################################

# plot gene exresssion using vlnplot: one cluster, multiple genes
gene_BoxPlot2 <- function(seurat_obj, genes) {
  require(Seurat)
  require(ggplot2)
  # Extract the UMAP coordinates for each cell and include information about genes to plot
  data_4_p <- FetchData(seurat_obj,
    vars = genes
  )

  data_4_p <- data_4_p %>%
    gather() %>%
    as.tibble() %>%
    mutate(key = fct_relevel(key, levels = genes))

  # Plot
  ggplot(data_4_p, aes(x = key, y = value)) +
    geom_boxplot(outlier.alpha = 0.1) +
    labs(y = "Expression Level") +
    cowplot::theme_cowplot() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      axis.title.x = element_blank()
    )
}

# plot feature distribution by boxplot of each cluster
feature_box <- function(x = NULL, seurat_obj = NULL, groupby = "seurat_clusters") {
  seurat_obj@meta.data %>%
    group_by(!!sym(groupby), !!sym(x)) %>%
    summarise() %>%
    ggplot(aes_string(groupby, x)) +
    stat_boxplot(geom = "errorbar", width = 0.4) +
    geom_boxplot(aes(fill = "blue"), outlier.alpha = 0.1, color = "#00008f") +
    scale_fill_manual(values = "#8f8fb1") +
    theme(
      axis.title.y = element_blank(),
      legend.position = "none"
    ) +
    coord_flip()
}

# Percent of VDJ of each cluster
percent_vdj <- function(seurat_obj, groupby = "seurat_clusters") {
  plot_grid(
    seurat_obj@meta.data %>%
      group_by(!!sym(groupby), VDJ) %>%
      summarise(n = n()) %>%
      mutate(Percent = n / sum(n) * 100) %>%
      ggplot(aes(x = !!sym(groupby), y = Percent, fill = VDJ)) +
      geom_bar(stat = "identity", color = "white") +
      scale_fill_manual(values = vdj_color_maps) +
      coord_flip() +
      theme(legend.position = "top") +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
      scale_y_continuous(expand = c(0, 0)),
    seurat_obj@meta.data %>%
      group_by(!!sym(groupby), IGHC) %>%
      summarise(n = n()) %>%
      mutate(Percent = n / sum(n) * 100) %>%
      ggplot(aes(x = !!sym(groupby), y = Percent, fill = IGHC)) +
      geom_bar(stat = "identity", color = "white") +
      scale_fill_manual(values = ig_color_maps) +
      coord_flip() +
      theme(legend.position = "top") +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
      scale_y_continuous(expand = c(0, 0)),
    ncol = 2
  )
}


# compositions of each cluster
# Origin2, Patient, cell num
cluster_stack <- function(seurat_obj = NULL,
                          group_by = NULL) {
  # check groupby
  if (is.null(group_by)) {
    stop("`group_by` is NULL.")
  }
  if (length(group_by) > 1) {
    stop("multiple `group_by` is currently not allowed.")
  }
  if (!(group_by %in% colnames(seurat_obj@meta.data))) {
    stop("`group_by` not found in seurat obj.")
  }
  cluster_byOrigin2 <- seurat_obj@meta.data %>%
    group_by(!!sym(group_by), Origin2_n) %>%
    summarise(n = n()) %>%
    mutate(Percent = n / sum(n) * 100)
  cluster_byPatient <- seurat_obj@meta.data %>%
    group_by(!!sym(group_by), Patient) %>%
    summarise(n = n()) %>%
    mutate(Percent = n / sum(n) * 100)
  cluster_cellNum <- seurat_obj@meta.data %>%
    group_by(!!sym(group_by)) %>%
    summarise(n = n())

  p1 <- ggplot(cluster_byOrigin2, aes(x = !!sym(group_by), y = Percent, fill = Origin2_n)) +
    geom_bar(stat = "identity", color = "white") +
    scale_fill_manual(values = origin2_color_maps) +
    labs(fill = NULL) +
    coord_flip() +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
    scale_y_continuous(expand = c(0, 0))
  p2 <- ggplot(cluster_byPatient, aes(x = !!sym(group_by), y = Percent, fill = Patient)) +
    geom_bar(stat = "identity", color = "white") +
    scale_fill_manual(values = patient_color_maps) +
    labs(fill = NULL) +
    coord_flip() +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank()
    )
  p3 <- ggplot(cluster_cellNum, aes(x = !!sym(group_by), y = log2(n), fill = "blue")) +
    geom_bar(stat = "identity", color = "#00008f") +
    scale_fill_manual(values = "#8f8fb1") +
    labs(fill = NULL) +
    coord_flip() +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none"
    )

  p1 + p2 + p3 +
    plot_layout(ncol = 3) & theme(plot.margin = margin(0))
}

# Tissue, Origin3, Origin
cluster_stack2 <- function(seurat_obj = NULL,
                           group_by = NULL) {
  # check groupby
  if (is.null(group_by)) {
    stop("`group_by` is NULL.")
  }
  if (length(group_by) > 1) {
    stop("multiple `group_by` is currently not allowed.")
  }
  if (!(group_by %in% colnames(seurat_obj@meta.data))) {
    stop("`group_by` not found in seurat obj.")
  }
  cluster_byTissue <- seurat_obj@meta.data %>%
    group_by(!!sym(group_by), Tissue) %>%
    summarise(n = n()) %>%
    mutate(Percent = n / sum(n) * 100)
  cluster_byOrigin3 <- seurat_obj@meta.data %>%
    group_by(!!sym(group_by), Origin3) %>%
    summarise(n = n()) %>%
    mutate(Percent = n / sum(n) * 100)
  cluster_byOrigin <- seurat_obj@meta.data %>%
    group_by(!!sym(group_by), Origin) %>%
    summarise(n = n()) %>%
    mutate(Percent = n / sum(n) * 100)

  p1 <- ggplot(cluster_byTissue, aes(x = !!sym(group_by), y = Percent, fill = Tissue)) +
    geom_bar(stat = "identity", color = "white") +
    scale_fill_manual(values = tissue_color_maps) +
    labs(fill = NULL) +
    coord_flip() +
    theme(legend.position = "right") +
    scale_y_continuous(expand = c(0, 0))
  p2 <- ggplot(cluster_byOrigin3, aes(x = !!sym(group_by), y = Percent, fill = Origin3)) +
    geom_bar(stat = "identity", color = "white") +
    scale_fill_manual(values = origin3_color_maps) +
    labs(fill = NULL) +
    coord_flip() +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      legend.position = "right",
      axis.text.y = element_blank(),
      axis.title.y = element_blank()
    )
  p3 <- ggplot(cluster_byOrigin, aes(x = !!sym(group_by), y = Percent, fill = Origin)) +
    geom_bar(stat = "identity", color = "white") +
    scale_fill_manual(values = origin_color_maps) +
    labs(fill = NULL) +
    coord_flip() +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      legend.position = "right",
      axis.text.y = element_blank(),
      axis.title.y = element_blank()
    )

  p1 + p2 + p3 +
    plot_layout(ncol = 3) & theme(plot.margin = margin(0))
}


#####################################################
# plot cell fraction of each sample
# bar plot
#####################################################

PlotSampleFraction <- function(
    seurat_obj = NULL,
    groupby = "seurat_clusters",
    color_maps = cluster_color_maps,
    pos_x = NULL) {
  require(patchwork)
  require(tidyverse)
  # get dat for plot
  df <- seurat_obj@meta.data %>%
    group_by(Source, !!sym(groupby)) %>%
    summarise(n = n()) %>%
    mutate(Percent = n / sum(n) * 100) %>%
    ungroup()

  # add space between patient
  if (groupby == "cellType") {
    tmp_df <- tibble(Source = letters[1:5], cellType = seurat_obj$cellType[1], n = 0, Percent = 0)
  } else if (groupby == "seurat_clusters") {
    tmp_df <- tibble(Source = letters[1:5], seurat_clusters = seurat_obj$seurat_clusters[1], n = 0, Percent = 0)
  } else if (groupby == "level_1") {
    tmp_df <- tibble(Source = letters[1:5], level_1 = seurat_obj$level_1[1], n = 0, Percent = 0)
  } else if (groupby == "level_2") {
    tmp_df <- tibble(Source = letters[1:5], level_2 = seurat_obj$level_2[1], n = 0, Percent = 0)
  } else if (groupby == "level_3") {
    tmp_df <- tibble(Source = letters[1:5], level_3 = seurat_obj$level_3[1], n = 0, Percent = 0)
  }
  df2 <- rbind(df, tmp_df)
  df2$Source <- factor(df2$Source, levels = c(
    "S0619.P1", "S0619.P2", "S0619.LN1", "S0619.LN2", "S0619.LN3", "S0619.LN4", "S0619.LN5", "S0619.LN6", "S0619.N1", "S0619.N2", "S0619.N3", "S0619.N4", "S0619.N5",
    "a",
    "S0730.P1", "S0730.P2", "S0730.LN1", "S0730.LN2", "S0730.LN3", "S0730.LN4", "S0730.N1", "S0730.N2", "S0730.N3", "S0730.N4", "S0730.N5",
    "b",
    "S0819.P1", "S0819.P2", "S0819.LN1", "S0819.LN2", "S0819.LN3", "S0819.LN4", "S0819.LN5", "S0819.N1", "S0819.N2", "S0819.N3", "S0819.N4", "S0819.N5",
    "c",
    "S0920.P1", "S0920.P2", "S0920.LN1", "S0920.LN2", "S0920.LN3", "S0920.LN4", "S0920.N1", "S0920.N2", "S0920.N3", "S0920.N4", "S0920.N5",
    "d",
    "S1125.P1", "S1125.P2", "S1125.LN1", "S1125.LN2", "S1125.LN3", "S1125.LN4", "S1125.LN5", "S1125.LN6", "S1125.N1", "S1125.N2", "S1125.N3", "S1125.N4", "S1125.N5",
    "e",
    "S1204.P1", "S1204.P2", "S1204.LN1", "S1204.LN2", "S1204.LN3", "S1204.LN5", "S1204.LN6", "S1204.LN7", "S1204.N1", "S1204.N2", "S1204.N3", "S1204.N4", "S1204.N5"
  ))

  # for labels and breaks
  df <- df %>%
    left_join(seurat_obj@misc$sam_info, by = "Source")

  # plot cell fraction
  p.bottom <- df2 %>%
    ggplot(aes(x = Source, y = Percent, fill = !!sym(groupby))) +
    geom_col(color = "white") +
    scale_fill_manual(values = color_maps) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(breaks = df$Source, labels = df$Origin) +
    labs(fill = NULL, x = NULL, y = "Proportion (%)") +
    theme_cowplot(font_size = 5) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

  # text annotation y location
  text_y <- (df %>% group_by(Source) %>% summarise(n = sum(n)) %>% pull(n) %>% max()) * .9
  # pos x
  if (is.null(pos_x)) {
    pos_x <- c(6, 20, 31, 45, 58, 72)
  }
  p.top <- df2 %>%
    ggplot(aes(x = Source, y = n, fill = !!sym(groupby))) +
    geom_col(color = "white") +
    scale_fill_manual(values = color_maps) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(breaks = df$Source, labels = df$Origin) +
    labs(fill = NULL, x = NULL, y = "Number") +
    theme_cowplot(font_size = 5) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    annotate("text", x = pos_x, y = text_y, label = c("S0619", "S0730", "S0819", "S0930", "S1125", "S1204"))

  p.top / p.bottom +
    plot_layout(ncol = 1, guides = "collect")
}


#####################################################
# plot pie chart of proportion of clusters in Tissue/Patient
#####################################################

# plot pie chart for tissue/metastasis/origin
PlotFrac <- function(
    origin = NULL,
    seurat_obj = NULL,
    plot_what = NULL,
    groupby = "seurat_clusters",
    color_maps = cluster_color_maps) {
  seurat_obj@meta.data %>%
    filter(!!sym(plot_what) == origin) %>%
    group_by(!!sym(groupby)) %>%
    summarise(n = n()) %>%
    mutate(Percent = n / sum(n) * 100) %>%
    ggplot(aes("", Percent, fill = !!sym(groupby))) +
    geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
    coord_polar("y") +
    geom_text_repel(aes(label = paste0(round(Percent, 1), "%")),
      position = position_stack(vjust = 0.5)
    ) +
    labs(x = NULL, y = NULL, fill = NULL, title = origin) +
    scale_fill_manual(values = color_maps) +
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}

# plot pie chart for patient
PlotFracPatient <- function(
    patient = NULL,
    seurat_obj = NULL,
    groupby = "seurat_clusters",
    color_maps = cluster_color_maps) {
  seurat_obj@meta.data %>%
    filter(Patient == patient) %>%
    group_by(!!sym(groupby)) %>%
    summarise(n = n()) %>%
    mutate(Percent = n / sum(n) * 100) %>%
    ggplot(aes("", Percent, fill = !!sym(groupby))) +
    geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
    coord_polar("y") +
    geom_text_repel(aes(label = paste0(round(Percent, 1), "%")),
      position = position_stack(vjust = 0.5)
    ) +
    labs(x = NULL, y = NULL, fill = NULL, title = patient) +
    scale_fill_manual(values = color_maps) +
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}

# plot pie chart for origin
PlotFracOrigin <- function(
    origin = NULL,
    seurat_obj = NULL,
    groupby = "seurat_clusters",
    color_maps = cluster_color_maps) {
  seurat_obj@meta.data %>%
    filter(Origin == origin) %>%
    group_by(!!sym(groupby)) %>%
    summarise(n = n()) %>%
    mutate(Percent = n / sum(n) * 100) %>%
    ggplot(aes("", Percent, fill = !!sym(groupby))) +
    geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
    coord_polar("y") +
    geom_text_repel(aes(label = paste0(round(Percent, 1), "%")),
      position = position_stack(vjust = 0.5)
    ) +
    labs(x = NULL, y = NULL, fill = NULL, title = origin) +
    scale_fill_manual(values = color_maps) +
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}

# plot pie chart for Metastasis
PlotFracMetastasis <- function(
    metastasis = NULL,
    seurat_obj = NULL,
    groupby = "seurat_clusters",
    color_maps = cluster_color_maps) {
  seurat_obj@meta.data %>%
    filter(Metastasis == metastasis) %>%
    group_by(!!sym(groupby)) %>%
    summarise(n = n()) %>%
    mutate(Percent = n / sum(n) * 100) %>%
    ggplot(aes("", Percent, fill = !!sym(groupby))) +
    geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
    coord_polar("y") +
    geom_text_repel(aes(label = paste0(round(Percent, 1), "%")),
      position = position_stack(vjust = 0.5)
    ) +
    labs(x = NULL, y = NULL, fill = NULL, title = metastasis) +
    scale_fill_manual(values = color_maps) +
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}

# plot pie chart for Drainage
PlotFracDrainage <- function(
    metastasis = NULL,
    seurat_obj = NULL,
    groupby = "seurat_clusters",
    color_maps = cluster_color_maps) {
  seurat_obj@meta.data %>%
    filter(Drainage == metastasis) %>%
    group_by(!!sym(groupby)) %>%
    summarise(n = n()) %>%
    mutate(Percent = n / sum(n) * 100) %>%
    ggplot(aes("", Percent, fill = !!sym(groupby))) +
    geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
    coord_polar("y") +
    geom_text_repel(aes(label = paste0(round(Percent, 1), "%")),
      position = position_stack(vjust = 0.5)
    ) +
    labs(x = NULL, y = NULL, fill = NULL, title = metastasis) +
    scale_fill_manual(values = color_maps) +
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}

# plot pie chart for Tissue
PlotFracTissue <- function(
    tissue = NULL,
    seurat_obj = NULL,
    groupby = "seurat_clusters",
    color_maps = cluster_color_maps) {
  seurat_obj@meta.data %>%
    filter(Tissue == tissue) %>%
    group_by(!!sym(groupby)) %>%
    summarise(n = n()) %>%
    mutate(Percent = n / sum(n) * 100) %>%
    ggplot(aes("", Percent, fill = !!sym(groupby))) +
    geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
    coord_polar("y") +
    geom_text_repel(aes(label = paste0(round(Percent, 1), "%")),
      position = position_stack(vjust = 0.5)
    ) +
    labs(x = NULL, y = NULL, fill = NULL, title = tissue) +
    scale_fill_manual(values = color_maps) +
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}

# plot pie chart for Source
PlotFracSource <- function(
    patient = NULL,
    seurat_obj = NULL,
    groupby = "seurat_clusters",
    color_maps = cluster_color_maps) {
  seurat_obj@meta.data %>%
    filter(Patient == patient) %>%
    group_by(Source, !!sym(groupby)) %>%
    summarise(n = n()) %>%
    mutate(Percent = n / sum(n) * 100) %>%
    ggplot(aes(Source, Percent, fill = !!sym(groupby))) +
    geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
    geom_text_repel(aes(label = paste0(round(Percent, 1), "%")),
      position = position_stack(vjust = 0.5)
    ) +
    labs(x = NULL, y = NULL, fill = NULL, title = patient) +
    scale_fill_manual(values = color_maps) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      axis.line = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}



#####################################################
# plot fraction changes of cells
#####################################################

PlotFracChange <- function(
    dat = NULL,
    x = NULL,
    title = NULL,
    comparisons = NULL,
    facet_by = NULL) {
  library(ggpubr)
  # cmp
  if (x == "Origin2_n") {
    comparisons <- list(
      c("prePBMC", "postPBMC"),
      c("nLN", "pLN"),
      c("Normal", "Adjacent"),
      c("Adjacent", "Tumor"),
      c("Normal", "Tumor"),
      c("pLN", "Tumor")
    )
    a <- dat %>%
      pull(Origin2_n) %>%
      unique()
    b <- sapply(comparisons, function(x) {
      length(intersect(x, a))
    })
    comparisons <- comparisons[b == 2]
  } else if (x == "Origin4") {
    comparisons <- list(
      c("prePBMC", "postPBMC"),
      c("nLN", "m_nLN"),
      c("m_nLN", "m_pLN"),
      c("nLN", "m_pLN"),
      c("Normal", "Adjacent"),
      c("Adjacent", "Tumor"),
      c("Normal", "Tumor"),
      c("m_pLN", "Tumor")
    )
    a <- dat %>%
      pull(Origin4) %>%
      unique()
    b <- sapply(comparisons, function(x) {
      length(intersect(x, a))
    })
    comparisons <- comparisons[b == 2]
  }
  p <- dat %>%
    ggplot(aes_string(x = x, y = "value", color = "Patient", group = x)) +
    geom_point(position = position_jitter(width = .1)) +
    stat_summary(
      fun.data = mean_sdl, fun.args = list(mult = 1),
      geom = "pointrange", color = "black",
      size = .4
    ) +
    stat_compare_means(
      comparisons = comparisons,
      method = "wilcox.test",
      size = 6 / (14 / 5)
    ) +
    #stat_n_text(size = 6 / (14 / 5)) +
    scale_color_manual(values = patient_color_maps) +
    labs(x = NULL, y = "Proportion (%)", title = title)

  if (is.null(facet_by)) {
    p +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6)
      )
  } else {
    p +
      facet_wrap(vars(!!sym(facet_by))) +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6)
      )
  }
}

#####################################################
# Dotplot for average expression of genes
#####################################################

# Dotplot for average expression of genes
AveDotPlot <- function(object, assay = NULL, features, cols = c(
                         "lightgrey",
                         "blue"
                       ),
                       col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
                       group.by = NULL, split.by = NULL, scale.by = "radius", scale.min = NA,
                       scale.max = NA) {
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  scale.func <- switch(EXPR = scale.by,
    size = scale_size,
    radius = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  data.features <- sapply(features, function(vars) {
    FetchData(object = object, vars = vars) %>%
      rowMeans()
  }) %>%
    as.data.frame()
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)
  } else {
    object[[group.by, drop = TRUE]]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  # print(data.features)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]]
    if (length(x = unique(x = splits)) > length(x = cols)) {
      stop("Not enought colors for the number of groups")
    }
    cols <- cols[1:length(x = unique(x = splits))]
    names(x = cols) <- unique(x = splits)
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(
      rep(x = id.levels, each = length(x = unique.splits)),
      "_", rep(x = unique(x = splits), times = length(x = id.levels))
    )
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident,
      1:(ncol(x = data.features) - 1),
      drop = FALSE
    ]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(
      X = data.use, MARGIN = 2, FUN = Seurat:::PercentAbove,
      threshold = 0
    )
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  avg.exp.scaled <- sapply(
    X = unique(x = data.plot$features.plot),
    FUN = function(x) {
      data.use <- data.plot[data.plot$features.plot ==
        x, "avg.exp"]
      data.use <- scale(x = data.use)
      data.use <- MinMax(
        data = data.use, min = col.min,
        max = col.max
      )
      return(data.use)
    }
  )
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (!is.null(x = split.by)) {
    avg.exp.scaled <- as.numeric(x = cut(
      x = avg.exp.scaled,
      breaks = 20
    ))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = rev(x = names(features))
  )
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (!is.null(x = split.by)) {
    splits.use <- vapply(
      X = strsplit(
        x = as.character(x = data.plot$id),
        split = "_"
      ), FUN = "[[", FUN.VALUE = character(length = 1L),
      2
    )
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled",
    no = "colors"
  )
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  # print(data.plot)
  plot <- ggplot(data = data.plot, mapping = aes_string(
    x = "features.plot",
    y = "id"
  )) +
    geom_point(mapping = aes_string(
      size = "pct.exp",
      color = color.by
    )) +
    scale.func(
      range = c(0, dot.scale),
      limits = c(scale.min, scale.max)
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    guides(size = guide_legend(title = "Percent Expressed")) +
    labs(x = "Features", y = ifelse(test = is.null(x = split.by),
      yes = "Identity", no = "Split Identity"
    )) +
    theme_cowplot()
  if (!is.null(x = split.by)) {
    plot <- plot + scale_color_identity()
  } else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  } else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (is.null(x = split.by)) {
    plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
  }
  return(plot)
}

# plot enrichment analysis results by barplot
vis_enrich_bar <- function(ek, eBP, show_n = 10) {
  plot_grid(
    ek %>%
      as.data.frame() %>%
      arrange(Count) %>%
      mutate(Description = factor(Description, levels = Description)) %>%
      top_n(show_n, wt = Count) %>%
      ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
      geom_bar(stat = "identity") +
      scale_fill_gradient(low = "red", high = "blue") +
      geom_hline(yintercept = 0) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(title = "KEGG") +
      coord_flip() +
      geom_text(size = 4, aes(y = 0.5, label = Description, hjust = 0)) +
      theme(
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank()
      ),
    eBP %>%
      as.data.frame() %>%
      arrange(Count) %>%
      mutate(Description = factor(Description, levels = Description)) %>%
      top_n(show_n, wt = Count) %>%
      ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
      geom_bar(stat = "identity") +
      scale_fill_gradient(low = "red", high = "blue") +
      geom_hline(yintercept = 0) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(title = "GO BP") +
      coord_flip() +
      geom_text(size = 4, aes(y = 0.5, label = Description, hjust = 0)) +
      theme(
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank()
      ),
    ncol = 2
  )
}

## get a matchSCore2-like plot for SingleR's matchReferences output
summary_ggplot <- function(data, ylab = "", xlab = "", geom.text.size = 6 / (14 / 5)) {
  my_df <- data.frame(data, check.names = FALSE, check.rows = FALSE)
  my_df.melt <- reshape2::melt(cbind(x = 1:nrow(my_df), my_df), id = "x")
  my_df.melt$x <- factor(my_df.melt$x)

  gg <- ggplot(my_df.melt, aes_string(x = "x", y = "variable", fill = "value")) +
    labs(x = xlab, y = ylab) +
    geom_tile(aes_string(fill = "value")) +
    scale_x_discrete(lab = rownames(my_df)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
    geom_text(aes(label = round(value, 2)), size = geom.text.size) +
    scale_fill_gradient(low = "white", high = "red", name = "Prob.") +
    coord_flip()

  return(gg)
}

#####################################################
# plot tissue enrich
#####################################################

# helper function
.plot_enrich <- function(data) {
  gg <- ggplot(data, aes_string(x = "x", y = "variable", fill = "value")) +
    geom_tile(aes_string(fill = "value")) +
    geom_text(aes(label = round(value, 1)), size = 6 / (14 / 5)) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "right") +
    scale_fill_gradient2(
      midpoint = 1,
      low = "blue", mid = "white", high = "red",
      name = expression(R[O / E])
    ) +
    theme_cowplot(font_size = 8) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 0),
      axis.title = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right",
      text = element_text(size = 6)
    )

  return(gg)
}

# tissue
plot_tissue_enrich <- function(data) {
  my_df <- data.frame(t(data), check.names = FALSE, check.rows = FALSE)
  my_df.melt <- reshape2::melt(cbind(x = rownames(my_df), my_df), id = "x")
  my_df.melt$x <- factor(my_df.melt$x, levels = c("PBMC", "LN", "Normal", "Adjacent", "Tumor"))

  return(.plot_enrich(my_df.melt))
}

# Origin
plot_origin_enrich <- function(data) {
  my_df <- data.frame(t(data), check.names = FALSE, check.rows = FALSE)
  my_df.melt <- reshape2::melt(cbind(x = rownames(my_df), my_df), id = "x")
  my_df.melt$x <- factor(my_df.melt$x, levels = c(
    "PBMC1", "PBMC2",
    "LN-LR", "LN-RR", "LN-UTP1", "LN-UTP2", "LN-ITP", "LN-LTP", "LN-LP", "LN-LP1", "LN-LP2", "LN-LG",
    "Normal", "Adjacent",
    "Tumor_core", "Tumor_invasion", "Tumor_middle"
  ))

  return(.plot_enrich(my_df.melt))
}

# metastasis
plot_metastasis_enrich <- function(data) {
  my_df <- data.frame(t(data), check.names = FALSE, check.rows = FALSE)
  my_df.melt <- reshape2::melt(cbind(x = rownames(my_df), my_df), id = "x")
  my_df.melt$x <- factor(my_df.melt$x, levels = c("N", "P"))

  return(.plot_enrich(my_df.melt))
}

# drainage
plot_drainage_enrich <- function(data) {
  my_df <- data.frame(t(data), check.names = FALSE, check.rows = FALSE)
  my_df.melt <- reshape2::melt(cbind(x = rownames(my_df), my_df), id = "x")
  my_df.melt$x <- factor(my_df.melt$x, levels = c("Up", "Down"))

  return(.plot_enrich(my_df.melt))
}

# Origin2
plot_origin2_enrich <- function(data) {
  my_df <- data.frame(t(data), check.names = FALSE, check.rows = FALSE)
  my_df.melt <- reshape2::melt(cbind(x = rownames(my_df), my_df), id = "x")
  my_df.melt$x <- factor(my_df.melt$x, levels = c(
    "prePBMC", "postPBMC",
    "nLN", "pLN",
    "Normal", "Adjacent",
    "Tumor"
  ))

  return(.plot_enrich(my_df.melt))
}

# Origin3
plot_origin3_enrich <- function(data) {
  my_df <- data.frame(t(data), check.names = FALSE, check.rows = FALSE)
  my_df.melt <- reshape2::melt(cbind(x = rownames(my_df), my_df), id = "x")
  my_df.melt$x <- factor(my_df.melt$x, levels = c(
    "PBMC1", "PBMC2",
    "LN_Up", "LN_Down",
    "Normal", "Adjacent",
    "Tumor"
  ))

  return(.plot_enrich(my_df.melt))
}


# PBMC1/PBMC2
plot_pbmc_enrich <- function(data) {
  my_df <- data.frame(t(data), check.names = FALSE, check.rows = FALSE)
  my_df.melt <- reshape2::melt(cbind(x = rownames(my_df), my_df), id = "x")
  my_df.melt$x <- factor(my_df.melt$x, levels = c("PBMC1", "PBMC2"))

  return(.plot_enrich(my_df.melt))
}


#####################################################
# cellphonedb interaction count
#####################################################
plot_heatmaps_cpdb <- function(col_cell_dir = "epi",
                               col_cell = "epi",
                               row_cell_dir = "cd4",
                               row_cell = "CD4",
                               ann_col = NA,
                               col_ann = NA,
                               cpdb_out_dir = "/niuyw-usb-disk/wangyx/cellphoneDB/scESCA/output") {
  require(dplyr)
  inter_count_file <- file.path(
    cpdb_out_dir,
    paste(col_cell_dir, row_cell_dir, sep = "_"),
    "heatmap_count.txt"
  )
  # file exist
  if (!file.exists(inter_count_file)) {
    stop("heatmap_count.txt not exist")
  }
  inter_count <- readr::read_delim(inter_count_file,
    delim = "\t",
    col_names = FALSE, skip = 1,
    progress = FALSE
  ) %>%
    column_to_rownames("X1") %>%
    as.data.frame()
  colnames(inter_count) <- rownames(inter_count)
  # print(inter_count)
  # select
  # col
  if (length(col_cell) == 1) {
    inter_count <- inter_count %>%
      rownames_to_column() %>%
      dplyr::select(rowname, starts_with(col_cell)) %>%
      column_to_rownames()
  } else {
    inter_count <- inter_count[, col_cell]
  }
  # row
  if (length(row_cell) == 1) {
    inter_count <- inter_count %>%
      rownames_to_column() %>%
      filter(str_detect(rowname, row_cell)) %>%
      column_to_rownames() %>%
      as.matrix()
  } else {
    inter_count <- inter_count[row_cell, ] %>%
      as.matrix()
  }

  # plot
  col1 <- "dodgerblue4"
  col2 <- "peachpuff"
  col3 <- "deeppink4"
  col.heatmap <- colorRampPalette(c(col1, col2, col3))(1000)
  pheatmap::pheatmap(inter_count,
    clustering_method = "average",
    color = col.heatmap,
    border_color = "white", fontsize_row = 11, fontsize_col = 11,
    annotation_col = ann_col, annotation_colors = col_ann, family = "Arial"
  )
}



#####################################################
# Diffusion map anslysis
#####################################################

# since Seurat::FeaturePlot works wrongly on DM
FeaturePlotDM <- function(seurat_obj = NULL, gene = NULL, dims = c(1, 2)) {
  dims <- paste0("DC_", dims)
  vars <- c(gene, dims)
  p <- FetchData(seurat_obj, vars = vars, slot = "data") %>%
    rownames_to_column("barcode") %>%
    ggplot(aes_string(x = dims[1], y = dims[2], color = gene)) +
    geom_point() +
    scale_color_gradient(low = "lightgrey", high = "blue") +
    labs(title = gene)
  return(p)
}

# traceplot of exp along DC
TracePlotDM <- function(seurat_obj = NULL, gene = NULL, dims = 1) {
  if (dims == 1) {
    df <- FetchData(seurat_obj, vars = c(gene, "DC_1"), slot = "data") %>%
      rownames_to_column("barcode") %>%
      arrange(DC_1)
  } else if (dims == 2) {
    df <- FetchData(seurat_obj, vars = c(gene, "DC_2"), slot = "data") %>%
      rownames_to_column("barcode") %>%
      arrange(DC_2)
  } else if (dims == 3) {
    df <- FetchData(seurat_obj, vars = c(gene, "DC_3"), slot = "data") %>%
      rownames_to_column("barcode") %>%
      arrange(DC_3)
  }

  a <- zoo::rollmean(df[, gene], dim(df)[1] / 20)
  b <- zoo::rollapply(df[, gene], width = dim(df)[1] / 20, FUN = sd)
  # print(a)
  tibble(value = a, n = 1:length(a), sd = b) %>%
    mutate(
      lower = value - sd,
      upper = value + sd
    ) %>%
    ggplot(aes(x = n, y = value)) +
    geom_line(color = "blue") +
    geom_ribbon(aes(ymin = lower, ymax = upper, color = NULL), alpha = .5, fill = "#dbd8ec") +
    labs(
      x = paste0("Cell along DC", dims),
      y = "Expression",
      title = gene
    )
}


#####################################################
# scatterplot3d
#####################################################

# from: http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r
#' Add grids to a scatterplot3d
#'
#' @description The goal of this function is to add grids on an existing
#'  plot created using the package scatterplot3d
#' @param x,y,z numeric vectors specifying the x, y, z coordinates of points.
#'  x can be a matrix or a data frame containing 3 columns corresponding to
#'  the x, y and z coordinates. In this case the arguments y and z are optional
#' @param grid specifies the facet(s) of the plot on which grids should be drawn.
#'  Possible values are the combination of "xy", "xz" or "yz".
#'  Example: grid = c("xy", "yz"). The default value is TRUE to add grids only on xy facet.
#' @param col.grid,lty.grid color and line type to be used for grids
#' @param lab a numerical vector of the form c(x, y, len).
#'  The values of x and y give the (approximate) number of tickmarks on the x and y axes.
#' @param lab.z the same as lab, but for z axis
#' @param scale.y of y axis related to x- and z axis
#' @param angle angle between x and y axis
#' @param "xlim, ylim, zlim" the x, y and z limits (min, max) of the plot.
#'
#' @note
#' Users who want to extend an existing scatterplot3d graphic with the
#'  function addgrids3d, should consider to set the arguments scale.y, angle, ...,
#'  to the value used in scatterplot3d.
#'
#' @author Alboukadel Kassambara \email{alboukadel.kassambara@@gmail.com}
#' @references http://www.sthda.com
#'
#' @example
#' library(scatterplot3d)
#' data(iris)
#' scatterplot3d(iris[, 1:3], pch = 16, grid=T, box=F)
#' addgrids3d(iris[, 1:3], grid = c("xy", "xz", "yz"))
addgrids3d <- function(x, y = NULL, z = NULL, grid = TRUE,
                       col.grid = "grey", lty.grid = par("lty"),
                       lab = par("lab"), lab.z = mean(lab[1:2]),
                       scale.y = 1, angle = 40,
                       xlim = NULL, ylim = NULL, zlim = NULL) {
  if (inherits(x, c("matrix", "data.frame"))) {
    x <- as.data.frame(x)
    y <- unlist(x[, 2])
    z <- unlist(x[, 3])
    x <- unlist(x[, 1])
  }

  p.lab <- par("lab")

  angle <- (angle %% 360) / 90
  yz.f <- scale.y * abs(if (angle < 1) angle else if (angle > 3) angle - 4 else 2 - angle)
  yx.f <- scale.y * (if (angle < 2) 1 - angle else angle - 3)


  # x axis range
  x.range <- range(x[is.finite(x)], xlim)
  x.prty <- pretty(x.range, n = lab[1], min.n = max(1, min(0.5 * lab[1], p.lab[1])))
  x.scal <- round(diff(x.prty[1:2]), digits = 12)
  x <- x / x.scal
  x.range <- range(x.prty) / x.scal
  x.max <- ceiling(x.range[2])
  x.min <- floor(x.range[1])
  if (!is.null(xlim)) {
    x.max <- max(x.max, ceiling(xlim[2] / x.scal))
    x.min <- min(x.min, floor(xlim[1] / x.scal))
  }
  x.range <- range(x.min, x.max)

  # y axis range
  y.range <- range(y[is.finite(y)], ylim)
  y.prty <- pretty(y.range, n = lab[2], min.n = max(1, min(0.5 * lab[2], p.lab[2])))
  y.scal <- round(diff(y.prty[1:2]), digits = 12)
  y.add <- min(y.prty)
  y <- (y - y.add) / y.scal
  y.max <- (max(y.prty) - y.add) / y.scal
  if (!is.null(ylim)) {
    y.max <- max(y.max, ceiling((ylim[2] - y.add) / y.scal))
  }

  # Z axis range
  z.range <- range(z[is.finite(z)], zlim)
  z.prty <- pretty(z.range, n = lab.z, min.n = max(1, min(0.5 * lab.z, p.lab[2])))
  z.scal <- round(diff(z.prty[1:2]), digits = 12)
  z <- z / z.scal
  z.range <- range(z.prty) / z.scal
  z.max <- ceiling(z.range[2])
  z.min <- floor(z.range[1])
  if (!is.null(zlim)) {
    z.max <- max(z.max, ceiling(zlim[2] / z.scal))
    z.min <- min(z.min, floor(zlim[1] / z.scal))
  }
  z.range <- range(z.min, z.max)

  # Add grid
  if ("xy" %in% grid || grid == TRUE) {
    i <- x.min:x.max
    segments(i, z.min, i + (yx.f * y.max), yz.f * y.max +
      z.min, col = col.grid, lty = lty.grid)
    i <- 0:y.max
    segments(x.min + (i * yx.f), i * yz.f + z.min, x.max +
      (i * yx.f), i * yz.f + z.min, col = col.grid, lty = lty.grid)
  }

  if ("xz" %in% grid) {
    i <- x.min:x.max
    segments(i + (yx.f * y.max), yz.f * y.max + z.min,
      i + (yx.f * y.max), yz.f * y.max + z.max,
      col = col.grid, lty = lty.grid
    )
    temp <- yx.f * y.max
    temp1 <- yz.f * y.max
    i <- z.min:z.max
    segments(x.min + temp, temp1 + i,
      x.max + temp, temp1 + i,
      col = col.grid, lty = lty.grid
    )
  }

  if ("yz" %in% grid) {
    i <- 0:y.max
    segments(x.min + (i * yx.f), i * yz.f + z.min,
      x.min + (i * yx.f), i * yz.f + z.max,
      col = col.grid, lty = lty.grid
    )
    temp <- yx.f * y.max
    temp1 <- yz.f * y.max
    i <- z.min:z.max
    segments(x.min + temp, temp1 + i,
      x.min, i,
      col = col.grid, lty = lty.grid
    )
  }
}
