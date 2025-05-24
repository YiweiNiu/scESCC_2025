

# modified from: https://github.com/thomas-neitmann/ggcharts
diverging_bar_chart <- function(data, x, y,
                                colors = c("darkred", "navyblue"),
                                line_size = 0.75,
                                text_size = 2.5,
                                point_size = 3) {
  x <- rlang::enquo(x)
  y <- rlang::enquo(y)

  data <- dplyr::mutate(
    .data = data,
    !!x := reorder(!!x, !!y)
  )

  limit <- max(dplyr::pull(data, !!y)) * 1.05

  plot <- ggplot(data, aes(!!x, !!y, fill = !!y, color = !!y)) +
    geom_bar(position = "dodge", stat="identity") +
    scale_color_gradient(low=colors[2], high=colors[1], guide=F) +
    scale_fill_gradient(low = adjustcolor(colors[2], 0.3),
                        high = adjustcolor(colors[1], 0.3),
                        guide = "none")

  plot +
    coord_flip() +
    geom_text(
      data = dplyr::filter(data, !!y >= 0),
      color = "black",
      size = text_size,
      aes(label = !!x, y = 0, hjust = "right"),
      nudge_y = -limit * .01
    ) +
    geom_text(
      data = dplyr::filter(data, !!y < 0),
      color = "black",
      size = text_size,
      aes(label = !!x, y = 0, hjust = "left"),
      nudge_y = limit * .013,
    ) +
    geom_hline(
      yintercept = 0,
      color = "black",
      size = .4
    ) +
    labs(x = NULL) +
    guides(y = "none")
}

# plot heatmap of glm
plot_heatmap_glm <- function(mat=NULL,
                             highlight_rows=NULL,
                             highlight_cols=NULL,
                             ...) {
  require(ComplexHeatmap)
  mat = pheatmap:::scale_rows(mat)
  x = ceiling(max(abs(mat)))
  col_fun = circlize::colorRamp2(c(-x, 0, x), c("blue", "white", "red"))

  # highligh rows
  if(!is.null(highlight_rows) & length(highlight_rows) != 0) {
    row_idx <- which(rownames(mat) %in% highlight_rows)
    fontsizes <- rep(6, nrow(mat))
    fontsizes[row_idx] <- 8
    fontcolors <- rep('black', nrow(mat))
    fontcolors[row_idx] <- 'red'

    # Create text annotation object for displaying row names
    rowAnno <- rowAnnotation(rows = anno_text(rownames(mat),
                                              gp = gpar(fontsize = fontsizes,
                                                        col = fontcolors)))
    show_row_names = FALSE
  } else {
    rowAnno = NULL
    show_row_names = TRUE
  }
  # highligh cols
  if(!is.null(highlight_cols) & length(highlight_cols) != 0) {
    col_idx <- which(colnames(mat) %in% highlight_cols)
    fontsizes <- rep(6, ncol(mat))
    fontsizes[col_idx] <- 8
    fontcolors <- rep('black', ncol(mat))
    fontcolors[col_idx] <- 'red'

    # Create text annotation object for displaying row names
    colAnno <- HeatmapAnnotation(cols = anno_text(colnames(mat),
                                                  gp = gpar(fontsize = fontsizes,
                                                            col = fontcolors),
                                                  rot = 70),
                                 which = "column")
    show_column_names = FALSE
  } else {
    colAnno = NULL
    show_column_names = TRUE
  }
  ht = Heatmap(mat, name = "GLM t",
               col = col_fun,
               # row name/anno
               show_row_names = show_row_names,
               row_names_gp = gpar(fontsize = 4),
               right_annotation = rowAnno,
               # col name/anno
               show_column_names = show_column_names,
               column_names_gp = gpar(fontsize = 4),
               column_names_rot = 70,
               bottom_annotation = colAnno,
               # legend
               heatmap_legend_param = list(direction = "vertical",
                                           labels_gp = gpar(fontsize = 4),
                                           title_gp = gpar(fontsize = 5),
                                           title_position = "topleft",
                                           grid_height = unit(2, "mm"),
                                           grid_width = unit(2, "mm")),
               ...
  )

  p <- plot_grid(grid.grabExpr(draw(ht, heatmap_legend_side = "right")))
  p
}



# compare cell scores tissue by tissue using glm
# two tissues, e.g. Tumor vs. Normal
tissue_compare <- function(seurat_obj, score_mat, tissues) {
  cells = intersect(colnames(seurat_obj), colnames(score_mat))
  score_mat = score_mat[, cells]
  my_glm <- function(x, seurat_obj, score_mat, tissues) {
    x = tibble(score = x,
               tissue = seurat_obj$Tissue[colnames(score_mat)],
               patient = seurat_obj$Patient[colnames(score_mat)]) %>%
      filter(tissue %in% tissues) %>%
      mutate(tissue = factor(tissue, levels = tissues)) %>%
      mutate(tissue = as.numeric(tissue))
    fit = glm(score ~ tissue + patient, family = gaussian, data = x)
    summary(fit)[["coefficients"]][2, c("t value", "Pr(>|t|)")]
  }
  a = apply(score_mat, 1, my_glm, seurat_obj = seurat_obj,
            score_mat = score_mat, tissues = tissues)

  # keep only significant ones
  df = tibble(signature = colnames(a), t = a[1, ], p = a[2, ]) %>%
    mutate(p_adj = p.adjust(p, "BH")) %>%
    filter(p_adj < 0.05)
  # plot
  # show only the top 10 for pos/neg
  d_4_p = df %>%
    filter(t > 0) %>%
    top_n(10, wt = t) %>%
    bind_rows(df %>%
                filter(t < 0) %>%
                top_n(10, wt = -(t)))

  return(list(df = df, d_4_p = d_4_p))
}


# compare each Origin2_n against other Origin2_n using glm
tissue_compare_each <- function(seurat_obj, score_mtx) {
  cells = intersect(colnames(seurat_obj), colnames(score_mtx))
  score_mtx = score_mtx[, cells]
  my_cmp2_helper <- function(x, seurat_obj, cells) {
    u = x[cells]
    a = tibble(score = u,
               cluster = 1,
               patient = seurat_obj$Patient[names(u)])
    v = x[!(names(x) %in% cells)]
    b = tibble(score = v,
               cluster = 0,
               patient = seurat_obj$Patient[names(v)]) %>%
      bind_rows(a) %>%
      mutate(cluster = as.factor(cluster))
    fit = glm(score ~ cluster + patient, family = gaussian, data = b)
    summary(fit)[["coefficients"]][2, c('t value', 'Pr(>|t|)')]
  }

  my_cmp2 <- function(cells, mtx, seurat_obj) {
    apply(mtx, 1, my_cmp2_helper, seurat_obj = seurat_obj, cells = cells)
  }

  seurat_obj$Origin2_n = droplevels(seurat_obj$Origin2_n)
  x = sapply(split(colnames(seurat_obj), seurat_obj$Origin2_n),
             my_cmp2, mtx = score_mtx, seurat_obj = seurat_obj)

  t = x[seq(1, dim(score_mtx)[1]*2-1, 2),]; rownames(t) = rownames(score_mtx)
  p = x[seq(2, dim(score_mtx)[1]*2, 2),]; rownames(p) = rownames(score_mtx)
  p_adj = apply(p, 2, p.adjust, method = 'BH')

  # neat
  p_l = p %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    gather('cluster', 'p', 2:(dim(p)[2]+1)) %>%
    mutate(id = paste(gene, cluster, sep = '-')) %>%
    dplyr::select(id, p)
  q_l = p_adj %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    gather('cluster', 'p_adj', 2:(dim(p)[2]+1)) %>%
    mutate(id = paste(gene, cluster, sep = '-')) %>%
    dplyr::select(id, p_adj)
  t_l = as_tibble(as.data.frame(t)) %>%
    mutate(gene = rownames(t)) %>%
    gather('cluster', 't', 1:dim(t)[2]) %>%
    mutate(id = paste(gene, cluster, sep = '-')) %>%
    left_join(p_l, by = 'id') %>%
    left_join(q_l, by = 'id') %>%
    dplyr::select(-id)

  return(list(df = t_l, t_mat = t))
}



# compare each tissue against other tissues using glm, using cell.meta other than seurat
# Origin4
tissue_compare_each2 <- function(score_mtx=NULL, meta=NULL) {
  # filter cells
  cells = intersect(rownames(meta), colnames(score_mtx))
  score_mtx = score_mtx[, cells]
  #
  my_cmp2_helper <- function(x, meta, cluster) {
    dat = data.frame(barcode = names(x),
                     score = x) %>%
      left_join(meta %>% rownames_to_column("barcode"), by = "barcode") %>%
      mutate(group = case_when(
        Origin4 == cluster ~ 1,
        TRUE ~ 0
      ),
      patient = as.factor(Patient)) %>%
      mutate(group = as.factor(group))
    fit = glm(score ~ group + patient, family = gaussian, data = dat)
    summary(fit)[["coefficients"]][2, c('t value', 'Pr(>|t|)')]
    #broom::tidy(fit)
  }

  my_cmp2 <- function(cluster, mtx, meta) {
    apply(mtx, 1, my_cmp2_helper, meta = meta, cluster = cluster)
  }

  meta$Origin4 = droplevels(meta$Origin4)
  x = sapply(levels(meta$Origin4),
             my_cmp2, mtx = score_mtx, meta = meta)

  t = x[seq(1, dim(score_mtx)[1]*2-1, 2),]; rownames(t) = rownames(score_mtx)
  p = x[seq(2, dim(score_mtx)[1]*2, 2),]; rownames(p) = rownames(score_mtx)
  p_adj = apply(p, 2, p.adjust, method = 'BH')
  # neat
  p_l = p %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    gather('cluster', 'p', 2:(dim(p)[2]+1)) %>%
    mutate(id = paste(gene, cluster, sep = '-')) %>%
    dplyr::select(id, p)
  q_l = p_adj %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    gather('cluster', 'p_adj', 2:(dim(p)[2]+1)) %>%
    mutate(id = paste(gene, cluster, sep = '-')) %>%
    dplyr::select(id, p_adj)
  t_l = as_tibble(as.data.frame(t)) %>%
    mutate(gene = rownames(t)) %>%
    gather('cluster', 't', 1:dim(t)[2]) %>%
    mutate(id = paste(gene, cluster, sep = '-')) %>%
    left_join(p_l, by = 'id') %>%
    left_join(q_l, by = 'id') %>%
    dplyr::select(-id)
  return(t_l)
}




# compare each cluster against other clusters using glm
cluster_compare <- function(seurat_obj, score_mtx) {
  cells = intersect(colnames(seurat_obj), colnames(score_mtx))
  score_mtx = score_mtx[, cells]
  my_cmp2_helper <- function(x, seurat_obj, cells) {
    u = x[cells]
    a = tibble(score = u,
               cluster = 1,
               patient = seurat_obj$Patient[names(u)])
    v = x[!(names(x) %in% cells)]
    b = tibble(score = v,
               cluster = 0,
               patient = seurat_obj$Patient[names(v)]) %>%
      bind_rows(a) %>%
      mutate(cluster = as.factor(cluster))
    fit = glm(score ~ cluster + patient, family = gaussian, data = b)
    summary(fit)[["coefficients"]][2, c('t value', 'Pr(>|t|)')]
  }

  my_cmp2 <- function(cells, mtx, seurat_obj) {
    apply(mtx, 1, my_cmp2_helper, seurat_obj = seurat_obj, cells = cells)
  }

  x = sapply(split(colnames(seurat_obj), seurat_obj$seurat_clusters),
             my_cmp2, mtx = score_mtx, seurat_obj = seurat_obj)

  t = x[seq(1, dim(score_mtx)[1]*2-1, 2),]; rownames(t) = rownames(score_mtx)
  p = x[seq(2, dim(score_mtx)[1]*2, 2),]; rownames(p) = rownames(score_mtx)
  p_adj = apply(p, 2, p.adjust, method = 'BH')

  # neat
  p_l = p %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    gather('cluster', 'p', 2:(dim(p)[2]+1)) %>%
    mutate(id = paste(gene, cluster, sep = '-')) %>%
    dplyr::select(id, p)
  q_l = p_adj %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    gather('cluster', 'p_adj', 2:(dim(p)[2]+1)) %>%
    mutate(id = paste(gene, cluster, sep = '-')) %>%
    dplyr::select(id, p_adj)
  t_l = as_tibble(as.data.frame(t)) %>%
    mutate(gene = rownames(t)) %>%
    gather('cluster', 't', 1:dim(t)[2]) %>%
    mutate(id = paste(gene, cluster, sep = '-')) %>%
    left_join(p_l, by = 'id') %>%
    left_join(q_l, by = 'id') %>%
    dplyr::select(-id)

  return(list(df = t_l, t_mat = t))
}


# compare each cluster against other clusters using glm, using cell.meta other than seurat
cluster_compare2 <- function(score_mtx=NULL, meta=NULL) {
  # filter cells
  cells = intersect(rownames(meta), colnames(score_mtx))
  score_mtx = score_mtx[, cells]
  #
  my_cmp2_helper <- function(x, meta, cluster) {
    dat = data.frame(barcode = names(x),
                     score = x) %>%
      left_join(meta %>% rownames_to_column("barcode"), by = "barcode") %>%
      mutate(group = case_when(
        seurat_clusters == cluster ~ 1,
        TRUE ~ 0
      ),
      patient = as.factor(Patient)) %>%
      mutate(group = as.factor(group))
    fit = glm(score ~ group + patient, family = gaussian, data = dat)
    summary(fit)[["coefficients"]][2, c('t value', 'Pr(>|t|)')]
    #broom::tidy(fit)
  }

  my_cmp2 <- function(cluster, mtx, meta) {
    apply(mtx, 1, my_cmp2_helper, meta = meta, cluster = cluster)
  }

  x = sapply(levels(meta$seurat_clusters),
             my_cmp2, mtx = score_mtx, meta = meta)

  t = x[seq(1, dim(score_mtx)[1]*2-1, 2),]; rownames(t) = rownames(score_mtx)
  p = x[seq(2, dim(score_mtx)[1]*2, 2),]; rownames(p) = rownames(score_mtx)
  p_adj = apply(p, 2, p.adjust, method = 'BH')
  # neat
  p_l = p %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    gather('cluster', 'p', 2:(dim(p)[2]+1)) %>%
    mutate(id = paste(gene, cluster, sep = '-')) %>%
    dplyr::select(id, p)
  q_l = p_adj %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    gather('cluster', 'p_adj', 2:(dim(p)[2]+1)) %>%
    mutate(id = paste(gene, cluster, sep = '-')) %>%
    dplyr::select(id, p_adj)
  t_l = as_tibble(as.data.frame(t)) %>%
    mutate(gene = rownames(t)) %>%
    gather('cluster', 't', 1:dim(t)[2]) %>%
    mutate(id = paste(gene, cluster, sep = '-')) %>%
    left_join(p_l, by = 'id') %>%
    left_join(q_l, by = 'id') %>%
    dplyr::select(-id)
  return(t_l)
}

# compare each cluster against other clusters using glm, using cell.meta other than seurat, no patient correction
cluster_compare3 <- function(score_mtx=NULL, meta=NULL) {
  # filter cells
  cells = intersect(rownames(meta), colnames(score_mtx))
  score_mtx = score_mtx[, cells]
  #
  my_cmp2_helper <- function(x, meta, cluster) {
    dat = data.frame(barcode = names(x),
                     score = x) %>%
      left_join(meta %>% rownames_to_column("barcode"), by = "barcode") %>%
      mutate(group = case_when(
        seurat_clusters == cluster ~ 1,
        TRUE ~ 0
      )) %>%
      mutate(group = as.factor(group))
    fit = glm(score ~ group, family = gaussian, data = dat)
    summary(fit)[["coefficients"]][2, c('t value', 'Pr(>|t|)')]
    #broom::tidy(fit)
  }

  my_cmp2 <- function(cluster, mtx, meta) {
    apply(mtx, 1, my_cmp2_helper, meta = meta, cluster = cluster)
  }

  x = sapply(levels(meta$seurat_clusters),
             my_cmp2, mtx = score_mtx, meta = meta)

  t = x[seq(1, dim(score_mtx)[1]*2-1, 2),]; rownames(t) = rownames(score_mtx)
  p = x[seq(2, dim(score_mtx)[1]*2, 2),]; rownames(p) = rownames(score_mtx)
  p_adj = apply(p, 2, p.adjust, method = 'BH')
  # neat
  p_l = p %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    gather('cluster', 'p', 2:(dim(p)[2]+1)) %>%
    mutate(id = paste(gene, cluster, sep = '-')) %>%
    dplyr::select(id, p)
  q_l = p_adj %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    gather('cluster', 'p_adj', 2:(dim(p)[2]+1)) %>%
    mutate(id = paste(gene, cluster, sep = '-')) %>%
    dplyr::select(id, p_adj)
  t_l = as_tibble(as.data.frame(t)) %>%
    mutate(gene = rownames(t)) %>%
    gather('cluster', 't', 1:dim(t)[2]) %>%
    mutate(id = paste(gene, cluster, sep = '-')) %>%
    left_join(p_l, by = 'id') %>%
    left_join(q_l, by = 'id') %>%
    dplyr::select(-id)
  return(t_l)
}


