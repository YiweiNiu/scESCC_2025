



#####################################################
# plot_cnv_violin
#####################################################

plot_cnv_violin <- function (sam=NULL, infercnv_dir=NULL, meta=NULL) {
  # read final infercnv obj
  x = readRDS(file.path(infercnv_dir, sam, 'run.final.infercnv_obj'))
  # x@expr.data - 1
  y = (x@expr.data - 1)^2
  z = as.data.frame(apply(y, 2, mean)) %>%
    rownames_to_column('barcode') %>%
    dplyr::rename(value = `apply(y, 2, mean)`) %>%
    left_join(meta, by = 'barcode') %>%
    mutate(cellType = factor(cellType,
                             levels = c('T cells', 'B cells', 'Myeloid',
                                        'Epithelia', 'Fibroblasts', 'Endothelia')))
  # cmp Epithelia with other non-immune cells
  my_comparisons <- list( c("Epithelia", "Fibroblasts"), c("Epithelia", "Endothelia") )
  gg = ggplot(z, aes(x=cellType, y=value, fill=cellType)) +
    geom_violin(trim=F, adjust=2) +
    geom_boxplot(width=.1, fill="black", outlier.colour=NA) +
    stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=1.5) +
    labs(x='', y='CNV level', title=sam, fill = NULL) +
    scale_fill_jco() +
    theme(legend.position = 'None',
          axis.text.x = element_text(angle = 60, hjust = 1))
  return(gg)
}

#####################################################
# plot CNV using complexheatmap
# works, improvement needed
#####################################################

plotCNV_complex = function(infercnv_obj=NULL,
                           obs_order = NULL,
                           cell_ann = NULL) {
  require(ComplexHeatmap)

  # expr
  exp = infercnv_obj@expr.data
  # ref mat
  ref_cell_lst = infercnv_obj@reference_grouped_cell_indices
  ref_mat_lst = lapply(ref_cell_lst, FUN = function(x){
    exp[,x]
  })
  ref_mat = do.call(cbind, ref_mat_lst)
  # obs mat
  obs_cell_lst = infercnv_obj@observation_grouped_cell_indices
  if (!is.null(obs_order)) {
    obs_mat_lst = lapply(obs_order, FUN = function(x){
      #cat(x)
      exp[,obs_cell_lst[[x]]]
    })
  } else {
    obs_mat_lst = lapply(obs_cell_lst, FUN = function(x){
      exp[,x]
    })
  }
  obs_mat = do.call(cbind, obs_mat_lst)
  # merge two mat
  mat_4_p = cbind(ref_mat, obs_mat)
  mat_4_p = t(mat_4_p)
  # color
  x = max(mat_4_p)
  y = min(mat_4_p)
  col_fun = circlize::colorRamp2(c(y, 1, x), c("#00008B", "white", "#8B0000"))
  # column anno
  gp = infercnv_obj@gene_order
  column_ha = HeatmapAnnotation(
    chr = gp$chr,
    col = list(chr = setNames(ggsci::pal_ucscgb()(nlevels(gp$chr)),
                              levels(gp$chr))),
    show_legend = F,
    show_annotation_name = F)

  # row anno
  ref_cluster = rep(names(ref_cell_lst), sapply(ref_cell_lst, length))
  if (!is.null(obs_order)) {
    obs_cluster = rep(obs_order, sapply(obs_order, FUN = function(x) {
      length(obs_cell_lst[[x]])
    }))
  } else {
    obs_cluster = rep(names(obs_cell_lst), sapply(obs_cell_lst, length))
  }
  cell_cluster = c(ref_cluster, obs_cluster)
  row_anno = data.frame(row.names = rownames(mat_4_p),
                        cluster = cell_cluster)
  if (!is.null(cell_ann)) {
    row_anno = cbind(row_anno, cell_ann)
  }
  row_ha = rowAnnotation(df = row_anno,
                         annotation_legend_param = list(
                           cluster = list(
                             at = unique(cell_cluster),
                             labels = unique(cell_cluster),
                             nrow = 2
                           )
                         ))
  # row split
  row_split = c(rep('ref', ncol(ref_mat)),
                rep('obs', ncol(obs_mat)))
  # plot
  ht1 = Heatmap(mat_4_p,
                col = col_fun,
                # not cluster columns and rows
                cluster_rows = FALSE, cluster_columns = FALSE,
                # name for color bar
                name = 'Inferred CNV',
                # column names and row names
                show_row_names = FALSE, show_column_names = FALSE,
                # column title and row title
                column_title = 'Chr',
                column_title_side = 'bottom',
                # column annotation
                bottom_annotation = column_ha,
                # row annotation
                right_annotation = row_ha,
                # row split
                row_split = row_split,
                row_gap = unit(2, "mm"),
                # legend
                heatmap_legend_param = list(direction = "horizontal")
  )
  draw(ht1,
       heatmap_legend_side = "bottom",
       annotation_legend_side = "bottom")
}



#####################################################
# plot CNV using complexheatmap 2
# add more color bars
# use the clustering results of infercnv
#####################################################

plotCNV_complex2 = function(infercnv_obj = NULL,
                            infercnv_dend = NULL,
                            srat_obj = NULL,
                            row_ann = c("Origin4")) {
  require(ComplexHeatmap)
  # expr
  exp = infercnv_obj@expr.data
  # ref mat
  ref_cell_lst = infercnv_obj@reference_grouped_cell_indices
  ref_mat_lst = lapply(ref_cell_lst, FUN = function(x){
    exp[,x]
  })
  ref_mat = do.call(cbind, ref_mat_lst)
  # obs mat
  obs_mat = exp[, labels(infercnv_dend)]
  # merge two mat
  mat_4_p = cbind(ref_mat, obs_mat)
  mat_4_p = t(mat_4_p)
  # 0.7 - 1.3
  mat_4_p[mat_4_p < 0.7] <- 0.7
  mat_4_p[mat_4_p > 1.3] <- 1.3
  # color
  x = max(mat_4_p)
  y = min(mat_4_p)
  col_fun = circlize::colorRamp2(c(y, 1, x), c("#00008B", "white", "#8B0000"))
  # column anno
  gp = infercnv_obj@gene_order
  column_ha = HeatmapAnnotation(
    chr = gp$chr,
    col = list(chr = setNames(ggsci::pal_ucscgb()(nlevels(gp$chr)),
                              levels(gp$chr))),
    show_legend = F,
    show_annotation_name = F)

  # row anno
  if (is.null(row_ann)) {
    row_ha = NULL
  } else {
    # row anno
    row_anno = FetchData(seurat, vars = row_ann)
    # order
    row_anno = row_anno[rownames(mat_4_p),]
    row_ha = rowAnnotation(df = row_anno,
                           col = list(
                             Origin2_n = c("PBMC1" = "#023fa5", "PBMC2" = "#7d87b9",
                                           "LN_N" = "#ff7f03", "LN_P" = "#EFC000FF",
                                           "Normal" = "#7AA6DCFF", "Adjacent" = "#868686FF",
                                           "Tumor" = "#CD534CFF")
                           ),
                           annotation_legend_param = list(
                             Origin2_n = list(
                               nrow = 4
                             ),
                             Source = list(
                               nrow = 4
                             )
                           ))
  }

  # row split
  row_split = factor(
    rep(c("ref", "obs"), c(ncol(ref_mat), ncol(obs_mat))),
    levels = c("ref", "obs")
  )
  # plot
  ht_list = Heatmap(mat_4_p,
                col = col_fun,
                # not cluster columns and rows
                cluster_rows = FALSE, cluster_columns = FALSE,
                # name for color bar
                name = 'ht1',
                # column names and row names
                show_row_names = FALSE, show_column_names = FALSE,
                # column title and row title
                column_title = 'Chr',
                column_title_side = 'bottom',
                # column annotation
                bottom_annotation = column_ha,
                # row annotation
                left_annotation = row_ha,
                # row split
                row_split = row_split,
                row_gap = unit(3, "mm"),
                # legend
                heatmap_legend_param = list(direction = "horizontal")
  )
  draw(ht_list,
       heatmap_legend_side = "bottom",
       annotation_legend_side = "bottom")
  # add vertical line decoration
  tot_gene = ncol(mat_4_p)
  agg_c = 0
  decorate_heatmap_body("ht1", {
    for (i in table(gp$chr)) {
      grid.lines(c(agg_c/tot_gene, agg_c/tot_gene), c(0, 1), gp = gpar(lty = 2, lwd = 1))
      agg_c = agg_c + i
    }
    #grid.lines(c(0.5, 0.5), c(0, 1), gp = gpar(lty = 2, lwd = 2))
  }, slice = 2)
}

#####################################################
# detect_malignant
# this was used in infercnv.03.Rmd
#####################################################

detect_malignant <- function(
  sam = NULL,
  sam_type = NULL,
  epi_meta = NULL,
  infercnv_dir = '/niuyw-usb-disk/Projects/scESCA/infercnv/200227_6samples'
) {
  x = readRDS(file.path(infercnv_dir, sam, 'run.final.infercnv_obj'))
  # epithelia
  epi_barcode = epi_meta %>% filter(Source == sam) %>% pull(barcode)
  x = x@expr.data[,epi_barcode]
  # -1, square
  y = (x - 1)^2
  # cnv sig, mean
  z = as.data.frame(apply(y, 2, mean)) %>%
    rownames_to_column('barcode') %>%
    dplyr::rename(cnv_sig = `apply(y, 2, mean)`)

  # classify
  if (sam_type == 'N') {
    q_10 = quantile(z$cnv_sig, 0.1)
    z = z %>%
      mutate(malignant = case_when(
        cnv_sig < q_10 ~ 'No', # non-malignant
        cnv_sig >= q_10 ~ 'Un' # unresolved
      ))
    # mean cnv of non-malignant
    mean_cnv = rowMeans(x[, z[z$malignant == 'No',]$barcode])
    # cnv R
    a = apply(x, 2, FUN = function(x) {cor(x, mean_cnv)})
    z$cnv_r = a
    z = z %>%
      mutate(malignant = case_when(
        malignant == "No" | cnv_r > 0.8 ~ 'No', # non-malignant
        TRUE  ~ "Un"
      ))
  } else {
    q_90 = quantile(z$cnv_sig, 0.9)
    z = z %>%
      mutate(malignant = case_when(
        cnv_sig > q_90 ~ 'Yes', # malignant
        cnv_sig <= q_90 ~ 'Un' # unresolved
      ))
    # mean cnv of non-malignant
    mean_cnv = rowMeans(x[, z[z$malignant == 'Yes',]$barcode])
    # cnv R
    a = apply(x, 2, FUN = function(x) {cor(x, mean_cnv)})
    z$cnv_r = a
    z = z %>%
      mutate(malignant = case_when(
        malignant == "Yes" | cnv_r > 0.8 ~ 'Yes', # malignant
        TRUE ~ "Un"
      ))
  }

  return(z)
}


#####################################################
# detect_malignant2
# this was used in infercnv.06.epi.Rmd
#####################################################

detect_malignant2 <- function(
  sam = NULL,
  sam_type = NULL,
  epi_meta = NULL,
  infercnv_dir = '/niuyw-usb-disk/Projects/scESCA/infercnv/200227_6samples'
) {
  x = readRDS(file.path(infercnv_dir, sam, 'run.final.infercnv_obj'))
  # epithelia
  epi_barcode = epi_meta %>% filter(Source == sam) %>% pull(barcode)
  x = x@expr.data[,epi_barcode]
  # -1, square
  y = (x - 1)^2
  # cnv sig, mean
  z = as.data.frame(apply(y, 2, mean)) %>%
    rownames_to_column('barcode') %>%
    dplyr::rename(cnv_sig = `apply(y, 2, mean)`)

  # classify
  if (sam_type == 'N') {
    q_10 = quantile(z$cnv_sig, 0.1)
    z = z %>%
      mutate(malignant = case_when(
        cnv_sig < q_10 ~ 'No', # non-malignant
        cnv_sig >= q_10 ~ 'Un' # unresolved
      ))
    # mean cnv of non-malignant
    mean_cnv = rowMeans(x[, z[z$malignant == 'No',]$barcode])
    # cnv R
    a = apply(x, 2, FUN = function(x) {cor(x, mean_cnv)})
    z$cnv_r = a
    z = z %>%
      mutate(malignant = case_when(
        malignant == "No" | cnv_r > 0.8 ~ 'No', # non-malignant
        TRUE  ~ "Un"
      ))
  } else {
    q_90 = quantile(z$cnv_sig, 0.9)
    z = z %>%
      mutate(malignant = case_when(
        cnv_sig > q_90 ~ 'Yes', # malignant
        cnv_sig <= q_90 ~ 'Un' # unresolved
      ))
    # mean cnv of non-malignant
    mean_cnv = rowMeans(x[, z[z$malignant == 'Yes',]$barcode])
    # cnv R
    a = apply(x, 2, FUN = function(x) {cor(x, mean_cnv)})
    z$cnv_r = a
    z = z %>%
      mutate(malignant = case_when(
        malignant == "Yes" | cnv_r > 0.8 ~ 'Yes', # malignant
        TRUE ~ "Un"
      ))
  }
  return(z)
}




