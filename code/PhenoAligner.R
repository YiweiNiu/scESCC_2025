

suppressPackageStartupMessages({
  library(tidyverse)
  library(BiocNeighbors)
  library(BiocParallel)
  library(dplyr)
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(Seurat)
  library(Matrix)
  library(Hmisc)
  library("future.apply")
})

# KNN
cross_knn <- function(srat_query = NULL, srat_ref = NULL, k = NULL) {
  # get cells
  coldata.query <- srat_query@meta.data %>%
    rownames_to_column("cell") %>%
    dplyr::select(cell, Origin4, Patient, cellType, cellType2, cellType3)
  coldata.ref <- srat_ref@meta.data %>%
    rownames_to_column("cell") %>%
    dplyr::select(cell, Origin4, Patient, cellType, cellType2, cellType3)
  # get PCA score of each cell
  pc.query <- srat_query@reductions$pca[[]]
  pc.ref <- srat_ref@reductions$pca[[]]
  cat("[INFO] Query KNN ...\n")
  N12 <- queryKNN(pc.ref, query = pc.query, k = k, get.distance = T, BPPARAM = MulticoreParam(workers = 8))
  N12$cellnames <- rownames(pc.query)
  knn.list <- list(
    pca.ref = pc.ref,
    pca.query = pc.query,
    queryToRef = N12,
    coldata.query = coldata.query,
    coldata.ref = coldata.ref
  )
  return(knn.list)
}

# cor
get_cor_future <- function(srat_query = NULL, srat_ref = NULL, knn.lst = NULL, cores = 8) {
  data.query <- GetAssayData(srat_query, slot = "data")
  data.ref <- GetAssayData(srat_ref, slot = "data")
  coldata.query <- knn.lst$coldata.query
  queryToRef_idx <- knn.lst$queryToRef$index
  plan(multisession, workers = cores) ## Run in parallel on local computer
  result <- future_lapply(1:ncol(data.query), function(i) {
    cell_query <- coldata.query[i, ]$cell
    y <- data.query[, cell_query]
    cell_ref <- queryToRef_idx[i, ]
    ref_mat <- data.ref[, cell_ref]
    cor.mat <- rcorr(as.matrix(y), as.matrix(ref_mat), type = "pearson")
    cor.result <- data.frame(cor.mat$r[1, 2:31])
    colnames(cor.result) <- "cor"
    cor.result$cell <- rownames(cor.result)
    cor.result <- cor.result[order(cor.result$cor, decreasing = T), ]
    pvalue <- data.frame(cor.mat$P[1, 2:31])
    colnames(pvalue) <- "pvalue"
    pvalue$cell <- rownames(pvalue)
    myresult <- plyr::join(cor.result, pvalue, by = "cell")
    return(myresult)
  })
  # Explicitly close multisession workers by switching plan
  plan(sequential)
  return(result)
}

# cor
get_cor <- function(srat_query = NULL, srat_ref = NULL, knn.lst = NULL) {
  data.query <- GetAssayData(srat_query, slot = "data")
  data.ref <- GetAssayData(srat_ref, slot = "data")
  coldata.query <- knn.lst$coldata.query
  queryToRef_idx <- knn.lst$queryToRef$index
  result <- list()
  for (i in 1:ncol(data.query)) {
    cell_query <- coldata.query[i, ]$cell
    y <- data.query[, cell_query]
    cell_ref <- queryToRef_idx[i, ]
    ref_mat <- data.ref[, cell_ref]
    cor.mat <- rcorr(as.matrix(y), as.matrix(ref_mat), type = "pearson")
    cor.result <- data.frame(cor.mat$r[1, 2:31])
    colnames(cor.result) <- "cor"
    cor.result$cell <- rownames(cor.result)
    cor.result <- cor.result[order(cor.result$cor, decreasing = T), ]
    pvalue <- data.frame(cor.mat$P[1, 2:31])
    colnames(pvalue) <- "pvalue"
    pvalue$cell <- rownames(pvalue)
    myresult <- plyr::join(cor.result, pvalue, by = "cell")
    result[[i]] <- myresult
  }
  #length(result)
  return(result)
}

# filter by cor
filter_by_cor <- function(cordata = NULL, knn.lst = NULL, threshold = .6) {
  #cordata <- readRDS(here::here("output", DOCNAME, t, paste0(cell, "_COR_MAT.rds")))
  cor_mat <- as.data.frame(do.call(rbind, cordata))
  cor_mat$BH <- p.adjust(cor_mat$pvalue, method = "BH", n = length(cor_mat$pvalue))
  cor_mat$origin <- rep(knn.list$coldata.query$cell, each = 30)
  cor_mat_keep <- cor_mat[which(cor_mat$cor > threshold), ]
  dim(cor_mat)
  dim(cor_mat_keep)
  length(unique(cor_mat_keep$origin))
  cor_mat_keep <- cor_mat_keep[which(cor_mat_keep$BH < 0.05), ]
  length(unique(cor_mat_keep$origin))

  phenoaligner <- cor_mat_keep %>%
    group_by(origin) %>%
    dplyr::mutate(
      cor_max = max(cor),
      cor_max2 = sort(cor, decreasing = TRUE)[2]
    ) %>%
    dplyr::filter(cor == max(cor))
  map_data <- knn.lst$coldata.ref
  phenoaligner <- left_join(phenoaligner, map_data, by = "cell")
  colnames(phenoaligner) <- c(
    "cor", "barcode", "pvalue", "BH", "origin", "cor_max", "cor_max2",
    "nn.tissue", "nn.patient", "nn.celltype", "nn.celltype2", "nn.celltype_sub"
  )

  map_data <- knn.lst$coldata.query
  colnames(map_data)[1] <- "origin"
  phenoaligner <- left_join(phenoaligner, map_data, by = "origin")
  return(list(
    cor_mat_keep = cor_mat_keep,
    phenoaligner = phenoaligner
  ))
}

filter_by_gap <- function(knn.lst = NULL, cor_mat_keep = NULL, phenoaligner=NULL, threshold = .01) {
  map_data <- knn.lst$coldata.ref
  different_tissue <- left_join(cor_mat_keep, map_data, by = "cell")
  colnames(different_tissue) <- c(
    "cor", "barcode", "pvalue", "BH", "origin",
    "nn.tissue", "nn.patient", "nn.celltype", "nn.celltype2", "nn.celltype_sub"
  )
  mph_tissue_max_cor <- different_tissue %>%
    group_by(origin, nn.tissue) %>%
    summarise(tissue_max = max(cor), cells = n())
  mph_tissue_max_diff <- mph_tissue_max_cor %>%
    group_by(origin) %>%
    dplyr::mutate(
      top_one = max(tissue_max),
      top_two = sort(tissue_max, decreasing = TRUE)[2]
    )
  mph_tissue_max_diff <- data.frame(mph_tissue_max_diff)
  mph_tissue_max_diff$gap <- mph_tissue_max_diff$top_one - mph_tissue_max_diff$top_two

  mph_tissue_max_diff <- mph_tissue_max_diff[!duplicated(mph_tissue_max_diff$origin), ]
  dim(mph_tissue_max_diff)
  # [1] 768   7
  mph_tissue_max_diff <- mph_tissue_max_diff[-which(mph_tissue_max_diff$gap < threshold), ]
  dim(mph_tissue_max_diff)
  ## [1] 440   7
  coldata.mt <- phenoaligner[which(phenoaligner$origin %in% mph_tissue_max_diff$origin), ]
  return(coldata.mt)
}

.shuffle_label <- function(cordata = NULL, knn.lst = NULL,
                           cor_threshold = .6, gap_threshold = .01,
                           times = 1000,
                           seed = 816) {
  #cordata <- readRDS(here::here("output", DOCNAME, t, paste0(cell, "_COR_MAT.rds")))
  fake.lst = list()
  set.seed(seed)
  for (i in 1:times) {
    knn.lst_fake = knn.lst
    # shuffle ref tissue labels
    knn.lst_fake$coldata.ref$Origin4 = sample(knn.lst_fake$coldata.ref$Origin4)
    cor_mat <- as.data.frame(do.call(rbind, cordata))
    cor_mat$BH <- p.adjust(cor_mat$pvalue, method = "BH", n = length(cor_mat$pvalue))
    cor_mat$origin <- rep(knn.list$coldata.query$cell, each = 30)
    cor_mat_keep <- cor_mat[which(cor_mat$cor > cor_threshold), ]
    #dim(cor_mat)
    #dim(cor_mat_keep)
    length(unique(cor_mat_keep$origin))
    cor_mat_keep <- cor_mat_keep[which(cor_mat_keep$BH < 0.05), ]
    length(unique(cor_mat_keep$origin))
    phenoaligner <- cor_mat_keep %>%
      group_by(origin) %>%
      dplyr::mutate(
        cor_max = max(cor),
        cor_max2 = sort(cor, decreasing = TRUE)[2]
      ) %>%
      dplyr::filter(cor == max(cor))
    map_data <- knn.lst_fake$coldata.ref
    phenoaligner <- left_join(phenoaligner, map_data, by = "cell")
    colnames(phenoaligner) <- c(
      "cor", "barcode", "pvalue", "BH", "origin", "cor_max", "cor_max2",
      "nn.tissue", "nn.patient", "nn.celltype", "nn.celltype2", "nn.celltype_sub"
    )
    map_data <- knn.lst_fake$coldata.query
    colnames(map_data)[1] <- "origin"
    phenoaligner <- left_join(phenoaligner, map_data, by = "origin")

    map_data <- knn.lst_fake$coldata.ref
    different_tissue <- left_join(cor_mat_keep, map_data, by = "cell")
    colnames(different_tissue) <- c(
      "cor", "barcode", "pvalue", "BH", "origin",
      "nn.tissue", "nn.patient", "nn.celltype", "nn.celltype2", "nn.celltype_sub"
    )
    mph_tissue_max_cor <- different_tissue %>%
      group_by(origin, nn.tissue) %>%
      summarise(tissue_max = max(cor), cells = n())
    mph_tissue_max_diff <- mph_tissue_max_cor %>%
      group_by(origin) %>%
      dplyr::mutate(
        top_one = max(tissue_max),
        top_two = sort(tissue_max, decreasing = TRUE)[2]
      )
    mph_tissue_max_diff <- data.frame(mph_tissue_max_diff)
    mph_tissue_max_diff$gap <- mph_tissue_max_diff$top_one - mph_tissue_max_diff$top_two

    mph_tissue_max_diff <- mph_tissue_max_diff[!duplicated(mph_tissue_max_diff$origin), ]
    #dim(mph_tissue_max_diff)
    # [1] 768   7
    mph_tissue_max_diff <- mph_tissue_max_diff[-which(mph_tissue_max_diff$gap < gap_threshold), ]
    #dim(mph_tissue_max_diff)
    ## [1] 440   7
    coldata.mt.fake <- phenoaligner[which(phenoaligner$origin %in% mph_tissue_max_diff$origin), ]
    fake.lst[[i]] = coldata.mt.fake
  }
  return(fake.lst)
}

.get_p_sub <- function(obs_mt = NULL, permutation_mt_lst = NULL, times = 1000) {
  obs_prop = obs_mt %>%
    group_by(nn.tissue, cellType3) %>%
    dplyr::summarise(cells = n()) %>%
    group_by(cellType3) %>%
    dplyr::mutate(frac = cells / sum(cells))
  permutation_lst = lapply(permutation_mt_lst, function(df) {
    permutation_prop =  df %>%
      group_by(nn.tissue, cellType3) %>%
      dplyr::summarise(cells = n()) %>%
      group_by(cellType3) %>%
      dplyr::mutate(frac = cells / sum(cells))
    left_join(obs_prop, permutation_prop, by = c('nn.tissue', 'cellType3')) %>%
      mutate(n = case_when(
        frac.y >= frac.x ~ 1,
        TRUE ~ 0
      )) %>%
      pull(n)
  })
  obs_prop$p = Reduce("+", permutation_lst)/times
  return(obs_prop)
}

.get_p_major <- function(obs_mt = NULL, permutation_mt_lst = NULL, times = 1000) {
  obs_prop = coldata.mt %>%
    group_by(nn.tissue, cellType) %>%
    dplyr::summarise(cells = n()) %>%
    group_by(cellType) %>%
    dplyr::mutate(frac = cells / sum(cells))
  permutation_lst = lapply(permutation_mt_lst, function(df) {
    permutation_prop =  df %>%
      group_by(nn.tissue, cellType) %>%
      dplyr::summarise(cells = n()) %>%
      group_by(cellType) %>%
      dplyr::mutate(frac = cells / sum(cells))
    left_join(obs_prop, permutation_prop, by = c('nn.tissue', 'cellType')) %>%
      mutate(n = case_when(
        frac.y >= frac.x ~ 1,
        TRUE ~ 0
      )) %>%
      pull(n)
  })
  Reduce("+", permutation_lst)
  obs_prop$p = Reduce("+", permutation_lst)/times
  return(obs_prop)
}


