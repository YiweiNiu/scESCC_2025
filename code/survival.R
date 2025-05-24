
# weighted z
weighted_z <- function(z = NULL, w = NULL) {
  bot = sqrt(sum(w^2))
  top = sum(z*w)
  return(top/bot)
}

# cox for overall survival
coxph_batch <- function(fpkm=NULL, clinical=NULL, markers_top=NULL) {

  # helper function for coxph
  my_cox <- function(i, fpkm=NULL, clinical=NULL, markers_top=NULL) {
    # ave exp of each sample
    exp = colMeans(fpkm[rownames(fpkm) %in% markers_top[[as.character(i)]],])
    # prep for survival
    surv.os = tibble(sampleID = names(exp),
                     value = as.numeric(exp)) %>%
      left_join(clinical, by = "sampleID") %>%
      filter(sample_type == "Primary Tumor") %>%
      dplyr::select(sampleID, OS.time, OS, value, gender, age_at_initial_pathologic_diagnosis, pathologic_stage) %>%
      mutate(pathologic_stage = as.numeric(as.factor(pathologic_stage))) %>%
      mutate(gender = as.numeric(as.factor(gender)))

    formula = as.formula("Surv(OS.time, OS) ~ value+gender+age_at_initial_pathologic_diagnosis+pathologic_stage")
    os.cox <- coxph(formula, data = surv.os)
    cox_res = summary(os.cox)
    return (c(
      z = cox_res$coefficients[1, 'z'],
      se = cox_res$coefficients[1, 'se(coef)'],
      p = cox_res$coefficients[1, 'Pr(>|z|)'],
      hr = cox_res$conf.int[1, 1],
      hr_lower = cox_res$conf.int[1, 3],
      hr_upper = cox_res$conf.int[1, 4]
    ))
  }

  res = sapply(as.numeric(names(markers_top)), FUN = my_cox, fpkm=fpkm, clinical=clinical, markers_top=markers_top)

  x = as.data.frame(t(res)) %>%
    mutate(cluster = factor(as.numeric(names(markers_top)), levels = names(markers_top))) %>%
    mutate(q = p.adjust(p, method = 'BH')) %>%
    dplyr::select(cluster, z, se, p, q, hr, hr_lower, hr_upper)

  # plot
  p2col <- function(p_vec) {
    q_vec = p.adjust(p_vec, method = 'BH')
    m = cbind(p_vec, q_vec)
    s <- function(x) {
      if (x[2] < 0.05) {'black'}
      else if (x[2] >= 0.05 & x[1] < 0.05) {'#36648b'}
      else {'#5cacee'}
    }
    apply(m, 1, s)
  }

  # z limit
  max_z = max(abs(x$z) + 0.5)
  # nominal p
  pos_p_lim = min(x[x$z > 0 & x$p < 0.05, 2])
  neg_p_lim = max(x[x$z < 0 & x$p < 0.05, 2])

  # FDR p
  pos_q_lim = min(x[x$z > 0 & x$q < 0.05, 2])
  neg_q_lim = min(x[x$z < 0 & x$q < 0.05, 2])

  p = ggplot(x, aes(x = cluster, y = z)) +
    geom_segment(aes(x = cluster, xend = cluster, y = 0, yend = z),
                 linetype = 'dashed', color = '#36648b', size = 0.8) +
    geom_hline(aes(yintercept = pos_p_lim), color = '#bebebe', linetype = 'dashed') +
    geom_hline(aes(yintercept = neg_p_lim), color = '#bebebe', linetype = 'dashed') +
    geom_hline(aes(yintercept = pos_q_lim), color = '#6a6a6a') +
    geom_hline(aes(yintercept = neg_q_lim), color = '#6a6a6a') +
    ylim(-max_z, max_z) +
    geom_point(size = 3, color = p2col(x$p)) +
    coord_flip() +
    theme_cowplot() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=1, fill = NA))
  return(list(x = x, p = p))
}

# TCGA: plot boxplot for exp in normal/tumor samples, for given cluster
plot_tcga_exp_box <- function(i, fpkm=NULL, clinical=NULL, markers_top=NULL) {
  exp = colMeans(fpkm[rownames(fpkm) %in% markers_top[[as.character(i)]], ])
  tb = tibble(sampleID = names(exp), value = as.numeric(exp)) %>%
    left_join(clinical, by = "sampleID") %>%
    mutate(sample_type = case_when(
      sample_type == 'Solid Tissue Normal' ~ 'Esophagus',
      sample_type == 'Primary Tumor' ~ 'ESCC'
    )) %>%
    mutate(sample_type = factor(sample_type,
                                levels = c('Esophagus', 'ESCC'))) %>%
    dplyr::select(sampleID, sample_type, value)

  p1 = ggplot(tb, aes(x = sample_type, y = value, fill = sample_type)) +
    geom_boxplot(outlier.alpha = .1, linetype = "dashed") +
    stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.alpha = .1) +
    stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = .2) +
    stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = .2) +
    labs(x = NULL, y = "FPKM", title = paste0("cluster ", i)) +
    stat_compare_means(comparisons = list(c("Esophagus", "ESCC"))) +
    scale_fill_manual(values = c('#7ccd7c', '#36648b')) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 60, hjust = 1))
  p1
}

# TCGA: plot boxplot for exp in normal/tumor samples, for given gene
plot_tcga_exp_box_single <- function(gene=NULL, fpkm=NULL, clinical=NULL) {
  if (!(gene %in% rownames(fpkm))) {
    cat("Gene not found")
    return (NULL)
  }
  exp = fpkm[gene,]
  tb = tibble(sampleID = names(exp), value = as.numeric(exp)) %>%
    left_join(clinical, by = "sampleID") %>%
    mutate(sample_type = case_when(
      sample_type == 'Solid Tissue Normal' ~ 'Esophagus',
      sample_type == 'Primary Tumor' ~ 'ESCC'
    )) %>%
    mutate(sample_type = factor(sample_type,
                                levels = c('Esophagus', 'ESCC'))) %>%
    dplyr::select(sampleID, sample_type, value)

  p1 = ggplot(tb, aes(x = sample_type, y = value, fill = sample_type)) +
    geom_boxplot(outlier.alpha = .1, linetype = "dashed") +
    stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.alpha = .1) +
    stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = .2) +
    stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = .2) +
    labs(x = NULL, y = "FPKM", title = gene) +
    stat_compare_means(comparisons = list(c("Esophagus", "ESCC"))) +
    scale_fill_manual(values = c('#7ccd7c', '#36648b')) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 60, hjust = 1))
  return(p1)
}

# TCGA: KM plot for OS/RFS, for average marker exp in given cluster
plot_tcga_km <- function(i, fpkm=NULL, clinical=NULL, markers_top=NULL) {
  exp = colMeans(fpkm[rownames(fpkm) %in% markers_top[[as.character(i)]], ])
  # OS
  surv.os = tibble(sampleID = names(exp), value = as.numeric(exp)) %>%
    left_join(clinical, by = "sampleID") %>%
    filter(sample_type == "Primary Tumor") %>%
    dplyr::select(sampleID, sample_type, value, OS.time, OS) %>%
    mutate(group = case_when(value > quantile(value, 0.5) ~ "High",
                             value < quantile(value, 0.5) ~ "Low",
                             TRUE ~ NA_character_))
  # RFS
  surv.rfs = tibble(sampleID = names(exp), value = as.numeric(exp)) %>%
    left_join(clinical, by = "sampleID") %>%
    filter(sample_type == "Primary Tumor") %>%
    dplyr::select(sampleID, sample_type, value, RFS.time, RFS) %>%
    mutate(group = case_when(value > quantile(value, 0.5) ~ "High",
                             value < quantile(value, 0.5) ~ "Low",
                             TRUE ~ NA_character_))

  fit.os = survfit(Surv(OS.time, OS) ~ group, data = surv.os)
  fit.rfs = survfit(Surv(RFS.time, RFS) ~ group, data = surv.rfs)
  ggsurv.os = ggsurvplot(fit = fit.os,
                         data = surv.os,
                         pval = TRUE,
                         xlab = "Time in days",
                         title = 'OS',
                         legend.title = '',
                         legend.labs = c('High', 'Low'),
                         ggtheme = theme(axis.text.x = element_text(angle = 60, hjust = 1)))

  ggsurv.rfs = ggsurvplot(fit = fit.rfs,
                          data = surv.rfs,
                          pval = TRUE,
                          xlab = "Time in days",
                          title = 'RFS',
                          legend.title = '',
                          legend.labs = c('High', 'Low'),
                          ggtheme = theme(axis.text.x = element_text(angle = 60, hjust = 1)))

  return(list(ggsurv.os = ggsurv.os, ggsurv.rfs = ggsurv.rfs))
}

# TCGA: KM plot for OS/RFS, for single gene
plot_tcga_km_single <- function(gene=NULL, fpkm=NULL, clinical=NULL) {
  if (!(gene %in% rownames(fpkm))) {
    cat("Gene not found")
    return (list(ggsurv.os = NULL, ggsurv.rfs = NULL))
  }
  exp = fpkm[gene,]
  if (length(exp[exp > 0]) < 50) {
    cat("Low exp")
    return (list(ggsurv.os = NULL, ggsurv.rfs = NULL))
  }
  # OS
  surv.os = tibble(sampleID = names(exp), value = as.numeric(exp)) %>%
    left_join(clinical, by = "sampleID") %>%
    filter(sample_type == "Primary Tumor") %>%
    dplyr::select(sampleID, sample_type, value, OS.time, OS) %>%
    mutate(group = case_when(value > quantile(value, 0.5) ~ "High",
                             value < quantile(value, 0.5) ~ "Low",
                             TRUE ~ NA_character_))
  # RFS
  surv.rfs = tibble(sampleID = names(exp), value = as.numeric(exp)) %>%
    left_join(clinical, by = "sampleID") %>%
    filter(sample_type == "Primary Tumor") %>%
    dplyr::select(sampleID, sample_type, value, RFS.time, RFS) %>%
    mutate(group = case_when(value > quantile(value, 0.5) ~ "High",
                             value < quantile(value, 0.5) ~ "Low",
                             TRUE ~ NA_character_))

  fit.os = survfit(Surv(OS.time, OS) ~ group, data = surv.os)
  fit.rfs = survfit(Surv(RFS.time, RFS) ~ group, data = surv.rfs)
  ggsurv.os = ggsurvplot(fit = fit.os,
                         data = surv.os,
                         pval = TRUE,
                         xlab = "Time in days",
                         title = 'OS',
                         legend.title = '',
                         legend.labs = c('High', 'Low'),
                         ggtheme = theme(axis.text.x = element_text(angle = 60, hjust = 1)))

  ggsurv.rfs = ggsurvplot(fit = fit.rfs,
                          data = surv.rfs,
                          pval = TRUE,
                          xlab = "Time in days",
                          title = 'RFS',
                          legend.title = '',
                          legend.labs = c('High', 'Low'),
                          ggtheme = theme(axis.text.x = element_text(angle = 60, hjust = 1)))

  return(list(ggsurv.os = ggsurv.os, ggsurv.rfs = ggsurv.rfs))
}

# TCGA: log-rank test for single gene
get_tcga_km_single <- function(gene=NULL, fpkm=NULL, clinical=NULL) {
  if (!(gene %in% rownames(fpkm))) {
    cat("Gene not found")
    return (c(NA, NA))
  }
  #print (gene)
  exp = fpkm[gene,]
  if (length(exp[exp > 0]) < 50) {
    cat("Low exp")
    return (c(NA, NA))
  }
  # OS
  surv.os = tibble(sampleID = names(exp), value = as.numeric(exp)) %>%
    left_join(clinical, by = "sampleID") %>%
    filter(sample_type == "Primary Tumor") %>%
    dplyr::select(sampleID, sample_type, value, OS.time, OS) %>%
    mutate(group = case_when(value > quantile(value, 0.5) ~ "High",
                             value < quantile(value, 0.5) ~ "Low",
                             TRUE ~ NA_character_))
  # RFS
  surv.rfs = tibble(sampleID = names(exp), value = as.numeric(exp)) %>%
    left_join(clinical, by = "sampleID") %>%
    filter(sample_type == "Primary Tumor") %>%
    dplyr::select(sampleID, sample_type, value, RFS.time, RFS) %>%
    mutate(group = case_when(value > quantile(value, 0.5) ~ "High",
                             value < quantile(value, 0.5) ~ "Low",
                             TRUE ~ NA_character_))

  fit.os = survfit(Surv(OS.time, OS) ~ group, data = surv.os)
  fit.rfs = survfit(Surv(RFS.time, RFS) ~ group, data = surv.rfs)

  diff.os = survdiff(Surv(OS.time, OS) ~ group, data = surv.os)
  diff.rfs = survdiff(Surv(RFS.time, RFS) ~ group, data = surv.rfs)

  p.os = pchisq(diff.os$chisq, length(diff.os$n)-1, lower.tail = FALSE)
  p.rfs = pchisq(diff.rfs$chisq, length(diff.rfs$n)-1, lower.tail = FALSE)
  return (c(p.os, p.rfs))
}

# Song: KM plot for OS, for average marker exp in given cluster
plot_song_km <- function(i, fpkm=NULL, clinical=NULL, markers_top=NULL) {
  exp = colMeans(fpkm[rownames(fpkm) %in% markers_top[[as.character(i)]], ])
  # OS
  surv.os = tibble(sampleID = names(exp), value = as.numeric(exp)) %>%
    left_join(clinical, by = "sampleID") %>%
    filter(sample_type == "Primary Tumor") %>%
    dplyr::select(sampleID, sample_type, value, OS.time, OS) %>%
    mutate(group = case_when(value > quantile(value, 0.5) ~ "High",
                             value < quantile(value, 0.5) ~ "Low",
                             TRUE ~ NA_character_))

  fit.os = survfit(Surv(OS.time, OS) ~ group, data = surv.os)
  ggsurv.os = ggsurvplot(fit = fit.os,
                         data = surv.os,
                         pval = TRUE,
                         xlab = "Time in days",
                         title = 'OS',
                         legend.title = '',
                         legend.labs = c('High', 'Low'),
                         ggtheme = theme(axis.text.x = element_text(angle = 60, hjust = 1)))

  return(ggsurv.os)
}

# Song: KM plot for OS, for single gene
plot_song_km_single <- function(gene=NULL, fpkm=NULL, clinical=NULL) {
  if (!(gene %in% rownames(fpkm))) {
    cat("Gene not found")
    return (NULL)
  }
  exp = fpkm[gene,]
  if (length(exp[exp == 0]) > 50) {
    cat("Low exp")
    return (NULL)
  }
  # OS
  surv.os = tibble(sampleID = names(exp), value = as.numeric(exp)) %>%
    left_join(clinical, by = "sampleID") %>%
    filter(sample_type == "Primary Tumor") %>%
    dplyr::select(sampleID, sample_type, value, OS.time, OS) %>%
    mutate(group = case_when(value > quantile(value, 0.5) ~ "High",
                             value < quantile(value, 0.5) ~ "Low",
                             TRUE ~ NA_character_))

  fit.os = survfit(Surv(OS.time, OS) ~ group, data = surv.os)
  ggsurv.os = ggsurvplot(fit = fit.os,
                         data = surv.os,
                         pval = TRUE,
                         xlab = "Time in days",
                         title = 'OS',
                         legend.title = '',
                         legend.labs = c('High', 'Low'),
                         ggtheme = theme(axis.text.x = element_text(angle = 60, hjust = 1)))

  return(ggsurv.os)
}

# Song: log-rank test for single gene
get_song_km_single <- function(gene=NULL, fpkm=NULL, clinical=NULL) {
  if (!(gene %in% rownames(fpkm))) {
    cat("Gene not found")
    return (NA)
  }
  exp = fpkm[gene,]
  if (length(exp[exp == 0]) > 50) {
    cat("Low exp")
    return (NA)
  }
  #print(gene)
  # OS
  surv.os = tibble(sampleID = names(exp), value = as.numeric(exp)) %>%
    left_join(clinical, by = "sampleID") %>%
    filter(sample_type == "Primary Tumor") %>%
    dplyr::select(sampleID, sample_type, value, OS.time, OS) %>%
    mutate(group = case_when(value > quantile(value, 0.5) ~ "High",
                             value < quantile(value, 0.5) ~ "Low",
                             TRUE ~ NA_character_))

  fit.os = survfit(Surv(OS.time, OS) ~ group, data = surv.os)
  diff.os = survdiff(Surv(OS.time, OS) ~ group, data = surv.os)

  p.os = pchisq(diff.os$chisq, length(diff.os$n)-1, lower.tail = FALSE)
  return (p.os)
}

