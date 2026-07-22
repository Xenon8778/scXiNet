library(dplyr)
library(ggplot2)
library(XICOR)
library(reshape2)
library(pROC)
library(caret)
library(energy)
library(minet)
library(scTenifoldNet)
library(cowplot)
library(ggh4x)
source('code/XICOR_mod.R')
source('code/scripts.R')

# ============================================================
# Helpers
# ============================================================
computeMI <- function(mat, estimator = "mi.empirical",
                      disc = "equalfreq", nbins = NULL) {
  if (is.null(nbins)) nbins <- max(2L, floor(sqrt(nrow(mat))))
  mim <- build.mim(mat, estimator = estimator, disc = disc, nbins = nbins)
  mim[mim < 0] <- 0
  diag(mim) <- 0
  mim
}

computeDcor <- function(mat) {
  p   <- ncol(mat)
  res <- matrix(0.0, p, p, dimnames = list(colnames(mat), colnames(mat)))
  for (j in seq_len(p - 1L))
    for (k in (j + 1L):p)
      res[j, k] <- res[k, j] <- dcor(mat[, j], mat[, k])
  res
}

# basename() returns "" on trailing-slash paths — strip it first
ds_basename <- function(path) basename(sub("/$", "", path))

# ============================================================
# Dataset folders
# ============================================================
ds_dirs <- c(
  DS4 = "SERGIO/data_sets/De-noised_100G_3T_300cPerT_dynamics_9_DS4/",
  DS5 = "SERGIO/data_sets/De-noised_100G_4T_300cPerT_dynamics_10_DS5/",
  DS6 = "SERGIO/data_sets/De-noised_100G_6T_300cPerT_dynamics_7_DS6/"
)

cMetric <- 'a95'
sym     <- TRUE
niter   <- 50

pal <- setNames(
  hcl.colors(7, palette = 'RdPu', rev = TRUE),
  c('Pearson', 'Spearman', 'pcNet', 'GRNBoost', 'Distance', 'MI', 'Chatterjee')
)
alg_levels <- names(pal)
met_levels <- c('FPR', 'Sensitivity(TPR)', 'Balanced Accuracy', 'AUROC')

# ============================================================
# Ground Truth loader (dataset-specific)
# ============================================================
LoadGT <- function(dir, symmetric = TRUE) {
  df  <- t(read.csv(paste0(dir, 'simulated_noNoise_T_0.csv'),
                    header = TRUE, row.names = 1))
  val <- read.csv(paste0(dir, 'gt_GRN.csv'), header = FALSE)
  el  <- as.matrix(apply(as.matrix(val[, 1:2]), 2, as.character))

  val_adj <- matrix(0, ncol(df), ncol(df),
                    dimnames = list(colnames(df), colnames(df)))
  val_adj[el] <- 1
  if (symmetric) val_adj <- pmax(val_adj, t(val_adj))
  diag(val_adj) <- 0

  val_lg <- apply(val_adj, 2, as.logical)
  rownames(val_lg) <- rownames(val_adj)
  val_lg
}

# ============================================================
# Per-timepoint processing
# ============================================================
process_timepoint <- function(i, dir, val_adj_logical, val_adj,
                              train_size, test_size) {
  message(sprintf("  T_%d", i))

  df <- t(read.csv(paste0(dir, 'simulated_noNoise_T_', i, '.csv'),
                   header = TRUE, row.names = 1))

  set.seed(i)
  sub_rows <- sample(rownames(df), train_size)
  df_train <- df[sub_rows, ]
  df_test  <- df[setdiff(rownames(df), sub_rows), ]

  # ---- Cutoff estimation ----
  cutoffs_Xi       <- CalcCut(df_train, niter = niter, nsamp = test_size,
                              algorithm = 'Xi',       val_mat = val_adj_logical, symmetric = sym)
  cutoffs_pcnet    <- CalcCut(df_train, niter = niter, nsamp = test_size,
                              algorithm = 'pcNet',    val_mat = val_adj_logical)
  cutoffs_Pearson  <- CalcCut(df_train, niter = niter, nsamp = test_size,
                              algorithm = 'Pearson',  val_mat = val_adj_logical)
  cutoffs_Spearman <- CalcCut(df_train, niter = niter, nsamp = test_size,
                              algorithm = 'Spearman', val_mat = val_adj_logical)
  cutoffs_dcor     <- CalcCut(df_train, niter = niter, nsamp = test_size,
                              algorithm = 'dcor',     val_mat = val_adj_logical)
  cutoffs_MI       <- CalcCut(df_train, niter = niter, nsamp = test_size,
                              algorithm = 'MI',       val_mat = val_adj_logical)

  # ---- GRN construction on test data ----
  # Scaling: all methods except MI are natively in [0,1] after abs().
  # MI is excluded from RMSE; AUROC uses raw scores (rank-invariant).
  cor_xi_test <- abs(xicor_mpar(df_test, nCores = 4))
  if (sym) cor_xi_test <- pmax(cor_xi_test, t(cor_xi_test))
  diag(cor_xi_test) <- 0
  cor_xi_filtered <- (cor_xi_test > cutoffs_Xi[cMetric]) * 1L

  cor_pcnet_test <- abs(as.matrix(pcNet(t(df_test))))
  diag(cor_pcnet_test) <- 0
  cor_pcnet_filtered <- (cor_pcnet_test > cutoffs_pcnet[cMetric]) * 1L

  cor_Spearman_test <- abs(cor(df_test, method = 'spearman'))
  diag(cor_Spearman_test) <- 0
  cor_Spearman_filtered <- (cor_Spearman_test > cutoffs_Spearman[cMetric]) * 1L

  cor_Pearson_test <- abs(cor(df_test, method = 'pearson'))
  diag(cor_Pearson_test) <- 0
  cor_Pearson_filtered <- (cor_Pearson_test > cutoffs_Pearson[cMetric]) * 1L

  cor_dcor_test <- computeDcor(df_test)
  diag(cor_dcor_test) <- 0
  cor_dcor_filtered <- (cor_dcor_test > cutoffs_dcor[cMetric]) * 1L

  cor_MI_test <- computeMI(df_test)
  diag(cor_MI_test) <- 0
  cor_MI_filtered <- (cor_MI_test > cutoffs_MI[cMetric]) * 1L

  # ---- Confusion matrices ----
  cm_xi       <- confusionMatrix(table(cor_xi_filtered,       val_adj), positive = "1")
  cm_pcNet    <- confusionMatrix(table(cor_pcnet_filtered,    val_adj), positive = "1")
  cm_spearman <- confusionMatrix(table(cor_Spearman_filtered, val_adj), positive = "1")
  cm_pearson  <- confusionMatrix(table(cor_Pearson_filtered,  val_adj), positive = "1")
  cm_dcor     <- confusionMatrix(table(cor_dcor_filtered,     val_adj), positive = "1")
  cm_MI       <- confusionMatrix(table(cor_MI_filtered,       val_adj), positive = "1")

  # ---- Metrics (RMSE excluded for MI; AUROC from raw scores) ----
  metrics_list <- list()
  add_metrics <- function(alg, cm, scores, cutoff, rmse = TRUE) {
    spec <- cm$byClass['Specificity']
    rows <- list(
      c(Algorithm = alg, Metric = 'Sensitivity(TPR)',  Value = unname(cm$byClass['Sensitivity'])),
      c(Algorithm = alg, Metric = 'Specificity(TNR)',  Value = unname(spec)),
      c(Algorithm = alg, Metric = 'FPR',               Value = unname(1 - spec)),
      c(Algorithm = alg, Metric = 'Balanced Accuracy', Value = unname(cm$byClass['Balanced Accuracy'])),
      c(Algorithm = alg, Metric = 'AUROC',
        Value = as.numeric(auc(roc(as.numeric(val_adj), as.numeric(scores), quiet = TRUE)))),
      c(Algorithm = alg, Metric = 'Cutoff', Value = unname(cutoff))
    )
    if (rmse) rows <- c(rows, list(
      c(Algorithm = alg, Metric = 'RMSE',
        Value = RMSE(pred = c(scores), obs = c(val_adj), na.rm = TRUE))
    ))
    metrics_list <<- c(metrics_list, rows)
  }

  add_metrics('Pearson',    cm_pearson,  cor_Pearson_test,  cutoffs_Pearson[cMetric])
  add_metrics('Spearman',   cm_spearman, cor_Spearman_test, cutoffs_Spearman[cMetric])
  add_metrics('Chatterjee', cm_xi,       cor_xi_test,       cutoffs_Xi[cMetric])
  add_metrics('pcNet',      cm_pcNet,    cor_pcnet_test,    cutoffs_pcnet[cMetric])
  add_metrics('Distance',   cm_dcor,     cor_dcor_test,     cutoffs_dcor[cMetric])
  add_metrics('MI',         cm_MI,       cor_MI_test,       cutoffs_MI[cMetric], rmse = FALSE)

  metrics_df <- as.data.frame(do.call(rbind, metrics_list), stringsAsFactors = FALSE)
  metrics_df$Value <- as.numeric(metrics_df$Value)
  metrics_df$Iter  <- i

  # ---- Score distributions for cutoff plot ----
  ut <- upper.tri(cor_xi_test)
  scores_df <- rbind(
    data.frame(Algorithm = 'Pearson',    Score = cor_Pearson_test[ut],  Cutoff = cutoffs_Pearson[cMetric]),
    data.frame(Algorithm = 'Spearman',   Score = cor_Spearman_test[ut], Cutoff = cutoffs_Spearman[cMetric]),
    data.frame(Algorithm = 'Chatterjee', Score = cor_xi_test[ut],       Cutoff = cutoffs_Xi[cMetric]),
    data.frame(Algorithm = 'pcNet',      Score = cor_pcnet_test[ut],    Cutoff = cutoffs_pcnet[cMetric]),
    data.frame(Algorithm = 'Distance',   Score = cor_dcor_test[ut],     Cutoff = cutoffs_dcor[cMetric]),
    data.frame(Algorithm = 'MI',         Score = cor_MI_test[ut],       Cutoff = cutoffs_MI[cMetric])
  )
  scores_df$Iter <- i

  list(metrics = metrics_df, scores = scores_df)
}

# ============================================================
# Main loop: iterate over DS folders
# ============================================================
all_metrics <- list()
all_scores  <- list()

for (ds_name in names(ds_dirs)) {
  dir <- ds_dirs[[ds_name]]
  message(sprintf("\n====== %s ======", ds_name))

  val_adj_logical <- LoadGT(dir, symmetric = sym)
  val_adj <- matrix(as.integer(val_adj_logical), nrow = nrow(val_adj_logical))

  df0        <- t(read.csv(paste0(dir, 'simulated_noNoise_T_0.csv'), header = TRUE, row.names = 1))
  train_size <- round(nrow(df0) * 3 / 5)
  test_size  <- nrow(df0) - train_size
  n_timepoints <- length(list.files(dir, pattern = 'simulated_noNoise_T_\\d+\\.csv')) - 1

  results <- lapply(0:n_timepoints, process_timepoint,
                    dir            = dir,
                    val_adj_logical = val_adj_logical,
                    val_adj        = val_adj,
                    train_size     = train_size,
                    test_size      = test_size)

  metrics_df <- do.call(rbind, lapply(results, `[[`, 'metrics'))
  scores_df  <- do.call(rbind, lapply(results, `[[`, 'scores'))
  metrics_df$DS <- ds_name
  scores_df$DS  <- ds_name

  res_tag <- paste0(ds_basename(dir), '_', cMetric, '_GTSym', sym)
  write.csv(metrics_df,
            paste0('Results/Benchmark/Cutoff_Selection_', res_tag, '_ver2.csv'),
            row.names = FALSE)

  all_metrics[[ds_name]] <- metrics_df
  all_scores[[ds_name]]  <- scores_df
}

metrics_all <- do.call(rbind, all_metrics)
scores_all  <- do.call(rbind, all_scores)

# ============================================================
# Merge with GRNBoost and plot — one figure per DS
# ============================================================
for (ds_name in names(ds_dirs)) {
  dir     <- ds_dirs[[ds_name]]
  res_tag <- paste0(ds_basename(dir), '_', cMetric, '_GTSym', sym)

  metrics_df  <- read.csv(paste0('Results/Benchmark/Cutoff_Selection_', res_tag, '_ver2.csv'))
  grn_path    <- paste0('Results/Benchmark/grnboost/Cutoff_Selection_',
                        ds_basename(dir), '_Symmetric.csv')
  grnboost_df <- read.csv(grn_path)
  colnames(grnboost_df) <- c('Iter', 'Algorithm', 'Sensitivity(TPR)', 'Specificity(TNR)',
                             'Balanced Accuracy', 'AUROC', 'Cutoff', 'Symmetric')
  grnboost_df <- grnboost_df %>%
    mutate(FPR = 1 - `Specificity(TNR)`) %>%
    reshape2::melt(id.vars = c('Iter', 'Algorithm')) %>%
    setNames(c('Iter', 'Algorithm', 'Metric', 'Value')) %>%
    select(Algorithm, Metric, Value, Iter)

  merged_df <- rbind(metrics_df %>% select(Algorithm, Metric, Value, Iter),
                     grnboost_df) %>%
    filter(Metric %in% met_levels) %>%
    mutate(Algorithm = factor(Algorithm, levels = alg_levels),
           Metric    = factor(Metric,    levels = met_levels),
           Value     = as.numeric(Value))

  write.csv(merged_df,
            paste0('Results/Cutoff_Selection_SERGIO_', res_tag, '_Merged_ver2.csv'),
            row.names = FALSE)

  # ---- Benchmark plot ----
  p_bench <- ggplot(merged_df, aes(x = Algorithm, y = Value, fill = Algorithm)) +
    geom_boxplot(width = 0.7, outlier.size = 0.5, linewidth = 0.5, colour = "black") +
    geom_jitter(width = 0.1, size = 0.5) +
    facet_wrap(~Metric, ncol = 4, scales = 'free') +
    labs(x = NULL, y = NULL, title = NULL) +
    theme_classic() +
    theme(legend.position = "none",
          text            = element_text(family = "Arial"),
          plot.title      = element_text(hjust = 0.5, face = 'bold'),
          axis.text.x     = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                         family = 'Arial', colour = 'black'),
          axis.text.y     = element_text(family = 'Arial', colour = 'black')) +
    scale_fill_manual(values = pal)+
    facetted_pos_scales(
      y = list(
        Metric == "AUROC" ~ scale_y_continuous(limits = c(0.6, 0.9))
      )
    )

  png(paste0('Results/Benchmark/Figures/Bench_SERGIO_', res_tag, '_ver2.png'),
      height = 3, width = 8, units = 'in', res = 300)
  print(p_bench)
  dev.off()

  # ---- RMSE subplot (MI excluded) ----
  rmse_algs <- c('Pearson', 'Spearman', 'Distance', 'Chatterjee')
  RMSE_df <- metrics_df %>%
    filter(Algorithm %in% rmse_algs, Metric == 'RMSE') %>%
    mutate(Algorithm = factor(Algorithm, levels = rmse_algs),
           Value     = as.numeric(Value))

  p_rmse <- ggplot(RMSE_df, aes(x = Algorithm, y = Value, fill = Algorithm)) +
    geom_boxplot(width = 0.7, outlier.size = 0.5, linewidth = 0.5, colour = "black") +
    geom_jitter(width = 0.1, size = 0.5) +
    facet_wrap(~Metric, scales = 'free') +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(legend.position = "none",
          text        = element_text(family = "Arial"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                     family = 'Arial', colour = 'black'),
          axis.text.y = element_text(family = 'Arial', colour = 'black')) +
    scale_fill_manual(values = pal[rmse_algs])

  png(paste0('Results/Benchmark/Figures/Bench_SERGIO_corrRMSE_', res_tag, '_ver2.png'),
      height = 3, width = 2, units = 'in', res = 300)
  print(p_rmse)
  dev.off()

  # ---- Score distribution plot ----
  sc_df <- scores_all %>% filter(DS == ds_name)

  # Merge GRNBoost score distributions if available
  grn_sc_path <- paste0('Results/Benchmark/grnboost/ScoreDist_',
                        ds_basename(dir), '_Symmetric.csv')
  if (file.exists(grn_sc_path)){
    gboost_df <- read.csv(grn_sc_path) %>%
      mutate(DS = ds_name)
    sc_df <- rbind(sc_df, gboost_df)}

  dist_alg_order <- c('Pearson', 'Spearman', 'pcNet', 'Distance', 'MI', 'GRNBoost', 'Chatterjee')
  cutoff_lines <- sc_df %>%
    filter(Algorithm %in% dist_alg_order) %>%
    group_by(Algorithm) %>%
    summarise(Cutoff = median(Cutoff), .groups = 'drop') %>%
    mutate(Algorithm = factor(Algorithm, levels = dist_alg_order))

  p_dist <- ggplot(sc_df %>% filter(Algorithm %in% dist_alg_order) %>%
                     mutate(Algorithm = factor(Algorithm, levels = dist_alg_order)),
                   aes(x = Score, fill = Algorithm)) +
    geom_density(alpha = 0.7, colour = NA) +
    geom_vline(data = cutoff_lines, aes(xintercept = Cutoff),
               linetype = 'dashed', linewidth = 0.6, colour = 'black') +
    facet_wrap(~Algorithm, ncol = 4, scales = 'free') +
    labs(x = 'Score', y = 'Density', title = ds_name) +
    theme_classic() +
    theme(legend.position  = "none",
          text             = element_text(family = "Arial"),
          plot.title       = element_text(hjust = 0.5, face = 'bold'),
          axis.text.x      = element_text(family = 'Arial', colour = 'black'),
          axis.text.y      = element_text(family = 'Arial', colour = 'black'),
          strip.background = element_blank(),
          strip.text       = element_text(family = 'Arial', face = 'bold', colour = 'black')) +
    scale_fill_manual(values = pal[dist_alg_order])

  png(paste0('Results/Benchmark/Figures/ScoreDist_SERGIO_', res_tag, '_ver2.png'),
      height = 4, width = 8, units = 'in', res = 300)
  print(p_dist)
  dev.off()
}

message("Done. All outputs written with '_ver2' suffix.")
