library(dplyr)
library(ggplot2)
library(XICOR)
library(reshape2)
library(pROC)
library(energy)
library(minet)
library(scTenifoldNet)
library(ggh4x)
source('code/XICOR_mod.R')
source('code/cutoffCalcScripts.R')   # requires updated CalcCut (FPR + FDR cutoffs)

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

# Distance correlation — upper triangle only, then mirrored
computeDcor <- function(mat, nCores = 1, cl = NULL) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  p   <- ncol(mat)

  rowFun <- function(j) {
    x <- mat[, j]
    vapply((j + 1L):p, function(k) energy::dcor2d(x, mat[, k], type = "V"), numeric(1))
  }

  if (nCores > 1 || !is.null(cl)) {
    own <- is.null(cl)
    if (own) cl <- parallel::makeCluster(nCores)
    on.exit(if (own) parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, c("mat", "p"), envir = environment())
    rows <- parallel::parLapplyLB(cl, seq_len(p - 1L), rowFun)  # LB: rows shrink with j
  } else {
    rows <- lapply(seq_len(p - 1L), rowFun)
  }

  res <- matrix(0.0, p, p, dimnames = list(colnames(mat), colnames(mat)))
  for (j in seq_len(p - 1L)) res[j, (j + 1L):p] <- rows[[j]]
  res <- res + t(res)
  sqrt(pmax(res, 0))   # zero-variance genes give NaN -> set to 0 in CalcCut
}

# basename() returns "" on trailing-slash paths — strip it first
ds_basename <- function(path) basename(sub("/$", "", path))

# ============================================================
# Parameters
# ============================================================
ds_dirs <- c(
  DS4 = "SERGIO/data_sets/De-noised_100G_3T_300cPerT_dynamics_9_DS4/",
  DS5 = "SERGIO/data_sets/De-noised_100G_4T_300cPerT_dynamics_10_DS5/",
  DS6 = "SERGIO/data_sets/De-noised_100G_6T_300cPerT_dynamics_7_DS6/"
)

niter     <- 50
alpha     <- 0.05
sym       <- TRUE
# Both cutoffs from CalcCut() are evaluated on every test network.
# CalcCut names its outputs from alpha (alpha = 0.2 -> a80, FDR80, FDR80_feasible),
# so the keys are derived here rather than hard-coded.
cut_sfx      <- as.character(round((1 - alpha) * 100))
cut_types    <- c(FPR = paste0('a', cut_sfx), FDR = paste0('FDR', cut_sfx))
feasible_key <- paste0('FDR', cut_sfx, '_feasible')
minach_key   <- paste0('FDR', cut_sfx, '_minAchievable')

alg_levels <- c('Pearson', 'Spearman', 'pcNet', 'GRNBoost', 'Distance', 'MI',
                'Chatterjee', 'ChatterjeeKNN')
xi_M <- 5   # Lin & Han (2023) revised xi: number of right nearest neighbours
n_cores    <- 4

met_levels <- c('FPR', 'FDR', 'Sensitivity(TPR)', 'Balanced Accuracy', 'RMSE', 'AUROC',
                'FDR_feasible', 'FDR_minAchievable')
pal <- setNames(
  hcl.colors(length(alg_levels) + 2, palette = 'RdPu', rev = TRUE)[seq_along(alg_levels)],
  alg_levels
)

# ============================================================
# Ground Truth loader (dataset-specific; gene order from T_0)
# ============================================================
LoadGT <- function(dir, symmetric = TRUE) {
  # Three-class GT: TRUE = direct edge, FALSE = no trek (true null),
  # NA = indirect / co-regulated (excluded from calibration and evaluation)
  genes <- colnames(t(read.csv(paste0(dir, 'simulated_noNoise_T_0.csv'),
                               header = TRUE, row.names = 1)))
  val <- read.csv(paste0(dir, 'gt_GRN.csv'), header = FALSE)
  el  <- data.frame(from = as.character(val[, 1]), to = as.character(val[, 2]))

  g <- igraph::graph_from_data_frame(el, directed = TRUE, vertices = genes)
  stopifnot(igraph::is_dag(g))
  R  <- is.finite(igraph::distances(g, v = genes, to = genes, mode = 'out')) * 1
  CA <- crossprod(R) > 0                       # common ancestor or ancestral relation

  direct <- igraph::as_adjacency_matrix(g, sparse = FALSE)[genes, genes] > 0
  if (symmetric) direct <- direct | t(direct)

  cls <- matrix(NA, length(genes), length(genes), dimnames = list(genes, genes))
  cls[!CA]    <- FALSE
  cls[direct] <- TRUE
  diag(cls)   <- NA
  ut <- upper.tri(cls)
  message(sprintf('  pos=%d  null=%d  excluded=%d', sum(cls[ut], na.rm = TRUE),
                  sum(!cls[ut], na.rm = TRUE), sum(is.na(cls[ut]))))
  cls
}

# ============================================================
# Evaluation helper (identical logic to CutoffSelect_Nonlinear ver3)
#
# Metrics over unique gene pairs (j < k) when GT is symmetric.
#   FPR = FP / (FP + TN)   targeted by 'a<..>'   (CutoffType = FPR)
#   FDR = FP / (FP + TP)   targeted by 'FDR<..>' (CutoffType = FDR)
#   FDR := NA when no pair is called (infeasible cutoff; see FDR_feasible)
# ============================================================
evaluate_method <- function(alg, scores, cuts, eval_mask, labels, rmse = TRUE) {
  s      <- as.numeric(scores[eval_mask])
  auroc  <- as.numeric(auc(roc(labels, s, quiet = TRUE)))
  rmse_v <- if (rmse) sqrt(mean((s - labels)^2, na.rm = TRUE)) else NA_real_

  do.call(rbind, lapply(names(cut_types), function(ct) {
    cutoff <- as.numeric(cuts[cut_types[[ct]]])
    pred   <- s > cutoff                      # Inf cutoff -> nothing called

    TP <- sum(pred & labels == 1L); FP <- sum(pred & labels == 0L)
    TN <- sum(!pred & labels == 0L); FN <- sum(!pred & labels == 1L)
    tpr <- TP / (TP + FN)
    tnr <- TN / (TN + FP)
    fdr <- if ((TP + FP) == 0) NA_real_ else FP / (TP + FP)

    out <- data.frame(
      Algorithm  = alg,
      CutoffType = ct,
      Metric = c('Sensitivity(TPR)', 'Specificity(TNR)', 'FPR', 'FDR',
                 'Balanced Accuracy', 'AUROC', 'Cutoff', 'nCalled'),
      Value  = c(tpr, tnr, 1 - tnr, fdr, (tpr + tnr) / 2, auroc, cutoff, TP + FP),
      stringsAsFactors = FALSE
    )
    if (ct == 'FDR') out <- rbind(out, data.frame(
      Algorithm = alg, CutoffType = ct,
      Metric = c('FDR_feasible', 'FDR_minAchievable'),
      Value  = unname(cuts[c(feasible_key, minach_key)])))
    if (rmse) out <- rbind(out, data.frame(
      Algorithm = alg, CutoffType = ct, Metric = 'RMSE', Value = rmse_v))
    out
  }))
}

# ============================================================
# Per-replicate (T_i) processing
# ============================================================
process_timepoint <- function(i, dir, val_adj_logical, eval_mask, labels,
                              train_size, test_size) {
  message(sprintf("  T_%d", i))

  df <- t(read.csv(paste0(dir, 'simulated_noNoise_T_', i, '.csv'),
                   header = TRUE, row.names = 1))
  df <- df[, colnames(val_adj_logical)]     # enforce GT gene order

  set.seed(i)
  sub_rows <- sample(rownames(df), train_size)
  df_train <- df[sub_rows, ]
  df_test  <- df[setdiff(rownames(df), sub_rows), ]

  # ---- 1. Cutoff estimation on training data (FPR and FDR) ----
  cc <- function(algo) CalcCut(df_train, alpha = alpha, niter = niter, nsamp = test_size,
                               algorithm = algo, val_mat = val_adj_logical, M = xi_M,
                               symmetric = sym)
  cuts <- list(
    Pearson    = cc('Pearson'),
    Spearman   = cc('Spearman'),
    Chatterjee = cc('Xi'),
    ChatterjeeKNN = cc('Xi_KNN'),   # CalcCut uses XIM_Mat(df, 10)
    pcNet      = cc('pcNet'),
    Distance   = cc('dcor'),
    MI         = cc('MI')
  )

  # ---- 2. GRN construction on test data ----
  sym_abs <- function(m) {
    m <- abs(unname(as.matrix(m)))
    if (sym) m <- pmax(m, t(m))
    diag(m) <- 0
    m
  }
  nets <- list(
    Pearson    = sym_abs(cor(df_test, method = 'pearson')),
    Spearman   = sym_abs(cor(df_test, method = 'spearman')),
    Chatterjee = sym_abs(xicor_mpar(df_test, nCores = 4)),
    ChatterjeeKNN = { set.seed(0); sym_abs(XIM_Mat(df_test, xi_M, verbose = FALSE)) },
    pcNet      = sym_abs(pcNet(t(df_test), verbose = FALSE)),
    Distance   = sym_abs(computeDcor(df_test, nCores = 4)),
    MI         = sym_abs(computeMI(df_test))
  )

  # ---- 3. Metrics at both cutoffs ----
  metrics_df <- do.call(rbind, lapply(names(nets), function(a)
    evaluate_method(a, nets[[a]], cuts[[a]], eval_mask, labels, rmse = a != 'MI')))
  metrics_df$Iter <- i

  # ---- 4. Score distributions (unique pairs) and cutoffs ----
  scores_df <- do.call(rbind, lapply(names(nets), function(a)
    data.frame(Algorithm = a, Score = nets[[a]][eval_mask], Iter = i)))
  cutoffs_df <- do.call(rbind, lapply(names(cuts), function(a)
    data.frame(Algorithm = a, CutoffType = names(cut_types),
               Cutoff = unname(cuts[[a]][cut_types]), Iter = i)))

  list(metrics = metrics_df, scores = scores_df, cutoffs = cutoffs_df)
}

# ============================================================
# Main loop over DS folders
# ============================================================
tag_of <- function(dir) paste0(ds_basename(dir), '_GTSym', sym)

for (ds_name in names(ds_dirs)) {
  dir <- ds_dirs[[ds_name]]
  message(sprintf("\n====== %s ======", ds_name))

  val_adj_logical <- LoadGT(dir, symmetric = sym)
  eval_mask <- (if (sym) upper.tri(val_adj_logical) else
    row(val_adj_logical) != col(val_adj_logical)) & !is.na(val_adj_logical)
  labels    <- as.integer(val_adj_logical[eval_mask])

  df0        <- t(read.csv(paste0(dir, 'simulated_noNoise_T_0.csv'), header = TRUE, row.names = 1))
  train_size <- round(nrow(df0) * 3 / 5)
  test_size  <- nrow(df0) - train_size
  tps <- sort(as.integer(sub('.*_T_(\\d+)\\.csv$', '\\1',
                             list.files(dir, pattern = '^simulated_noNoise_T_\\d+\\.csv$'))))

  results <- lapply(tps, process_timepoint,
                    dir = dir, val_adj_logical = val_adj_logical,
                    eval_mask = eval_mask, labels = labels,
                    train_size = train_size, test_size = test_size)

  res_tag <- tag_of(dir)
  for (part in c('metrics', 'scores', 'cutoffs')) {
    out <- do.call(rbind, lapply(results, `[[`, part))
    fn  <- c(metrics = 'Cutoff_Selection_', scores = 'ScoreDist_', cutoffs = 'Cutoffs_')[part]
    write.csv(out, paste0('Results/Benchmark/', fn, 'SERGIO_', res_tag, '_ver3.csv'),
              row.names = FALSE)
  }
}

# ============================================================
# Quick look: R methods only (before GRNBoost merge), shown on screen
# ============================================================
quick_metrics <- c('FPR', 'FDR', 'Sensitivity(TPR)', 'Balanced Accuracy', 'AUROC',
                   'FDR_feasible', 'FDR_minAchievable')
quick_df <- do.call(rbind, lapply(names(ds_dirs), function(ds_name) {
  read.csv(paste0('Results/Benchmark/Cutoff_Selection_SERGIO_',
                  tag_of(ds_dirs[[ds_name]]), '_ver3.csv')) %>%
    mutate(DS = ds_name)
})) %>%
  filter(Metric %in% quick_metrics) %>%
  mutate(Algorithm = factor(Algorithm, levels = alg_levels),
         Metric    = factor(Metric,    levels = quick_metrics),
         DS        = factor(DS,        levels = names(ds_dirs)))

for (ct in names(cut_types)) {
  p_quick <- ggplot(quick_df %>% filter(CutoffType == ct),
                    aes(x = Algorithm, y = Value, fill = Algorithm)) +
    geom_boxplot(width = 0.7, outlier.size = 0.5, linewidth = 0.4, colour = 'black') +
    geom_jitter(width = 0.1, height = 0, size = 0.4) +
    { if (ct == 'FPR') geom_hline(data = data.frame(Metric = factor('FPR', levels = quick_metrics)),
                                  aes(yintercept = alpha), linetype = 'dashed', colour = 'red')
      else geom_hline(data = data.frame(Metric = factor('FDR', levels = quick_metrics)),
                      aes(yintercept = alpha), linetype = 'dashed', colour = 'red') } +
    facet_grid(DS ~ Metric, scales = 'free_y') +
    scale_fill_manual(values = pal, drop = TRUE) +
    labs(x = NULL, y = NULL,
         title = sprintf('R methods, %s-controlled cutoff (alpha = %.2f; red = alpha)', ct, alpha)) +
    theme_classic() +
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = 'black'),
          axis.text.y = element_text(colour = 'black'),
          strip.background = element_blank(),
          strip.text = element_text(face = 'bold'))
  print(p_quick)
}

# Median per DS x method x cutoff type, for a quick numeric check
quick_df %>%
  group_by(DS, CutoffType, Algorithm, Metric) %>%
  summarise(Median = round(median(Value, na.rm = TRUE), 3), .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = Metric, values_from = Median) %>%
  arrange(DS, CutoffType, Algorithm) %>%
  print(n = Inf)

# ============================================================
# Merge with GRNBoost (xicor_benchmark_SERGIO_ver3.py) and plot per DS
# ============================================================
bench_theme <- theme_classic() +
  theme(
    legend.position = "none",
    text            = element_text(family = "Arial"),
    axis.text.x     = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   family = 'Arial', colour = 'black'),
    axis.text.y     = element_text(family = 'Arial', colour = 'black'),
    plot.title      = element_text(family = 'Arial', face = 'bold', size = 10)
  )

bench_box <- function(d, ncol) {
  ggplot(d, aes(x = Algorithm, y = Value, fill = Algorithm)) +
    geom_boxplot(width = 0.7, outlier.size = 0.5, linewidth = 0.5, colour = "black") +
    geom_jitter(width = 0.1, height = 0, size = 0.5) +
    facet_wrap(~Metric, ncol = ncol, scales = 'free') +
    labs(x = NULL, y = NULL) +
    bench_theme +
    scale_fill_manual(values = pal)
}

dist_alg_order <- c('Pearson', 'Spearman', 'pcNet', 'Distance', 'MI', 'GRNBoost', 'Chatterjee',
                    'ChatterjeeKNN')

for (ds_name in names(ds_dirs)) {
  dir     <- ds_dirs[[ds_name]]
  res_tag <- tag_of(dir)
  rpath   <- function(fn) paste0('Results/Benchmark/', fn, 'SERGIO_', res_tag, '_ver3.csv')
  gpath   <- function(fn) paste0('Results/Benchmark/grnboost/', fn,
                                 ds_basename(dir), '_Symmetric_ver3.csv')

  metrics_df  <- read.csv(rpath('Cutoff_Selection_'))
  grnboost_df <- read.csv(gpath('Cutoff_Selection_'), check.names = FALSE) %>%
    select(Iter, Algorithm, CutoffType, `Sensitivity(TPR)`, `Specificity(TNR)`, FPR, FDR,
           `Balanced Accuracy`, AUROC, Cutoff, nCalled, FDR_feasible, FDR_minAchievable) %>%
    reshape2::melt(id.vars = c('Iter', 'Algorithm', 'CutoffType'),
                   variable.name = 'Metric', value.name = 'Value') %>%
    mutate(Metric = as.character(Metric)) %>%
    filter(!(Metric %in% c('FDR_feasible', 'FDR_minAchievable') & CutoffType == 'FPR')) %>%
    select(Algorithm, CutoffType, Metric, Value, Iter)

  merged_all <- rbind(metrics_df, grnboost_df) %>% mutate(DS = ds_name)
  write.csv(merged_all,
            paste0('Results/Cutoff_Selection_SERGIO_', res_tag, '_Merged_ver3.csv'),
            row.names = FALSE)

  merged_df <- merged_all %>%
    filter(Metric %in% met_levels) %>%
    mutate(Algorithm = factor(Algorithm, levels = alg_levels),
           Metric    = factor(Metric,    levels = met_levels),
           Value     = as.numeric(Value))

  # ---- Benchmark plots: one per cutoff type ----
  for (ct in names(cut_types)) {
    d  <- merged_df %>% filter(CutoffType == ct, Metric != 'RMSE') %>% droplevels()
    nm <- nlevels(d$Metric)
    p  <- bench_box(d, ncol = nm) +
      ggtitle(sprintf('%s — %s-controlled cutoff (alpha = %.2f)', ds_name, ct, alpha)) +
      facetted_pos_scales(y = list(
        Metric == "AUROC" ~ scale_y_continuous(limits = c(0.6, 1), oob = scales::squish)))
    png(paste0('Results/Benchmark/Figures/Bench_SERGIO_', res_tag, '_', ct, 'cut_ver3.png'),
        height = 3.3, width = 1.9 * nm, units = 'in', res = 300)
    print(p)
    dev.off()
  }

  # ---- RMSE subplot (cutoff-independent; [0,1]-bounded methods) ----
  rmse_algs <- c('Pearson', 'Spearman', 'Distance', 'Chatterjee', 'ChatterjeeKNN')
  RMSE_df <- merged_df %>%
    filter(CutoffType == 'FPR', Algorithm %in% rmse_algs, Metric == 'RMSE') %>%
    mutate(Algorithm = factor(Algorithm, levels = rmse_algs))

  p_rmse <- bench_box(RMSE_df, ncol = 1) + scale_fill_manual(values = pal[rmse_algs])
  png(paste0('Results/Benchmark/Figures/Bench_SERGIO_corrRMSE_', res_tag, '_ver3.png'),
      height = 3, width = 2, units = 'in', res = 300)
  print(p_rmse)
  dev.off()

  # ---- Score distributions with both cutoffs ----
  # Dashed = FPR cutoff, dotted = FDR cutoff (median over replicates); Inf not drawn.
  scores_df  <- read.csv(rpath('ScoreDist_'))
  cutoffs_df <- read.csv(rpath('Cutoffs_'))
  if (file.exists(gpath('ScoreDist_'))) scores_df  <- rbind(scores_df,  read.csv(gpath('ScoreDist_')))
  if (file.exists(gpath('Cutoffs_')))   cutoffs_df <- rbind(cutoffs_df, read.csv(gpath('Cutoffs_')))

  cutoff_lines <- cutoffs_df %>%
    group_by(Algorithm, CutoffType) %>%
    summarise(Cutoff = median(Cutoff), .groups = 'drop') %>%
    filter(is.finite(Cutoff), Algorithm %in% dist_alg_order) %>%
    mutate(Algorithm  = factor(Algorithm, levels = dist_alg_order),
           CutoffType = factor(CutoffType, levels = names(cut_types)))

  p_dist <- ggplot(scores_df %>% filter(Algorithm %in% dist_alg_order) %>%
                     mutate(Algorithm = factor(Algorithm, levels = dist_alg_order)),
                   aes(x = Score, fill = Algorithm)) +
    geom_density(alpha = 0.7, colour = 'black') +
    geom_vline(data = cutoff_lines,
               aes(xintercept = Cutoff, linetype = CutoffType),
               linewidth = 0.6, colour = 'black') +
    scale_linetype_manual(values = c(FPR = 'dashed', FDR = 'dotted'), name = 'Cutoff') +
    facet_wrap(~Algorithm, ncol = 4, scales = 'free') +
    labs(x = 'Score', y = 'Density', title = ds_name) +
    guides(fill = 'none') +
    theme_classic() +
    theme(legend.position  = c(0.92, 0.2),
          text             = element_text(family = "Arial"),
          plot.title       = element_text(hjust = 0.5, face = 'bold'),
          axis.text.x      = element_text(family = 'Arial', colour = 'black'),
          axis.text.y      = element_text(family = 'Arial', colour = 'black'),
          strip.background = element_blank(),
          strip.text       = element_text(family = 'Arial', face = 'bold', colour = 'black')) +
    scale_fill_manual(values = pal[dist_alg_order])

  png(paste0('Results/Benchmark/Figures/ScoreDist_SERGIO_', res_tag, '_ver3.png'),
      height = 4, width = 8, units = 'in', res = 300)
  print(p_dist)
  dev.off()
}

message("Done. All outputs written with '_ver3' suffix.")
