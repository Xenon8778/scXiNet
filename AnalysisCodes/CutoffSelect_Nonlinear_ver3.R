library(dplyr)
library(ggplot2)
library(XICOR)
library(ggpubr)
library(reshape2)
library(pROC)
library(caret)
library(energy)
library(minet)        # MI via build.mim
library(ggsignif)
library(scTenifoldNet)
library(cowplot)
source('code/XICOR_mod.R')
source('code/cutoffCalcScripts.R')   # requires updated CalcCut (FPR + FDR cutoffs)

# ============================================================
# Helpers
# ============================================================

# Mutual Information matrix (raw MIM — no ARACNE pruning)
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

# ============================================================
# Ground Truth
# ============================================================
dir    <- "Simulation_scRNAseq/"
nCells <- 1000

LoadGT <- function(symmetric = TRUE, edgelist = FALSE) {
  # Three-class GT: TRUE = direct edge, FALSE = no trek (true null),
  # NA = indirect / co-regulated (excluded from calibration and evaluation)
  genes <- colnames(t(read.csv(paste0(dir, 'counts_1_', nCells, '.csv'),
                               header = TRUE, row.names = 1)))
  val <- read.csv(paste0(dir, 'gt_GRN.csv'), header = TRUE, row.names = 1)
  el  <- as.matrix(apply(as.matrix(val[, 1:3]), 2, as.character))
  if (edgelist) return(el)

  g <- igraph::graph_from_data_frame(data.frame(from = el[, 1], to = el[, 2]),
                                     directed = TRUE, vertices = genes)
  if (!igraph::is_dag(g)) warning("GT contains cycles; trek labelling uses reachability as is.")
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

val_adj_logical     <- LoadGT(symmetric = TRUE)
validation_edgelist <- LoadGT(symmetric = TRUE, edgelist = TRUE)

# ---- Scatter plots of representative interactions ----
df_plot <- t(read.csv(paste0(dir, 'counts_1_', nCells, '.csv'),
                      header = TRUE, row.names = 1))
make_scatter <- function(genes, title) {
  d <- as.data.frame(df_plot[, genes])
  ggplot(d, aes_string(x = genes[1], y = genes[2])) +
    geom_point(size = 0.6) +
    geom_smooth(method = 'gam', se = FALSE, colour = 'magenta3') +
    theme_classic() + ggtitle(title)
}
p1 <- make_scatter(validation_edgelist[1, 1:2], 'Linear')
p2 <- make_scatter(validation_edgelist[6, 1:2], 'Parabolic')
p3 <- make_scatter(validation_edgelist[7, 1:2], 'Exponential')
p4 <- make_scatter(validation_edgelist[8, 1:2], 'Sinusoidal')

png('figures/scRNAseq_nonlinear_simulation_scatter_ver3.png',
    height = 3, width = 12, units = 'in', res = 300)
plot_grid(p1, p2, p3, p4, ncol = 4)
dev.off()

# ============================================================
# Benchmarking parameters
# ============================================================
niter      <- 50
train_size <- round(nrow(df_plot) * 3 / 5)
test_size  <- nrow(df_plot) - train_size
alpha      <- 0.05
# Both cutoffs from CalcCut() are evaluated on every test network.
# CalcCut names outputs from alpha (0.05 -> a95, FDR95, ...), so keys are derived.
cut_sfx      <- as.character(round((1 - alpha) * 100))
cut_types    <- c(FPR = paste0('a', cut_sfx), FDR = paste0('FDR', cut_sfx))
feasible_key <- paste0('FDR', cut_sfx, '_feasible')
minach_key   <- paste0('FDR', cut_sfx, '_minAchievable')
sym        <- TRUE
xi_M       <- 5    # revised xi (Lin & Han 2023): right nearest neighbours; must match CalcCut's Xi_KNN
n_cores    <- 4

val_adj_logical <- LoadGT(symmetric = sym)

# ============================================================
# Evaluation helper
#
# All metrics are computed over unique gene pairs (j < k) when the
# ground truth is symmetric, matching Methods 3.5.1 / 3.9. Using the
# full p x p matrix double-counts every pair and adds p diagonal
# entries as guaranteed true negatives, which deflates FPR and RMSE.
#
#   FPR = FP / (FP + TN)   targeted by the 'a<..>'   cutoff (CutoffType = FPR)
#   FDR = FP / (FP + TP)   targeted by the 'FDR<..>' cutoff (CutoffType = FDR)
#   Both metrics are reported under both cutoffs.
#   Indirect / co-regulated pairs (NA in GT) are excluded from all metrics.
#   FDR := NA when no pair is called (infeasible cutoff; see FDR_feasible)
# ============================================================
eval_mask <- (if (sym) upper.tri(val_adj_logical) else
  row(val_adj_logical) != col(val_adj_logical)) & !is.na(val_adj_logical)
labels    <- as.integer(val_adj_logical[eval_mask])

evaluate_method <- function(alg, scores, cuts, rmse = TRUE) {
  # cuts: full CalcCut() output, c(a<..>, FDR<..>, FDR<..>_feasible, FDR<..>_minAchievable)
  s     <- as.numeric(scores[eval_mask])
  auroc <- as.numeric(auc(roc(labels, s, quiet = TRUE)))
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
# Per-dataset processing
# ============================================================
process_dataset <- function(i) {
  message(sprintf("\n=== Dataset %d ===", i))

  df <- t(read.csv(paste0(dir, 'counts_', i, '_', nCells, '.csv'),
                   header = TRUE, row.names = 1))

  # Train / test split (seed per dataset for independent splits)
  set.seed(i)
  sub_rows <- sample(rownames(df), train_size)
  df_train <- df[sub_rows, ]
  df_test  <- df[setdiff(rownames(df), sub_rows), ]

  # ----------------------------------------------------------
  # 1. Cutoff estimation on training data (FPR and FDR cutoffs)
  # ----------------------------------------------------------
  cc <- function(algo) CalcCut(df_train, alpha = alpha, niter = niter, nsamp = test_size,
                               algorithm = algo, val_mat = val_adj_logical, nCores = n_cores,
                               symmetric = sym, M = xi_M)
  cuts <- list(
    Pearson    = cc('Pearson'),
    Spearman   = cc('Spearman'),
    Chatterjee = cc('Xi'),
    ChatterjeeKNN = cc('Xi_KNN'),
    pcNet      = cc('pcNet'),
    Distance   = cc('dcor'),
    MI         = cc('MI')
  )

  # ----------------------------------------------------------
  # 2. GRN construction on test data
  #   No scaling: bounded methods are in [0,1] after abs().
  #   MI is unbounded -> excluded from RMSE; AUROC is rank-based.
  # ----------------------------------------------------------
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

  # ----------------------------------------------------------
  # 3. Metrics at both cutoffs
  # ----------------------------------------------------------
  metrics_df <- do.call(rbind, lapply(names(nets), function(a)
    evaluate_method(a, nets[[a]], cuts[[a]], rmse = a != 'MI')))
  metrics_df$Iter <- i

  # ----------------------------------------------------------
  # 4. Score distributions (unique pairs) and cutoffs, stored separately
  # ----------------------------------------------------------
  scores_df <- do.call(rbind, lapply(names(nets), function(a)
    data.frame(Algorithm = a, Score = nets[[a]][eval_mask], Iter = i)))
  cutoffs_df <- do.call(rbind, lapply(names(cuts), function(a)
    data.frame(Algorithm = a, CutoffType = names(cut_types),
               Cutoff = unname(cuts[[a]][cut_types]), Iter = i)))

  list(metrics = metrics_df, scores = scores_df, cutoffs = cutoffs_df)
}

# ============================================================
# Run over all 10 datasets
# ============================================================
results    <- lapply(1:10, process_dataset)
metrics_df <- do.call(rbind, lapply(results, `[[`, 'metrics'))
scores_df  <- do.call(rbind, lapply(results, `[[`, 'scores'))
cutoffs_df <- do.call(rbind, lapply(results, `[[`, 'cutoffs'))

# ============================================================
# Quick look: R methods only (before GRNBoost merge), shown on screen
# ============================================================
quick_levels  <- c('Pearson', 'Spearman', 'pcNet', 'Distance', 'MI', 'Chatterjee', 'ChatterjeeKNN')
quick_metrics <- c('FPR', 'FDR', 'Sensitivity(TPR)', 'Balanced Accuracy', 'AUROC',
                   'FDR_feasible', 'FDR_minAchievable')
quick_pal <- setNames(hcl.colors(length(quick_levels) + 2, palette = 'RdPu',
                                 rev = TRUE)[seq_along(quick_levels)], quick_levels)
quick_df <- metrics_df %>%
  filter(Metric %in% quick_metrics) %>%
  mutate(Algorithm = factor(Algorithm, levels = quick_levels),
         Metric    = factor(Metric,    levels = quick_metrics),
         Value     = as.numeric(Value))

for (ct in names(cut_types)) {
  target <- data.frame(Metric = factor(ct, levels = quick_metrics), y = alpha)  # FPR or FDR line
  p_quick <- ggplot(quick_df %>% filter(CutoffType == ct),
                    aes(x = Algorithm, y = Value, fill = Algorithm)) +
    geom_boxplot(width = 0.7, outlier.size = 0.5, linewidth = 0.4, colour = 'black') +
    geom_jitter(width = 0.1, height = 0, size = 0.4) +
    geom_hline(data = target, aes(yintercept = y), linetype = 'dashed', colour = 'red') +
    facet_wrap(~Metric, nrow = 1, scales = 'free_y') +
    scale_fill_manual(values = quick_pal) +
    labs(x = NULL, y = NULL,
         title = sprintf('Nonlinear sim, R methods, %s-controlled cutoff (alpha = %.2f; red = alpha)',
                         ct, alpha)) +
    theme_classic() +
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = 'black'),
          axis.text.y = element_text(colour = 'black'),
          strip.background = element_blank(),
          strip.text = element_text(face = 'bold'))
  print(p_quick)
}

# Medians per cutoff type x method, for a quick numeric check
quick_df %>%
  group_by(CutoffType, Algorithm, Metric) %>%
  summarise(Median = round(median(Value, na.rm = TRUE), 3), .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = Metric, values_from = Median) %>%
  arrange(CutoffType, Algorithm) %>%
  print(n = Inf)

# ---- Export raw results ----
res_tag <- paste0(tail(unlist(strsplit(dir, '/')), 1), '_GTSym', sym, '_', nCells)
write.csv(metrics_df, paste0('Results/Benchmark/Cutoff_Selection_', res_tag, '_ver3.csv'),
          row.names = FALSE)
write.csv(cutoffs_df, paste0('Results/Benchmark/Cutoffs_', res_tag, '_ver3.csv'),
          row.names = FALSE)

# ============================================================
# Merge with GRNBoost (CutoffSelect_GRNBoost_ver3.py)
# Same pair set, same two cutoffs, FDR computed natively.
# ============================================================
metrics_df <- read.csv(paste0('Results/Benchmark/Cutoff_Selection_', res_tag, '_ver3.csv'))
alg_levels <- c('Pearson', 'Spearman', 'pcNet', 'GRNBoost',
                'Distance', 'MI', 'Chatterjee', 'ChatterjeeKNN')
met_levels <- c('FPR', 'FDR', 'Sensitivity(TPR)', 'Balanced Accuracy',
                'RMSE', 'AUROC', 'FDR_feasible', 'FDR_minAchievable')

grn_tag <- paste0('Simulation_scRNAseq_', nCells, '_Symmetric_ver3.csv')
grnboost_df <- read.csv(paste0('Results/Benchmark/grnboost/Cutoff_Selection_', grn_tag),
                        check.names = FALSE) %>%
  select(Iter, Algorithm, CutoffType, `Sensitivity(TPR)`, `Specificity(TNR)`, FPR, FDR,
         `Balanced Accuracy`, AUROC, Cutoff, nCalled, FDR_feasible, FDR_minAchievable) %>%
  reshape2::melt(id.vars = c('Iter', 'Algorithm', 'CutoffType'),
                 variable.name = 'Metric', value.name = 'Value') %>%
  mutate(Metric = as.character(Metric)) %>%
  filter(!(Metric %in% c('FDR_feasible', 'FDR_minAchievable') & CutoffType == 'FPR')) %>%
  select(Algorithm, CutoffType, Metric, Value, Iter)

merged_all <- rbind(metrics_df, grnboost_df)
write.csv(merged_all,
          'Results/Cutoff_Selection_Simulation_scRNAseq_Merged_ver3.csv',
          row.names = FALSE)

merged_df <- merged_all %>%
  filter(Metric %in% met_levels) %>%
  mutate(
    Algorithm = factor(Algorithm, levels = alg_levels),
    Metric    = factor(Metric,    levels = met_levels),
    Value     = as.numeric(Value)
  )

# Named colour palette (one colour per algorithm)
pal <- setNames(
  hcl.colors(length(alg_levels)+2, palette = 'RdPu', rev = TRUE)[1:(length(alg_levels))],
  alg_levels
)

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
    geom_boxplot(width = 0.7, outlier.size = 0.5, linewidth = 0.5, color = "black") +
    geom_jitter(width = 0.1, height = 0, size = 0.5) +
    facet_wrap(~Metric, ncol = ncol, scales = 'free') +
    labs(x = NULL, y = NULL) +
    bench_theme +
    scale_fill_manual(values = pal)
}

# ---- Benchmark plots: one per cutoff type ----
#   FPR cutoff: FPR, FDR, TPR, Balanced Accuracy, AUROC
#   FDR cutoff: same + FDR_feasible, FDR_minAchievable
#   (AUROC/RMSE are cutoff-independent; shown in both for reference)
for (ct in names(cut_types)) {
  d <- merged_df %>%
    filter(CutoffType == ct, Metric != 'RMSE',
           !(ct == 'FPR' & Metric %in% c('FDR_feasible', 'FDR_minAchievable'))) %>%
    droplevels()
  nm <- length(unique(d$Metric))
  p  <- bench_box(d, ncol = nm) +
    ggtitle(sprintf('%s-controlled cutoff (alpha = %.2f)', ct, alpha))
  png(paste0('Results/Benchmark/Figures/Bench_Simulation_scRNAseq_nonlinear_',
             nCells, '_', ct, 'cut_ver3.png'),
      height = 3.3, width = 1.9 * nm, units = 'in', res = 300)
  print(p)
  dev.off()
}

# ---- RMSE subplot (cutoff-independent): [0,1]-bounded methods only ----
rmse_algs <- c('Pearson', 'Spearman', 'Distance', 'Chatterjee', 'ChatterjeeKNN')
RMSE_df <- merged_df %>%
  filter(CutoffType == 'FPR', Algorithm %in% rmse_algs, Metric == 'RMSE') %>%
  mutate(Algorithm = factor(Algorithm, levels = rmse_algs))

p_rmse <- bench_box(RMSE_df, ncol = 1) + scale_fill_manual(values = pal[rmse_algs])
png(paste0('Results/Benchmark/Figures/Bench_Simulation_scRNAseq_nonlinear_corrRMSE_',
           nCells, '_ver3.png'),
    height = 3, width = 2, units = 'in', res = 300)
print(p_rmse)
dev.off()

# ============================================================
# Score distributions with both cutoffs
# Dashed = FPR cutoff, dotted = FDR cutoff (median across datasets).
# An FDR cutoff of Inf (infeasible) is not drawn.
# ============================================================
alg_order <- c('Pearson', 'Spearman', 'pcNet', 'Distance', 'MI', 'GRNBoost', 'Chatterjee',
               'ChatterjeeKNN')

grn_scores_path  <- paste0('Results/Benchmark/grnboost/ScoreDist_',  grn_tag)
grn_cutoffs_path <- paste0('Results/Benchmark/grnboost/Cutoffs_',    grn_tag)
if (file.exists(grn_scores_path))  scores_df  <- rbind(scores_df,  read.csv(grn_scores_path))
if (file.exists(grn_cutoffs_path)) cutoffs_df <- rbind(cutoffs_df, read.csv(grn_cutoffs_path))

cutoff_lines <- cutoffs_df %>%
  group_by(Algorithm, CutoffType) %>%
  summarise(Cutoff = median(Cutoff), .groups = 'drop') %>%
  filter(is.finite(Cutoff), Algorithm %in% alg_order) %>%
  mutate(Algorithm  = factor(Algorithm, levels = alg_order),
         CutoffType = factor(CutoffType, levels = names(cut_types)))

scores_plot_df <- scores_df %>%
  filter(Algorithm %in% alg_order) %>%
  mutate(Algorithm = factor(Algorithm, levels = alg_order))

p_dist <- ggplot(scores_plot_df, aes(x = Score, fill = Algorithm)) +
  geom_density(alpha = 0.7, colour = 'black') +
  geom_vline(data = cutoff_lines,
             aes(xintercept = Cutoff, linetype = CutoffType),
             linewidth = 0.6, colour = 'black') +
  scale_linetype_manual(values = c(FPR = 'dashed', FDR = 'dotted'), name = 'Cutoff') +
  facet_wrap(~Algorithm, ncol = 4, scales = 'free') +
  labs(x = 'Score', y = 'Density') +
  guides(fill = 'none') +
  theme_classic() +
  theme(
    legend.position  = c(0.92, 0.2),
    text             = element_text(family = "Arial"),
    axis.text.x      = element_text(family = 'Arial', colour = 'black'),
    axis.text.y      = element_text(family = 'Arial', colour = 'black'),
    strip.background = element_blank(),
    strip.text       = element_text(family = 'Arial', face = 'bold', colour = 'black')
  ) +
  scale_fill_manual(values = pal[alg_order])

png(paste0('Results/Benchmark/Figures/ScoreDist_Simulation_scRNAseq_nonlinear_',
           nCells, '_ver3.png'),
    height = 4, width = 9, units = 'in', res = 300)
print(p_dist)
dev.off()

message("Done. All outputs written with '_ver3' suffix.")
