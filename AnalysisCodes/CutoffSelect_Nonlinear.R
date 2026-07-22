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
source('code/scripts.R')

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
computeDcor <- function(mat) {
  p   <- ncol(mat)
  res <- matrix(0.0, p, p, dimnames = list(colnames(mat), colnames(mat)))
  for (j in seq_len(p - 1L)) {
    for (k in (j + 1L):p) {
      res[j, k] <- res[k, j] <- dcor(mat[, j], mat[, k])
    }
  }
  res
}

# ============================================================
# Ground Truth
# ============================================================
dir    <- "Simulation_scRNAseq/"
nCells <- 1000

LoadGT <- function(symmetric = TRUE, edgelist = FALSE) {
  df  <- t(read.csv(paste0(dir, 'counts_1_', nCells, '.csv'),
                    header = TRUE, row.names = 1))
  val <- read.csv(paste0(dir, 'gt_GRN.csv'), header = TRUE, row.names = 1)
  el  <- as.matrix(apply(as.matrix(val[, 1:3]), 2, as.character))

  val_adj <- matrix(0, ncol(df), ncol(df),
                    dimnames = list(colnames(df), colnames(df)))
  val_adj[el[, 1:2]] <- 1
  if (symmetric) val_adj <- pmax(val_adj, t(val_adj))
  diag(val_adj) <- 0

  val_lg <- apply(val_adj, 2, as.logical)
  rownames(val_lg) <- rownames(val_adj)
  if (edgelist) el else val_lg
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

png('figures/scRNAseq_nonlinear_simulation_scatter_ver2.png',
    height = 3, width = 12, units = 'in', res = 300)
plot_grid(p1, p2, p3, p4, ncol = 4)
dev.off()

# ============================================================
# Benchmarking parameters
# ============================================================
niter      <- 50
train_size <- round(nrow(df_plot) * 3 / 5)
test_size  <- nrow(df_plot) - train_size
cMetric    <- 'a95'
sym        <- TRUE

val_adj_logical <- LoadGT(symmetric = sym)
val_adj         <- matrix(as.integer(val_adj_logical),
                          nrow = nrow(val_adj_logical))

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
  # 1. Cutoff estimation on training data
  # ----------------------------------------------------------
  cutoffs_Xi       <- CalcCut(df_train, niter = niter, nsamp = test_size,
                              algorithm = 'Xi',       val_mat = val_adj_logical,
                              symmetric = sym)
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

  # ----------------------------------------------------------
  # 2. GRN construction on test data
  #
  # Scaling notes:
  #   - No scaling applied. All bounded methods (Pearson, Spearman,
  #     Xi, dcor, pcNet) are already in [0,1] after abs(); raw scores are
  #     used directly for both RMSE and AUROC.
  #   - MI scores are unbounded (nats) and excluded from RMSE entirely.
  #     AUROC still uses raw MI scores (rank-invariant, scale-free).
  # ----------------------------------------------------------

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

  # ----------------------------------------------------------
  # 3. Confusion matrices
  # ----------------------------------------------------------
  cm_xi       <- confusionMatrix(table(cor_xi_filtered,       val_adj), positive = "1")
  cm_pcNet    <- confusionMatrix(table(cor_pcnet_filtered,    val_adj), positive = "1")
  cm_spearman <- confusionMatrix(table(cor_Spearman_filtered, val_adj), positive = "1")
  cm_pearson  <- confusionMatrix(table(cor_Pearson_filtered,  val_adj), positive = "1")
  cm_dcor     <- confusionMatrix(table(cor_dcor_filtered,     val_adj), positive = "1")
  cm_MI       <- confusionMatrix(table(cor_MI_filtered,       val_adj), positive = "1")

  # ----------------------------------------------------------
  # 4. Collect performance metrics
  #
  #   RMSE: only computed for methods with scores in [0,1].
  #         MI is excluded — raw MI scores are unbounded (nats)
  #         and not meaningfully comparable to val_adj in {0,1}.
  #   AUROC: raw scores for all methods (rank-invariant, scale-free).
  # ----------------------------------------------------------
  metrics_list <- list()

  add_metrics <- function(alg, cm, scores, cutoff, rmse = TRUE) {
    spec <- cm$byClass['Specificity']
    rows <- list(
      c(Algorithm = alg, Metric = 'Sensitivity(TPR)',  Value = unname(cm$byClass['Sensitivity'])),
      c(Algorithm = alg, Metric = 'Specificity(TNR)',  Value = unname(spec)),
      c(Algorithm = alg, Metric = 'FPR',               Value = unname(1 - spec)),
      c(Algorithm = alg, Metric = 'Balanced Accuracy', Value = unname(cm$byClass['Balanced Accuracy'])),
      c(Algorithm = alg, Metric = 'AUROC',
        Value = as.numeric(auc(roc(as.numeric(val_adj), as.numeric(scores),
                                   quiet = TRUE)))),
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
  add_metrics('MI',         cm_MI,       cor_MI_test,       cutoffs_MI[cMetric],  rmse = FALSE)

  metrics_df <- as.data.frame(do.call(rbind, metrics_list), stringsAsFactors = FALSE)
  metrics_df$Value <- as.numeric(metrics_df$Value)
  metrics_df$Iter  <- i

  # ----------------------------------------------------------
  # 5. Score distributions + cutoffs (upper triangle only, no diagonal)
  # ----------------------------------------------------------
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
# Run over all 10 datasets
# ============================================================
results    <- lapply(1:10, process_dataset)
metrics_df <- do.call(rbind, lapply(results, `[[`, 'metrics'))
scores_df  <- do.call(rbind, lapply(results, `[[`, 'scores'))

# ---- Export raw results ----
res_tag <- paste0(tail(unlist(strsplit(dir, '/')), 1),
                  '_', cMetric, '_GTSym', sym, '_', nCells)
write.csv(metrics_df,
          paste0('Results/Benchmark/Cutoff_Selection_', res_tag, '_ver2.csv'),
          row.names = FALSE)

# ============================================================
# Merge with GRNBoost (Python) and plot performance
# ============================================================
alg_levels <- c('Pearson', 'Spearman', 'pcNet', 'GRNBoost',
                'Distance', 'MI', 'Chatterjee')
met_levels <- c('FPR', 'Sensitivity(TPR)', 'Balanced Accuracy', 'RMSE', 'AUROC')

metrics_df <- read.csv(paste0('Results/Benchmark/Cutoff_Selection_', res_tag, '_ver2.csv'))
grnboost_df <- read.csv(paste0('Results/Benchmark/grnboost/',
                               'Cutoff_Selection_Simulation_scRNAseq_',
                               nCells, '_Symmetric.csv'))
colnames(grnboost_df) <- c('Iter', 'Algorithm', 'Sensitivity(TPR)', 'Specificity(TNR)',
                           'Balanced Accuracy', 'AUROC', 'Cutoff', 'Symmetric')
grnboost_df <- grnboost_df %>%
  mutate(FPR = 1 - `Specificity(TNR)`) %>%
  reshape2::melt(id.vars = c('Iter', 'Algorithm')) %>%
  setNames(c('Iter', 'Algorithm', 'Metric', 'Value')) %>%
  mutate(Algorithm = if_else(Algorithm == "pcNet", "pcNetpy", Algorithm)) %>%
  select(Algorithm, Metric, Value, Iter)

merged_df <- rbind(metrics_df, grnboost_df) %>%
  filter(Algorithm != 'pcNetpy', Metric %in% met_levels) %>%
  mutate(
    Algorithm = factor(Algorithm, levels = alg_levels),
    Metric    = factor(Metric,    levels = met_levels),
    Value     = as.numeric(Value)
  )

write.csv(merged_df,
          'Results/Cutoff_Selection_Simulation_scRNAseq_Merged_ver2.csv',
          row.names = FALSE)

# Named colour palette (one colour per algorithm)
pal <- setNames(
  hcl.colors(length(alg_levels), palette = 'RdPu', rev = TRUE),
  alg_levels
)

# ---- Main benchmark plot (RMSE excluded — separate subplot) ----
merged_df <- read.csv('Results/Cutoff_Selection_Simulation_scRNAseq_Merged_ver2.csv')
merged_df <- merged_df %>%
  mutate(
    Algorithm = factor(Algorithm, levels = alg_levels),
    Metric    = factor(Metric,    levels = met_levels),
    Value     = as.numeric(Value)
  )

p_bench <- ggplot(merged_df %>% filter(Metric != 'RMSE'),
                  aes(x = Algorithm, y = Value, fill = Algorithm)) +
  geom_boxplot(width = 0.7, outlier.size = 0.5, linewidth = 0.5, color = "black") +
  geom_jitter(width = 0.1, size = 0.5) +
  facet_wrap(~Metric, ncol = 5, scales = 'free') +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    legend.position = "none",
    text            = element_text(family = "Arial"),
    axis.text.x     = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   family = 'Arial', colour = 'black'),
    axis.text.y     = element_text(family = 'Arial', colour = 'black')
  ) +
  scale_fill_manual(values = pal)

png(paste0('Results/Benchmark/Figures/Bench_Simulation_scRNAseq_nonlinear_',
           nCells, '_ver2.png'),
    height = 3, width = 8, units = 'in', res = 300)
print(p_bench)
dev.off()

# ---- RMSE subplot: only methods with scores in [0,1]; MI excluded ----
rmse_algs <- c('Pearson', 'Spearman', 'Distance', 'Chatterjee')
RMSE_df <- merged_df %>%
  filter(Algorithm %in% rmse_algs, Metric == 'RMSE') %>%
  mutate(Algorithm = factor(Algorithm, levels = rmse_algs),
         Value     = as.numeric(Value))

p_rmse <- ggplot(RMSE_df, aes(x = Algorithm, y = Value, fill = Algorithm)) +
  geom_boxplot(width = 0.7, outlier.size = 0.5, linewidth = 0.5, color = "black") +
  geom_jitter(width = 0.1, size = 0.5) +
  facet_wrap(~Metric, scales = 'free') +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    legend.position = "none",
    text            = element_text(family = "Arial"),
    axis.text.x     = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   family = 'Arial', colour = 'black'),
    axis.text.y     = element_text(family = 'Arial', colour = 'black')
  ) +
  scale_fill_manual(values = pal[rmse_algs])

png(paste0('Results/Benchmark/Figures/Bench_Simulation_scRNAseq_nonlinear_corrRMSE_',
           nCells, '_ver2.png'),
    height = 3, width = 2, units = 'in', res = 300)
print(p_rmse)
dev.off()

# ============================================================
# Score distribution plot with cutoff lines
# One density panel per algorithm; vertical dashed line = cutoff.
# Cutoff is the median across all 10 dataset iterations.
# ============================================================
alg_order <- c('Pearson', 'Spearman', 'pcNet', 'Distance', 'MI', 'GRNBoost', 'Chatterjee')

# Merge GRNBoost score distributions from Python output if available
grn_scores_path <- paste0('Results/Benchmark/grnboost/ScoreDist_Simulation_scRNAseq_',
                          nCells, '_Symmetric.csv')
if (file.exists(grn_scores_path)) {
  grn_scores_df <- read.csv(grn_scores_path)
  scores_df     <- rbind(scores_df, grn_scores_df)
}

cutoff_lines <- scores_df %>%
  group_by(Algorithm) %>%
  summarise(Cutoff = median(Cutoff), .groups = 'drop') %>%
  mutate(Algorithm = factor(Algorithm, levels = alg_order))

scores_plot_df <- scores_df %>%
  filter(Algorithm %in% alg_order) %>%
  mutate(Algorithm = factor(Algorithm, levels = alg_order))

p_dist <- ggplot(scores_plot_df, aes(x = Score, fill = Algorithm)) +
  geom_density(alpha = 0.7, colour = 'black') +
  geom_vline(data = cutoff_lines,
             aes(xintercept = Cutoff),
             linetype = 'dashed', linewidth = 0.6, colour = 'black') +
  facet_wrap(~Algorithm, ncol = 4, scales = 'free') +
  labs(x = 'Score', y = 'Density') +
  theme_classic() +
  theme(
    legend.position = "none",
    text            = element_text(family = "Arial"),
    axis.text.x     = element_text(family = 'Arial', colour = 'black'),
    axis.text.y     = element_text(family = 'Arial', colour = 'black'),
    strip.background = element_blank(),
    strip.text       = element_text(family = 'Arial', face = 'bold', colour = 'black')
  ) +
  scale_fill_manual(values = pal[alg_order])

png(paste0('Results/Benchmark/Figures/ScoreDist_Simulation_scRNAseq_nonlinear_',
           nCells, '_ver2.png'),
    height = 3, width = 8, units = 'in', res = 300)
print(p_dist)
dev.off()

message("Done. All outputs written with '_ver2' suffix.")
