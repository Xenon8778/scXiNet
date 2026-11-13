# ============================================================
# SweepM_ChatterjeeKNN.R
# Parameter sweep for M (number of right nearest neighbours) in the
# revised xi (Lin & Han 2023; XIM_Mat) on the nonlinear simulation and on
# SERGIO DS4/DS5/DS6.
#
# Same data handling as CutoffSelect_Nonlinear_ver3.R and
# CutoffSelect_SERGIO_noAnc_ver3.R:
#   - the same ground truths (3-class for nonlinear, 2-class noAnc for SERGIO)
#   - the same 60/40 train/test split (set.seed(i))
#   - the same CalcCut logic (FPR + FDR cutoffs), with M exposed as a parameter
#   - the same test-network construction and metrics over unique pairs
# With niter = 50 and all replicates, the M = 10 row reproduces the
# benchmark's ChatterjeeKNN numbers (a sanity check).
#
# Speed knobs: niter_sweep, n_rep, do_cutoffs (FALSE = AUROC only, no CalcCut).
# ============================================================
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(pROC); library(tidyr)
})
source('code/cutoffCalcScripts.R')   # XIM_Mat (via XIMCOR_mod.R), fdpCurve

# ---------------- Parameters ----------------
M_grid      <- c(1, 3, 5, 10, 20, 50)
alpha       <- 0.05
niter_sweep <- 20      # CalcCut iterations per M (benchmark uses 50)
n_rep       <- 5       # replicates per dataset (first n_rep; Inf = all)
do_cutoffs  <- FALSE    # FALSE -> AUROC only (much faster)
n_cores     <- 4
out_dir     <- 'Results/Benchmark/SweepM/'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------- Ground-truth loaders (copied from the benchmarks) ----------------
LoadGT_nonlin <- function(dir, nCells) {
  genes <- colnames(t(read.csv(paste0(dir, 'counts_1_', nCells, '.csv'),
                               header = TRUE, row.names = 1)))
  val <- read.csv(paste0(dir, 'gt_GRN.csv'), header = TRUE, row.names = 1)
  el  <- as.matrix(apply(as.matrix(val[, 1:3]), 2, as.character))
  g   <- igraph::graph_from_data_frame(data.frame(from = el[, 1], to = el[, 2]),
                                       directed = TRUE, vertices = genes)
  R   <- is.finite(igraph::distances(g, v = genes, to = genes, mode = 'out')) * 1
  CA  <- crossprod(R) > 0
  direct <- igraph::as_adjacency_matrix(g, sparse = FALSE)[genes, genes] > 0
  direct <- direct | t(direct)
  cls <- matrix(NA, length(genes), length(genes), dimnames = list(genes, genes))
  cls[!CA] <- FALSE; cls[direct] <- TRUE; diag(cls) <- NA
  cls
}

LoadGT_sergio <- function(dir) {
  genes <- colnames(t(read.csv(paste0(dir, 'simulated_noNoise_T_0.csv'),
                               header = TRUE, row.names = 1)))
  val <- read.csv(paste0(dir, 'gt_GRN.csv'), header = FALSE)
  g   <- igraph::graph_from_data_frame(
    data.frame(from = as.character(val[, 1]), to = as.character(val[, 2])),
    directed = TRUE, vertices = genes)
  direct <- igraph::as_adjacency_matrix(g, sparse = FALSE)[genes, genes] > 0
  cls <- direct | t(direct)
  dimnames(cls) <- list(genes, genes); diag(cls) <- NA
  cls
}

sergio_files <- function(dir) {
  tps <- sort(as.integer(sub('.*_T_(\\d+)\\.csv$', '\\1',
                             list.files(dir, pattern = '^simulated_noNoise_T_\\d+\\.csv$'))))
  setNames(paste0(dir, 'simulated_noNoise_T_', tps, '.csv'), tps)
}

nl_dir <- 'Simulation_scRNAseq/'; nl_cells <- 1000
datasets <- list(
  Nonlinear = list(gt = function() LoadGT_nonlin(nl_dir, nl_cells),
                   files = setNames(paste0(nl_dir, 'counts_', 1:10, '_', nl_cells, '.csv'), 1:10)),
  DS4 = list(dir = 'SERGIO/data_sets/De-noised_100G_3T_300cPerT_dynamics_9_DS4/'),
  DS5 = list(dir = 'SERGIO/data_sets/De-noised_100G_4T_300cPerT_dynamics_10_DS5/'),
  DS6 = list(dir = 'SERGIO/data_sets/De-noised_100G_6T_300cPerT_dynamics_7_DS6/')
)
for (nm in c('DS4', 'DS5', 'DS6')) local({
  d <- datasets[[nm]]$dir
  datasets[[nm]]$gt    <<- function() LoadGT_sergio(d)
  datasets[[nm]]$files <<- sergio_files(d)
})

# ---------------- Core ----------------
xim_scores <- function(df, M, mask) {
  set.seed(0)                                      # tie-breaking, as in benchmarks
  r <- abs(unname(as.matrix(XIM_Mat(df, M, nCores = n_cores, verbose = FALSE))))
  r <- pmax(r, t(r)); diag(r) <- 0
  s <- r[mask]; s[is.na(s)] <- 0
  s
}

# CalcCut('Xi_KNN') with M exposed (identical logic otherwise)
calib_M <- function(x, M, mask, lab_lgl, nsamp, niter, ngrid = 2000) {
  replace <- nsamp >= nrow(x)
  S  <- matrix(NA_real_, niter, length(lab_lgl)); cv <- numeric(niter)
  for (b in seq_len(niter)) {
    set.seed(b)
    s <- xim_scores(x[sample.int(nrow(x), nsamp, replace = replace), , drop = FALSE], M, mask)
    S[b, ] <- s
    cv[b]  <- quantile(s[!lab_lgl], 1 - alpha, names = FALSE)
  }
  grid <- sort(unique(quantile(S, seq(0, 1, length.out = ngrid), names = FALSE)))
  cur  <- lapply(seq_len(niter), function(b) fdpCurve(S[b, ], lab_lgl, grid))
  FDRb <- rowMeans(sapply(cur, `[[`, 'fdp')); Rb <- rowMeans(sapply(cur, `[[`, 'R'))
  ok   <- FDRb <= alpha & Rb >= 1
  c(FPR = median(cv), FDR = if (any(ok)) min(grid[ok]) else Inf,
    FDR_feasible = as.numeric(any(ok)), FDR_minAchievable = min(FDRb[Rb >= 1]))
}

eval_at <- function(s, labels, cutoff) {
  pred <- s > cutoff
  TP <- sum(pred & labels == 1); FP <- sum(pred & labels == 0)
  TN <- sum(!pred & labels == 0); FN <- sum(!pred & labels == 1)
  tpr <- TP / (TP + FN); tnr <- TN / (TN + FP)
  c(FPR = 1 - tnr, FDR = if (TP + FP == 0) NA_real_ else FP / (TP + FP),
    TPR = tpr, BA = (tpr + tnr) / 2)
}

run_dataset <- function(ds_name) {
  spec   <- datasets[[ds_name]]
  gt     <- spec$gt()
  mask   <- upper.tri(gt) & !is.na(gt)
  labels <- as.integer(gt[mask])
  files  <- head(spec$files, n_rep)
  message(sprintf('\n=== %s: %d reps, pos=%d null=%d ===', ds_name, length(files),
                  sum(labels), sum(!labels)))

  out <- list()
  for (i_chr in names(files)) {
    i  <- as.integer(i_chr)
    df <- t(read.csv(files[[i_chr]], header = TRUE, row.names = 1))[, colnames(gt)]
    train_size <- round(nrow(df) * 3 / 5); test_size <- nrow(df) - train_size
    set.seed(i)
    sub_rows <- sample(rownames(df), train_size)
    df_train <- df[sub_rows, ]; df_test <- df[setdiff(rownames(df), sub_rows), ]

    for (M in M_grid) {
      t0    <- Sys.time()
      s     <- xim_scores(df_test, M, mask)
      auroc <- as.numeric(auc(roc(labels, s, quiet = TRUE)))
      rmse  <- sqrt(mean((s - labels)^2))
      row   <- data.frame(DS = ds_name, Rep = i, M = M, Metric = c('AUROC', 'RMSE'),
                          CutoffType = 'none', Value = c(auroc, rmse))
      if (do_cutoffs) {
        cuts <- calib_M(df_train, M, mask, labels == 1L, test_size, niter_sweep)
        for (ct in c('FPR', 'FDR')) {
          m <- eval_at(s, labels, cuts[[ct]])
          row <- rbind(row, data.frame(DS = ds_name, Rep = i, M = M, Metric = names(m),
                                       CutoffType = ct, Value = unname(m)))
        }
        row <- rbind(row, data.frame(DS = ds_name, Rep = i, M = M,
                                     Metric = c('Cutoff', 'Cutoff', 'FDR_feasible', 'FDR_minAchievable'),
                                     CutoffType = c('FPR', 'FDR', 'FDR', 'FDR'),
                                     Value = unname(cuts[c('FPR', 'FDR', 'FDR_feasible', 'FDR_minAchievable')])))
      }
      row <- rbind(row, data.frame(DS = ds_name, Rep = i, M = M, Metric = 'Runtime_s',
                                   CutoffType = 'none',
                                   Value = as.numeric(difftime(Sys.time(), t0, units = 'secs'))))
      out[[length(out) + 1]] <- row
      message(sprintf('  rep %d  M=%2d  AUROC=%.3f', i, M, auroc))
    }
  }
  res <- do.call(rbind, out)
  write.csv(res, paste0(out_dir, 'SweepM_', ds_name, '.csv'), row.names = FALSE)  # checkpoint
  res
}

# ---------------- Run ----------------
res <- do.call(rbind, lapply(names(datasets), run_dataset))
res$DS <- factor(res$DS, levels = names(datasets))
write.csv(res, paste0(out_dir, 'SweepM_all.csv'), row.names = FALSE)

# ---------------- Summary + choice of M ----------------
key <- res %>%
  mutate(Key = ifelse(CutoffType == 'none', Metric, paste0(Metric, '@', CutoffType))) %>%
  group_by(DS, M, Key) %>%
  summarise(Median = median(Value, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Key, values_from = Median) %>%
  arrange(DS, M)
print(key, n = Inf, width = Inf)
write.csv(key, paste0(out_dir, 'SweepM_summary.csv'), row.names = FALSE)

# Primary criterion: median BA at the FPR cutoff (AUROC if do_cutoffs = FALSE).
# Ties go to the smaller M.
crit <- if (do_cutoffs) 'BA@FPR' else 'AUROC'
best <- key %>% group_by(DS) %>%
  arrange(desc(.data[[crit]]), M, .by_group = TRUE) %>% slice(1) %>%
  select(DS, M, all_of(crit))
cat('\nBest M per dataset (', crit, '):\n'); print(best)

# One M for everything: the lowest mean rank across datasets
overall <- key %>% group_by(DS) %>% mutate(rank = rank(-.data[[crit]], ties.method = 'min')) %>%
  group_by(M) %>% summarise(mean_rank = mean(rank), mean_crit = mean(.data[[crit]]), .groups = 'drop') %>%
  arrange(mean_rank, M)
cat('\nOverall (mean rank across datasets):\n'); print(overall)

# ---------------- Plot ----------------
plot_metrics <- if (do_cutoffs)
  c('AUROC', 'RMSE', 'BA@FPR', 'TPR@FPR', 'FPR@FPR', 'FDR@FDR', 'TPR@FDR', 'Runtime_s') else
    c('AUROC', 'RMSE')
pd <- res %>%
  mutate(Key = ifelse(CutoffType == 'none', Metric, paste0(Metric, '@', CutoffType))) %>%
  filter(Key %in% plot_metrics) %>%
  mutate(Key = factor(Key, levels = plot_metrics)) %>%
  group_by(DS, M, Key) %>%
  summarise(med = median(Value, na.rm = TRUE), lo = quantile(Value, .25, na.rm = TRUE),
            hi = quantile(Value, .75, na.rm = TRUE), .groups = 'drop')

if (do_cutoffs){
  hl <- data.frame(Key = factor(intersect(c('FPR@FPR', 'FDR@FDR'), plot_metrics), levels = plot_metrics),
                 y = alpha)
}

p <- ggplot(pd, aes(M, med, colour = DS, fill = DS)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line() + geom_point(size = 1.2) +
  geom_vline(xintercept = 5, linetype = 'dotted', colour = 'grey40') +   # current default
  facet_wrap(~Key, scales = 'free_y', ncol = 3) +
  scale_x_continuous(breaks = M_grid) +
  labs(x = 'M (right nearest neighbours)', y = 'Median (IQR band)',
       fill = 'Dataset', color='Dataset',
       title = sprintf('ChatterjeeKNN M sweep',
                       alpha, niter_sweep, n_rep)) +
  theme_classic() + theme(strip.background = element_blank(), strip.text = element_text(face = 'bold'))
if (do_cutoffs){p +
    geom_hline(data = hl, aes(yintercept = y), linetype = 'dashed', colour = 'red')
}
print(p)
ggsave(paste0(out_dir, 'SweepM_ChatterjeeKNN.png'), p,
       width = 5, height = 3, dpi = 300)

message('Done. Outputs in ', out_dir)
