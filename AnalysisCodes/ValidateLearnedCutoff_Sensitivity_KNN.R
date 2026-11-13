# ============================================================
# Validation of the Algorithm-2 cutoff on the sensitivity simulations (known truth).
# Per setting:
#   reference : rep 2, n' = 7500 cells; null = nonmonotonic == FALSE pairs
#   learn     : LearnCutoff() at each target n (subsample n of n' cells, B = 50, median of quantiles)
#   held-out  : reps 1, 3, 4, 5 at n = 1000, 2000, 3750, cutoff applied unchanged
#   evaluate  : realised FPR / FDR / TPR; compared with the baseline formula and cross-setting transfer
#   baseline  : FPR-criterion fit C(n) = c1 * n^-c2 from CutoffSelect_Directed_multiple_KNN.R
#               (same criterion as Algorithm 2 here, so the comparison is like-for-like)
# ============================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(XICOR)
source('code/XICOR_mod.R')
source('code/XIMCOR_mod.R')
source('code/learnCutoffCurve.R')

xi_alg      <- 'Xi_KNN'                    # 'Xi_KNN' | 'Xi'
xi_M        <- 5
stopifnot(xi_alg %in% c('Xi', 'Xi_KNN'))
alpha       <- 0.05
B_learn     <- 50
n_cores     <- 8
tag         <- paste0(xi_alg, '_ver3')
ref_rep     <- 2;  ref_n <- 7500
target_n    <- c(1000, 2000, 3750)
target_reps <- c(1,3,4,5)
matrix_n    <- 2000

sim_dir <- 'Simulation_Sensitivity'
out_dir <- file.path('Results/Benchmark/Sensitivity_Alg2', tag)
dir.create(file.path(out_dir, 'Figures'), showWarnings = FALSE, recursive = TRUE)
base_fit_file <- paste0('Results/Benchmark/Cutoff_function_fits_FPR_FDR_', tag, '.csv')
base_f <- NULL
if (file.exists(base_fit_file)) {
  ref <- read.csv(base_fit_file)[, c('Criterion', 'c1', 'c2')]
  ref$Criterion <- recode(as.character(ref$Criterion), a95 = 'FPR')
  base_f <- ref[ref$Criterion == 'FPR', ]   # FPR fit: Algorithm 2 here also uses the FPR criterion
  stopifnot(nrow(base_f) == 1)
} else message('No baseline formula file; the baseline comparison is skipped.')

#---------------------
# Helpers
#---------------------
LoadGT <- function(d) {
  val <- read.csv(file.path(d, 'gt_GRN.csv'), header = TRUE, row.names = 1)
  el  <- as.matrix(apply(as.matrix(val[, 1:3]), 2, as.character))
  colnames(el) <- c('source', 'target', 'nonmonotonic')
  el
}
CalcDelta <- function(df, edgelist) {
  set.seed(0)
  delta <- vector('numeric', length = nrow(edgelist))
  for (m in 1:nrow(edgelist)) {
    pair   <- as.matrix(df[, edgelist[m, 1:2]])
    cormat <- as.matrix(if (xi_alg == 'Xi_KNN') XIM_Mat(pair, xi_M, verbose = FALSE) else
      xicor_mpar(pair, nCores = 1))
    delta[m] <- abs(cormat[2, 1] - cormat[1, 2])
  }
  delta[!is.finite(delta)] <- 0
  delta
}
perf <- function(delta, nonmono, cutoff) {
  rej <- delta > cutoff; nul <- nonmono == 'FALSE'
  TP <- sum(rej & !nul); FP <- sum(rej & nul); FN <- sum(!rej & !nul); TN <- sum(!rej & nul)
  c(TPR = TP / (TP + FN), FPR = FP / (FP + TN), FDR = if (TP + FP > 0) FP / (TP + FP) else NA_real_)
}

message('Statistic: ', xi_alg, ' | outputs -> ', out_dir)

cl <- parallel::makeCluster(n_cores)
parallel::clusterEvalQ(cl, {
  library(XICOR)
  source('code/XICOR_mod.R'); source('code/XIMCOR_mod.R')
})

#---------------------
# Learn (Algorithm 2) on the reference at each target n; held-out deltas. Restartable.
#---------------------
learn_setting <- function(setting) {
  f <- file.path(out_dir, paste0('learned_', setting, '.rds'))
  if (file.exists(f)) return(readRDS(f))
  message('== learn: ', setting)
  d  <- file.path(sim_dir, setting)
  el <- LoadGT(d)
  X  <- t(readRDS(file.path(d, sprintf('counts_%d_%d.rds', ref_rep, ref_n))))      # cells x genes
  pr <- data.frame(g1 = el[, 'source'], g2 = el[, 'target'])
  out <- lapply(target_n, function(n)
    LearnCutoff(X, pr, null = el[, 'nonmonotonic'] == 'FALSE', n = n, B = B_learn,
                alpha = alpha, xi_alg = xi_alg, M = xi_M, cl = cl))
  names(out) <- target_n
  saveRDS(out, f); out
}
heldout_setting <- function(setting) {
  f <- file.path(out_dir, paste0('targets_', setting, '.rds'))
  if (file.exists(f)) return(readRDS(f))
  message('== held-out deltas: ', setting)
  d  <- file.path(sim_dir, setting)
  el <- LoadGT(d)
  out <- list()
  for (n in target_n) for (r in target_reps) {
    X <- t(readRDS(file.path(d, sprintf('counts_%d_%d.rds', r, n))))
    out[[length(out) + 1]] <- list(n = n, rep = r, delta = CalcDelta(X, el), nonmono = el[, 'nonmonotonic'])
  }
  saveRDS(out, f); out
}

setting_names <- list.dirs(sim_dir, recursive = FALSE, full.names = FALSE)
setting_names <- c('baseline', setdiff(setting_names, 'baseline'))
parallel::clusterExport(cl, c('CalcDelta', 'xi_alg', 'xi_M'))
tryCatch({
  learned <- setNames(lapply(setting_names, learn_setting),   setting_names)   # learned[[setting]][[as.character(n)]]$cutoff
  targets <- setNames(lapply(setting_names, heldout_setting), setting_names)
}, finally = parallel::stopCluster(cl))
cut_of  <- function(s, n) learned[[s]][[as.character(n)]]$cutoff

cuts_tab <- expand.grid(Setting = setting_names, n = target_n, stringsAsFactors = FALSE)
cuts_tab$Cutoff <- mapply(cut_of, cuts_tab$Setting, cuts_tab$n)
write.csv(cuts_tab, file.path(out_dir, paste0('Alg2_cutoffs_', tag, '.csv')), row.names = FALSE)
print(cuts_tab)

#---------------------
# (1) Held-out realised error: Algorithm 2 vs baseline formula
#---------------------
rows <- list()
for (s in setting_names) for (t in targets[[s]]) {
  cuts <- c(`Algorithm 2` = cut_of(s, t$n))
  if (!is.null(base_f)) cuts <- c(cuts, `Baseline formula` = base_f$c1 * t$n^(-base_f$c2))
  for (src in names(cuts)) {
    m <- perf(t$delta, t$nonmono, cuts[[src]])
    rows[[length(rows) + 1]] <- data.frame(Setting = s, nTarget = t$n, rep = t$rep, Source = src,
                                           Cutoff = unname(cuts[[src]]),
                                           TPR = m[['TPR']], FPR = m[['FPR']], FDR = m[['FDR']])
  }
}
metrics <- bind_rows(rows)
metrics$Setting <- factor(metrics$Setting, levels = setting_names)
metrics$Source  <- factor(metrics$Source, levels = intersect(c('Baseline formula', 'Algorithm 2'),
                                                             unique(metrics$Source)))
write.csv(metrics, file.path(out_dir, paste0('Alg2_heldout_metrics_', tag, '.csv')), row.names = FALSE)

summ <- metrics %>% filter(nTarget == matrix_n) %>%
  group_by(Setting, Source) %>%
  summarise(across(c(Cutoff, TPR, FPR, FDR), ~ mean(.x, na.rm = TRUE)), .groups = 'drop')
write.csv(summ, file.path(out_dir, paste0('Alg2_summary_n', matrix_n, '_', tag, '.csv')), row.names = FALSE)
print(as.data.frame(summ))

metric_plot <- metrics %>% filter(Setting != 'dropout_after')
mk_box <- function(metric, ttl) {
  ggplot(metric_plot, aes(Setting, .data[[metric]], fill = Source)) +
    geom_boxplot(outliers = FALSE, linewidth = 0.3, position = position_dodge(0.8), width = 0.7) +
    geom_jitter(size = 1, alpha = 0.6, position = position_dodge(0.8)) +
    geom_hline(yintercept = alpha, linetype = 2, colour = 'grey30') +
    scale_fill_manual(values = c('Baseline formula' = '#F4A6A6', 'Algorithm 2' = '#9EC3E6'), name = NULL) +
    coord_flip() + facet_wrap(~ nTarget, nrow = 1) +
    labs(x = NULL, y = ttl, subtitle = 'Columns: held-out dataset size n; dashed = 0.05') +
    theme_classic() +
    theme(axis.text = element_text(colour = 'black'),
          strip.background = element_blank(), legend.position = 'bottom',
          panel.grid.major.y = element_line(colour = 'grey70', linetype = 'dotted', linewidth = 0.4))
}
for (mt in c('FPR', 'FDR', 'TPR')) {
  png(file.path(out_dir, 'Figures', sprintf('Alg2_heldout_%s_%s.png', mt, tag)),
      height = 6, width = 9, units = 'in', res = 300)
  print(mk_box(mt, sprintf('Realised %s on held-out data (cutoff applied unchanged)', mt))); dev.off()
}

#---------------------
# (2) Cross-setting transfer: cutoff learned on setting A applied to setting B (n = matrix_n)
#---------------------
rows <- list()
for (a in setting_names) for (b in setting_names) for (t in targets[[b]]) if (t$n == matrix_n) {
  m <- perf(t$delta, t$nonmono, cut_of(a, matrix_n))
  rows[[length(rows) + 1]] <- data.frame(LearnedOn = a, AppliedTo = b, rep = t$rep, FPR = m[['FPR']])
}
tm <- bind_rows(rows) %>% group_by(LearnedOn, AppliedTo) %>%
  summarise(FPR = mean(FPR), .groups = 'drop') %>%
  mutate(LearnedOn = factor(LearnedOn, levels = rev(setting_names)),
         AppliedTo = factor(AppliedTo, levels = setting_names))
write.csv(tm, file.path(out_dir, paste0('Alg2_transfer_matrix_n', matrix_n, '_', tag, '.csv')), row.names = FALSE)
tm_plot <- tm %>% filter(AppliedTo != 'dropout_after', LearnedOn != 'dropout_after')
p3 <- ggplot(tm_plot, aes(AppliedTo, LearnedOn, fill = FPR)) +
  geom_tile(colour = 'white') +
  geom_tile(data = ~ dplyr::filter(.x, as.character(LearnedOn) == as.character(AppliedTo)),
            fill = NA, colour = 'black', linewidth = 0.8) +
  geom_text(aes(label = sprintf('%.2f', FPR)), size = 2.4) +
  scale_fill_gradient2(low = '#9EC3E6', mid = 'white', high = '#D6604D', midpoint = alpha,
                       limits = c(0, 0.6), oob = scales::squish, name = 'Realised FPR') +
  labs(x = 'Applied to (held-out data of this setting)', y = 'Cutoff learned on',
       title = sprintf('Cross-setting transfer of the Algorithm 2 cutoff (n = %d)', matrix_n)) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                          axis.text = element_text(colour = 'black'),
                          legend.position = 'bottom')
png(file.path(out_dir, 'Figures', paste0('Alg2_transfer_matrix_', tag, '.png')),
    height = 8, width = 7, units = 'in', res = 300)
print(p3); dev.off()
