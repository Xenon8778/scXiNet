library(XICOR)
library(ggplot2)
library(dplyr)
library(pROC)
library(patchwork)
source('code/realDataSetup.R')   # gene panel, STRING GT, sampleCells, computeXi, alpha, nCells

niter   <- 50
n_cores <- 8

# ---- Method: ChatterjeeKNN (revised xi, Lin & Han 2023) ----
xi_alg <- 'Xi_KNN'   # CalcCut algorithm key ('Xi' for the original coefficient)
xi_M   <- 5         # right nearest neighbours; must match CalcCut's Xi_KNN
tag    <- paste0(xi_alg, '_ver3')

# Network on a cell x gene matrix: |xi|, symmetrised by max, NA -> 0, zero diagonal
computeNet <- function(df) {
  m <- if (xi_alg == 'Xi_KNN') { set.seed(0); XIM_Mat(df, xi_M, nCores = n_cores,
                                                      verbose = TRUE) } else
    xicor_mpar(df, nCores = n_cores)
  m <- abs(as.matrix(m)); dimnames(m) <- list(colnames(df), colnames(df))
  m[is.na(m)] <- 0
  m <- pmax(m, t(m))
  diag(m) <- 0
  m
}

# Same gene panel in every dataset -> one GT matrix
gt <- LoadGT_STRING(commonGenes)
eval_mask <- upper.tri(gt) & !is.na(gt)
labels    <- as.integer(gt[eval_mask])

#---------------------
# Helpers
#---------------------
# Metrics at both cutoffs over unique, non-excluded pairs (as in ver3)
evaluate_cuts <- function(dataset, scores, cuts) {
  s     <- as.numeric(scores[eval_mask])
  auroc <- as.numeric(auc(roc(labels, s, quiet = TRUE)))
  do.call(rbind, lapply(names(cut_types), function(ct) {
    cutoff <- as.numeric(cuts[cut_types[[ct]]])
    pred <- s > cutoff
    TP <- sum(pred & labels == 1L); FP <- sum(pred & labels == 0L)
    TN <- sum(!pred & labels == 0L); FN <- sum(!pred & labels == 1L)
    tpr <- TP / (TP + FN); tnr <- TN / (TN + FP)
    fdr <- if ((TP + FP) == 0) NA_real_ else FP / (TP + FP)
    data.frame(Dataset = dataset, CutoffType = ct,
               Metric = c('Sensitivity(TPR)', 'Specificity(TNR)', 'FPR', 'FDR',
                          'Balanced Accuracy', 'AUROC', 'Cutoff', 'nCalled'),
               Value  = c(tpr, tnr, 1 - tnr, fdr, (tpr + tnr) / 2, auroc, cutoff, TP + FP))
  }))
}

# FPR / FDR / TPR / F1 as a function of cutoff (vectorised; replaces caret loop)
evalCurve <- function(scores, cutoff_seq = seq(0.01, 0.99, 0.01)) {
  s <- as.numeric(scores[eval_mask])
  do.call(rbind, lapply(cutoff_seq, function(c) {
    pred <- s > c
    TP <- sum(pred & labels == 1L); FP <- sum(pred & labels == 0L)
    FN <- sum(!pred & labels == 1L); TN <- sum(!pred & labels == 0L)
    prec <- if ((TP + FP) == 0) NA_real_ else TP / (TP + FP)
    rec  <- TP / (TP + FN)
    data.frame(cutoff = c, FPR = FP / (FP + TN), FDR = 1 - prec, TPR = rec,
               F1 = if (is.na(prec) || (prec + rec) == 0) NA_real_ else 2 * prec * rec / (prec + rec))
  }))
}

cut_lines <- function(cuts) {
  d <- data.frame(CutoffType = names(cut_types), Cutoff = as.numeric(cuts[cut_types]))
  d[is.finite(d$Cutoff), ]               # infeasible FDR cutoff (Inf) not drawn
}

plotCurve <- function(curve, cuts, metric, title, type = 'FPR') {
  vl <- subset(cut_lines(cuts), CutoffType == type)
  ggplot(curve, aes(cutoff, .data[[metric]])) +
    geom_line(linewidth = 1, colour = 'black') +
    geom_hline(yintercept = alpha, colour = 'firebrick', linewidth = 0.7) +
    geom_vline(data = vl, aes(xintercept = Cutoff),
               colour = 'firebrick', linewidth = 0.7, linetype = 'dashed') +
    theme_classic() + labs(x = 'Cutoff', y = metric, title = title)
}

#---------------------
# Merged Harmony-corrected T-cell object
# (Harmony only alters the embedding; the 'data' layer used for Xi is untouched)
#---------------------
soPBMC <- readRDS('data/PBMC/PBMC_Tcells_merged.rds')

getCounts <- function(dataset_id, n = nCells) {
  so <- subset(soPBMC, subset = orig.dataset == dataset_id)
  sampleCells(so[commonGenes, ], n)   # n = NULL -> all cells
}

#---------------------
# Reference dataset (GSE239592): learn FPR + FDR cutoffs (Algorithm 1, bootstrap n = n')
#---------------------
# Pool = ALL reference cells (n' > n), so each iteration subsamples nCells WITHOUT
# replacement. Pre-sampling nCells first (n' = n) forces a bootstrap with duplicated
# cells, which inflates null xi and the cutoff.
dfRefAll <- getCounts(refDataset, n = NULL)
stopifnot(nrow(dfRefAll) > nCells)
dfCount1 <- getCounts(refDataset)    # fixed nCells sample, for the in-sample curve only

# Align GT to the data's gene order (data columns are not reordered; see sampleCells)
gt <- gt[colnames(dfCount1), colnames(dfCount1)]
eval_mask <- upper.tri(gt) & !is.na(gt)
labels    <- as.integer(gt[eval_mask])

cutoffs_Xi <- CalcCut(dfRefAll, alpha = alpha, niter = niter, nsamp = nCells, M = xi_M,
                      algorithm = xi_alg, val_mat = gt, nCores = n_cores, symmetric = TRUE,
                      seed = 42)
print(cutoffs_Xi)
if (!isTRUE(as.logical(cutoffs_Xi[feasible_key])))
  warning(sprintf('FDR target %.2f infeasible on reference; min achievable = %.3f',
                  alpha, cutoffs_Xi[minach_key]))

#---------------------
# Apply reference cutoffs unmodified to all datasets (reference = in-sample check)
#---------------------
allDatasets <- c(refDataset, 'GSM8331607', 'GSM7873657', 'GSM7873658')
res <- lapply(allDatasets, function(ds) {
  message('== ', ds)
  dfc <- if (ds == refDataset) dfCount1 else getCounts(ds)
  stopifnot(identical(colnames(dfc), rownames(gt)))
  net <- computeNet(dfc)
  list(net = net, metrics = evaluate_cuts(ds, net, cutoffs_Xi), curve = evalCurve(net))
})
names(res) <- allDatasets

metrics_df <- do.call(rbind, lapply(res, `[[`, 'metrics'))
metrics_df %>%
  filter(Metric %in% c('FPR', 'FDR', 'Sensitivity(TPR)', 'nCalled')) %>%
  tidyr::pivot_wider(names_from = Metric, values_from = Value) %>%
  print(n = Inf)

write.csv(metrics_df, paste0('Results/RealData_CutoffEval_', tag, '.csv'), row.names = FALSE)
write.csv(data.frame(Key = names(cutoffs_Xi), Value = unname(cutoffs_Xi)),
          paste0('Results/RealData_Cutoffs_GSE239592_', tag, '.csv'), row.names = FALSE)
saveRDS(lapply(res, `[[`, 'net'), paste0('Results/RealData_XiNets_', tag, '.rds'))  # for dependency/active-regulation tests

#---------------------
# Figures
#---------------------
for (m in c('FPR', 'FDR')) {
  pl <- lapply(allDatasets, function(ds) plotCurve(res[[ds]]$curve, cutoffs_Xi, m, ds))
  png(paste0('figures/RealData_', m, '_', tag, '.png'), height = 3, width = 10, units = 'in', res = 300)
  print(wrap_plots(pl, ncol = 4, guides = 'collect') & xlim(0, 0.3))
  dev.off()
}

f1_plots <- lapply(allDatasets, function(ds)
  ggplot(res[[ds]]$curve, aes(cutoff, F1)) +
    geom_line(linewidth = 1) +
    geom_vline(data = cut_lines(cutoffs_Xi), aes(xintercept = Cutoff, linetype = CutoffType),
               colour = 'firebrick', linewidth = 0.7) +
    scale_linetype_manual(values = c(FPR = 'dashed', FDR = 'dotted'), name = 'Cutoff') +
    theme_classic() + labs(x = 'Cutoff', title = ds))
png(paste0('figures/RealData_F1Score_', tag, '.png'), height = 3, width = 12, units = 'in', res = 300)
print(wrap_plots(f1_plots, ncol = 4, guides = 'collect') & xlim(0, 0.8))
dev.off()

# Realised FPR and FDR on held-out datasets, under each cutoff
heldout_df <- metrics_df %>%
  filter(Metric %in% c('FPR', 'FDR')) %>%
  mutate(Metric = factor(Metric, levels = c('FPR', 'FDR')))

errBar <- ggplot(heldout_df[heldout_df$Metric == 'FPR' & heldout_df$CutoffType == 'FPR',], aes(x = Dataset, y = Value, fill = Dataset)) +
  geom_col(colour = 'black', position = position_dodge(0.9), width = 0.8) +
  geom_hline(yintercept = alpha, colour = 'firebrick', linetype = 'dashed', linewidth = 0.7) +
  scale_fill_manual(values = scales::hue_pal()(4)[1:4]) +
  theme_classic() + ylim(0, 0.4) +
  labs(x = NULL, y = "FPR") +
  theme(strip.background = element_blank(), strip.text = element_text(face = 'bold'),
        legend.position = 'none',
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
errBar

png(paste0('figures/RealData_ErrorRates_heldout_', tag, '.png'), height = 3, width = 3, units = 'in', res = 300)
errBar
dev.off()
