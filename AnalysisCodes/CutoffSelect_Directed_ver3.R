library(dplyr)
library(ggplot2)
library(XICOR)
library(ggpubr)
library(reshape2)
library(pROC)
library(caret)
library(energy)
library(ggsignif)
library(scTenifoldNet)
library(cowplot)
library(ggh4x)
source('code/XICOR_mod.R')
source('code/cutoffCalcScripts.R')
source('code/XIMCOR_mod.R')       # XIM_Mat: revised xi (Lin & Han 2023)

dir = "Simulation_Directed/"

df = read.csv(paste0(dir,'counts_1.csv'), header = T, row.names = 1)
df = t(df)  # cells x genes

LoadGT <- function(edgelist = TRUE){
  val = read.csv(paste0(dir,'gt_GRN.csv'), header = T, row.names = 1)
  validation_edgelist = as.matrix(apply(as.matrix(val[,1:3]), 2, as.character))
  if (edgelist) return(validation_edgelist)
  genes <- rownames(df)
  val_adj = matrix(0, nrow = length(genes), ncol = length(genes),
                   dimnames = list(genes, genes))
  val_adj[validation_edgelist[,1:2]] = 1
  diag(val_adj) = 0
  apply(val_adj, 2, as.logical)
}

# (Scatter-plot block unchanged; omitted here for brevity)

########################################
# Cutoff Selection and Benchmarking
########################################
niter      = 50
train_size = nrow(df)*3/5
test_size  = nrow(df) - train_size
alpha      = 0.05
cut_types  = c('a95', 'FDR95')   # FPR-controlled and FDR-controlled cutoffs, both evaluated
cut_labels = c(a95 = 'FPR-controlled', FDR95 = 'FDR-controlled')
xi_M       = 5                          # right nearest neighbours for ChatterjeeKNN
methods    = c(Chatterjee = 'Xi', ChatterjeeKNN = 'Xi_KNN')

validation_edgelist <- LoadGT(edgelist = TRUE)

# Signed asymmetry xi(TF->Target) - xi(Target->TF) for every pair in `edgelist`.
# xicor_mpar and XIM_Mat: r[i,j] = xi(column_i -> column_j); columns are (TF, Target).
DirStat <- function(x, edgelist, method = 'Xi'){
  xi_fun <- switch(method,
                   Xi     = function(m) xicor_mpar(m),
                   Xi_KNN = function(m) XIM_Mat(m, xi_M, verbose = FALSE))
  vapply(seq_len(nrow(edgelist)), function(j){
    r <- xi_fun(x[, edgelist[j, 1:2]])
    r[1,2] - r[2,1]
  }, numeric(1))
}

# Empirical FDR cutoff: largest rejection set (ranked by |stat|) whose
# false-discovery proportion is <= alpha. Returned value is used as `stat > cut`.
FDRCut <- function(stat, is_null, alpha = 0.05){
  ok <- !is.na(stat)
  stat <- stat[ok]; is_null <- is_null[ok]
  o   <- order(stat, decreasing = TRUE)
  s   <- stat[o]
  fdp <- cumsum(is_null[o]) / seq_along(s)
  k   <- which(fdp <= alpha)
  if (length(k) == 0) return(max(s))            # nothing rejectable at this alpha
  k <- max(k)                                    # step-up: largest admissible set
  nxt <- if (k < length(s)) s[k + 1] else 0
  (s[k] + nxt) / 2                               # midpoint -> exactly top-k rejected
}

# Algorithm 2 extended: per resample, compute |Delta| for ALL labelled pairs once,
# then derive (i) FPR cutoff from nulls only, (ii) FDR cutoff from nulls + alternatives.
CalcCutDirected <- function(x, niter = 100, nsamp = 1500,
                            edgelist = validation_edgelist, alpha = 0.05,
                            method = 'Xi'){
  is_null <- edgelist[, 'nonmonotonic'] == 'FALSE'
  cutA <- cutF <- numeric(niter)
  pb <- txtProgressBar(min = 0, max = niter, style = 3, width = 40)
  for (b in seq_len(niter)){
    set.seed(b)
    df_sub <- x[sample(rownames(x), size = nsamp, replace = FALSE), ]
    stat <- abs(DirStat(df_sub, edgelist, method))
    cutA[b] <- quantile(stat[is_null], 1 - alpha, na.rm = TRUE)
    cutF[b] <- FDRCut(stat, is_null, alpha)
    setTxtProgressBar(pb, b)
  }
  close(pb)
  c(a95 = median(cutA, na.rm = TRUE), FDR95 = median(cutF, na.rm = TRUE))
}

process_dataset <- function(i) {
  df <- t(read.csv(paste0(dir, 'counts_', i, '.csv'), header = TRUE, row.names = 1))

  set.seed(0)
  sub_rows <- sample(rownames(df), train_size)
  df_train <- df[sub_rows, ]
  df_test  <- df[!(rownames(df) %in% sub_rows), ]

  truth    <- validation_edgelist[, 'nonmonotonic']
  dir_rows <- truth == 'TRUE'

  eval_method <- function(alg){
    cutoffs_Xi <- CalcCutDirected(df_train, niter = niter, nsamp = test_size,
                                  edgelist = validation_edgelist, alpha = alpha,
                                  method = methods[[alg]])
    # Test statistics computed once, then scored under each cutoff type
    set.seed(0)
    d        <- DirStat(df_test, validation_edgelist, methods[[alg]])
    direction_accuracy <- mean(d[dir_rows] > 0, na.rm = TRUE)   # cutoff-independent

    eval_cut <- function(type){
      cut      <- cutoffs_Xi[[type]]
      rejected <- abs(d) > cut
      cm <- confusionMatrix(table(factor(rejected, levels = c(FALSE, TRUE)),
                                  factor(truth,    levels = c('FALSE', 'TRUE'))),
                            positive = "TRUE")
      # Realised FDR (false discovery proportion) on the held-out test split
      n_rej <- sum(rejected, na.rm = TRUE)
      n_fp  <- sum(rejected & !dir_rows, na.rm = TRUE)
      fdp   <- if (n_rej > 0) n_fp / n_rej else 0
      dir_acc_det <- mean(d[dir_rows & rejected] > 0, na.rm = TRUE)

      data.frame(
        Algorithm  = alg,
        CutoffType = type,
        Metric = c('Sensitivity(TPR)', 'Specificity(TNR)', 'FPR', 'FDR',
                   'Balanced Accuracy', 'DirectionAccuracy', 'DirectionAccuracyDetected', 'Cutoff'),
        Value = as.numeric(c(cm$byClass['Sensitivity'], cm$byClass['Specificity'],
                             1 - cm$byClass['Specificity'], fdp,
                             cm$byClass['Balanced Accuracy'],
                             direction_accuracy, dir_acc_det, cut)),
        Iter = i)
    }

    do.call(rbind, lapply(cut_types, eval_cut))
  }

  list(metrics = do.call(rbind, lapply(names(methods), eval_method)))
}

results    <- lapply(1:10, process_dataset)
metrics_df <- do.call(rbind, lapply(results, `[[`, 'metrics'))

out_csv <- paste0('Results/Benchmark/Cutoff_Selection_', tail(unlist(strsplit(dir,'/')),1),
                  '_FPR_FDR_Directed_ver3.csv')
write.csv(metrics_df, out_csv)

#------
# Plotting
#------
metrics_df <- read.csv(out_csv, row.names = 1)

grnboost_df <- read.csv('Results/Benchmark/Cutoff_Selection_Simulation_scRNAseq_Directed_GRNBoost_FDR_ver2.csv',
                        row.names = NULL, check.names = FALSE)
colnames(grnboost_df) <- c('Iter','Algorithm','CutoffType','Sensitivity(TPR)','Specificity(TNR)','FPR','FDR',
                           'Balanced Accuracy','DirectionAccuracy','DirectionAccuracyDetected','Cutoff')
grnboost_df <- reshape2::melt(grnboost_df, id.vars = c('Iter','Algorithm','CutoffType'),
                              variable.name = 'Metric', value.name = 'Value')
grnboost_df <- grnboost_df[c('Algorithm','CutoffType','Metric','Value','Iter')]

merged_df <- rbind(metrics_df, grnboost_df) %>%
  filter(!Metric %in% c('Specificity(TNR)', 'Cutoff'))
merged_df$Metric <- factor(merged_df$Metric,
                           levels = c('FPR','FDR','Sensitivity(TPR)','Balanced Accuracy',
                                      'DirectionAccuracy','DirectionAccuracyDetected'))
merged_df$Algorithm  <- factor(merged_df$Algorithm, levels = c('GRNBoost','Chatterjee','ChatterjeeKNN'))
merged_df$CutoffType <- factor(cut_labels[as.character(merged_df$CutoffType)], levels = cut_labels)
merged_df$Value <- as.numeric(merged_df$Value)

write.csv(merged_df, 'Results/Cutoff_Selection_Simulation_scRNAseq_Directed_Merged_FPR_FDR_ver3.csv')

plot_df <- merged_df %>% filter(Metric != 'DirectionAccuracy')%>%
  droplevels()
plot_df$Metric <- recode_factor(plot_df$Metric, "DirectionAccuracyDetected" = "DirectionAccuracy\n(Detected)")
plot_df$Metric <- factor(plot_df$Metric, levels = c('FPR','FDR','Sensitivity(TPR)','Balanced Accuracy',
                                                    'DirectionAccuracy\n(Detected)'))

p1 <- ggplot(plot_df, aes(x = Algorithm, y = Value, fill = Algorithm)) +
  geom_boxplot(width = 0.5, outliers = F) +
  geom_jitter(width = 0.1, size = 1, height = 0) +
  facet_grid2(CutoffType ~ Metric, scales = 'free_y', independent = 'y') +
  theme_classic() +
  scale_fill_manual(values = hcl.colors(7, palette = 'RdPu', rev = T)[c(2, 4, 6)]) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none") +
  facetted_pos_scales(y = list(
    Metric == "DirectionAccuracy\n(Detected)" ~ scale_y_continuous(limits = c(0.8, 1)),
    Metric %in% c("FPR", "FDR") ~ scale_y_continuous(limits = c(0, 0.25))))
p1

png('Results/Benchmark/Figures/Bench_Simulation_scRNAseq_nonlinear_Directed_FPR_FDR_ver3.png',
    height = 5, width = 9, units = 'in', res = 300)
p1
dev.off()
