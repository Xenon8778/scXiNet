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

#-----------------
# Loading data
#-----------------

# Load necessary scripts containing functions and dependencies
source('code/scripts.R')

# Define directory containing dataset files
dir = "Simulation_Directed/"

# -------------------------------
# Loading Ground Truth Data
# -------------------------------
## Load the simulated gene expression data (without noise), used only for the
## illustrative scatter plots below.
df = read.csv(paste0(dir,'counts_1.csv'), header = T, row.names = 1)
df = t(df)  # cells x genes

## Returns the ground-truth edgelist (TF, Target, directed) for the simulation.
## "directed" == TRUE  -> parabolic/sinusoidal (active regulation, TF -> Target)
## "directed" == FALSE -> linear/exponential   (dependent but no active regulation; true nulls)
LoadGT <- function(edgelist = TRUE){
  val = read.csv(paste0(dir,'gt_GRN.csv'), header = T, row.names = 1)
  validation_edgelist = as.matrix(apply(as.matrix(val[,1:3]), 2, as.character))

  if (edgelist){
    return(validation_edgelist)
  }

  # Adjacency-matrix form, only built when explicitly requested (not used
  # downstream in this script; kept for compatibility with other analyses).
  genes <- rownames(df)
  val_adj = matrix(0, nrow = length(genes), ncol = length(genes),
                   dimnames = list(genes, genes))
  val_adj[validation_edgelist[,1:2]] = 1
  diag(val_adj) = 0
  apply(val_adj, 2, as.logical)
}

# Create scatter plots for different types of gene expression relationships
p1 = ggplot(df[,c(398,421)], aes(x = g398, y=g421))+
  geom_point() + geom_smooth(method = 'gam', se = F, colour = 'magenta3') + theme_classic()+
  ggtitle('Linear')
p2 = ggplot(df[,c(324,456)], aes(x = g324, y=g456))+
  geom_point() + geom_smooth(method = 'gam', se = F, colour = 'magenta3') + theme_classic()+
  ggtitle('Parabolic')
p3 = ggplot(df[,c(167,269)], aes(x = g167, y=g269))+
  geom_point() + geom_smooth(method = 'gam', se = F, colour = 'magenta3') + theme_classic()+
  ggtitle('Exponential')
p4 = ggplot(df[,c(129,409)], aes(x=g129, y=g409))+
  geom_point() + geom_smooth(method = 'gam', se = F, colour = 'magenta3') + theme_classic()+
  ggtitle('Sinusoidal')
p5 = ggplot(df[,c(418,200)], aes(x=g418, y=g200))+
  geom_point() + geom_smooth(method = 'loess', se = F, colour = 'magenta3') + theme_classic()+
  ggtitle('Random')
plot_grid(p1,p2,p3,p4,p5, ncol = 5)

# Save the plot grid as a PNG file
png('figures/scRNAseq_directed_simulation_scatter_ver2.png', height = 3, width = 12, units = 'in', res = 300)
plot_grid(p1,p2,p3,p4,p5, ncol = 5)
dev.off()

########################################
# Cutoff Selection and Benchmarking
########################################

# Set parameters for benchmarking
niter = 50                        # Number of iterations for cutoff calculation
train_size = nrow(df)*3/5         # Number of samples used for training
test_size = nrow(df)-train_size   # Number of samples used for testing
cMetric = 'a95'                   # Metric used to select cutoff ('a95' or 'FDR95')

# -------------------------------
# Ground truth edgelist (shared across all datasets/iterations)
# -------------------------------
validation_edgelist <- LoadGT(edgelist = TRUE)
# True nulls for cutoff estimation (Algorithm 2, N0): dependent pairs with NO
# active regulation, i.e. linear + exponential interactions.
null_edgelist <- validation_edgelist %>% as.data.frame() %>%
  filter(directed == 'FALSE') %>% as.matrix()

# Cutoff estimation for testing active regulation (Algorithm 2), based on the
# 0.95 quantile of |xi(TF->Target) - xi(Target->TF)| among true-null pairs.
CalcCutDirected <- function(x,
                            niter = 100,
                            nsamp = 1500,
                            null_edgelist = null_edgelist){

  cut95 <- c()
  pb <- txtProgressBar(min = 0, max = niter, initial = 0, style = 3, width = 40)

  for (b in 1:niter){
    set.seed(b)
    df_sub <- x[sample(rownames(x), size = nsamp, replace = FALSE), ]

    null_stat <- vector('numeric', length = nrow(null_edgelist))
    for (j in 1:nrow(null_edgelist)){
      cormat <- xicor_mpar(df_sub[, null_edgelist[j,1:2]])
      null_stat[j] <- abs(cormat[2,1] - cormat[1,2])
    }

    cut95 <- append(cut95, quantile(null_stat, 0.95, na.rm = TRUE))
    setTxtProgressBar(pb, b)
  }
  close(pb)

  c("a95" = median(cut95, na.rm = TRUE))
}

# Process a single simulated dataset: estimate the cutoff on the training
# split, then evaluate (1) detection of active regulation and (2) direction
# accuracy on the held-out test split.
process_dataset <- function(i) {
  df <- read.csv(paste0(dir, 'counts_',i,'.csv'), header = TRUE, row.names = 1)
  df <- t(df)

  # Train-Test Split
  set.seed(0)
  sub_rows <- sample(rownames(df), train_size)
  df_train <- df[sub_rows, ]
  df_test <- df[!(rownames(df) %in% sub_rows), ]

  # Cut-off Calculation
  cutoffs_Xi <- CalcCutDirected(df_train, niter = 50, nsamp = test_size,
                                null_edgelist = null_edgelist)

  Dependent_links <- vector('list', length = nrow(validation_edgelist))
  for (k in 1:nrow(validation_edgelist)){
    cormat <- xicor_mpar(df_test[,validation_edgelist[k,1:2]])
    # xicor_mpar returns r[i,j] = xi(column_i -> column_j) (order by i, predict j).
    # Columns here are (TF, Target), so r[1,2] = xi(TF->Target) = Forward,
    # and r[2,1] = xi(Target->TF) = Backward.
    Dependent_links[[k]] <- c(cormat[1,2],cormat[2,1])
  }
  Dependent_links <- as.data.frame(do.call('rbind',Dependent_links))
  colnames(Dependent_links) <- c('Forward','Backward')  # Forward = xi(TF -> Target)
  Dependent_links <- Dependent_links %>% mutate(diff = Forward - Backward)
  Dependent_links <- Dependent_links %>% mutate(absdiff = abs(diff))
  Dependent_links$directed <- validation_edgelist[,'directed']
  Dependent_links <- Dependent_links %>% mutate(rejected = absdiff > cutoffs_Xi)

  # Confusion matrix: detection of active regulation (procedure 8)
  cm_xi <- confusionMatrix(table(Dependent_links$rejected, Dependent_links$directed), positive = "TRUE")

  # Direction accuracy: among pairs truly under active regulation
  # (directed == TRUE), does Forward (TF -> Target) exceed Backward
  # (Target -> TF)? NOTE: verify this orientation against xicor_mpar's
  # indexing convention in code/XICOR_mod.R before trusting this number.
  # Reported two ways:
  #   - unconditional: over ALL truly-directed pairs, including ones whose
  #     |diff| never cleared the cutoff (noise-dominated, coin-flip sign).
  #   - detected-only: restricted to truly-directed pairs that were also
  #     flagged as actively regulated (rejected == TRUE), matching how the
  #     paper's real-data section evaluates direction (only on already-
  #     confirmed links).
  dir_rows <- Dependent_links$directed == 'TRUE'
  direction_accuracy <- mean(Dependent_links$diff[dir_rows] > 0, na.rm = TRUE)

  dir_detected_rows <- dir_rows & Dependent_links$rejected
  direction_accuracy_detected <- mean(Dependent_links$diff[dir_detected_rows] > 0, na.rm = TRUE)

  # Extract Performance Metrics
  metrics_list <- list()
  add_metric <- function(algorithm, metric, value) {
    metrics_list <<- append(metrics_list, list(c(Algorithm = algorithm, Metric = metric, Value = as.numeric(value))))
  }

  add_metric('Chatterjee', 'Sensitivity(TPR)', cm_xi$byClass['Sensitivity'])
  add_metric('Chatterjee', 'Specificity(TNR)', cm_xi$byClass['Specificity'])
  add_metric('Chatterjee', 'FPR', 1-cm_xi$byClass['Specificity'])
  add_metric('Chatterjee', 'Balanced Accuracy', cm_xi$byClass['Balanced Accuracy'])
  add_metric('Chatterjee', 'DirectionAccuracy', direction_accuracy)
  add_metric('Chatterjee', 'DirectionAccuracyDetected', direction_accuracy_detected)
  add_metric('Chatterjee', 'Cutoff', cutoffs_Xi[cMetric])

  metrics_df <- do.call("rbind", metrics_list)
  metrics_df <- transform(metrics_df, Value = as.numeric(Value))
  metrics_df$Iter <- i
  list(metrics = metrics_df, confusion_matrix = cm_xi)
}

# Process each dataset
results <- lapply(1:10, process_dataset)

# Combine results
metrics_df <- do.call(rbind, lapply(results, function(x) x$metrics))

# Exporting the metrics data frame to a CSV file
write.csv(metrics_df,paste0('Results/Benchmark/Cutoff_Selection_',tail(unlist(strsplit(dir,'/')),1),'_',cMetric,'_Directed_ver2.csv'))

#------
# Plotting
#------

# Reading the previously saved CSV file back into R and adding GRNBoost data from python
metrics_df <- read.csv(paste0('Results/Benchmark/Cutoff_Selection_',tail(unlist(strsplit(dir,'/')),1),'_',cMetric,'_Directed_ver2.csv'), row.names = 1)
p1 = ggplot(metrics_df,aes(x = Algorithm, y = Value, fill = Algorithm)) +
  geom_boxplot(width = 0.5, outliers = F)+
  geom_jitter(aes(x = Algorithm, y = Value), width = 0.1, size = 1)+
  facet_wrap(~Metric, ncol = 4, scales = 'free')+
  theme_classic()+
  scale_fill_manual(values = hcl.colors(6,palette = 'Zissou'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = 'bold'))
p1

# Read GRNBoost results (produced by xicor_benchmark_directedv2_ver2.py)
grnboost_df <- read.csv(paste0('Results/Benchmark/Cutoff_Selection_Simulation_scRNAseq_Directed_GRNBoost_ver2.csv'),
                        row.names = NULL)
colnames(grnboost_df) <-  c('Iter','Algorithm','Sensitivity(TPR)','Specificity(TNR)','FPR','Balanced Accuracy','DirectionAccuracy','DirectionAccuracyDetected','Cutoff')
grnboost_df <- grnboost_df%>%
  mutate(FPR = 1-`Specificity(TNR)`)
grnboost_df <- reshape2::melt(grnboost_df, id.vars = c('Iter','Algorithm'))
colnames(grnboost_df) <-  c('Iter','Algorithm','Metric','Value')

grnboost_df <- grnboost_df[c('Algorithm','Metric','Value','Iter')]
merged_df <- rbind(metrics_df,grnboost_df) %>%
  filter(Metric != 'Specificity(TNR)') %>%
  filter(Metric != 'Cutoff')
merged_df$Metric = factor(merged_df$Metric,
                          levels = c('FPR','Specificity(TNR)','Sensitivity(TPR)',
                                     'Balanced Accuracy','DirectionAccuracy',
                                     'DirectionAccuracyDetected','Cutoff'))
merged_df$Algorithm = factor(merged_df$Algorithm,
                             levels = c('GRNBoost','Chatterjee'))
merged_df$Value <- as.numeric(merged_df$Value)

# Save results
write.csv(merged_df,'Results/Cutoff_Selection_Simulation_scRNAseq_Directed_Merged_ver2.csv')

# Plot results
plot_df <- merged_df %>%
  filter(Metric != 'DirectionAccuracy') %>%
  filter(Metric != 'Cutoff')
p1 = ggplot(plot_df,aes(x = Algorithm, y = Value, fill = Algorithm)) +
  geom_boxplot(width = 0.5, outliers = F)+
  geom_jitter(aes(x = Algorithm, y = Value), width = 0.1, size = 1)+
  facet_wrap(~Metric, ncol = 4, scales = 'free')+
  theme_classic()+
  scale_fill_manual(values = hcl.colors(7, palette = 'RdPu', rev = T)[2:6][c(1,5)])+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = 'bold'))+
  facetted_pos_scales(
    y = list(
      Metric == "DirectionAccuracyDetected" ~ scale_y_continuous(limits = c(0.8, 1)),
      Metric == "FPR" ~ scale_y_continuous(limits = c(0, 0.2))
    )
  )
p1

# Saving the plot as a PNG file with specified dimensions and resolution
png(paste0('Results/Benchmark/Figures/Bench_Simulation_scRNAseq_nonlinear_Directed_ver2.png'),
    height = 3, width = 9, units = 'in', res = 300)
p1
dev.off()
