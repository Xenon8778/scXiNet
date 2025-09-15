library(dplyr)
library(ggplot2)
library(XICOR)
library(ggpubr)
library(reshape2)
library(pROC)
library(caret)
library(ggsignif)
library(scTenifoldNet)
library(energy)
source('code/XICOR_mod.R')

#-----------------
# Loading data
#-----------------

# Load necessary scripts containing functions and dependencies
source('code/scripts.R')

# Define directory containing dataset files
dir = "SERGIO/data_sets/De-noised_100G_6T_300cPerT_dynamics_7_DS6/"

# -------------------------------
# Loading Ground Truth Data
# -------------------------------
## Load the simulated gene expression data (without noise)
df = read.csv(paste0(dir,'simulated_noNoise_T_0.csv'), header = T, row.names = 1) # Test load
# Transpose the dataframe to have genes as columns and cells as rows
df = t(df)

LoadGT <- function(symmetric = TRUE){
  ## Load the simulated gene expression data (without noise)
  df = read.csv(paste0(dir,'simulated_noNoise_T_0.csv'),
                header = T, row.names = 1) # Test load
  # Transpose the dataframe to have genes as columns and cells as rows
  df = t(df)

  # Load the ground truth gene regulatory network (GRN) as an edgelist
  val = read.csv(paste0(dir,'gt_GRN.csv'), header = F)

  # Convert the first two columns of the ground truth file into a character matrix representing gene interactions (regulatory relationships)
  validation_edgelist = as.matrix(apply(as.matrix(val[,1:2]),2, function(x) as.character(x)))

  # Initialize an adjacency matrix of zeros with dimensions equal to the number of genes
  val_adj = matrix(0,nrow = ncol(df), ncol = ncol(df))
  colnames(val_adj) <- colnames(df)
  rownames(val_adj) <- colnames(df)

  # Populate the adjacency matrix with ones where regulatory interactions exist
  val_adj[validation_edgelist] = 1
  if (symmetric){
    val_adj <- pmax(val_adj, t(val_adj)) # Make Symmetric
  }

  # Ensure that diagonal elements remain zero (no self-regulation)
  diag(val_adj) = 0

  # Convert the adjacency matrix into a logical matrix (TRUE = edge exists, FALSE = no edge)
  val_adj_logical = apply(val_adj,2,  as.logical)
  return(val_adj_logical)
}

val_adj_logical <- LoadGT(symmetric = T)
table(val_adj_logical)

########################################
# Cutoff Selection and Benchmarking
########################################

# Set parameters for benchmarking
niter = 50                       # Number of iterations for cutoff calculation
cms = list()                      # List to store confusion matrices
train_size = nrow(df)*3/5         # Number of samples used for training
test_size = nrow(df)-train_size   # Number of samples used for testing
cMetric = 'a95'                   # Metric used to select cutoff ('a95' or 'FDR95')
sym = TRUE                       # Metric used to select whether GT should be symmetric.

# -------------------------------
# Benchmarking Loop
# -------------------------------
# Load Ground Truth
val_adj_logical <- LoadGT(symmetric = sym)
val_adj <- matrix(as.integer(val_adj_logical), nrow = nrow(val_adj_logical))
table(val_adj_logical)

# Function to process a single dataset
process_dataset <- function(i) {
  # Load simulated gene expression data
  df <- read.csv(paste0(dir, 'simulated_noNoise_T_',i,'.csv'), header = TRUE, row.names = 1)
  df <- t(df)

  # Train-Test Split
  set.seed(0)
  sub_rows <- sample(rownames(df), train_size)
  df_train <- df[sub_rows, ]
  df_test <- df[!(rownames(df) %in% sub_rows), ]

  # Cut-off Calculation
  cutoffs_Xi <- CalcCut(df_train, niter = niter, nsamp = test_size, algorithm = 'Xi', val_mat = val_adj_logical)
  cutoffs_pcnet <- CalcCut(df_train, niter = niter, nsamp = test_size, algorithm = 'pcNet', val_mat = val_adj_logical)
  cutoffs_Pearson <- CalcCut(df_train, niter = niter, nsamp = test_size, algorithm = 'Pearson', val_mat = val_adj_logical)
  cutoffs_Spearman <- CalcCut(df_train, niter = niter, nsamp = test_size, algorithm = 'Spearman', val_mat = val_adj_logical)

  # Correlation Matrices and Apply Cutoffs
  cor_xi_test <- abs(xicor_mpar(df_test, nCores = 4))
  cor_xi_test <- pmax(cor_xi_test, t(cor_xi_test)) # Max Value Symmetric
  diag(cor_xi_test) <- 0
  cor_xi_scaled <- range01(cor_xi_test)
  cor_xi_filtered <- apply(cor_xi_test > cutoffs_Xi[cMetric], c(1, 2), function(x) as.integer(x))

  cor_pcnet_test <- abs(as.matrix(pcNet(t(df_test))))
  diag(cor_pcnet_test) <- 0
  cor_pcnet_scaled <- range01(cor_pcnet_test)
  cor_pcnet_filtered <- apply(cor_pcnet_test > cutoffs_pcnet[cMetric], c(1, 2), function(x) as.integer(x))

  cor_Spearman_test <- abs(cor(df_test, method = 'spearman'))
  diag(cor_Spearman_test) <- 0
  cor_Spearman_scaled <- range01(cor_Spearman_test)
  cor_Spearman_filtered <- apply(cor_Spearman_test > cutoffs_Spearman[cMetric], c(1, 2), function(x) as.integer(x))

  cor_Pearson_test <- abs(cor(df_test, method = 'pearson'))
  diag(cor_Pearson_test) <- 0
  cor_Pearson_scaled <- range01(cor_Pearson_test)
  cor_Pearson_filtered <- apply(cor_Pearson_test > cutoffs_Pearson[cMetric], c(1, 2), function(x) as.integer(x))

  # Confusion Matrices
  cm_xi <- confusionMatrix(table(cor_xi_filtered, val_adj), positive = "1")
  cm_pcNet <- confusionMatrix(table(cor_pcnet_filtered, val_adj), positive = "1")
  cm_spearman <- confusionMatrix(table(cor_Spearman_filtered, val_adj), positive = "1")
  cm_pearson <- confusionMatrix(table(cor_Pearson_filtered, val_adj), positive = "1")

  # Store confusion matrices
  cms <- append(cms, list(list(Xi = cm_xi, Spearman = cm_spearman, Pearson = cm_pearson)))

  # Extract Performance Metrics
  metrics_list <- list()
  add_metric <- function(algorithm, metric, value) {
    metrics_list <<- append(metrics_list, list(c(Algorithm = algorithm, Metric = metric, Value = as.numeric(value))))
  }

  # Add Cutoffs to metrics_list
  add_metric('Chatterjee', 'Cutoff', cutoffs_Xi[cMetric])
  add_metric('pcNet', 'Cutoff', cutoffs_pcnet[cMetric])
  add_metric('Pearson', 'Cutoff', cutoffs_Pearson[cMetric])
  add_metric('Spearman', 'Cutoff', cutoffs_Spearman[cMetric])

  add_metric('Spearman', 'Sensitivity(TPR)', cm_spearman$byClass['Sensitivity'])
  add_metric('Spearman', 'Specificity(TNR)', cm_spearman$byClass['Specificity'])
  add_metric('Spearman', 'FPR', 1-cm_spearman$byClass['Specificity'])
  add_metric('Spearman', 'Balanced Accuracy', cm_spearman$byClass['Balanced Accuracy'])
  add_metric('Spearman', 'RMSE', RMSE(pred = c(cor_Spearman_scaled), obs = c(val_adj), na.rm = TRUE))
  auroc_Spearman <- auc(roc(as.numeric(val_adj), as.numeric(cor_Spearman_scaled)))

  add_metric('Pearson', 'Sensitivity(TPR)', cm_pearson$byClass['Sensitivity'])
  add_metric('Pearson', 'Specificity(TNR)', cm_pearson$byClass['Specificity'])
  add_metric('Pearson', 'FPR', 1-cm_pearson$byClass['Specificity'])
  add_metric('Pearson', 'Balanced Accuracy', cm_pearson$byClass['Balanced Accuracy'])
  add_metric('Pearson', 'RMSE', RMSE(pred = c(cor_Pearson_scaled), obs = c(val_adj), na.rm = TRUE))
  auroc_Pearson <- auc(roc(as.numeric(val_adj), as.numeric(cor_Pearson_scaled)))

  add_metric('Chatterjee', 'Sensitivity(TPR)', cm_xi$byClass['Sensitivity'])
  add_metric('Chatterjee', 'Specificity(TNR)', cm_xi$byClass['Specificity'])
  add_metric('Chatterjee', 'FPR', 1-cm_xi$byClass['Specificity'])
  add_metric('Chatterjee', 'Balanced Accuracy', cm_xi$byClass['Balanced Accuracy'])
  add_metric('Chatterjee', 'RMSE', RMSE(pred = c(cor_xi_scaled), obs = c(val_adj), na.rm = TRUE))
  auroc_xi <- auc(roc(as.numeric(val_adj), as.numeric(cor_xi_scaled)))

  add_metric('pcNet', 'Sensitivity(TPR)', cm_pcNet$byClass['Sensitivity'])
  add_metric('pcNet', 'Specificity(TNR)', cm_pcNet$byClass['Specificity'])
  add_metric('pcNet', 'FPR', 1-cm_pcNet$byClass['Specificity'])
  add_metric('pcNet', 'Balanced Accuracy', cm_pcNet$byClass['Balanced Accuracy'])
  add_metric('pcNet', 'RMSE', RMSE(pred = c(cor_pcnet_scaled), obs = c(val_adj), na.rm = TRUE))
  auroc_pcnet <- auc(roc(as.numeric(val_adj), as.numeric(cor_pcnet_scaled)))

  # Combine Metrics into Dataframe
  metrics_df <- do.call("rbind", metrics_list)
  metrics_df <- transform(metrics_df, Value = as.numeric(Value))
  metrics_df$Iter <- i
  return(list(metrics = metrics_df, confusion_matrices = cms))
}

# Use sapply to process each dataset
results <- lapply(0:14, process_dataset)

# Combine results
metrics_df <- do.call(rbind, lapply(results, function(x) x$metrics))

# Exporting the metrics data frame to a CSV file
write.csv(metrics_df,paste0('Results/Benchmark/Cutoff_Selection_',tail(unlist(strsplit(dir,'/')),1),'_',cMetric,'_GTSym',sym,'.csv'))

# -------------------------------
# Reading the previously saved CSV file back into R and adding GRNBoost data from python
# -------------------------------
metrics_df <- read.csv(paste0('Results/Benchmark/Cutoff_Selection_',tail(unlist(strsplit(dir,'/')),1),'_a95_GTSym',sym,'.csv'), row.names = 1)
p1 = ggplot(metrics_df,aes(x = Algorithm, y = Value, fill = Algorithm)) +
  geom_boxplot(width = 0.5, position = position_dodge(0.9), outlier.size=0, size=0.5)+
  geom_jitter(aes(x = Algorithm, y = Value), width = 0.1, size = 1, m=)+
  facet_wrap(~Metric, ncol = 6, scales = 'free')+
  theme_classic()+
  scale_fill_manual(values = hcl.colors(6,palette = 'Zissou'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = 'bold'))
p1

grnboost_df <- read.csv(paste0('Results/Benchmark/grnboost/Cutoff_Selection_',tail(unlist(strsplit(dir,'/')),1),'_Symmetric.csv'))
colnames(grnboost_df) <-  c('Iter','Algorithm','Sensitivity(TPR)','Specificity(TNR)','Balanced Accuracy','RMSE')
grnboost_df <- grnboost_df%>%
  mutate(FPR = 1-`Specificity(TNR)`)
grnboost_df <- reshape2::melt(grnboost_df, id.vars = c('Iter','Algorithm'))
colnames(grnboost_df) <-  c('Iter','Algorithm','Metric','Value')

grnboost_df <- grnboost_df %>%
  mutate(Algorithm = replace(Algorithm, Algorithm == "pcNet", "pcNetpy"))
grnboost_df <- grnboost_df[c('Algorithm','Metric','Value','Iter')]
merged_df <- rbind(metrics_df,grnboost_df) %>%
  filter(Algorithm != 'pcNetpy') %>%
  filter(Metric != 'MCC') %>%
  filter(Metric != 'Specificity(TNR)')%>%
  filter(Metric != 'Cutoff')
merged_df$Metric = factor(merged_df$Metric,
                          levels = c('FPR','Specificity(TNR)','Sensitivity(TPR)',
                                     'Balanced Accuracy','RMSE'))
merged_df$Algorithm = factor(merged_df$Algorithm,
                             levels = c('Pearson','Spearman','pcNet', 'GRNBoost',
                                        'Chatterjee'))
# Save results
write.csv(merged_df,paste0('Results/Cutoff_Selection_SERGIO_',tail(unlist(strsplit(dir,'/')),1),'_Merged.csv'))

# Creating a plot for the benchmarking results
p1 = ggplot(merged_df,aes(x = Algorithm, y = Value, fill = Algorithm)) +
  geom_boxplot(width = 0.5, position = position_dodge(0.9), outlier.size=0.5, size=0.5)+
  geom_jitter(aes(x = Algorithm, y = Value), width = 0.1, size = 0.5, m=)+
  facet_wrap(~Metric, ncol = 4, scales = 'free')+
  labs(x=NULL,y=NULL)+
  theme_classic()+ theme(legend.position="none")+
  scale_fill_manual(values = hcl.colors(7, palette = 'RdPu', rev = T)[2:6])+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, family = 'Arial', colour = 'black'),
        axis.text.y = element_text(family = 'Arial', colour = 'black'))
p1
# Saving the plot as a PNG file with specified dimensions and resolution
png(paste0('Results/Benchmark/Figures/Bench_',tail(unlist(strsplit(dir,'/')),1),'.png'),
   height = 3, width = 7, units = 'in', res = 300)
p1
dev.off()

