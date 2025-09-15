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
source('code/XICOR_mod.R')

#-----------------
# Loading data
#-----------------

# Load necessary scripts containing functions and dependencies
source('code/scripts.R')

# Define directory containing dataset files
dir = "Simulation_scRNAseq/"

# -------------------------------
# Loading Ground Truth Data
# -------------------------------
nCells <- 1000
## Load the simulated gene expression data (without noise)
df = read.csv(paste0(dir,'counts_1_',nCells,'.csv'), header = T, row.names = 1) # Test load
# Transpose the dataframe to have genes as columns and cells as rows
df = t(df)

LoadGT <- function(symmetric = TRUE, edgelist = FALSE){
  ## Load the simulated gene expression data (without noise)
  df = read.csv(paste0(dir,'counts_1_',nCells,'.csv'),
                header = T, row.names = 1) # Test load
  # Transpose the dataframe to have genes as columns and cells as rows
  df = t(df)

  # Load the ground truth gene regulatory network (GRN) as an edgelist
  val = read.csv(paste0(dir,'gt_GRN.csv'), header = T, row.names = 1)

  # Convert the first two columns of the ground truth file into a character matrix representing gene interactions (regulatory relationships)
  validation_edgelist = as.matrix(apply(as.matrix(val[,1:3]),2, function(x) as.character(x)))

  # Initialize an adjacency matrix of zeros with dimensions equal to the number of genes
  val_adj = matrix(0,nrow = ncol(df), ncol = ncol(df))
  colnames(val_adj) <- colnames(df)
  rownames(val_adj) <- colnames(df)

  # Populate the adjacency matrix with ones where regulatory interactions exist
  val_adj[validation_edgelist[,1:2]] = 1
  if (symmetric){
    val_adj <- pmax(val_adj, t(val_adj)) # Make Symmetric
  }

  # Ensure that diagonal elements remain zero (no self-regulation)
  diag(val_adj) = 0

  # Convert the adjacency matrix into a logical matrix (TRUE = edge exists, FALSE = no edge)
  val_adj_logical = apply(val_adj,2,  as.logical)
  rownames(val_adj_logical) <- rownames(val_adj)
  if (edgelist){
    return(validation_edgelist)
  } else {
    return(val_adj_logical)
  }
}


val_adj_logical <- LoadGT(symmetric = T)
table(val_adj_logical)
validation_edgelist <- LoadGT(symmetric = T,edgelist = T)


#df <- df / rowSums(df) * 10000
p1 = ggplot(df[,validation_edgelist[1,1:2]], aes(x = g14, y=g99))+
  geom_point() + geom_smooth(method = 'gam', se = F, colour = 'firebrick2') + theme_classic()+
  ggtitle('Linear')
p2 = ggplot(df[,validation_edgelist[6,1:2]], aes(x = g87, y=g50))+
  geom_point() + geom_smooth(method = 'gam', se = F, colour = 'firebrick2') + theme_classic()+
  ggtitle('Parabolic')
p2
p3 = ggplot(df[,validation_edgelist[7,1:2]], aes(x = g43, y=g65))+
  geom_point() + geom_smooth(method = 'gam', se = F, colour = 'firebrick2') + theme_classic()+
  ggtitle('Exponential')
p3
p4 = ggplot(df[,validation_edgelist[8,1:2]], aes(x=g100, y=g11))+
  geom_point() +
  geom_smooth(method = 'gam', se = F, colour = 'firebrick2') + theme_classic()+
  ggtitle('Sinusoidal')
p5 = ggplot(df[,c(39,3)], aes(x=g39, y=g3))+
  geom_point() + theme_classic()+
  ggtitle('Random')
plot_grid(p1,p2,p3,p4,p5, ncol = 3)

png('figures/scRNAseq_nonlinear_simulation_scatter.png', height = 3, width = 12, units = 'in', res = 300)
plot_grid(p1,p2,p3,p4, ncol = 4)
dev.off()

########################################
# Cutoff Selection and Benchmarking
########################################

# Set parameters for benchmarking
niter = 50                        # Number of iterations for cutoff calculation
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
  df <- read.csv(paste0(dir, 'counts_',i,'_',nCells,'.csv'), header = TRUE, row.names = 1)
  df <- t(df)

  # Train-Test Split
  set.seed(0)
  sub_rows <- sample(rownames(df), train_size)
  df_train <- df[sub_rows, ]
  df_test <- df[!(rownames(df) %in% sub_rows), ]

  # Cut-off Calculation
  cutoffs_Xi <- CalcCut(df_train, niter = niter, nsamp = test_size, algorithm = 'Xi', val_mat = val_adj_logical, symmetric = sym)
  cutoffs_pcnet <- CalcCut(df_train, niter = niter, nsamp = test_size, algorithm = 'pcNet', val_mat = val_adj_logical)
  cutoffs_Pearson <- CalcCut(df_train, niter = niter, nsamp = test_size, algorithm = 'Pearson', val_mat = val_adj_logical)
  cutoffs_Spearman <- CalcCut(df_train, niter = niter, nsamp = test_size, algorithm = 'Spearman', val_mat = val_adj_logical)
  # cutoffs_dcor <- CalcCut(df_train, niter = niter, nsamp = test_size, algorithm = 'dcor', val_mat = val_adj_logical)

  # Correlation Matrices and Apply Cutoffs
  cor_xi_test <- abs(xicor_mpar(df_test, nCores = 4))
  if (sym){
    cor_xi_test <- pmax(cor_xi_test, t(cor_xi_test)) # Max Value Symmetric
  }
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

  # Distance Correlation
  # cor_dcor_test <- matrix(0, nrow = ncol(df_test), ncol = ncol(df_test))
  # for (j in 1:ncol(df_test)) {
  #  for (k in 1:j){
  #   cor_dcor_test[j, k] <- dcor(df_test[, j], df_test[, k])
  #  }
  # }
  # cor_dcor_test <- pmax(cor_dcor_test, t(cor_dcor_test)) # Make Symmetric
  # diag(cor_dcor_test) <- 0
  # cor_dcor_scaled <- range01(cor_dcor_test)
  # cor_dcor_filtered <- apply(cor_dcor_test > cutoffs_dcor[cMetric], c(1, 2), function(x) as.integer(x))

  # Confusion Matrices
  cm_xi <- confusionMatrix(table(cor_xi_filtered, val_adj), positive = "1")
  cm_pcNet <- confusionMatrix(table(cor_pcnet_filtered, val_adj), positive = "1")
  cm_spearman <- confusionMatrix(table(cor_Spearman_filtered, val_adj), positive = "1")
  cm_pearson <- confusionMatrix(table(cor_Pearson_filtered, val_adj), positive = "1")
  # cm_dcor <- confusionMatrix(table(cor_dcor_filtered, val_adj), positive = "1")

  # Store confusion matrices
  cms <- append(cms, list(list(Xi = cm_xi, Spearman = cm_spearman,
                               # dcor = cm_dcor,
                               pcNet = cm_pcNet,
                               Pearson = cm_pearson
  )))

  # Extract Performance Metrics
  metrics_list <- list()
  add_metric <- function(algorithm, metric, value) {
    metrics_list <<- append(metrics_list, list(c(Algorithm = algorithm, Metric = metric, Value = as.numeric(value))))
  }

  add_metric('Spearman', 'Sensitivity(TPR)', cm_spearman$byClass['Sensitivity'])
  add_metric('Spearman', 'Specificity(TNR)', cm_spearman$byClass['Specificity'])
  add_metric('Spearman', 'FPR', 1-cm_spearman$byClass['Specificity'])
  add_metric('Spearman', 'Balanced Accuracy', cm_spearman$byClass['Balanced Accuracy'])
  add_metric('Spearman', 'RMSE', RMSE(pred = c(cor_Spearman_scaled), obs = c(val_adj), na.rm = TRUE))
  auroc_Spearman <- auc(roc(as.numeric(val_adj), as.numeric(cor_Spearman_scaled)))
  add_metric('Spearman', 'AUROC', auroc_Spearman)
  add_metric('Spearman', 'Cutoff', cutoffs_Spearman[cMetric])

  add_metric('Pearson', 'Sensitivity(TPR)', cm_pearson$byClass['Sensitivity'])
  add_metric('Pearson', 'Specificity(TNR)', cm_pearson$byClass['Specificity'])
  add_metric('Pearson', 'FPR', 1-cm_pearson$byClass['Specificity'])
  add_metric('Pearson', 'Balanced Accuracy', cm_pearson$byClass['Balanced Accuracy'])
  add_metric('Pearson', 'RMSE', RMSE(pred = c(cor_Pearson_scaled), obs = c(val_adj), na.rm = TRUE))
  auroc_Pearson <- auc(roc(as.numeric(val_adj), as.numeric(cor_Pearson_scaled)))
  add_metric('Pearson', 'AUROC', auroc_Pearson)
  add_metric('Pearson', 'Cutoff', cutoffs_Pearson[cMetric])

  add_metric('Chatterjee', 'Sensitivity(TPR)', cm_xi$byClass['Sensitivity'])
  add_metric('Chatterjee', 'Specificity(TNR)', cm_xi$byClass['Specificity'])
  add_metric('Chatterjee', 'FPR', 1-cm_xi$byClass['Specificity'])
  add_metric('Chatterjee', 'Balanced Accuracy', cm_xi$byClass['Balanced Accuracy'])
  add_metric('Chatterjee', 'RMSE', RMSE(pred = c(cor_xi_scaled), obs = c(val_adj), na.rm = TRUE))
  auroc_xi <- auc(roc(as.numeric(val_adj), as.numeric(cor_xi_scaled)))
  add_metric('Chatterjee', 'AUROC', auroc_xi)
  add_metric('Chatterjee', 'Cutoff', cutoffs_Xi[cMetric])

  add_metric('pcNet', 'Sensitivity(TPR)', cm_pcNet$byClass['Sensitivity'])
  add_metric('pcNet', 'Specificity(TNR)', cm_pcNet$byClass['Specificity'])
  add_metric('pcNet', 'FPR', 1-cm_pcNet$byClass['Specificity'])
  add_metric('pcNet', 'Balanced Accuracy', cm_pcNet$byClass['Balanced Accuracy'])
  add_metric('pcNet', 'RMSE', RMSE(pred = c(cor_pcnet_scaled), obs = c(val_adj), na.rm = TRUE))
  auroc_pcnet <- auc(roc(as.numeric(val_adj), as.numeric(cor_pcnet_scaled)))
  add_metric('pcNet', 'AUROC', auroc_pcnet)
  add_metric('pcNet', 'Cutoff', cutoffs_pcnet[cMetric])

  # add_metric('Distance', 'Sensitivity(TPR)', cm_dcor$byClass['Sensitivity'])
  # add_metric('Distance', 'Specificity(TNR)', cm_dcor$byClass['Specificity'])
  # add_metric('Distance', 'FPR', 1-cm_dcor$byClass['Specificity'])
  # add_metric('Distance', 'Balanced Accuracy', cm_dcor$byClass['Balanced Accuracy'])
  # add_metric('Distance', 'RMSE', RMSE(pred = c(cor_dcor_scaled), obs = c(val_adj), na.rm = TRUE))
  # auroc_dcor <- auc(roc(as.numeric(val_adj), as.numeric(cor_dcor_scaled)))

  # Combine Metrics into Dataframe
  metrics_df <- do.call("rbind", metrics_list)
  metrics_df <- transform(metrics_df, Value = as.numeric(Value))
  metrics_df$Iter <- i
  return(list(metrics = metrics_df, confusion_matrices = cms))
}


# Use sapply to process each dataset
results <- lapply(1:10, process_dataset)

# Combine results
metrics_df <- do.call(rbind, lapply(results, function(x) x$metrics))

# Exporting the metrics data frame to a CSV file
write.csv(metrics_df,paste0('Results/Benchmark/Cutoff_Selection_',tail(unlist(strsplit(dir,'/')),1),'_',cMetric,'_GTSym',sym,'_',nCells,'.csv'))

#------
# Plotting
#------
sym = TRUE   # Metric used to select whether GT should be symmetric.

# Reading the previously saved CSV file back into R and adding GRNBoost data from python
metrics_df <- read.csv(paste0('Results/Benchmark/Cutoff_Selection_',tail(unlist(strsplit(dir,'/')),1),'_a95_GTSym',sym,'_',nCells,'.csv'), row.names = 1)
p1 = ggplot(metrics_df,aes(x = Algorithm, y = Value, fill = Algorithm)) +
  geom_boxplot(width = 0.5, outliers = F)+
  geom_jitter(aes(x = Algorithm, y = Value), width = 0.1, size = 1, m=)+
  facet_wrap(~Metric, ncol = 4, scales = 'free')+
  theme_classic()+
  scale_fill_manual(values = hcl.colors(6,palette = 'Zissou'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = 'bold'))
p1

grnboost_df <- read.csv(paste0('Results/Benchmark/grnboost/Cutoff_Selection_Simulation_scRNAseq_',nCells,'_Symmetric.csv'))
colnames(grnboost_df) <-  c('Iter','Algorithm','Sensitivity(TPR)','Specificity(TNR)','Balanced Accuracy','RMSE','AUROC','Cutoff','Symmetric')
grnboost_df <- grnboost_df%>%
  mutate(FPR = 1-`Specificity(TNR)`)
grnboost_df <- reshape2::melt(grnboost_df, id.vars = c('Iter','Algorithm'))
colnames(grnboost_df) <-  c('Iter','Algorithm','Metric','Value')

grnboost_df <- grnboost_df %>%
  mutate(Algorithm = replace(Algorithm, Algorithm == "pcNet", "pcNetpy"))
grnboost_df <- grnboost_df[c('Algorithm','Metric','Value','Iter')]
merged_df <- rbind(metrics_df,grnboost_df) %>%
  filter(Algorithm != 'pcNetpy') %>%
  filter(Metric != 'Specificity(TNR)') %>%
  filter(Metric != 'Symmetric') %>%
  filter(Metric != 'Cutoff')
merged_df$Metric = factor(merged_df$Metric,
                          levels = c('FPR','Specificity(TNR)','Sensitivity(TPR)',
                                     'Balanced Accuracy','RMSE','AUROC','Cutoff'))
merged_df$Algorithm = factor(merged_df$Algorithm,
                          levels = c('Pearson','Spearman','pcNet','GRNBoost',
                                     'Chatterjee'))
merged_df$Value <- as.numeric(merged_df$Value)

# Save results
write.csv(merged_df,'Results/Cutoff_Selection_Simulation_scRNAseq_Merged.csv')

p1 = ggplot(merged_df,aes(x = Algorithm, y = Value, fill = Algorithm)) +
  geom_boxplot(width = 0.5, position = position_dodge(0.9), outlier.size=0.5, size=0.5)+
  geom_jitter(aes(x = Algorithm, y = Value), width = 0.1, size = 0.5, m=)+
  facet_wrap(~Metric, ncol = 6, scales = 'free')+
  labs(x=NULL,y=NULL)+
  theme_classic()+ theme(legend.position="none")+
  theme(text=element_text(family="Arial"), axis.text.x = element_text(family="Arial", colour = 'black'))+
  scale_fill_manual(values = hcl.colors(7, palette = 'RdPu', rev = T)[2:6])+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, family = 'Arial', colour = 'black'),
        axis.text.y = element_text(family = 'Arial', colour = 'black'))
p1

# Saving the plot as a PNG file with specified dimensions and resolution
png(paste0('Results/Benchmark/Figures/Bench_Simulation_scRNAseq_nonlinear_',nCells,'.png'),
    height = 3, width = 9, units = 'in', res = 300)
p1
dev.off()

