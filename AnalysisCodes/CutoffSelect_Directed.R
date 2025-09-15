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
dir = "Simulation_Directed/"

# -------------------------------
# Loading Ground Truth Data
# -------------------------------
## Load the simulated gene expression data (without noise)
df = read.csv(paste0(dir,'counts_1.csv'), header = T, row.names = 1) # Test load
# Transpose the dataframe to have genes as columns and cells as rows
df = t(df)

LoadGT <- function(symmetric = TRUE, edgelist = FALSE){
  ## Load the simulated gene expression data (without noise)
  df = read.csv(paste0(dir,'counts_1.csv'),
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

val_adj_logical <- LoadGT(symmetric = F)
val_adj <- matrix(as.integer(val_adj_logical), nrow = nrow(val_adj_logical))
table(val_adj_logical)
validation_edgelist <- LoadGT(symmetric = F, edgelist = T)

# Create scatter plots for different types of gene expression relationships
p1 = ggplot(df[,c(398,421)], aes(x = g398, y=g421))+
  geom_point() + geom_smooth(method = 'gam', se = F, colour = 'firebrick2') + theme_classic()+
  ggtitle('Linear')
p2 = ggplot(df[,c(324,456)], aes(x = g324, y=g456))+
  geom_point() + geom_smooth(method = 'gam', se = F, colour = 'firebrick2') + theme_classic()+
  ggtitle('Parabolic')
p3 = ggplot(df[,c(167,269)], aes(x = g167, y=g269))+
  geom_point() + geom_smooth(method = 'gam', se = F, colour = 'firebrick2') + theme_classic()+
  ggtitle('Exponential')
p3
p4 = ggplot(df[,c(129,409)], aes(x=g129, y=g409))+
  geom_point() + geom_smooth(method = 'gam', se = F, colour = 'firebrick2') + theme_classic()+
  ggtitle('Sinusoidal')
p5 = ggplot(df[,c(418,200)], aes(x=g418, y=g200))+
  geom_point() + geom_smooth(method = 'loess', se = F, colour = 'firebrick2') + theme_classic()+
  ggtitle('Random')
plot_grid(p1,p2,p3,p4,p5, ncol = 3)

# Save the plot grid as a PNG file
png('figures/scRNAseq_directed_simulation_scatter.png', height = 6, width = 9, units = 'in', res = 300)
plot_grid(p1,p2,p3,p4,p5, ncol = 3)
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
sym = FALSE                       # Metric used to select whether GT should be symmetric.

# -------------------------------
# Benchmarking Loop
# -------------------------------
# Load Ground Truth
val_adj_logical <- LoadGT(symmetric = sym)
validation_edgelist <- LoadGT(symmetric = T, edgelist = T)
linear_edgelist <- validation_edgelist %>% as.data.frame() %>% filter(directed == 'FALSE') %>%
  as.matrix()
val_adj_logical[linear_edgelist[,c(2,1)]] = TRUE
val_adj <- matrix(as.integer(val_adj_logical), nrow = nrow(val_adj_logical))
table(val_adj_logical)

# Define a function to calculate a cutoff for directed relationships
CalcCutDirected <- function(x,
                    niter = 100,
                    nsamp = 1500,
                    algorithm = 'Xi',
                    seed = 0,
                    linear_edgelist = linear_edgelist){

  cut95 <- c()
  FDRCut95 <- c()

  pb <- txtProgressBar(min = 0, max = niter, initial = 0, style = 3, width = 40)

  PR_val_adj <- val_adj
  diag(PR_val_adj) <- NA  # Remove diagonal values from ground truth

  for (i in 1:niter){
    # Subsampling
    set.seed(i)
    df <- x[sample(rownames(x), size = nsamp, replace = FALSE), ]

    linear_stat <- vector('numeric', length = nrow(linear_edgelist))
    for (i in 1:nrow(linear_edgelist)){
      cormat <- xicor_mpar(df[,linear_edgelist[i,1:2]])
      linear_stat[i] <- abs(cormat[2,1] - cormat[1,2])
    }

    # Compute 0.95 quantiles and store results
    cut95 <- append(cut95, quantile(linear_stat,0.95))

    setTxtProgressBar(pb, i)
  }
  close(pb)

  # Compute median cutoff values
  median_cutoff <- c("a95" = median(cut95, na.rm = TRUE))

  return(median_cutoff)
}

# Function to process a single dataset
process_dataset <- function(i) {
  # Load simulated gene expression data
  df <- read.csv(paste0(dir, 'counts_',i,'.csv'), header = TRUE, row.names = 1)
  df <- t(df)

  # Train-Test Split
  set.seed(0)
  sub_rows <- sample(rownames(df), train_size)
  df_train <- df[sub_rows, ]
  df_test <- df[!(rownames(df) %in% sub_rows), ]

  # Cut-off Calculation
  cutoffs_Xi <- CalcCutDirected(df_train, niter = 50, nsamp = test_size,
                                linear_edgelist = linear_edgelist)
  cutoffs_Xi

  Dependent_links <- vector('list', length = nrow(validation_edgelist))
  for (k in 1:nrow(validation_edgelist)){
    cormat <- xicor_mpar(df_test[,validation_edgelist[k,1:2]])
    Dependent_links[[k]] <- c(cormat[2,1],cormat[1,2])
  }
  Dependent_links <- as.data.frame(do.call('rbind',Dependent_links))
  colnames(Dependent_links) <- c('Forward','Backward')
  Dependent_links <- Dependent_links %>% mutate(diff = Forward - Backward)
  Dependent_links <- Dependent_links %>% mutate(absdiff = abs(diff))
  Dependent_links$directed <- validation_edgelist[,'directed']
  Dependent_links <- Dependent_links %>% mutate(rejected = absdiff > cutoffs_Xi)

  # Confusion Matrices
  cm_xi <- confusionMatrix(table(Dependent_links$rejected, Dependent_links$directed), positive = "TRUE")
  cm_xi

  # Store confusion matrices
  cms <- append(cms, list(list(Xi = cm_xi
  )))

  # Extract Performance Metrics
  metrics_list <- list()
  add_metric <- function(algorithm, metric, value) {
    metrics_list <<- append(metrics_list, list(c(Algorithm = algorithm, Metric = metric, Value = as.numeric(value))))
  }

  add_metric('Chatterjee', 'Sensitivity(TPR)', cm_xi$byClass['Sensitivity'])
  add_metric('Chatterjee', 'Specificity(TNR)', cm_xi$byClass['Specificity'])
  add_metric('Chatterjee', 'FPR', 1-cm_xi$byClass['Specificity'])
  add_metric('Chatterjee', 'Balanced Accuracy', cm_xi$byClass['Balanced Accuracy'])
  #add_metric('Chatterjee', 'RMSE', RMSE(pred = c(cor_xi_scaled), obs = c(val_adj), na.rm = TRUE))
  #auroc_xi <- auc(roc(as.numeric(val_adj), as.numeric(cor_xi_scaled)))
  #add_metric('Chatterjee', 'AUROC', auroc_xi)
  add_metric('Chatterjee', 'Cutoff', cutoffs_Xi[cMetric])

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
write.csv(metrics_df,paste0('Results/Benchmark/Cutoff_Selection_',tail(unlist(strsplit(dir,'/')),1),'_',cMetric,'_Directed.csv'))

#------
# Plotting
#------

# Reading the previously saved CSV file back into R and adding GRNBoost data from python
metrics_df <- read.csv(paste0('Results/Benchmark/Cutoff_Selection_',tail(unlist(strsplit(dir,'/')),1),'_',cMetric,'_Directed.csv'), row.names = 1)
p1 = ggplot(metrics_df,aes(x = Algorithm, y = Value, fill = Algorithm)) +
  geom_boxplot(width = 0.5, outliers = F)+
  geom_jitter(aes(x = Algorithm, y = Value), width = 0.1, size = 1, m=)+
  facet_wrap(~Metric, ncol = 4, scales = 'free')+
  theme_classic()+
  scale_fill_manual(values = hcl.colors(6,palette = 'Zissou'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = 'bold'))
p1

# Read GRNBoost results
grnboost_df <- read.csv(paste0('Results/Benchmark/Cutoff_Selection_Simulation_scRNAseq_Directed_GRNBoost.csv'))
colnames(grnboost_df) <-  c('Iter','Algorithm','Sensitivity(TPR)','Specificity(TNR)','FPR','Balanced Accuracy','Cutoff')
grnboost_df <- grnboost_df%>%
  mutate(FPR = 1-`Specificity(TNR)`)
grnboost_df <- reshape2::melt(grnboost_df, id.vars = c('Iter','Algorithm'))
colnames(grnboost_df) <-  c('Iter','Algorithm','Metric','Value')

grnboost_df <- grnboost_df %>%
  mutate(Algorithm = replace(Algorithm, Algorithm == "pcNet", "pcNetpy"))
grnboost_df <- grnboost_df[c('Algorithm','Metric','Value','Iter')]
merged_df <- rbind(metrics_df,grnboost_df) %>%
  filter(Algorithm != 'pcNetpy') %>%
  filter(Metric != 'Symmetric') %>%
  filter(Metric != 'Cutoff') %>%
  filter(Metric != 'Specificity(TNR)')
merged_df$Metric = factor(merged_df$Metric,
                          levels = c('FPR','Specificity(TNR)','Sensitivity(TPR)',
                                     'Balanced Accuracy','RMSE','AUROC','Cutoff'))
merged_df$Algorithm = factor(merged_df$Algorithm,
                             levels = c('Pearson','Spearman','pcNet','GRNBoost',
                                        'Chatterjee'))
merged_df$Value <- as.numeric(merged_df$Value)

# Save results
write.csv(merged_df,'Results/Cutoff_Selection_Simulation_scRNAseq_Directed_Merged.csv')

# Plot results
p1 = ggplot(merged_df,aes(x = Algorithm, y = Value, fill = Algorithm)) +
  geom_boxplot(width = 0.5, outliers = F)+
  geom_jitter(aes(x = Algorithm, y = Value), width = 0.1, size = 1, m=)+
  facet_wrap(~Metric, ncol = 4, scales = 'free')+
  theme_classic()+
  scale_fill_manual(values = hcl.colors(7, palette = 'RdPu', rev = T)[2:6][c(1,5)])+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = 'bold'))
p1

# Saving the plot as a PNG file with specified dimensions and resolution
png(paste0('Results/Benchmark/Figures/Bench_Simulation_scRNAseq_nonlinear_Directed.png'),
    height = 3, width = 6, units = 'in', res = 300)
p1
dev.off()
