library(dplyr)
library(XICOR)
library(scTenifoldNet)
library(energy)
library(caret)
library(pROC)
library(PRROC)

# Load modified Xi Correlation functions
source('code/XIMCOR_mod.R')
source('code/XICOR_mod.R')

##############
# Function: range01
# Description: Normalizes a numeric vector to a range of [0,1].
# Arguments:
#   x - Numeric vector.
# Returns:
#   Normalized numeric vector within [0,1] range.
##############
range01 <- function(x){
  (x - min(x)) / (max(x) - min(x))
}

##############
# Function: CalcCut
# Description: Computes statistical cutoffs for correlation values based on subsampling and different correlation algorithms.
# Arguments:
#   x        - Data frame containing numerical variables.
#   niter    - Number of iterations for subsampling (default: 100).
#   nsamp    - Number of samples to draw in each iteration (default: 1500).
#   algorithm - Correlation method to use ('Xi', 'Xi_KNN', 'Spearman', or 'Pearson') (default: 'Xi').
#   seed     - Random seed for reproducibility (default: 0).
#   val_mat  - Ground Truth Matrix of values to be ignored (set as NA).
# Returns:
#   Named vector of median cutoff values at 90%, 95%, and 99% quantiles.
##############
CalcCut <- function(x,
                    alpha = 0.05,
                    niter = 100,
                    nsamp = 1500,
                    algorithm = 'Xi',
                    seed = 0,
                    symmetric = TRUE,
                    nCores = 1,
                    val_mat){

  cut95 <- c()
  FDRCut95 <- c()

  pb <- txtProgressBar(min = 0, max = niter, initial = 0, style = 3, width = 20)

  PR_val_adj <- val_mat
  diag(PR_val_adj) <- NA  # Remove diagonal values from ground truth

  for (i in 1:niter){
    # Subsampling
    set.seed(i)
    df <- x[sample(rownames(x), size = nsamp, replace = FALSE), ]

    # Compute correlation matrix based on the selected algorithm
    if (algorithm == 'Xi'){
      set.seed(seed)
      res <- xicor_mpar(df, nCores = nCores)  # Compute Xi Correlation
      if(symmetric){res = pmax(res,t(res))} # Max Value Symmetric
    } else if (algorithm == 'Xi_KNN'){
      set.seed(seed)
      res <- XIM_Mat(df, 10, verbose = FALSE)  # Compute Xi_KNN Correlation
    } else if (algorithm == 'Spearman'){
      res <- cor(df, method = 'spearman')  # Compute Spearman's Rank Correlation
    } else if (algorithm == 'Pearson'){
      res <- cor(df, method = 'pearson')  # Compute Pearson's Correlation
    } else if (algorithm == 'pcNet'){
      res <- as.matrix(pcNet(t(df), verbose = FALSE))  # Compute PC regression based network
    } else if (algorithm == 'dcor'){
      res <- matrix(0, nrow = ncol(df), ncol = ncol(df))
      for (j in 1:ncol(df)) {
        for (k in 1:j){
          res[j, k] <- dcor(df[, j], df[, k])
        }
      }
      res <- pmax(res, t(res)) # Make Symmetric
    } else {
      stop("Invalid algorithm specified.")
    }

    res <- abs(res)
    diag(res) <- NA  # Remove diagonal values

    # PR Curve to estimate cutoff associated to FDR
    pr <- pr.curve(scores.class0 = na.omit(c(res)),
                   scores.class1 = na.omit(c(PR_val_adj)), curve=T)
    pr.curve <- data.frame(pr$curve)
    colnames(pr.curve) <- c('Recall', 'Precision', 'Val')
    Rc <- pr.curve[which.min(abs((1-alpha)-pr.curve$Precision)),]$Val # FDR = 1 - Precision
    FDRCut95 <- append(FDRCut95, Rc)

    res[val_mat] <- NA  # Remove True Positive values

    # Compute alpha quantiles and store results
    cut95 <- append(cut95, quantile(res, 1-alpha, na.rm = TRUE))

    setTxtProgressBar(pb, i)
  }
  close(pb)

  # Compute median cutoff values
  mean_cutoff <- c(median(cut95, na.rm = TRUE),
                   median(FDRCut95, na.rm = TRUE))
  names(mean_cutoff) <- c(paste0("a",100-alpha*100),paste0("FDR",100-alpha*100))

  return(mean_cutoff)
}
