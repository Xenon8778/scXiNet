library(XICOR)
library(Seurat)
library(ggplot2)
library(pbapply)
library(patchwork)
source('code/XICOR_mod.R')
source('code/scripts.R')

commonGenes <- read.csv('data/300CommonGenes.csv', row.names = 1)$x

#---------------------
# Load STRING database
#---------------------
STRINGedglist <- read.table('data/STRING/9606.protein.links.v12.0.txt.gz',
                            sep = ' ', header = T)
ProteinInfo <- read.csv('data/STRING/9606.protein.info.v12.0.txt.gz',
                        sep = '\t', header = T)
ProteinMap <- c(ProteinInfo$preferred_name)
names(ProteinMap) <- ProteinInfo$X.string_protein_id

STRINGedglist$protein1 <- ProteinMap[STRINGedglist$protein1]
STRINGedglist$protein2 <- ProteinMap[STRINGedglist$protein2]

#---------------------
# Load Data 1
#---------------------
# Load scRNA-seq data 1
so1 <- readRDS('data/PBMC/GSE239592_Tcells.rds')

soFilt1 <- so1[commonGenes,]
soFilt1

## Filter weak STRING links
STRINGedglist_sub <- STRINGedglist[STRINGedglist$protein1 %in% rownames(soFilt1) &
                                     STRINGedglist$protein2 %in% rownames(soFilt1) &
                                     STRINGedglist$combined_score > 400,]

STRINGadj <- matrix(0, nrow(soFilt1), nrow(soFilt1))
colnames(STRINGadj) <- rownames(soFilt1)
rownames(STRINGadj) <- rownames(soFilt1)
STRINGadj[as.matrix(STRINGedglist_sub[,1:2])] <- 1 # Fill adjacency matrix
STRINGadj[1:5,1:5]
dim(STRINGadj)
STRINGadj_logical = apply(STRINGadj,2,  as.logical)
STRINGadj_logical[1:5,1:5]

## Prep for cutoff selection
nCells <- 1800
dfCount <- as.data.frame(t(as.matrix(GetAssayData(soFilt1, layer='data'))))
set.seed(0)
dfCount <- dfCount[sample(rownames(dfCount), nCells),]
dim(dfCount)
dfCount[1:5,1:5]

train_size = as.integer(nrow(dfCount)*0.60)         # Number of samples used for training
test_size = nrow(dfCount)-train_size    # Number of samples used for testing

### Train-Test Split
set.seed(0)
sub_rows <- sample(rownames(dfCount), train_size)
df_train <- dfCount[sub_rows, ]
df_test <- dfCount[!(rownames(dfCount) %in% sub_rows), ]

cutoffs_Xi <- CalcCut(df_train, niter=50, nsamp=test_size, algorithm='Xi', alpha = 0.05,
                      val_mat=STRINGadj_logical, symmetric=TRUE, nCores = 4)
cutoffs_Xi

# Correlation Matrices and Apply Cutoffs
cor_xi_test <- abs(xicor_mpar(df_test, nCores = 8))
cor_xi_test <- pmax(cor_xi_test, t(cor_xi_test)) # Max Value Symmetric
diag(cor_xi_test) <- 0
cor_xi_filtered <- apply(cor_xi_test > cutoffs_Xi['a95'], c(1, 2), function(x) as.integer(x))

# Confusion Matrices
cm_xi <- confusionMatrix(table(cor_xi_filtered, STRINGadj), positive = "1")
cm_xi

getPrecision <- function(corData, c) {
  corData[corData < c] <- 0
  corData[corData > 0] <- 1
  all_class <- union(corData, STRINGadj)
  newtable <- table(factor(corData, all_class), factor(STRINGadj, all_class))
  cmI <- confusionMatrix(newtable, positive = "1")
  if(is.na(cmI$byClass['Precision'])){
    return(c(cmI$byClass,cmI$overall))
  } else {
    return(c(cmI$byClass,cmI$overall))
  }
}

getPrecision(corData=cor_xi_test, c=cutoffs_Xi['a95'])
Precisions <- as.data.frame(t(pbsapply(seq(0.01, 0.99, 0.01), FUN = getPrecision,
                                       corData=cor_xi_test, USE.NAMES = T)))
Precisions$cutoff = seq(0.01, 0.99, 0.01)


p1 <- ggplot(Precisions, aes(x=cutoff))+
  geom_line(aes(y=1-Specificity), na.rm = F, linewidth = 1, color='black')+
  geom_vline(xintercept = cutoffs_Xi['a95'], color='firebrick', linetype='dashed', linewidth = 0.7)+
  geom_hline(yintercept = 0.05, color='firebrick', linewidth = 0.7)+
  theme_classic()+xlab('Cutoff')+ylab('FPR')+
  ggtitle('GSE239592')
p1
q1 <- ggplot(Precisions, aes(x=cutoff))+
  geom_line(aes(y=F1), na.rm = F, linewidth = 1, color='black')+
  geom_vline(xintercept = cutoffs_Xi['a95'], color='firebrick', linewidth = 0.7)+
  theme_classic()+xlab('Cutoff')+
  ggtitle('GSE239592')

#---------------------
# Load Data 3
#---------------------
so3 <- readRDS('data/PBMC/GSM8331607_Tcells.rds')
so3
soFilt <- so3[commonGenes,]
soFilt

## Filter weak STRING links
STRINGedglist_sub <- STRINGedglist[STRINGedglist$protein1 %in% rownames(soFilt) &
                                     STRINGedglist$protein2 %in% rownames(soFilt) &
                                     STRINGedglist$combined_score > 400,]

STRINGadj <- matrix(0, nrow(soFilt), nrow(soFilt))
colnames(STRINGadj) <- rownames(soFilt)
rownames(STRINGadj) <- rownames(soFilt)
STRINGadj[as.matrix(STRINGedglist_sub[,1:2])] <- 1 # Fill adjacency matrix
STRINGadj[1:5,1:5]
dim(STRINGadj)
STRINGadj_logical = apply(STRINGadj,2,  as.logical)
STRINGadj_logical[1:5,1:5]

## Prep for cutoff selection
# nCells <- 2000
dfCount <- as.data.frame(t(as.matrix(GetAssayData(soFilt, layer='data'))))
set.seed(0)
dfCount <- dfCount[sample(rownames(dfCount), nCells),]
dim(dfCount)
dfCount[1:5,1:5]

train_size = as.integer(nrow(dfCount)*0.60)         # Number of samples used for training
test_size = nrow(dfCount)-train_size    # Number of samples used for testing

### Train-Test Split
set.seed(0)
sub_rows <- sample(rownames(dfCount), train_size)
df_train <- dfCount[sub_rows, ]
df_test <- dfCount[!(rownames(dfCount) %in% sub_rows), ]

# Correlation Matrices and Apply Cutoffs
cor_xi_test <- abs(xicor_mpar(df_test, nCores = 8))
cor_xi_test <- pmax(cor_xi_test, t(cor_xi_test)) # Max Value Symmetric
diag(cor_xi_test) <- 0
cor_xi_filtered <- apply(cor_xi_test > cutoffs_Xi['a95'], c(1, 2), function(x) as.integer(x))

# Confusion Matrices
cm_xi <- confusionMatrix(table(cor_xi_filtered, STRINGadj), positive = "1")
cm_xi

getPrecision <- function(corData, c) {
  corData[corData < c] <- 0
  corData[corData > 0] <- 1
  all_class <- union(corData, STRINGadj)
  newtable <- table(factor(corData, all_class), factor(STRINGadj, all_class))
  cmI <- confusionMatrix(newtable, positive = "1")
  if(is.na(cmI$byClass['Precision'])){
    return(c(cmI$byClass,cmI$overall))
  } else {
    return(c(cmI$byClass,cmI$overall))
  }
}

getPrecision(corData=cor_xi_test, c=cutoffs_Xi['a95'])
Precisions <- as.data.frame(t(pbsapply(seq(0.01, 0.99, 0.01), FUN = getPrecision,
                                       corData=cor_xi_test, USE.NAMES = T)))
Precisions$cutoff = seq(0.01, 0.99, 0.01)


p2 <- ggplot(Precisions, aes(x=cutoff))+
  geom_line(aes(y=1-Specificity), na.rm = F, linewidth = 1, color='black')+
  geom_vline(xintercept = cutoffs_Xi['a95'], color='firebrick', linetype='dashed', linewidth = 0.7)+
  geom_hline(yintercept = 0.05, color='firebrick', linewidth = 0.7)+
  theme_classic()+xlab('Cutoff')+ylab('FPR')+
  ggtitle('GSM8331607')
q2 <- ggplot(Precisions, aes(x=cutoff))+
  geom_line(aes(y=F1), na.rm = F, linewidth = 1, color='black')+
  geom_vline(xintercept = cutoffs_Xi['a95'], color='firebrick', linewidth = 0.7)+
  theme_classic()+xlab('Cutoff')+
  ggtitle('GSM8331607')

#---------------------
# Load Data 4
#---------------------
so4 <- readRDS('data/PBMC/GSM7873657_Tcells.rds')
so4
soFilt <- so4[commonGenes,]
soFilt

## Filter weak STRING links
STRINGedglist_sub <- STRINGedglist[STRINGedglist$protein1 %in% rownames(soFilt) &
                                     STRINGedglist$protein2 %in% rownames(soFilt) &
                                     STRINGedglist$combined_score > 400,]

STRINGadj <- matrix(0, nrow(soFilt), nrow(soFilt))
colnames(STRINGadj) <- rownames(soFilt)
rownames(STRINGadj) <- rownames(soFilt)
STRINGadj[as.matrix(STRINGedglist_sub[,1:2])] <- 1 # Fill adjacency matrix
STRINGadj[1:5,1:5]
dim(STRINGadj)
STRINGadj_logical = apply(STRINGadj,2,  as.logical)
STRINGadj_logical[1:5,1:5]

## Prep for cutoff selection
# nCells <- 2000
dfCount <- as.data.frame(t(as.matrix(GetAssayData(soFilt, layer='data'))))
set.seed(0)
dfCount <- dfCount[sample(rownames(dfCount), nCells),]
dim(dfCount)
dfCount[1:5,1:5]

train_size = as.integer(nrow(dfCount)*0.60)         # Number of samples used for training
test_size = nrow(dfCount)-train_size    # Number of samples used for testing

### Train-Test Split
set.seed(0)
sub_rows <- sample(rownames(dfCount), train_size)
df_train <- dfCount[sub_rows, ]
df_test <- dfCount[!(rownames(dfCount) %in% sub_rows), ]

# Correlation Matrices and Apply Cutoffs
cor_xi_test <- abs(xicor_mpar(df_test, nCores = 8))
cor_xi_test <- pmax(cor_xi_test, t(cor_xi_test)) # Max Value Symmetric
diag(cor_xi_test) <- 0
cor_xi_filtered <- apply(cor_xi_test > cutoffs_Xi['a95'], c(1, 2), function(x) as.integer(x))

# Confusion Matrices
cm_xi <- confusionMatrix(table(cor_xi_filtered, STRINGadj), positive = "1")
cm_xi

getPrecision <- function(corData, c) {
  corData[corData < c] <- 0
  corData[corData > 0] <- 1
  all_class <- union(corData, STRINGadj)
  newtable <- table(factor(corData, all_class), factor(STRINGadj, all_class))
  cmI <- confusionMatrix(newtable, positive = "1")
  if(is.na(cmI$byClass['Precision'])){
    return(c(cmI$byClass,cmI$overall))
  } else {
    return(c(cmI$byClass,cmI$overall))
  }
}

getPrecision(corData=cor_xi_test, c=cutoffs_Xi['a95'])
Precisions <- as.data.frame(t(pbsapply(seq(0.01, 0.99, 0.01), FUN = getPrecision,
                                       corData=cor_xi_test, USE.NAMES = T)))
Precisions$cutoff = seq(0.01, 0.99, 0.01)


p3 <- ggplot(Precisions, aes(x=cutoff))+
  geom_line(aes(y=1-Specificity), na.rm = F, linewidth = 1, color='black')+
  geom_vline(xintercept = cutoffs_Xi['a95'], color='firebrick', linetype='dashed', linewidth = 0.7)+
  geom_hline(yintercept = 0.05, color='firebrick', linewidth = 0.7)+
  theme_classic()+xlab('Cutoff')+ylab('FPR')+
  ggtitle('GSM7873657')
q3 <- ggplot(Precisions, aes(x=cutoff))+
  geom_line(aes(y=F1), na.rm = F, linewidth = 1, color='black')+
  geom_vline(xintercept = cutoffs_Xi['a95'], color='firebrick', linewidth = 0.7)+
  theme_classic()+xlab('Cutoff')+
  ggtitle('GSM7873657')


#---------------------
# Load Data 5
#---------------------
so5 <- readRDS('data/PBMC/GSM7873658_Tcells.rds')
so5
soFilt <- so5[commonGenes,]
soFilt

## Filter weak STRING links
STRINGedglist_sub <- STRINGedglist[STRINGedglist$protein1 %in% rownames(soFilt) &
                                     STRINGedglist$protein2 %in% rownames(soFilt) &
                                     STRINGedglist$combined_score > 400,]

STRINGadj <- matrix(0, nrow(soFilt), nrow(soFilt))
colnames(STRINGadj) <- rownames(soFilt)
rownames(STRINGadj) <- rownames(soFilt)
STRINGadj[as.matrix(STRINGedglist_sub[,1:2])] <- 1 # Fill adjacency matrix
STRINGadj[1:5,1:5]
dim(STRINGadj)
STRINGadj_logical = apply(STRINGadj,2,  as.logical)
STRINGadj_logical[1:5,1:5]

## Prep for cutoff selection
# nCells <- 2000
dfCount <- as.data.frame(t(as.matrix(GetAssayData(soFilt, layer='data'))))
set.seed(0)
dfCount <- dfCount[sample(rownames(dfCount), nCells),]
dim(dfCount)
dfCount[1:5,1:5]

train_size = as.integer(nrow(dfCount)*0.60)         # Number of samples used for training
test_size = nrow(dfCount)-train_size    # Number of samples used for testing

### Train-Test Split
set.seed(0)
sub_rows <- sample(rownames(dfCount), train_size)
df_train <- dfCount[sub_rows, ]
df_test <- dfCount[!(rownames(dfCount) %in% sub_rows), ]

# Correlation Matrices and Apply Cutoffs
cor_xi_test <- abs(xicor_mpar(df_test, nCores = 8))
cor_xi_test <- pmax(cor_xi_test, t(cor_xi_test)) # Max Value Symmetric
diag(cor_xi_test) <- 0
cor_xi_filtered <- apply(cor_xi_test > cutoffs_Xi['a95'], c(1, 2), function(x) as.integer(x))

# Confusion Matrices
cm_xi <- confusionMatrix(table(cor_xi_filtered, STRINGadj), positive = "1")
cm_xi

getPrecision <- function(corData, c) {
  corData[corData < c] <- 0
  corData[corData > 0] <- 1
  all_class <- union(corData, STRINGadj)
  newtable <- table(factor(corData, all_class), factor(STRINGadj, all_class))
  cmI <- confusionMatrix(newtable, positive = "1")
  if(is.na(cmI$byClass['Precision'])){
    return(c(cmI$byClass,cmI$overall))
  } else {
    return(c(cmI$byClass,cmI$overall))
  }
}

getPrecision(corData=cor_xi_test, c=cutoffs_Xi['a95'])
Precisions <- as.data.frame(t(pbsapply(seq(0.01, 0.99, 0.01), FUN = getPrecision,
                                       corData=cor_xi_test, USE.NAMES = T)))
Precisions$cutoff = seq(0.01, 0.99, 0.01)

p4 <- ggplot(Precisions, aes(x=cutoff))+
  geom_line(aes(y=1-Specificity), na.rm = F, linewidth = 1, color='black')+
  geom_vline(xintercept = cutoffs_Xi['a95'], color='firebrick', linetype='dashed', linewidth = 0.7)+
  geom_hline(yintercept = 0.05, color='firebrick', linewidth = 0.7)+
  theme_classic()+xlab('Cutoff')+ylab('FPR')+
  ggtitle('GSM7873658')
q4 <- ggplot(Precisions, aes(x=cutoff))+
  geom_line(aes(y=F1), na.rm = F, linewidth = 1, color='black')+
  geom_vline(xintercept = cutoffs_Xi['a95'], color='firebrick', linewidth = 0.7)+
  theme_classic()+xlab('Cutoff')+
  ggtitle('GSM7873658')




#---------------------
# Comparison plots
#---------------------

png(paste0('figures/RealData_FPR.png'),height = 3, width = 12, units = 'in', res = 300)
p1+p2+p3+p4+plot_layout(ncol = 4)&
  ylim(0,0.4)&
  xlim(0,0.5)
dev.off()

png(paste0('figures/RealData_F1Score.png'),height = 3, width = 12, units = 'in', res = 300)
q1+q2+q3+q4+plot_layout(ncol = 4)&
  xlim(0,0.7)
dev.off()

#---------------------
# Load Data 8
#---------------------
so8 <- readRDS('data/Intestine/GSM8683884_Intestine.rds')
so8
soFilt <- so8[commonGenes,]
soFilt

## Filter weak STRING links
STRINGedglist_sub <- STRINGedglist[STRINGedglist$protein1 %in% rownames(soFilt) &
                                     STRINGedglist$protein2 %in% rownames(soFilt) &
                                     STRINGedglist$combined_score > 400,]

STRINGadj <- matrix(0, nrow(soFilt), nrow(soFilt))
colnames(STRINGadj) <- rownames(soFilt)
rownames(STRINGadj) <- rownames(soFilt)
STRINGadj[as.matrix(STRINGedglist_sub[,1:2])] <- 1 # Fill adjacency matrix
STRINGadj[1:5,1:5]
dim(STRINGadj)
STRINGadj_logical = apply(STRINGadj,2,  as.logical)
STRINGadj_logical[1:5,1:5]

## Prep for cutoff selection
# nCells <- 2000
dfCount <- as.data.frame(t(as.matrix(GetAssayData(soFilt, layer='data'))))
set.seed(0)
dfCount <- dfCount[sample(rownames(dfCount), nCells),]
dim(dfCount)
dfCount[1:5,1:5]

train_size = as.integer(nrow(dfCount)*0.60)         # Number of samples used for training
test_size = nrow(dfCount)-train_size    # Number of samples used for testing

### Train-Test Split
set.seed(0)
sub_rows <- sample(rownames(dfCount), train_size)
df_train <- dfCount[sub_rows, ]
df_test <- dfCount[!(rownames(dfCount) %in% sub_rows), ]

cutoffs_Xi_8 <- CalcCut(df_train, niter=50, nsamp=test_size, algorithm='Xi', alpha = 0.05,
                      val_mat=STRINGadj_logical, symmetric=TRUE, nCores = 8)
cutoffs_Xi_8

# Correlation Matrices and Apply Cutoffs
cor_xi_test <- abs(xicor_mpar(df_test, nCores = 8))
cor_xi_test <- pmax(cor_xi_test, t(cor_xi_test)) # Max Value Symmetric
diag(cor_xi_test) <- 0
cor_xi_filtered <- apply(cor_xi_test > cutoffs_Xi['a95'], c(1, 2), function(x) as.integer(x))

# Confusion Matrices
cm_xi <- confusionMatrix(table(cor_xi_filtered, STRINGadj), positive = "1")
cm_xi

getPrecision <- function(corData, c) {
  corData[corData < c] <- 0
  corData[corData > 0] <- 1
  all_class <- union(corData, STRINGadj)
  newtable <- table(factor(corData, all_class), factor(STRINGadj, all_class))
  cmI <- confusionMatrix(newtable, positive = "1")
  if(is.na(cmI$byClass['Precision'])){
    return(c(cmI$byClass,cmI$overall))
  } else {
    return(c(cmI$byClass,cmI$overall))
  }
}

getPrecision(corData=cor_xi_test, c=cutoffs_Xi['a95'])
Precisions <- as.data.frame(t(pbsapply(seq(0.01, 0.99, 0.01), FUN = getPrecision,
                                       corData=cor_xi_test, USE.NAMES = T)))
Precisions$cutoff = seq(0.01, 0.99, 0.01)


p7 <- ggplot(Precisions, aes(x=cutoff))+
  geom_line(aes(y=1-Specificity), na.rm = F, linewidth = 1, color='black')+
  geom_vline(xintercept = cutoffs_Xi['a95'], color='firebrick', linetype='dashed', linewidth = 0.7)+
  geom_vline(xintercept = cutoffs_Xi_8['a95'], color='forestgreen', linetype='dashed', linewidth = 0.7)+
  geom_hline(yintercept = 0.05, color='firebrick', linewidth = 0.7)+
  theme_classic()+xlab('Cutoff')+ylab('FPR')+
  ggtitle('GSM8683884')

q7 <- ggplot(Precisions, aes(x=cutoff))+
  geom_line(aes(y=F1), na.rm = F, linewidth = 1, color='black')+
  geom_vline(xintercept = cutoffs_Xi['a95'], color='firebrick', linewidth = 0.7)+
  theme_classic()+xlab('Cutoff')+
  ggtitle('GSM8683884')


#---------------------
# Load Data 9
#---------------------
so9 <- readRDS('data/Intestine/GSM8683885_Intestine.rds')
so9

soFilt <- so9[commonGenes,]
soFilt

## Filter weak STRING links
STRINGedglist_sub <- STRINGedglist[STRINGedglist$protein1 %in% rownames(soFilt) &
                                     STRINGedglist$protein2 %in% rownames(soFilt) &
                                     STRINGedglist$combined_score > 400,]

STRINGadj <- matrix(0, nrow(soFilt), nrow(soFilt))
colnames(STRINGadj) <- rownames(soFilt)
rownames(STRINGadj) <- rownames(soFilt)
STRINGadj[as.matrix(STRINGedglist_sub[,1:2])] <- 1 # Fill adjacency matrix
STRINGadj[1:5,1:5]
dim(STRINGadj)
STRINGadj_logical = apply(STRINGadj,2,  as.logical)
STRINGadj_logical[1:5,1:5]

## Prep for cutoff selection
# nCells <- 2000
dfCount <- as.data.frame(t(as.matrix(GetAssayData(soFilt, layer='data'))))
set.seed(0)
dfCount <- dfCount[sample(rownames(dfCount), nCells),]
dim(dfCount)
dfCount[1:5,1:5]

train_size = as.integer(nrow(dfCount)*0.60)         # Number of samples used for training
test_size = nrow(dfCount)-train_size    # Number of samples used for testing

### Train-Test Split
set.seed(0)
sub_rows <- sample(rownames(dfCount), train_size)
df_train <- dfCount[sub_rows, ]
df_test <- dfCount[!(rownames(dfCount) %in% sub_rows), ]

# Correlation Matrices and Apply Cutoffs
cor_xi_test <- abs(xicor_mpar(df_test, nCores = 8))
cor_xi_test <- pmax(cor_xi_test, t(cor_xi_test)) # Max Value Symmetric
diag(cor_xi_test) <- 0
cor_xi_filtered <- apply(cor_xi_test > cutoffs_Xi['a95'], c(1, 2), function(x) as.integer(x))

# Confusion Matrices
cm_xi <- confusionMatrix(table(cor_xi_filtered, STRINGadj), positive = "1")
cm_xi

getPrecision <- function(corData, c) {
  corData[corData < c] <- 0
  corData[corData > 0] <- 1
  all_class <- union(corData, STRINGadj)
  newtable <- table(factor(corData, all_class), factor(STRINGadj, all_class))
  cmI <- confusionMatrix(newtable, positive = "1")
  if(is.na(cmI$byClass['Precision'])){
    return(c(cmI$byClass,cmI$overall))
  } else {
    return(c(cmI$byClass,cmI$overall))
  }
}

getPrecision(corData=cor_xi_test, c=cutoffs_Xi['a95'])
Precisions <- as.data.frame(t(pbsapply(seq(0.01, 0.99, 0.01), FUN = getPrecision,
                                       corData=cor_xi_test, USE.NAMES = T)))
Precisions$cutoff = seq(0.01, 0.99, 0.01)


p8 <- ggplot(Precisions, aes(x=cutoff))+
  geom_line(aes(y=1-Specificity), na.rm = F, linewidth = 1, color='black')+
  geom_vline(xintercept = cutoffs_Xi['a95'], color='firebrick', linetype='dashed', linewidth = 0.7)+
  geom_vline(xintercept = cutoffs_Xi_8['a95'], color='forestgreen', linetype='dashed', linewidth = 0.7)+
  geom_hline(yintercept = 0.05, color='firebrick', linewidth = 0.7)+
  theme_classic()+xlab('Cutoff')+ylab('FPR')+
  ggtitle('GSM8683885')
p8
q8 <- ggplot(Precisions, aes(x=cutoff))+
  geom_line(aes(y=F1), na.rm = F, linewidth = 1, color='black')+
  geom_vline(xintercept = cutoffs_Xi['a95'], color='firebrick', linewidth = 0.7)+
  theme_classic()+xlab('Cutoff')+
  ggtitle('GSM8683885')


#---------------------
# Comparison plots different celltype
#---------------------

png(paste0('figures/RealData_FPR_diffCT.png'),height = 3, width = 6, units = 'in', res = 300)
p7+p8&
  ylim(0,0.3)&
  xlim(0,0.3)
dev.off()

png(paste0('figures/RealData_F1Score_diffCT.png'),height = 3, width = 6, units = 'in', res = 300)
q7+q8&
  xlim(0,0.4)&
  geom_vline(xintercept = cutoffs_Xi_8['a95'], color='forestgreen', linetype='dashed', linewidth = 0.7)
dev.off()
