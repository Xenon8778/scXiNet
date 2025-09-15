source('code/XICOR_mod.R')
source('code/scripts.R')
library(igraph)
library(ggraph)
library(cowplot)

#---------------------
# Load TFLink database
#---------------------
TcellCommonGenes <- read.csv('data/TcellCommonGenes.csv', row.names = 1)$x

TFLink <- read.csv('data/TFTargetDBs/TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv.gz', sep = '\t')
TFLinkTFTarget <- TFLink[c('Name.TF','Name.Target')] %>%
  mutate('TF.Target' = paste0(Name.TF,"_",Name.Target))
TFLinkTFTarget$database <- 'TFLink'
head(TFLinkTFTarget)

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
so1
# Identify highly variable gene by dispersion
# so1 <- so1[rownames(so1)[!startsWith(rownames(so1),'RPL') &
#                            !startsWith(rownames(so1),'RPS') &
#                            !startsWith(rownames(so1),'MT-')],]
# so1 <- FindVariableFeatures(so1, nfeatures = 2000,
#                                   selection.method = 'dispersion')

# soFilt1 <- so1[TcellCommonGenes,]
# soFilt1
soFilt1 <- so1

## Prep for cutoff selection
nCells <- 4000
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

## Filter weak STRING links
STRINGedglist_sub <- STRINGedglist[STRINGedglist$protein1 %in% rownames(soFilt1) &
                                     STRINGedglist$protein2 %in% rownames(soFilt1) &
                                     STRINGedglist$combined_score > 400,]
STRINGedglist_sub <- STRINGedglist_sub %>%
  mutate('Gene.Gene' = paste0(protein1,'_',protein2))
STRINGadj <- matrix(0, nrow(soFilt1), nrow(soFilt1))
colnames(STRINGadj) <- rownames(soFilt1)
rownames(STRINGadj) <- rownames(soFilt1)
STRINGadj[as.matrix(STRINGedglist_sub[,1:2])] <- 1 # Fill adjacency matrix
STRINGadj[1:5,1:5]
dim(STRINGadj)
STRINGadj_logical = apply(STRINGadj,2,  as.logical)


### Cutoff Selection, FPR = 0.05
cutoffs_Xi <- CalcCut(df_train, niter=50, nsamp=test_size, algorithm='Xi', alpha = 0.05,
                      val_mat=STRINGadj_logical, symmetric=TRUE, nCores = 8)
cutoffs_Xi
# cutoffs_Xi <- c('a95' = 0.06117073)

# Correlation Matrices and Apply Cutoffs
correlation_matrices_list <- list()
n_iters <- 10
for (i in 1:n_iters) {
  set.seed(i)
  current_cor_matrix <- abs(xicor_mpar(df_test, nCores = 8))
  correlation_matrices_list[[i]] <- current_cor_matrix
  cat(paste0("Completed repetition ", i, " of ", n_iters, "\n"))
}
summed_cor_matrix <- Reduce('+', correlation_matrices_list)
cor_xi_test <- summed_cor_matrix / n_iters
table(is.na(correlation_matrices_list[[2]]))

# Filter non-dependent link
# cor_xi_test <- pmax(cor_xi_test, t(cor_xi_test)) # Max Value Symmetric
diag(cor_xi_test) <- 0
cor_xi_filtered <- cor_xi_test
cor_xi_filtered[pmax(cor_xi_test, t(cor_xi_test)) < cutoffs_Xi['a95']] <- 0
cor_xi_filtered[1:10,1:10]
table(cor_xi_filtered>0)


ActRegCut <- 0.455*(test_size)**-0.215
cor_xi_diff <- cor_xi_filtered - t(cor_xi_filtered)
cor_xi_diff[1:10,1:10]


# generating edgelist
xi_edge = melt(cor_xi_diff,value.name = 'diffCor')
xi_edge$Var1 <- as.character(xi_edge$Var1)
xi_edge$Var2 <- as.character(xi_edge$Var2)
xi_edge = xi_edge %>% arrange(-diffCor) %>% filter(abs(diffCor) > ActRegCut)
xi_edge_graph = xi_edge %>% arrange(-diffCor) %>% filter(diffCor > ActRegCut)%>%
  mutate("pairs" = paste0(Var1,"_",Var2))
xi_edge_graph <- xi_edge_graph %>%
  mutate('Validation' = case_when(
    pairs %in% STRINGedglist_sub$Gene.Gene & pairs %in% TFLinkTFTarget$TF.Target ~ "Both",
    pairs %in% STRINGedglist_sub$Gene.Gene ~ "STRING",
    pairs %in% TFLinkTFTarget$TF.Target ~ "TF Gene",
    TRUE ~ "None"
  )) %>%
  mutate(Validation = factor(Validation, levels = c('None','STRING','TF Gene','Both'), ordered = TRUE)) %>%
  arrange(desc(Validation))

head(xi_edge_graph,10)
table(xi_edge_graph$Validation)

# Get known TF target links
xi_edgeSub <- xi_edge %>% mutate("pairs" = paste0(Var1,"_",Var2)) %>%
  filter(pairs %in% TFLinkTFTarget$TF.Target)
xi_edgeSub

xi_edgeSub_graph <- xi_edgeSub %>%
  mutate(Var1 = ifelse(diffCor > 0, xi_edgeSub$Var1, stringr::str_split(pairs, "_", simplify = T)[, 2])) %>%
  mutate(Var2 = ifelse(diffCor > 0, xi_edgeSub$Var2, stringr::str_split(pairs, "_", simplify = T)[, 1]))
xi_edgeSub_graph

signs <- table(sign(xi_edgeSub$diffCor))
paste0("Correct TF identification = ", round(signs['1']/sum(signs)*100, 2),'%')

# plotting
networkObj = graph_from_edgelist(as.matrix(xi_edge_graph[,1:2]), directed = T)
networkObj <- set_edge_attr(networkObj, name = "Validation", value = xi_edge_graph$Validation)
edge_attr(networkObj)

set.seed(0)
layout <- layout_with_fr(networkObj)

p1 <- ggraph(networkObj, layout = layout) +
  geom_edge_link(aes(color=`Validation`, group = `Validation`),
                 arrow = arrow(length = unit(0.05, 'in'), type = 'closed'),
                 edge_width = 0.5,
                 end_cap = circle(0.04, 'in'), start_cap = circle(0.02, 'in'))+
  geom_node_point(size = 1.5,
                  colour = alpha('magenta4', 0.7),
                  shape = 16) +
  geom_node_text(aes(label = name),
                 colour = 'black',
                 size = 3,
                 repel = T, max.overlaps = 3) +
  scale_edge_colour_manual(values = c(alpha('black',0.2),'indianred1','gold','darkorchid1'))+
  theme_graph()
p1
png('figures/Network_GSE239592_AllGene.png', height = 6, width = 8, res=300, unit='in')
p1
dev.off()

# TF-Gene
networkObj = graph_from_edgelist(as.matrix(xi_edgeSub_graph[,1:2]), directed = T)
networkObj <- set_edge_attr(networkObj, name = "Correct Direction", value = ifelse(xi_edgeSub_graph$diffCor > 0, 'Yes', 'No'))
edge_attr(networkObj)

set.seed(0)
layout <- layout_with_fr(networkObj)

p1 <- ggraph(networkObj, layout = layout) +
  geom_edge_link(aes(color=`Correct Direction`),
                 arrow = arrow(length = unit(0.05, 'in'), type = 'open'),
                 #colour = alpha('darkgrey',0.6),
                 edge_width = 0.5,
                 end_cap = circle(0.1, 'in'), start_cap = circle(0.05, 'in'))+
  geom_node_point(size = 2,
                  colour = alpha('magenta4', 0.7),
                  shape = 16) +
  geom_node_text(aes(label = name),
                 colour = 'black',
                 size = 3,
                 repel = T, max.overlaps = 3) +
  scale_edge_colour_manual(values = c('No' = alpha('black',0.5),
                                      'Yes'='indianred3'))+
  theme_graph()
p1
png('figures/Network_GSE239592_AllGene_TFGene.png', height = 5, width = 5, res=300, unit='in')
p1
dev.off()

