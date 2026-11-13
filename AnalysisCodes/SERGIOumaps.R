# Load necessary libraries
library(umap)
library(stats)
library(ggplot2)
library(cowplot)
set.seed(42) # For reproducibility


dirList = list.dirs("SERGIO/data_sets/")[-1]
dirList <- dirList[lengths(sapply(dirList,strsplit,"_"))>7]
dirList

PCAplots <- list()

for (i in 3:length(dirList)-1){
  df = read.csv(paste0(dirList[i],'/simulated_noNoise_T_1.csv'), header = T, row.names = 1) # Test load

  # 1. Perform PCA
  pca_result <- prcomp(t(df), scale. = TRUE) # scale. = TRUE for scaling the data
  pca_scores <- pca_result$x

  df_pca <- data.frame(pca_scores[,1:2])
  colnames(df_pca) <- c('pca1', 'pca2')

  # Plot the PCA results
  p1 <- ggplot(df_pca, aes(x = pca1, y = pca2)) +
    geom_point(size = 1) +
    labs(title = NULL,
         x = "PC 1",
         y = "PC 2") +
    theme_classic(base_family = 'arial')
  p1

  PCAplots[[i-1]] <- p1
  fname <- paste0('figures/SERGIO/PCA_',tail(unlist(strsplit(dirList[i],'_')),1),'.png')
  png(fname, width = 3, height = 3, res = 300, units = 'in')
  print(p1)
  dev.off()
}

p2 <- plot_grid(plotlist = PCAplots, nrow = 1)
png('figures/SERGIO/PCA_All.png', width = 9, height = 3, res = 300, units = 'in')
p2
dev.off()


# Sergio Topologies ####
library(igraph)
library(ggraph)
library(dplyr)
library(tidyr)
library(cowplot)

dir.create('figures/SERGIO', showWarnings = FALSE, recursive = TRUE)

read_sergio_grn <- function(dir) {
  f <- list.files(dir, pattern = '^Interaction_cID.*\\.txt$', full.names = TRUE)[1]
  if (is.na(f)) stop('No Interaction_cID file in ', dir)
  ln <- readLines(f)
  ln <- ln[nzchar(trimws(ln))]
  do.call(rbind, lapply(ln, function(l) {
    v  <- as.numeric(strsplit(l, ',')[[1]])
    nr <- as.integer(v[2])
    data.frame(from = as.character(as.integer(v[3:(2 + nr)])),
               to   = as.character(as.integer(v[1])),
               K    = v[(3 + nr):(2 + 2 * nr)])
  }))
}

read_num_table <- function(f) {             # comma / tab / space delimited
  ln <- trimws(gsub('\r', '', readLines(f, warn = FALSE)))
  ln <- ln[nzchar(ln)]
  rows <- lapply(ln, function(l)
    suppressWarnings(as.numeric(strsplit(l, '[,\t ]+')[[1]])))
  rows <- rows[!sapply(rows, function(r) all(is.na(r)))]
  do.call(rbind, rows)
}

dsDirs  <- dirList[3:length(dirList)-1]       # same indexing as PCA loop
dsNames <- sapply(dsDirs, function(d) tail(unlist(strsplit(d, '_')), 1))

grns <- lapply(dsDirs, function(d) {
  gt <- read_sergio_grn(d)
  gt <- gt[gt$from != gt$to, ]
  gt[!duplicated(gt[, c('from', 'to')]), ]
})
names(grns) <- dsNames

# ---- Verify shared topology ----
key  <- function(x) sort(paste(x$from, x$to, sign(x$K)))
same <- sapply(grns[-1], function(x) identical(key(x), key(grns[[1]])))
cat('Topology + edge signs identical to', dsNames[1], ':\n'); print(same)

# ---- Shared GRN ----
gt <- grns[[1]]
gt$Sign <- ifelse(gt$K > 0, 'Activation', 'Repression')
g  <- graph_from_data_frame(gt, directed = TRUE)

indeg  <- degree(g, mode = 'in')
outdeg <- degree(g, mode = 'out')
V(g)$outdeg <- outdeg
V(g)$role   <- ifelse(outdeg > 0, 'Regulator', 'Target')
V(g)$MR     <- indeg == 0
V(g)$lab    <- paste0('g', V(g)$name)

depth <- NA
if (is_dag(g)) {
  lvl <- setNames(numeric(vcount(g)), V(g)$name)
  for (v in as_ids(topo_sort(g))) {
    nb <- names(neighbors(g, v, mode = 'in'))
    if (length(nb)) lvl[v] <- max(lvl[nb]) + 1
  }
  depth <- max(lvl)
}

ug   <- as_undirected(g, mode = 'collapse')
comm <- cluster_louvain(ug)
pl   <- fit_power_law(degree(g, mode = 'all'), implementation = 'plfit',
                      p.value = TRUE, p.precision = 0.01)
hubThr <- mean(outdeg) + 2*sd(outdeg)

topoTab <- data.frame(
  Datasets         = paste(dsNames, collapse = '/'),
  Genes            = vcount(g),
  Edges            = ecount(g),
  Density          = edge_density(g),
  MeanDegree       = mean(degree(g, mode = 'all')),
  MaxOutDegree     = max(outdeg),
  MaxInDegree      = max(indeg),
  Regulators       = sum(outdeg > 0),
  MasterRegulators = sum(indeg == 0 & outdeg > 0),
  Targets          = sum(outdeg == 0),
  Hubs_outdeg2SD   = sum(outdeg > hubThr),
  IsDAG            = is_dag(g),
  RegDepth         = depth,
  WeakComponents   = igraph::components(g, mode = 'weak')$no,
  Transitivity     = transitivity(ug, type = 'global'),
  Modularity       = modularity(comm),
  Communities      = length(comm),
  MeanPathLength   = mean_distance(g, directed = TRUE),
  PowerLaw_alpha   = pl$alpha,
  PowerLaw_KS_p    = pl$KS.p,
  Activations      = sum(gt$K > 0),
  Repressions      = sum(gt$K < 0)
)
write.csv(topoTab, 'figures/SERGIO/SERGIO_GT_topology.csv', row.names = FALSE)
print(t(topoTab))


# fix layout once and reuse it for both plots
set.seed(42); lay <- layout_with_kk(g)

p_grn <- ggraph(g, layout = lay) +
  geom_edge_link(aes(edge_colour = Sign),
                 arrow = arrow(length = unit(1.2,'mm'), type = 'closed'),
                 end_cap = circle(1.5,'mm'),
                 edge_width = 0.25, edge_alpha = 0.6) +
  scale_edge_colour_manual(values = c(Activation = 'grey40',
                                      Repression = '#E53935'),
                           name = 'Edge') +
  geom_node_point(aes(size = outdeg, colour = role)) +
  geom_node_label(aes(label = lab, filter = MR), size = 3, repel = TRUE,
                  label.padding = unit(0.12, 'lines'), fill = 'white') +
  scale_size_continuous(range = c(0.8, 4), name = 'Out-degree') +
  scale_colour_manual(values = c('Regulator' = '#4A148C',
                                 'Target'    = '#F48FB1'),
                      name = 'Role') +
  labs(title = sprintf('%d genes, %d edges (%d activating, %d repressing)',
                          vcount(g), ecount(g), sum(gt$K > 0), sum(gt$K < 0))) +
  theme_graph(base_family = 'sans') +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'plain'),
        legend.text   = element_text(size = 9),)
ggsave('figures/SERGIO/GRN_shared.png', p_grn, width = 6, height = 5, dpi = 300, bg = 'white')

# 1-trek pairs from ground truth: non-adjacent, joined by chain (i->k->j) or fork (i<-k->j)
# ---- Three-class labelling, identical to load_ground_truth(symmetric=TRUE) ----
A  <- as_adjacency_matrix(g, sparse = FALSE)            # regulator -> target
Rm <- (is.finite(distances(g, mode = 'out'))) * 1       # i reaches j, incl. self (R = I | A | closure)
CA <- (t(Rm) %*% Rm) > 0                                # common ancestor or ancestral relation
direct <- (A + t(A)) > 0

cls <- matrix(-1L, nrow(A), ncol(A), dimnames = dimnames(A))
cls[!CA]   <- 0L
cls[direct] <- 1L
diag(cls)  <- -1L

ut <- upper.tri(cls)
cat(sprintf('pos=%d null=%d excluded=%d\n', sum(cls[ut] == 1), sum(cls[ut] == 0), sum(cls[ut] == -1)))

# ---- Excluded pairs (cls == -1), typed ----
excl <- (cls == -1L) & ut
idx  <- which(excl, arr.ind = TRUE)
anc  <- (Rm > 0) | t(Rm > 0)                            # one reaches the other (any path length)
tk <- data.frame(x = lay[idx[, 1], 1], y = lay[idx[, 1], 2],
                 xend = lay[idx[, 2], 1], yend = lay[idx[, 2], 2],
                 Trek = ifelse(anc[idx], 'Ancestral (chain)', 'Common ancestor (fork)'))

p_trek <- ggraph(g, layout = lay) +
  geom_segment(data = tk, aes(x = x, y = y, xend = xend, yend = yend, colour = Trek),
               linewidth = 0.25, alpha = 0.5) +                 # excluded trek pairs (bottom)
  geom_edge_link(colour = 'black', edge_width = 0.3) +          # ground truth (on top)
  geom_node_point(aes(size = outdeg), colour = '#4A148C') +
  scale_size_continuous(range = c(0.6, 3.5), name = 'Out-degree') +
  scale_colour_manual(values = c('Ancestral (chain)' = '#FB8C00',
                                 'Common ancestor (fork)' = '#1E88E5'), name = 'Excluded pair') +
  labs(title = sprintf('%d excluded trek pairs vs %d true edges', nrow(tk), ecount(g))) +
  theme_graph(base_family = 'arial') +
  theme(plot.title = element_text(hjust = .5, face='plain', size = 12))
ggsave('figures/SERGIO/GRN_trek_pairs.png', p_trek, width = 6, height = 5, dpi = 300, bg = 'white')

# ---- Per-gene excluded partners vs in-degree ----
tp <- data.frame(indeg = degree(g, mode = 'in'), trekdeg = rowSums(cls == -1L))
p_hub <- ggplot(tp, aes(indeg, trekdeg)) +
  geom_jitter(width = .15, height = 0, size = 1.2, colour = '#4A148C', alpha = .6) +
  labs(x = 'In-degree', y = 'Number of excluded trek partners') +
  theme_classic(base_size = 11, base_family = 'arial')
p_hub
ggsave('figures/SERGIO/Trek_degree.png', p_hub, width = 3, height = 3, dpi = 300, bg = 'white')


# Key number for the text: share of the "null" set that is actually dependent
n_all <- choose(vcount(g), 2); n_null <- n_all - ecount(g)
c(pairs = n_all, true_edges = ecount(g), trek = nrow(tk), trek_frac_of_null = nrow(tk) / n_null)

# Degree distribution (CCDF, log-log)
dd <- data.frame(k = degree(g, mode = 'all')) %>%
  dplyr::count(k) %>% arrange(k) %>%
  mutate(ccdf = rev(cumsum(rev(n))) / sum(n))
p_deg <- ggplot(dd, aes(k, ccdf)) +
  geom_point(size = 1.5) +
  # scale_x_log10() +
  # scale_y_log10() +
  labs(x = 'Degree (k)', y = 'Fraction of genes with degree ≥ k')+
  theme_classic(base_family = 'arial') +
  theme(axis.text = element_text(size = 9))
ggsave('figures/SERGIO/Degree_shared.png', p_deg, width = 3, height = 3, dpi = 300, bg = 'white')

# ---- Differentiation graphs (bMat) + MR production rates (Regs) ----
bPlots <- list(); rPlots <- list(); nCT <- c()
for (k in seq_along(dsDirs)) {
  d <- dsDirs[k]; ds <- dsNames[k]

  bm <- read_num_table(list.files(d, '^bMat', full.names = TRUE)[1])
  ct <- paste0('C', seq_len(nrow(bm)))
  dimnames(bm) <- list(ct, ct)
  gb <- graph_from_adjacency_matrix(bm, mode = 'directed', weighted = TRUE, diag = FALSE)
  V(gb)$type <- ifelse(degree(gb, mode = 'in') == 0, 'Progenitor',
                       ifelse(degree(gb, mode = 'out') == 0, 'Terminal', 'Transitional'))

  bPlots[[ds]] <- ggraph(gb, layout = 'tree',
                         root = which(degree(gb, mode = 'in') == 0), mode = 'out') +
    geom_edge_link(aes(label = round(weight, 2)),
                   arrow = arrow(length = unit(2, 'mm'), type = 'closed'),
                   start_cap = circle(5, 'mm'), end_cap = circle(5, 'mm'),
                   angle_calc = 'along', label_dodge = unit(2.5, 'mm'),
                   label_size = 3, edge_colour = 'grey40') +
    geom_node_point(aes(colour = type), size = 12) +
    geom_node_text(aes(label = name), colour = 'white', size = 3.5, fontface = 'bold') +
    scale_colour_manual(values = c(Progenitor = '#4A148C', Transitional = '#F48FB1',
                                   Terminal = '#D81B60'), name = 'Cell type') +
    scale_x_continuous(expand = expansion(0.15)) +
    scale_y_continuous(expand = expansion(0.1)) +
    labs(title = paste0(nrow(bm), ' cell types')) +
    theme_graph(base_family = 'arial') +
    theme(plot.title = element_text(hjust = 0.5, face = 'plain', size = 14))

  rg <- read_num_table(list.files(d, '^Regs_cID', full.names = TRUE)[1])
  df <- data.frame(MR = paste0('g', as.integer(rg[, 1])), rg[, -1, drop = FALSE])
  colnames(df)[-1] <- paste0('C', seq_len(ncol(rg) - 1))
  long <- pivot_longer(df, -MR, names_to = 'CellType', values_to = 'Rate') %>%
    mutate(CellType = factor(CellType, levels = paste0('C', seq_len(ncol(rg) - 1))))

  rPlots[[ds]] <- ggplot(long, aes(CellType, MR, fill = Rate)) +
    geom_tile(colour = 'white') +
    geom_text(aes(label = round(Rate, 2)), size = 3) +
    scale_fill_gradient(low = 'white', high = 'mediumpurple1', name = 'Production\nrate') +
    labs(title = NULL, x = 'Cell type', y = 'Master regulator') +
    theme_classic(base_family = 'arial') +
    coord_fixed() +
    theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 11))
  nCT[ds] <- ncol(rg) - 1

  ggsave(paste0('figures/SERGIO/bMat_', ds, '.png'), bPlots[[ds]],
         width = 4, height = 4, dpi = 300, bg = 'white')
  ggsave(paste0('figures/SERGIO/MRrates_', ds, '.png'), rPlots[[ds]],
         width = 4, height = 3.5, dpi = 300, bg = 'white')
}

# ---- Combined figure ----
p_top <- plot_grid(p_grn, p_deg, nrow = 1, rel_widths = c(1.5, 1),
                   labels = c('a', 'b'), label_size = 18)
p_b   <- plot_grid(plotlist = bPlots, nrow = 1)
p_r   <- plot_grid(plotlist = rPlots, nrow = 1, align = 'h', axis = 'tb',
                   rel_widths = nCT + 3)

p_all <- plot_grid(p_top, p_b, p_r, ncol = 1,
                   rel_heights = c(1, 1, 1),
                   labels = c('', 'c', 'd'), label_size = 18)
ggsave('figures/SERGIO/SERGIO_topology_trajectory.png', p_all,
       width = 12, height = 12, dpi = 300, bg = 'white')
