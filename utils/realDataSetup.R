# ============================================================
# Shared setup for real-data cutoff selection and network construction.
# Sourced by RealDataCutoffSelection_ver3.R and RealDataNetwork_GSM8331607_ver3.R
# so the gene panel, ground truth and Xi computation are identical in both.
# ============================================================
library(Seurat)
library(igraph)
library(dplyr)
source('code/XICOR_mod.R')
source('code/cutoffCalcScripts.R')   # updated CalcCut (FPR + FDR cutoffs, NA-aware val_mat)

nCells     <- 1500
alpha      <- 0.05
refDataset <- 'GSE239592'

# CalcCut names outputs from alpha (0.05 -> a95, FDR95, ...)
cut_sfx      <- as.character(round((1 - alpha) * 100))
cut_types    <- c(FPR = paste0('a', cut_sfx), FDR = paste0('FDR', cut_sfx))
feasible_key <- paste0('FDR', cut_sfx, '_feasible')
minach_key   <- paste0('FDR', cut_sfx, '_minAchievable')

# Ground-truth parameters (STRING analogue of direct / no-trek / indirect)
hi_conf          <- 900   # combined_score > hi_conf           -> direct edge (TRUE)
lo_conf          <- NULL  # NULL = off (2-hop exclusion only). If set, STRING links with
# score > lo_conf are excluded from the null (0 = any link).
null_min_dist    <- 3     # hi-conf shortest path >= this / none -> null eligible; dist 2 -> NA
drop_unannotated <- FALSE # TRUE: genes with no hi-conf STRING edge -> all pairs NA
# ---- Gene panel (one definition for both scripts) ----
drop_ribo_mito <- FALSE    # RPL/RPS/MT-/unannotated ENSG; set FALSE to reproduce the 1000-gene panel
set.seed(0)
commonGenes <- read.csv('data/TcellCommonGenes.csv', row.names = 1)$x
if (drop_ribo_mito) commonGenes <- commonGenes[!grepl('^(RPL|RPS|MT-|ENSG)', commonGenes)]
write.csv(commonGenes, 'data/TcellCommonGenes_panel.csv', row.names = FALSE)
message('Gene panel: ', length(commonGenes), ' genes')

# ---- STRING ----
STRINGedglist <- read.table('data/STRING/9606.protein.links.v12.0.txt.gz',
                            sep = ' ', header = TRUE)
ProteinInfo <- read.csv('data/STRING/9606.protein.info.v12.0.txt.gz',
                        sep = '\t', header = TRUE)
ProteinMap <- setNames(ProteinInfo$preferred_name, ProteinInfo$X.string_protein_id)
STRINGedglist$protein1 <- ProteinMap[STRINGedglist$protein1]
STRINGedglist$protein2 <- ProteinMap[STRINGedglist$protein2]
STRINGedglist <- STRINGedglist[!is.na(STRINGedglist$protein1) & !is.na(STRINGedglist$protein2), ]

# High-confidence graph on the FULL proteome (indirect paths via out-of-panel genes count)
g_hi <- STRINGedglist[STRINGedglist$combined_score > hi_conf, 1:2] |>
  graph_from_data_frame(directed = FALSE) |>
  simplify()

# Three-class ground truth (mirrors LoadGT in CutoffSelect_Nonlinear_ver3.R)
#   TRUE  = direct hi-conf STRING edge
#   FALSE = no hi-conf path shorter than null_min_dist hops (incl. genes with no hi-conf edge)
#   NA    = indirect (2 hops via hi-conf STRING); optionally lo_conf evidence
LoadGT_STRING <- function(genes) {
  p   <- length(genes)
  dn  <- list(genes, genes)
  inG <- genes[genes %in% V(g_hi)$name]

  D <- matrix(Inf, p, p, dimnames = dn)
  if (length(inG)) D[inG, inG] <- distances(g_hi, v = inG, to = inG)

  ev <- matrix(FALSE, p, p, dimnames = dn)
  if (!is.null(lo_conf)) {
    lo <- STRINGedglist[STRINGedglist$combined_score > lo_conf &
                          STRINGedglist$protein1 %in% genes &
                          STRINGedglist$protein2 %in% genes, 1:2]
    ev[as.matrix(lo)] <- TRUE
    ev <- ev | t(ev)
  }

  cls <- matrix(NA, p, p, dimnames = dn)
  cls[D >= null_min_dist & !ev] <- FALSE
  cls[D == 1]                   <- TRUE
  if (drop_unannotated) {
    un <- !(genes %in% inG)
    cls[un, ] <- NA; cls[, un] <- NA
  }
  diag(cls) <- NA

  ut <- upper.tri(cls)
  message(sprintf('  pos=%d  null=%d  excluded=%d  (unannotated genes=%d)',
                  sum(cls[ut], na.rm = TRUE), sum(!cls[ut], na.rm = TRUE),
                  sum(is.na(cls[ut])), sum(!(genes %in% inG))))
  cls
}

# ---- Cells and Xi ----
# Column order is left as in the Seurat object (as in the original scripts): reordering
# columns changes the random tie-break stream in xicor_mpar. nCells = NULL -> all cells.
sampleCells <- function(soFilt, nCells, seed = 0) {
  dfCount <- as.data.frame(t(as.matrix(GetAssayData(soFilt, layer = 'data'))))
  if (is.null(nCells)) return(dfCount)
  set.seed(seed)
  dfCount[sample(rownames(dfCount), nCells), ]
}

# Directed |xi| matrix, entry [j, k] = xi(j -> k).
# n_tiebreak > 1 averages over random tie-breaks. CalcCut / Algorithm 2 calibrate on
# single draws, and averaging shrinks the null spread, so n_tiebreak > 1 makes the
# transferred cutoffs conservative (not anti-conservative).
computeXi <- function(dfCount, symmetric = TRUE, n_tiebreak = 1, nCores = 8) {
  m <- Reduce(`+`, lapply(seq_len(n_tiebreak), function(i) {
    set.seed(i); abs(xicor_mpar(dfCount, nCores = nCores))
  })) / n_tiebreak
  m[is.na(m)] <- 0
  if (symmetric) m <- pmax(m, t(m))
  diag(m) <- 0
  m
}
