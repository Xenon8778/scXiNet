library(dplyr)
library(XICOR)
library(scTenifoldNet)
library(energy)
library(caret)
library(pROC)
library(minet)        # Mutual information via build.mim

# Load modified Xi Correlation functions
source('code/XIMCOR_mod.R')
source('code/XICOR_mod.R')

# ============================================================
# range01: min-max normalisation to [0, 1]
# ============================================================
range01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(x * 0)   # constant input → all zeros
  (x - rng[1]) / diff(rng)
}

# ============================================================
# computeMI
# Builds a pairwise mutual-information matrix using minet's
# build.mim (raw MIM — no ARACNE pruning).
# Returns a symmetric, non-negative p x p matrix.
#
# Arguments:
#   mat       - n x p matrix (cells x genes)
#   estimator - MI estimator passed to build.mim (default: "mi.empirical")
#   disc      - Discretisation method             (default: "equalfreq")
#   nbins     - Number of bins; defaults to floor(sqrt(n))
# ============================================================
computeMI <- function(mat,
                      estimator = "mi.empirical",
                      disc      = "equalfreq",
                      nbins     = NULL) {
  if (is.null(nbins)) nbins <- max(2L, floor(sqrt(nrow(mat))))
  mim <- build.mim(mat, estimator = estimator, disc = disc, nbins = nbins)
  mim[mim < 0] <- 0
  diag(mim)    <- 0
  mim
}

# ============================================================
# computeDcor
# Exact distance-correlation matrix using the O(n log n) algorithm of
# Huo & Szekely (2016), via energy::dcor2d (V-statistic).
# dcor2d returns dCor^2, so sqrt() gives the same value as energy::dcor().
# Computes each pair once (upper triangle) and mirrors the result.
# Returns a symmetric p x p matrix (diagonal = 0).
#
# Arguments:
#   mat    - n x p matrix (cells x genes)
#   nCores - cores for parallel computation over genes (fork-based;
#            use 1 on Windows)
# ============================================================
computeDcor <- function(mat, nCores = 1, cl = NULL) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  p   <- ncol(mat)

  rowFun <- function(j) {
    x <- mat[, j]
    vapply((j + 1L):p, function(k) energy::dcor2d(x, mat[, k], type = "V"), numeric(1))
  }

  if (nCores > 1 || !is.null(cl)) {
    own <- is.null(cl)
    if (own) cl <- parallel::makeCluster(nCores)
    on.exit(if (own) parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, c("mat", "p"), envir = environment())
    rows <- parallel::parLapplyLB(cl, seq_len(p - 1L), rowFun)  # LB: rows shrink with j
  } else {
    rows <- lapply(seq_len(p - 1L), rowFun)
  }

  res <- matrix(0.0, p, p, dimnames = list(colnames(mat), colnames(mat)))
  for (j in seq_len(p - 1L)) res[j, (j + 1L):p] <- rows[[j]]
  res <- res + t(res)
  sqrt(pmax(res, 0))   # zero-variance genes give NaN -> set to 0 in CalcCut
}

# ============================================================
# fdrCutoff
# Smallest score cutoff whose call set (score > cutoff) has realised
# FDP = FP / (FP + TP) <= alpha, i.e. the largest rejection set with
# FDP <= alpha. Evaluated only at tie-block ends so ties are never split.
# Returns Inf if no non-empty call set reaches FDP <= alpha.
# ============================================================
fdrCutoff <- function(s, lab, alpha = 0.05) {
  o    <- order(s, decreasing = TRUE)
  s    <- s[o]; lab <- lab[o]
  tp   <- cumsum(lab); fp <- cumsum(!lab)
  fdp  <- fp / (tp + fp)
  ends <- !duplicated(s, fromLast = TRUE)
  ok   <- which(ends & fdp <= alpha)
  if (!length(ok)) return(Inf)
  k <- max(ok)
  if (k == length(s)) return(-Inf)       # every pair called
  (s[k] + s[k + 1]) / 2                  # strictly between called / uncalled
}

# FDP at each grid cutoff (score > C called); FDP := 0 when nothing called
fdpCurve <- function(s, lab, grid) {
  fp  <- cumsum(!lab[order(s, decreasing = TRUE)])
  R   <- length(s) - findInterval(grid, sort(s))     # n called at each C
  out <- numeric(length(grid))
  out[R > 0] <- fp[R[R > 0]] / R[R > 0]
  list(fdp = out, R = R)
}

# ============================================================
# CalcCut  (Algorithm 1)
# Estimates dependency-score cutoffs from repeated subsampling of a
# training set with known ground truth.
#
#   a<..>   : (1-alpha) quantile of null-pair scores  -> controls FPR
#   FDR<..> : cutoff giving realised FDP <= alpha     -> controls FDR
#
# Subsampling follows Algorithm 1: without replacement if nsamp < n',
# bootstrap (with replacement) otherwise.
#
# Pairs: unique pairs j < k when `symmetric` (scores symmetrised by pmax
# for asymmetric methods: Xi, pcNet); all off-diagonal pairs otherwise.
#
# Arguments:
#   x         - n x p data matrix (cells x genes) — training data.
#   alpha     - Target FPR / FDR level (default: 0.05).
#   niter     - Number of subsampling iterations (default: 100).
#   nsamp     - Subsample size per iteration (default: 1500).
#   algorithm - One of: 'Xi', 'Xi_KNN', 'Spearman', 'Pearson',
#               'pcNet', 'dcor', 'MI'.
#   seed      - RNG seed for stochastic algorithms (default: 0).
#   symmetric - Symmetrise asymmetric scores / evaluate j<k (default: TRUE).
#   nCores    - Cores for parallel Xi computation (default: 1).
#   val_mat   - Logical p x p ground-truth matrix (TRUE = true edge).
#
# Returns:
#   Named numeric vector, e.g. for alpha = 0.05:
#   c(a95 = <FPR cutoff>, FDR95 = <FDR cutoff>,
#     FDR95_feasible = <fraction of iterations where FDR <= alpha was reachable>)
# ============================================================

CalcCut <- function(x, alpha = 0.05, niter = 100, nsamp = 1500, algorithm = 'Xi',
                    seed = 0, symmetric = TRUE, nCores = 1, val_mat, ngrid = 2000,
                    M = 10) {   # M: right nearest neighbours for Xi_KNN (default keeps old behaviour)

  algorithm <- match.arg(algorithm, c('Xi','Xi_KNN','Spearman','Pearson','pcNet','dcor','MI'))
  val_mat <- unname(as.matrix(val_mat)); storage.mode(val_mat) <- 'logical'
  if (symmetric && !isSymmetric(val_mat)) stop("val_mat not symmetric.")

  mask   <- (if (symmetric) upper.tri(val_mat) else row(val_mat) != col(val_mat)) & !is.na(val_mat)
  labels <- val_mat[mask]
  if (!any(labels) || all(labels)) stop("val_mat needs both true and null pairs.")

  replace <- nsamp >= nrow(x)
  S       <- matrix(NA_real_, niter, length(labels))
  cut_vec <- numeric(niter)

  pb <- txtProgressBar(min = 0, max = niter, style = 3, width = 30); on.exit(close(pb))
  for (i in seq_len(niter)) {
    set.seed(i)
    df  <- x[sample.int(nrow(x), nsamp, replace = replace), , drop = FALSE]
    res <- switch(algorithm,
                  Xi       = { set.seed(seed); xicor_mpar(df, nCores = nCores) },
                  Xi_KNN   = { set.seed(seed); XIM_Mat(df, M, nCores = nCores, verbose = FALSE) },
                  Spearman = cor(df, method = 'spearman'),
                  Pearson  = cor(df, method = 'pearson'),
                  pcNet    = as.matrix(pcNet(t(df), verbose = FALSE)),
                  dcor     = computeDcor(df, nCores = nCores),
                  MI       = computeMI(df))
    res <- abs(unname(as.matrix(res)))
    if (symmetric) res <- pmax(res, t(res))
    s <- res[mask]; s[is.na(s)] <- 0

    S[i, ]     <- s
    cut_vec[i] <- quantile(s[!labels], 1 - alpha, names = FALSE)   # FPR (unchanged)
    setTxtProgressBar(pb, i)
  }

  # ---- FDR: average FDP curve over iterations, then threshold ----
  grid  <- sort(unique(quantile(S, seq(0, 1, length.out = ngrid), names = FALSE)))
  curvs <- lapply(seq_len(niter), function(i) fdpCurve(S[i, ], labels, grid))
  FDRb  <- rowMeans(sapply(curvs, `[[`, 'fdp'))
  Rb    <- rowMeans(sapply(curvs, `[[`, 'R'))
  ok    <- FDRb <= alpha & Rb >= 1
  fdr_cut <- if (any(ok)) min(grid[ok]) else Inf          # smallest C = largest call set

  sfx <- as.character(round((1 - alpha) * 100))
  c(setNames(median(cut_vec), paste0('a', sfx)),
    setNames(fdr_cut,         paste0('FDR', sfx)),
    setNames(as.numeric(any(ok)), paste0('FDR', sfx, '_feasible')),
    setNames(min(FDRb[Rb >= 1]),  paste0('FDR', sfx, '_minAchievable')))
}
