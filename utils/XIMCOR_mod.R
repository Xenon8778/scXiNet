# Revised Chatterjee rank correlation xi_{n,M} (Lin & Han, 2023) — optimized.
#
#   xi_{n,M} = -2 + 6 * S / ((n+1) * (n*M + M*(M+1)/4)),
#   S = sum_{m=1..M} [ sum_{i<=n-m} min(r_(i), r_(i+m)) + sum_{i>n-m} r_(i) ]
# where r_(i) is the y-rank of the cell with the i-th smallest x
# (the m-th right nearest neighbour of the last m cells is the cell itself).
#
# Same outputs and orientation as the original functions; ties are broken at
# random, as before. Speed-ups:
#   * ranks are computed once per gene, not once per pair;
#   * for each x-gene, all p y-genes are scored at once as a matrix, using
#     min(a, b) = (a + b - |a - b|) / 2 with prefix sums (faster than pmin);
#   * the "tail" term depends only on M and the permutation, so it is precomputed
#     with one weighted column sum.
# XIM_Mat: O(n log n * p + n * M * p^2) arithmetic, no R-level pair loop.

.xim_norm <- function(n, M) (n + 1) * (n * M + M * (M + 1) / 4)

# Core: Y is an n x k matrix of y-ranks already ordered by x. Returns S per column.
.xim_S <- function(Y, M) {
  n <- nrow(Y)
  # tail: sum_{m=1..M} sum of last m entries = sum_i w_i * Y[i, ], w_i = #{m : i > n-m}
  w <- pmax(0, M - (n - seq_len(n)))
  S <- colSums(Y * w)
  # sum min(a, b) = (sum a + sum b - sum |a - b|) / 2  (abs/subtract beat pmin)
  tot <- colSums(Y)
  csum <- apply(Y, 2, cumsum)    # prefix sums: sums of head/tail blocks in O(1)
  for (m in seq_len(M)) {
    head <- csum[n - m, ]                          # sum Y[1:(n-m), ]
    tail <- tot - csum[m, ]                        # sum Y[(m+1):n, ]
    S <- S + (head + tail -
                colSums(abs(Y[(m + 1):n, , drop = FALSE] - Y[seq_len(n - m), , drop = FALSE]))) / 2
  }
  S
}

# Single pair
XIMcalculate <- function(xvec, yvec, M) {
  n  <- length(xvec)
  yr <- rank(yvec, ties.method = "random")[order(rank(xvec, ties.method = "random"))]
  -2 + 6 * .xim_S(matrix(yr), M) / .xim_norm(n, M)
}

# All pairs: out[i, j] = xi_{n,M}(x_i -> x_j), diagonal 0 (same as original).
# nCores > 1 splits the x-genes across a PSOCK cluster.
XIM_Mat <- function(xmat, M, verbose = TRUE, nCores = 1) {
  xmat <- as.matrix(xmat)
  n <- nrow(xmat); p <- ncol(xmat)
  R   <- apply(xmat, 2, rank, ties.method = "random")   # y-ranks, once per gene
  ord <- apply(xmat, 2, function(v) order(rank(v, ties.method = "random")))
  storage.mode(R) <- "double"
  nrm <- .xim_norm(n, M)

  row_fun <- function(i) {
    r <- -2 + 6 * .xim_S(R[ord[, i], , drop = FALSE], M) / nrm
    r[i] <- 0
    r
  }

  if (nCores > 1) {
    cl <- parallel::makeCluster(nCores); on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, c("R", "ord", "nrm", "M", ".xim_S"), envir = environment())
    chunks <- split(seq_len(p), cut(seq_len(p), nCores, labels = FALSE))
    out <- do.call(rbind, parallel::parLapply(cl, chunks, function(ii)
      t(vapply(ii, row_fun, numeric(p)))))
  } else {
    if (verbose) { pb <- txtProgressBar(min = 0, max = p, style = 3); on.exit(close(pb), add = TRUE) }
    out <- matrix(0, p, p)
    for (i in seq_len(p)) {
      out[i, ] <- row_fun(i)
      if (verbose) setTxtProgressBar(pb, i)
    }
  }
  dimnames(out) <- list(colnames(xmat), colnames(xmat))
  out
}

# Test statistic xi^{+-}_{n,M}
XIMstat <- function(xvec, yvec, M) {
  n  <- length(xvec)
  yr <- rank(yvec, ties.method = "random")[order(rank(xvec, ties.method = "random"))]
  S  <- .xim_S(cbind(yr, n + 1 - yr), M)
  -2 + 6 * max(S) / .xim_norm(n, M)
}

# Null simulation of xi^{+-}_{n,M}: B random permutations, scored in one matrix pass
XIMsim <- function(n, M, B) {
  Y <- replicate(B, sample.int(n))                          # n x B
  S <- pmax(.xim_S(Y, M), .xim_S(n + 1 - Y, M))
  -2 + 6 * S / .xim_norm(n, M)
}

# Simulation-based test from precomputed statistic and null draws
XIMtestT <- function(XIMstat, XIMsim, alpha = 0.05) {
  B <- length(XIMsim)
  (1 + sum(XIMstat <= XIMsim)) / (1 + B) <= alpha
}

# Simulation-based test directly (original referenced an undefined `n`)
XIMtest <- function(xvec, yvec, M, B, alpha = 0.05) {
  XIMtestT(XIMstat(xvec, yvec, M), XIMsim(length(xvec), M, B), alpha)
}
