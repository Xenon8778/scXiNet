# ============================================================
# LearnCutoff(): Algorithm 2, as in the manuscript. Nothing else.
#   repeat B times:
#     subsample n cells (without replacement if n < n', bootstrap otherwise)
#     delta = |xi(g1 -> g2) - xi(g2 -> g1)| for every known-null pair (N0)
#     C_b   = (1 - alpha) quantile of the null deltas
#   cutoff = median(C_1, ..., C_B)
# Requires known null pairs (`null` = TRUE): simulation ground truth or a curated reference.
# ============================================================

LearnCutoff <- function(X, pairs, null, n,
                        B      = 50,
                        alpha  = 0.05,
                        xi_alg = 'Xi',   # 'Xi_KNN' (XIM_Mat) | 'Xi' (xicor_mpar)
                        M      = 10,
                        seed   = 1,
                        cl     = NULL) {
  #' X     : cell x gene matrix (reference dataset, n' cells)
  #' pairs : data.frame(g1, g2) of known dependent pairs
  #' null  : logical, one per row of pairs; TRUE = no active regulation (N0)
  #' n     : sample size the cutoff is for (size of the dataset it will be applied to)
  X <- as.matrix(X)
  stopifnot(is.logical(null), length(null) == nrow(pairs), any(null))
  P <- cbind(as.character(pairs$g1[null]), as.character(pairs$g2[null]))
  stopifnot(all(P %in% colnames(X)))
  X <- X[, unique(c(P)), drop = FALSE]

  one <- function(b) {
    Xr <- get('.LC_X', envir = globalenv())
    set.seed(seed + b)
    S <- Xr[sample.int(nrow(Xr), n, replace = n >= nrow(Xr)), , drop = FALSE]
    d <- apply(P, 1, function(p) {
      m <- as.matrix(if (xi_alg == 'Xi_KNN') XIM_Mat(S[, p], M, verbose = FALSE) else
        xicor_mpar(S[, p], nCores = 1))
      abs(m[2, 1] - m[1, 2])
    })
    d[!is.finite(d)] <- 0              # constant gene -> no asymmetry, as in CalcDelta
    unname(quantile(d, 1 - alpha))
  }
  env <- new.env(parent = globalenv())  # ship only small settings with each task
  env$P <- P; env$n <- n; env$alpha <- alpha; env$xi_alg <- xi_alg; env$M <- M; env$seed <- seed
  environment(one) <- env

  put <- function(d) { assign('.LC_X', d, envir = globalenv()); NULL }
  environment(put) <- globalenv()
  if (is.null(cl)) {
    put(X); on.exit(rm('.LC_X', envir = globalenv()), add = TRUE)
  } else {
    parallel::clusterCall(cl, put, X)
    on.exit(parallel::clusterEvalQ(cl, rm('.LC_X', envir = globalenv())), add = TRUE)
  }
  Cb <- unlist(pbapply::pblapply(seq_len(B), one, cl = cl))

  list(cutoff = median(Cb), per_iter = Cb,
       settings = list(n = n, nRef = nrow(X), B = B, alpha = alpha, xi_alg = xi_alg, M = M,
                       nNull = nrow(P)))
}
