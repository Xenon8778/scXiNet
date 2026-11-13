library(foreach)
library(doParallel)

xicor_mpar <- function(x, factor = FALSE, nCores = 1, seed = NULL) {
  if (is.data.frame(x)) {
    if (factor) x[] <- lapply(x, function(c) if (is.numeric(c)) c else as.integer(factor(c)))
    x <- as.matrix(x)
  }
  stopifnot(is.matrix(x), is.numeric(x) || is.logical(x), ncol(x) > 0, nrow(x) > 1)
  n <- nrow(x); p <- ncol(x)
  if (!is.null(seed)) set.seed(seed)

  ord <- apply(x, 2, function(v) order(rank(v, ties.method = "random")))  # n x p
  R   <- apply(x, 2, rank, ties.method = "max")                            # r_i (integer-valued)
  L   <- apply(-x, 2, rank, ties.method = "max")                           # l_i
  den <- 2 * colSums(L * (n - L)) / n          # xi = 1 - sum|dr| / den[j]

  row_fun <- function(i) {
    Ri <- R[ord[, i], , drop = FALSE]
    1 - colSums(abs(Ri[-1, , drop = FALSE] - Ri[-n, , drop = FALSE])) / den
  }

  r <- if (nCores > 1) {
    cl <- parallel::makeCluster(nCores); on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl, c("R", "ord", "den", "n"), envir = environment())
    chunks <- split(seq_len(p), cut(seq_len(p), nCores, labels = FALSE))
    do.call(rbind, parallel::parLapply(cl, chunks, function(ii) t(vapply(ii, row_fun, numeric(p)))))
  } else t(vapply(seq_len(p), row_fun, numeric(p)))

  dimnames(r) <- list(colnames(x), colnames(x))
  r   # r[i, j] = xi(i -> j)
}
