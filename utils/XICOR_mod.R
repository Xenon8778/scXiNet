library(foreach)
library(doParallel)

# Function to calculate the Xi correlation matrix with parallel support
xicor_mpar <- function(x, factor = FALSE, nCores = 1) {
  # x: Input matrix or data frame.
  # factor: A logical value. If TRUE, converts factor variables to integers before calculation.
  # nCores: The number of CPU cores to use for parallel processing. Defaults to 1 (no parallelism).

  # --- Input validation and data preparation ---

  # If 'factor' is TRUE, convert non-numeric columns (if they exist) to integers.
  # This is useful for categorical data, treating factor levels as ordered numbers.
  if (factor) {
    if (!is.numeric(x)) x <- apply(x, 2, function(col) as.numeric(factor(col)))
  }

  # Ensure input is a numeric matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !(is.numeric(x) || is.logical(x))) {
    stop("Input 'x' must be a numeric matrix or data frame.")
  }
  stopifnot(is.atomic(x))
  ncx <- ncol(x)
  if (ncx == 0) stop("'x' has no columns.")
  n <- nrow(x)
  if (n == 0) stop("'x' has no rows.")

  # Precompute ranks and orderings for all columns of x
  ranks_x <- apply(x, 2, rank, ties.method = "random")
  orderings_x <- apply(ranks_x, 2, order)

  # Precompute ranks and orderings for all columns of y (which is also x here)
  ranks_y <- apply(x, 2, rank, ties.method = "max")
  ranks_y_neg <- apply(-x, 2, rank, ties.method = "max")

  calculate_xicor_pair <- function(i, j, x, ranks_x, orderings_x, ranks_y, ranks_y_neg, n) {
    ord_i <- orderings_x[, i]
    fr <- ranks_y[, j] / n
    gr <- ranks_y_neg[, j] / n
    fr_ordered <- fr[ord_i]
    A1 <- sum(abs(fr_ordered[1:(n - 1)] - fr_ordered[2:n])) / (2 * n)
    CU <- mean(gr * (1 - gr))
    return(1 - A1 / CU)
  }

  if (nCores == 1) {
    r <- matrix(0, nrow = ncx, ncol = ncx)
    for (i in seq_len(ncx)) {
      for (j in seq_len(ncx)) {
        r[i, j] <- calculate_xicor_pair(i, j, x, ranks_x, orderings_x, ranks_y, ranks_y_neg, n)
      }
    }
  } else {
    # Initialize multicore
    cl <- makeCluster(nCores)
    registerDoParallel(cl)

    indices <- expand.grid(i = seq_len(ncx), j = seq_len(ncx))
    results <- parApply(cl, indices, 1, function(row_indices) {
      i <- row_indices["i"]
      j <- row_indices["j"]
      calculate_xicor_pair(i, j, x, ranks_x, orderings_x, ranks_y, ranks_y_neg, n)
    })

    r <- matrix(results, nrow = ncx, ncol = ncx, byrow = FALSE)

    stopCluster(cl)
  }

  rownames(r) <- colnames(x)
  colnames(r) <- colnames(x)
  return(r)
}
