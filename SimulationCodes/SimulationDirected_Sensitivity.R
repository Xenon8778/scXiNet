# ============================================================
# Sensitivity simulations for the directionality (asymmetric-association) cutoff.
# One-at-a-time sweeps around the baseline of SimulationDirected.R (Scenario 2):
#   noise SD, dropout rate, expression level, dispersion, count family, dropout stage.
# Output: Simulation_Sensitivity/<setting>/counts_<rep>_<nCells>.rds  (genes x cells)
#         Simulation_Sensitivity/<setting>/gt_GRN.csv                 (source, target, nonmonotonic)
#         Simulation_Sensitivity/settings.csv
# Restartable: existing files are skipped.
#
# Differences from SimulationDirected.R (deliberate; the BASELINE is regenerated with this
# simulator, so compare settings against 'baseline' here, not against the paper's files):
#  1. Vectorised (rpois / rnbinom) instead of per-entry rbinom(1, >=1000 trials, lambda/trials),
#     which approximates Poisson. Same model, different RNG stream.
#  2. The original calls set.seed(0) before choosing gene pairs and never reseeds, so every
#     replicate drew the SAME interaction-noise vectors. Here the stream is reseeded per replicate.
#  3. Options for count family ('poisson' / 'nb') and dropout stage ('before' = original,
#     'after' = independent dropout also on the interaction targets).
# ============================================================
library(dplyr)

out_dir        <- 'Simulation_Sensitivity'
nGrid          <- c(1000, 2000, 3750, 5000, 7500)   # test size = 0.4 n -> 400 ... 3000 (3750 -> 1500 = real-data nCells)
nReps          <- 5
n_genes        <- 500
n_interactions <- 250

simulate_scrna <- function(n_genes, n_cells, normal_sd, mean_expression, dispersion, dropout_rate,
                           count_dist = 'poisson', nb_size = 2, dropout_stage = 'before',
                           n_interactions, seed) {
  stopifnot(n_interactions <= n_genes / 2)
  set.seed(seed)

  rcount <- function(lambda) {          # vectorised count draw, keeps dim
    lambda <- pmax(lambda, 0)
    v <- if (count_dist == 'nb') rnbinom(length(lambda), size = nb_size, mu = lambda) else
      rpois(length(lambda), lambda)
    if (!is.null(dim(lambda))) dim(v) <- dim(lambda)
    v
  }
  drop <- function(m) {
    m[matrix(runif(length(m)) < dropout_rate, nrow(m), ncol(m))] <- 0L
    m
  }

  gene_means   <- rgamma(n_genes, shape = 1 / dispersion, scale = mean_expression * dispersion)
  size_factors <- rgamma(n_cells, shape = 1, scale = 1)
  counts <- rcount(outer(gene_means, size_factors))
  rownames(counts) <- paste0('g', seq_len(n_genes))
  colnames(counts) <- paste0('c', seq_len(n_cells))
  if (dropout_stage == 'before') counts <- drop(counts)

  # Gene pairs: fixed across seeds/settings (as in the original, set.seed(0))
  set.seed(0)
  gpairs <- matrix(sample(rownames(counts), size = 2 * n_interactions, replace = FALSE), ncol = 2)
  set.seed(seed + 1e5)                  # independent noise stream per replicate

  X <- t(counts)                        # cells x genes
  for (i in seq_len(n_interactions)) {
    s  <- X[, gpairs[i, 1]]
    mx <- max(s, 1)
    y <- switch(as.character(i %% 4),
                '1' = abs(s * 2 + rnorm(n_cells, sd = normal_sd)),                       # linear
                '2' = { y <- -0.02 * s^2 + s + 5; y[y < 0] <- 0; y + abs(rnorm(n_cells, sd = normal_sd)) },  # parabola
                '3' = { y <- exp(s / mx * 5) + rnorm(n_cells, sd = normal_sd); y[y < 0] <- 0; y },           # exponential
                '0' = sin(s / (mx / 15)) * 5 + 10 + rnorm(n_cells, sd = normal_sd))      # sinusoid
    X[, gpairs[i, 2]] <- rcount(y)
  }
  counts <- t(X)
  if (dropout_stage == 'after') counts <- drop(counts)

  gt <- data.frame(source = gpairs[, 1], target = gpairs[, 2])
  gt$nonmonotonic <- !((seq_len(nrow(gt)) %% 4) %in% c(1, 3))   # linear/exp = FALSE (null), parabola/sin = TRUE
  list(counts = counts, gpairs = gt)
}

#---------------------
# Settings: baseline = SimulationDirected.R Scenario 2
#---------------------
base <- list(normal_sd = 0.5, dropout_rate = 0.3, mean_expression = 15, dispersion = 0.5,
             count_dist = 'poisson', dropout_stage = 'before')
vary <- function(par, values)
  setNames(lapply(values, function(v) modifyList(base, setNames(list(v), par))),
           paste0(par, '_', values))
settings <- c(list(baseline = base),
              vary('normal_sd',        c(0.2, 1, 1.5)),
              vary('dropout_rate',     c(0, 0.1, 0.2, 0.5)),
              vary('mean_expression',  c(10, 20, 30)),     # low = sparse counts
              vary('dispersion',       c(0.2, 1, 1.5)),       # spread of gene means
              list(count_nb      = modifyList(base, list(count_dist = 'nb'))))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(do.call(rbind, lapply(names(settings), function(s)
  data.frame(setting = s, as.data.frame(settings[[s]])))), file.path(out_dir, 'settings.csv'),
  row.names = FALSE)

for (s in names(settings)) {
  message('Setting: ', s)
  d <- file.path(out_dir, s); dir.create(d, showWarnings = FALSE)
  for (n in nGrid) for (i in seq_len(nReps)) {
    f <- file.path(d, sprintf('counts_%d_%d.rds', i, n))
    if (file.exists(f)) next
    sim <- do.call(simulate_scrna, c(list(n_genes = n_genes, n_cells = n,
                                          n_interactions = n_interactions, seed = i), settings[[s]]))
    saveRDS(sim$counts, f)
    if (!file.exists(file.path(d, 'gt_GRN.csv'))) write.csv(sim$gpairs, file.path(d, 'gt_GRN.csv'))
  }
}
message('Done.')
