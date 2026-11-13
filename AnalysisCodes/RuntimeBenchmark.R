# Benchmark runtime of 8 correlation/network measures on simulated scRNA-seq data
# Run interactively in RStudio

library(minet)
library(energy)
library(scTenifoldNet)
library(dplyr)
library(ggplot2)
library(foreach)
library(doParallel)

source('code/XICOR_mod.R')   # defines xicor_mpar()
source('code/XIMCOR_mod.R')  # defines XIM_Mat() (ChatterjeeKNN, Lin & Han 2023)

## ---- Simulate NB scRNA-seq data ----

simulate_scRNAseq_nb <- function(n_genes, n_cells,
                                 mean_expr = 15, dispersion = 0.5,
                                 dropout_rate = 0.3, seed = 1) {
  set.seed(seed)

  gene_means   <- rgamma(n_genes, shape = 1 / dispersion, scale = mean_expr * dispersion)
  size_factors <- rgamma(n_cells, shape = 1, scale = 1)
  mu_matrix    <- outer(size_factors, gene_means)

  counts <- matrix(
    rnbinom(n_cells * n_genes, size = 1 / dispersion, mu = as.vector(mu_matrix)),
    nrow = n_cells, ncol = n_genes
  )

  dropout_mask <- matrix(runif(n_cells * n_genes) < dropout_rate, nrow = n_cells)
  counts[dropout_mask] <- 0

  rownames(counts) <- paste0("Cell", seq_len(n_cells))
  colnames(counts) <- paste0("Gene", seq_len(n_genes))
  counts
}

## ---- Correlation/network measures ----

run_pearson  <- function(X) cor(X, method = "pearson")
run_spearman <- function(X) cor(X, method = "spearman")

run_minet <- function(X) {
  minet::build.mim(as.data.frame(X), estimator = "mi.empirical", disc = "equalwidth")
}

run_pcnet <- function(X, nComp = 3) {
  scTenifoldNet::pcNet(t(X), nComp = nComp, scaleScores = TRUE, symmetric = TRUE, q = 0)
}

# Same formula as computeDcor() in cutoffCalcScripts.R: O(n log n) V-statistic; dcor2d returns dCor^2
dcor_pair <- function(x, y) sqrt(energy::dcor2d(x, y, type = "V"))

run_chatterjee <- function(X, nCores = 1) {
  r <- xicor_mpar(X, nCores = nCores)
  pmax(r, t(r))  # symmetrize: max(xi(i->j), xi(j->i))
}

run_chatterjee_knn <- function(X, M = 5, nCores = 1) {
  r <- as.matrix(XIM_Mat(X, M, nCores = nCores, verbose = FALSE))
  pmax(r, t(r))  # same symmetrization as run_chatterjee
}

# dcor has no whole-matrix implementation in R, so at large gene counts
# we time a random subsample of pairs and extrapolate to the full pair count.
# With nCores > 1 the sampled pairs are split over a PSOCK cluster; cluster start-up and data
# export happen before the clock starts, so only the pairwise computation is timed.
time_pairwise_method <- function(X, pair_fun, max_pairs = 2000, seed = 42, nCores = 1) {
  p <- ncol(X)
  n_total_pairs <- p * (p - 1) / 2
  all_pairs <- combn(p, 2)

  if (n_total_pairs <= max_pairs) {
    idx <- seq_len(n_total_pairs)
    extrapolated <- FALSE
  } else {
    set.seed(seed)
    idx <- sample(seq_len(n_total_pairs), max_pairs)
    extrapolated <- TRUE
  }

  if (nCores > 1) {
    cl <- parallel::makeCluster(nCores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, c("X", "all_pairs", "pair_fun"), envir = environment())
    t0 <- Sys.time()
    parallel::parSapply(cl, idx, function(k) pair_fun(X[, all_pairs[1, k]], X[, all_pairs[2, k]]))
  } else {
    t0 <- Sys.time()
    for (k in idx) {
      pair_fun(X[, all_pairs[1, k]], X[, all_pairs[2, k]])
    }
  }
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (extrapolated) {
    runtime <- (elapsed / length(idx)) * n_total_pairs
    mode <- "extrapolated"
  } else {
    runtime <- elapsed
    mode <- "measured"
  }
  list(runtime_sec = runtime, mode = mode)
}

## ---- Benchmark driver ----

benchmark_grid <- function(gene_grid, cell_grid, max_pairwise_ops = 2000, nCores = 1, seed = 1) {
  results <- list()

  for (n_genes in gene_grid) {
    for (n_cells in cell_grid) {

      message(sprintf("%d genes x %d cells", n_genes, n_cells))
      X <- simulate_scRNAseq_nb(n_genes, n_cells, seed = seed)

      t_pearson    <- system.time(run_pearson(X))["elapsed"]
      t_spearman   <- system.time(run_spearman(X))["elapsed"]
      t_minet      <- system.time(run_minet(X))["elapsed"]
      t_pcnet      <- system.time(run_pcnet(X))["elapsed"]
      t_chatterjee <- system.time(run_chatterjee(X, nCores = nCores))["elapsed"]
      t_chat_knn   <- system.time(run_chatterjee_knn(X, nCores = nCores))["elapsed"]
      dcor_out     <- time_pairwise_method(X, dcor_pair, max_pairwise_ops, seed, nCores = nCores)

      results[[length(results) + 1]] <- data.frame(
        n_genes, n_cells,
        method = c("Pearson", "Spearman", "MI (minet)", "pcNet",
                   "Chatterjee's xi", "ChatterjeeKNN", "Distance correlation"),
        runtime_sec = c(t_pearson, t_spearman, t_minet, t_pcnet,
                        t_chatterjee, t_chat_knn, dcor_out$runtime_sec),
        mode = c("measured", "measured", "measured", "measured",
                 "measured", "measured", dcor_out$mode),
        nCores = nCores
      )
    }
  }

  do.call(rbind, results)
}

## ---- Run ----

gene_grid <- c(50, 100, 1000)
cell_grid <- c(100, 1000, 3000)

bench_results <- benchmark_grid(gene_grid, cell_grid, max_pairwise_ops = 2000, nCores = 1, seed = 1)
bench_results$method <- factor(bench_results$method)
bench_results$method <- recode_factor(bench_results$method,
                                      "MI (minet)" = "MI",
                                      "Chatterjee's xi" = 'Chatterjee',
                                      "Distance correlation" = 'Distance')
View(bench_results)

dir.create("Results", showWarnings = FALSE)
write.csv(bench_results, "Results/correlation_runtime_benchmark.csv", row.names = FALSE)

## ---- Run with 8 cores ----
# Multi-core where the implementation supports it: Chatterjee (xicor_mpar), ChatterjeeKNN (XIM_Mat),
# Distance (pairs split over cores). Pearson, Spearman, MI (minet) and pcNet have no multi-core
# option here and run single-threaded; they are marked '*' in the plot.
par_cores   <- 8
bench_par <- benchmark_grid(gene_grid, cell_grid, max_pairwise_ops = 2000, nCores = par_cores, seed = 1)
bench_par$method <- recode_factor(factor(bench_par$method),
                                  "MI (minet)" = "MI",
                                  "Chatterjee's xi" = 'Chatterjee',
                                  "Distance correlation" = 'Distance')
write.csv(bench_par, "Results/correlation_runtime_benchmark_8cores.csv", row.names = FALSE)

## ---- Plot ----
method_levels <- c("Pearson", "Spearman","pcNet", "MI", "Distance", "Chatterjee", "ChatterjeeKNN")
serial_only   <- c("Pearson", "Spearman", "pcNet", "MI")   # no multi-core option

plotRuntime <- function(csv_file, mark_serial = FALSE) {
  d <- read.csv(csv_file)
  d$n_genes <- factor(d$n_genes, levels = gene_grid)
  d$n_cells <- as.character(d$n_cells)
  d$method  <- factor(d$method, levels = method_levels)
  lab <- if (mark_serial) function(x) ifelse(x %in% serial_only, paste0(x, '*'), x) else waiver()

  ggplot(d, aes(method, runtime_sec, color = n_cells, group = n_cells, shape = n_cells)) +
    geom_point(size = 2) +
    geom_line() +
    facet_wrap(~ n_genes, scales = "free_y") +
    # scale_y_log10() +
    scale_x_discrete(labels = lab) +
    labs(x = NULL, y = "Runtime (sec)", color = "Number of cells", shape = "Number of cells") +
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
    scale_color_manual(values = c('pink2','magenta3','mediumpurple'))
}

runtime_plot <- plotRuntime("Results/correlation_runtime_benchmark.csv")
runtime_plot
ggsave("figures/correlation_runtime_benchmark.png", runtime_plot, width = 8, height = 3, dpi = 300)

# 8 cores; * = single-threaded (no multi-core option)
runtime_plot_par <- plotRuntime("Results/correlation_runtime_benchmark_8cores.csv", mark_serial = TRUE)
runtime_plot_par
ggsave("figures/correlation_runtime_benchmark_8cores.png", runtime_plot_par, width = 8, height = 3, dpi = 300)


## ---- Parallelisation advantage: parallelisable methods only, 1 vs 8 cores ----
par_methods <- c("Distance", "Chatterjee", "ChatterjeeKNN")

d_par <- rbind(read.csv("Results/correlation_runtime_benchmark.csv"),
               read.csv("Results/correlation_runtime_benchmark_8cores.csv")) %>%
  filter(method %in% par_methods) %>%
  mutate(method  = factor(method, levels = par_methods),
         nCores  = factor(nCores),
         n_genes = factor(paste(n_genes, 'genes'), levels = paste(gene_grid, 'genes')),
         n_cells = factor(paste(n_cells, 'cells'), levels = paste(cell_grid, 'cells')))

par_plot <- ggplot(d_par, aes(nCores, runtime_sec, colour = method, group = method)) +
  geom_point(size = 2) +
  geom_line() +
  facet_grid(n_genes ~ n_cells, scales = "free_y") +   # same y-range for all methods within a panel
  labs(x = "Cores", y = "Runtime (sec)", colour = NULL) +
  theme_classic() +
  scale_y_log10() +
  theme(strip.background = element_blank(), legend.position = 'bottom') +
  scale_colour_manual(values = c(Distance = 'grey50', Chatterjee = 'magenta3', ChatterjeeKNN = 'mediumpurple'))

par_plot
ggsave("figures/correlation_runtime_parallel_1vs8.png", par_plot, width = 6, height = 5, dpi = 300)
