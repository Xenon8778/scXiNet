# Load necessary libraries
library(XICOR)
library(dplyr)
library(igraph)
library(ggraph)
library(ggplot2)

# Title: Single-Cell RNA Sequencing Data Simulation with Gene Interactions

simulate_scrna_with_interactions <- function(n_genes = 100, n_cells = 500, normal_sd = 0.5,
                                             mean_expression = 5, dispersion = 0.5,
                                             dropout_rate = 0.2, n_interactions = 100, seed = 0) {

  #'  Simulates single-cell RNA sequencing (scRNA-seq) data with gene interactions.
  #'  This function approximates Poisson-distributed counts using Bernoulli trials "rbinom"
  #'  and incorporates various types of gene-gene interactions.
  #'
  #'  n_genes: Numeric. The number of genes to simulate.
  #'  n_cells: Numeric. The number of cells to simulate.
  #'  normal_sd Numeric. The standard deviation of the normal distribution used
  #'    to add noise to interaction effects.
  #'  mean_expression: Numeric. The average expression level across genes before
  #'    applying cell-specific size factors.
  #'  dispersion: Numeric. The overdispersion parameter for simulating gene means
  #'    from a Gamma distribution. A smaller value indicates more dispersion.
  #'  dropout_rate: Numeric. The probability of a count being set to zero (dropout).
  #'  n_interactions: Numeric. The number of gene-gene interactions to simulate.
  #'  seed: Numeric. A seed for random number generation to ensure reproducibility.
  #'
  #' Returns a list containing:
  #'   counts: A matrix of simulated scRNA-seq counts, with genes as rows and cells as columns.
  #'   gpairs: A data frame of simulated gene interaction pairs.

  # Set the seed for reproducibility of random number generation
  set.seed(seed)

  # Check if the number of requested interactions is feasible given the number of genes
  if (n_interactions > n_genes / 2) {
    stop("Too many interactions requested. Either increase the number of genes being simulated or lower the number of interactions requested.")
  }

  # Simulate gene-specific mean expression levels
  # These are drawn from a Gamma distribution, allowing for variability and a skew towards lower expression.
  gene_means <- rgamma(n_genes, shape = 1 / dispersion, scale = mean_expression * dispersion)

  # Simulate cell-specific size factors
  # These account for differences in sequencing depth or total RNA content across cells.
  size_factors <- rgamma(n_cells, shape = 1, scale = 1)

  # Initialize a matrix to store the simulated counts
  counts <- matrix(0, nrow = n_genes, ncol = n_cells)

  # Simulate counts for each gene and cell
  for (i in 1:n_genes) {
    for (j in 1:n_cells) {
      # Calculate the true mean expression for the current gene-cell pair
      true_mean <- gene_means[i] * size_factors[j]

      # Approximate Poisson-like counts with overdispersion using rbinom (Bernoulli trials).
      # This method is used to simulate count data, which often exhibit overdispersion
      # (variance greater than the mean), characteristic of Negative Binomial distributions.
      # While it mentions "Negative binomial approximated with Poisson" and then "rbinom",
      # the direct simulation is through rbinom to approximate a Poisson where
      # lambda = true_mean.
      lambda <- true_mean
      # Determine the number of Bernoulli trials. This ensures enough trials to approximate Poisson,
      # especially for small lambda values, and scales with lambda for larger values.
      trials <- max(1000, round(lambda * 10))
      # Calculate the probability of 'success' for each trial.
      # prob = lambda / trials ensures that trials * prob approximates lambda.
      prob <- lambda / trials

      # Simulate one count value using rbinom.
      counts[i, j] <- rbinom(1, trials, prob)
    }
  }

  # Add zero inflation (dropout) to the simulated data
  # A random mask is created, and counts corresponding to 'TRUE' in the mask are set to zero.
  dropout_mask <- matrix(runif(n_genes * n_cells) < dropout_rate, nrow = n_genes, ncol = n_cells)
  counts[dropout_mask] <- 0

  # Randomly shuffle cell columns to break any implicit ordering from simulation.
  counts <- counts[, sample(ncol(counts))]

  # Assign row and column names for clarity
  rownames(counts) <- paste0('g', 1:n_genes)
  colnames(counts) <- paste0('c', 1:n_cells)

  # Simulate gene-pair interactions
  # Reset seed for consistent gene pair selection across simulations if needed, though
  # the main function seed handles overall reproducibility.
  set.seed(0)
  # Randomly select 'n_interactions' unique gene pairs for interaction.
  gpairs <- matrix(sample(x = rownames(counts), size = 2 * n_interactions, replace = FALSE), ncol = 2)
  gpairs <- as.data.frame(gpairs) # Convert to data frame

  # Convert counts matrix to a data frame for easier manipulation with dplyr-like operations
  # Transposing ensures cells are rows and genes are columns, which is common for data frames
  # when operating on gene expression.
  scdata <- as.data.frame(t(counts))

  # Apply various interaction types to selected gene pairs
  for (i in 1:n_interactions) {
    g_source <- gpairs[i, 1] # Source gene
    g_target <- gpairs[i, 2] # Target gene

    # Apply different interaction models based on the loop index (i) modulo 4
    if (i %% 4 == 1) {
      # Linear interaction: target = 2 * source + noise
      scdata[g_target] <- abs(scdata[g_source] * 2 + rnorm(n_cells, sd = normal_sd))
    } else if (i %% 4 == 2) {
      # Parabolic interaction: target = -0.02 * source^2 + 1 * source + 5 + noise
      # Ensures no negative values are introduced by the parabolic function.
      scdata[g_target] <- -0.02 * scdata[g_source]^2 + 1 * scdata[g_source] + 5
      scdata[g_target][scdata[g_target] < 0] <- 0 # Clip negative values to 0
      scdata[g_target] <- scdata[g_target] + abs(rnorm(n_cells, sd = normal_sd)) # Add absolute normal noise
    } else if (i %% 4 == 3) {
      # Exponential interaction: target = exp(source / 15) + noise
      scdata[g_target] <- exp((scdata[g_source]/max(scdata[g_source])*5)) + rnorm(n_cells, sd = normal_sd)
      scdata[g_target][scdata[g_target] < 0] <- 0
    } else if (i %% 4 == 0) {
      # Sinusoidal interaction: target = sin(source / (max(source)/15)) * 5 + 10 + noise
      # Normalizes the source gene expression for the sine wave period.
      scdata[g_target] <- sin(scdata[g_source] / max(scdata[g_source] / 15)) * 5 + 10 + rnorm(n_cells, sd = normal_sd)
    }

    # Convert the modified (continuous) interaction results back to count data.
    # This simulates the effect of the interaction on the underlying true expression,
    # which is then sampled as counts. Negative values are implicitly handled by pmax for trials.
    lambda <- pmax(scdata[[g_target]], 0)
    trials <- pmax(1000, round(lambda * 10)) # Ensure at least 1000 trials
    prob <- lambda / trials

    scdata[g_target] <- rbinom(n_cells, trials, prob)
  }

  # Transpose the data frame back to a matrix with genes as rows and cells as columns
  # and ensure it's a matrix for consistency with the initial 'counts' object.
  counts_with_interactions <- t(as.matrix(scdata))

  # Return the simulated counts and the ground truth gene interaction pairs
  return(list(counts = counts_with_interactions, gpairs = gpairs))
}

#----
# Simulation Scenarios
#----

# Scenario 1: Single simulation run for a fixed number of cells (1000 cells)
# This block simulates data for 10 iterations with 1000 cells and saves the results.

message("Starting single simulation run for 1000 cells...")
nCells <- 1000
# Create a directory to store the simulation results if it doesn't exist.
if (!dir.exists('Simulation_scRNAseq')) {
  dir.create('Simulation_scRNAseq')
}

for (i in 1:10) {
  # Call the simulation function
  simulated_data <- simulate_scrna_with_interactions(n_genes = 100, n_cells = nCells, dropout_rate = 0.3,
                                                     mean_expression = 15, normal_sd = 0.2,
                                                     n_interactions = 50, seed = i)
  counts_result <- simulated_data$counts
  message(paste0("Completed iteration ", i, " for ", nCells, " cells."))
  # Save the simulated count data to a CSV file.
  write.csv(counts_result, paste0('Simulation_scRNAseq/counts_', i, '_', nCells, '.csv'))
}

# Process and save the ground truth gene regulatory network (GRN) from the last simulation run.
# Note: gpairs_result will be the same across all 'i' iterations because the 'seed' for gpairs
#       selection within simulate_scrna_with_interactions is hardcoded to 0.
gpairs_result <- simulated_data$gpairs
nonmonotonic <- vector('logical', length = nrow(gpairs_result))
# Determine if an interaction is "nonmonotonic" based on its type (modulo 4).
# Linear (i %% 4 == 1) are considered monotonic for this specific labeling,
# while others are considered nonmonotonic due to the nature of their functional form.
for (i in 1:nrow(gpairs_result)) {
  if (i %% 4 == 1) {
    nonmonotonic[i] = FALSE # Linear interactions might be treated as monotonic in some contexts
  } else {
    nonmonotonic[i] = TRUE  # Parabolic, exponential, sinusoidal are inherently nonmonotonic
  }
}
gpairs_result$nonmonotonic <- nonmonotonic # Add a 'nonmonotonic' column to the gene pairs data frame

# Save the ground truth gene regulatory network (interaction pairs with directionality).
write.csv(gpairs_result, 'Simulation_scRNAseq/gt_GRN.csv')

message("Single simulation run completed and saved.")

#----
# Network topology
#----
# Ground-truth GRN network built from gt_GRN.csv
# Columns: source, target, i, type, nonmonotonic
edges <- read.csv("Simulation_scRNAseq/gt_GRN.csv", stringsAsFactors = FALSE)

types <- c("Linear", "Parabolic", "Exponential", "Sinusoidal")
edges$type <- types[(edges$X - 1) %% 4 + 1]
gt_graph <- graph_from_data_frame(edges[, c("V1", "V2", "type", "nonmonotonic")])

type_cols <- c(Linear = "#4C78A8", Exponential = "#54A24B",
               Parabolic = "#E45756", Sinusoidal = "#B279A2")

message("Nodes: ", vcount(gt_graph), " | Edges: ", ecount(gt_graph),
        " | Components: ", components(gt_graph)$no)
print(table(edges$type))

set.seed(1)
p <- ggraph(gt_graph, layout = "fr") +
  geom_edge_link(aes(edge_colour = type), edge_width = 0.6,
                 arrow = arrow(length = unit(0.06, "in"), type = "closed"),
                 start_cap = circle(0, "in"), end_cap = circle(0.13, "in")) +
  geom_node_point(size = 8, shape = 21, fill = "white", colour = "black") +
  geom_node_text(aes(label = name), size = 3, family = "Arial") +
  scale_edge_colour_manual(values = type_cols, name = "Interaction") +
  theme_void(base_family = "Arial") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))

# Save plot
ggsave("figures/Simulation_scRNAseq_gt_GRN_network.png", p, width = 8, height = 8,
       dpi = 300, bg = "white", device = ragg::agg_png)


# Further notes:
# The simulated data can be used for benchmarking gene regulatory network inference methods,
# especially in the context of single-cell RNA sequencing data.
# The 'gt_GRN.csv' serves as the ground truth against which inferred networks can be compared.
