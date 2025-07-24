# Load necessary libraries
library(XICOR)
library(dplyr)

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
      scdata[g_target] <- exp((scdata[g_source]) / 15) + rnorm(n_cells, sd = normal_sd)
    } else if (i %% 4 == 0) {
      # Sinusoidal interaction: target = sin(source / (max(source)/15)) * 5 + 10 + noise
      # Normalizes the source gene expression for the sine wave period.
      scdata[g_target] <- sin(scdata[g_source] / max(scdata[g_source] / 15)) * 5 + 10 + rnorm(n_cells, sd = normal_sd)
    }

    # Convert the modified (continuous) interaction results back to count data.
    # This simulates the effect of the interaction on the underlying true expression,
    # which is then sampled as counts. Negative values are implicitly handled by pmax for trials.
    lambda <- scdata[[g_target]]
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
directed <- vector('logical', length = nrow(gpairs_result))
# Determine if an interaction is "directed" based on its type (modulo 4).
# Linear (i %% 4 == 1) are considered undirected for this specific labeling,
# while others are considered directed due to the nature of their functional form.
for (i in 1:nrow(gpairs_result)) {
  if (i %% 4 == 1) {
    directed[i] = FALSE # Linear interactions might be treated as undirected in some contexts
  } else {
    directed[i] = TRUE  # Parabolic, exponential, sinusoidal are inherently directed
  }
}
gpairs_result$directed <- directed # Add a 'directed' column to the gene pairs data frame

# Save the ground truth gene regulatory network (interaction pairs with directionality).
write.csv(gpairs_result, 'Simulation_scRNAseq/gt_GRN.csv')

message("Single simulation run completed and saved.")


# Scenario 2: Multiple simulation runs across varying cell numbers
# This block performs simulations for a range of cell numbers (500 to 10000, in steps of 500)
# and for 10 iterations each, saving the results to a different directory.

message("\nStarting multiple simulation runs across varying cell numbers...")
# Define the sequence of cell numbers to simulate.
nCells_range <- seq(500, 10000, 500)

# Create a directory to store the multiple simulation results if it doesn't exist.
if (!dir.exists('Simulation_scRNAseq_multiple')) {
  dir.create('Simulation_scRNAseq_multiple')
}

for (k in nCells_range) {
  message(paste0("Simulating for nCells = ", k))
  for (i in 1:10) {
    # Call the simulation function with varying cell numbers and iteration-specific seeds.
    simulated_data <- simulate_scrna_with_interactions(n_genes = 100, n_cells = k, dropout_rate = 0.3,
                                                       mean_expression = 15, normal_sd = 0.2,
                                                       n_interactions = 50, seed = i)
    counts_result <- simulated_data$counts

    # Save the simulated count data.
    write.csv(counts_result, paste0('Simulation_scRNAseq_multiple/counts_', i, '_', k, '.csv'))
  }
}

# Process and save the ground truth gene regulatory network (GRN) for the multiple simulation scenario.
# This will be identical to the previous 'gt_GRN.csv' since the gene pairs and their
# assigned directionality are determined by the fixed seed within the function.
gpairs_result <- simulated_data$gpairs # Takes the gpairs from the very last simulation (k=10000, i=10)
directed <- vector('logical', length = nrow(gpairs_result))
for (i in 1:nrow(gpairs_result)) {
  if (i %% 4 == 1) {
    directed[i] = FALSE
  } else {
    directed[i] = TRUE
  }
}
gpairs_result$directed <- directed

# Save the ground truth gene regulatory network.
write.csv(gpairs_result, 'Simulation_scRNAseq_multiple/gt_GRN.csv')

message("Multiple simulation runs completed and saved.")

# Further notes:
# The simulated data can be used for benchmarking gene regulatory network inference methods,
# especially in the context of single-cell RNA sequencing data.
# The 'gt_GRN.csv' serves as the ground truth against which inferred networks can be compared.
# The varying cell numbers allow for assessing method performance with different sample sizes.
