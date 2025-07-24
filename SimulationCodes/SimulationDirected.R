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

# Scenario 1: Fixed number of cells (1000 cells), multiple iterations
# This block performs 10 simulation runs with a consistent set of parameters
# (500 genes, 1000 cells, 250 interactions) and saves the output.

message("Starting simulation for directed networks with fixed cell count (1000 cells)...")
# Create a directory to store the simulation results if it doesn't already exist.
if (!dir.exists('Simulation_Directed')) {
  dir.create('Simulation_Directed')
}

for (i in 1:10) {
  # Call the simulation function with specific parameters for this scenario.
  simulated_data <- simulate_scrna_with_interactions(n_genes = 500, n_cells = 1000, dropout_rate = 0.3,
                                                     mean_expression = 15, normal_sd = 0.5,
                                                     n_interactions = 250, seed = i)
  counts_result <- simulated_data$counts
  message(paste0("  Completed iteration ", i))
  # Save the simulated count data to a CSV file.
  write.csv(counts_result, paste0('Simulation_Directed/counts_', i, '.csv'))
}

# Extract and process the ground truth gene regulatory network (GRN).
# Note: `gpairs_result` will hold the `gpairs` from the *last* simulation iteration (i=10)
# because the seed for `gpairs` generation within the function is fixed (set.seed(0)).
gpairs_result <- simulated_data$gpairs
directed <- vector('logical', length = nrow(gpairs_result))
# Assign a 'directed' flag based on the interaction type.
# Linear interactions (i %% 4 == 1) are flagged as FALSE (potentially undirected),
# while others are TRUE (implying a directed functional relationship).
for (i in 1:nrow(gpairs_result)) {
  if (i %% 4 == 1 | i %% 4 == 3) {
    directed[i] = FALSE # Linear interactions can be symmetric, hence 'undirected'
  } else {
    directed[i] = TRUE  # Non-linear functions imply a clear input-output direction
  }
}
gpairs_result$directed <- directed # Add the 'directed' column to the ground truth pairs.

# Save the ground truth gene interaction pairs to a CSV file.
write.csv(gpairs_result, 'Simulation_Directed/gt_GRN.csv')

message("Fixed cell count simulation completed and ground truth GRN saved.")

# --- Scenario 2: Multiple Simulations across Varying Cell Numbers ---
# This block runs simulations for a range of cell counts, allowing for
# evaluating method performance at different sample sizes.

message("\nStarting multiple simulations for directed networks across varying cell numbers...")
# Define the range of cell numbers to simulate, from 500 to 10000 in steps of 500.
nCells_range <- seq(500, 10000, 500)

# Create a directory for these simulation results.
if (!dir.exists('Simulation_Directed_multiple')) {
  dir.create('Simulation_Directed_multiple')
}

for (k in nCells_range) {
  message(paste0("  Simulating for nCells = ", k))
  for (i in 1:10) {
    # Call the simulation function for each combination of cell count and iteration.
    simulated_data <- simulate_scrna_with_interactions(n_genes = 500, n_cells = k, dropout_rate = 0.3,
                                                       mean_expression = 15, normal_sd = 0.5,
                                                       n_interactions = 250, seed = i)
    counts_result <- simulated_data$counts

    # Save the simulated count data, including the cell count in the filename.
    write.csv(counts_result, paste0('Simulation_Directed_multiple/counts_', i, '_', k, '.csv'))
  }
}

# Extract and process the ground truth GRN for this scenario.
# Similar to the previous block, this `gpairs_result` will be consistent due to `set.seed(0)`
# for interaction pair selection within the function.
gpairs_result <- simulated_data$gpairs
directed <- vector('logical', length = nrow(gpairs_result))
for (i in 1:nrow(gpairs_result)) {
  if (i %% 4 == 1 | i %% 4 == 3) {
    directed[i] = FALSE
  } else {
    directed[i] = TRUE
  }
}
gpairs_result$directed <- directed
# Save the ground truth GRN for the multiple simulation scenario.
write.csv(gpairs_result, 'Simulation_Directed_multiple/gt_GRN.csv')

message("Multiple simulations across varying cell numbers completed and ground truth GRN saved.")

# Concluding remarks:
# This R script provides a robust framework for simulating scRNA-seq data with
# controlled gene-gene interactions. The simulated data and the associated
# ground truth interaction networks (`gt_GRN.csv` files) are valuable resources
# for developing and benchmarking gene regulatory network inference algorithms,
# particularly those designed for single-cell data and aiming to infer directed relationships.
