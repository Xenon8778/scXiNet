library(XICOR)
library(dplyr)

# Simulation
# Model - Negative binomial or Poisson (approximated)

# Simulate scRNA-seq data using rbinom (Bernoulli trials) to approximate Poisson-distributed counts.

# Parameters
n_genes <- 100         # Number of genes (reduced for faster simulation)
n_cells <- 500         # Number of cells
mean_expression <- 5   # Average expression level (lambda for Poisson approximation)
dispersion <- 0.5      # Overdispersion parameter (higher values mean more variability in gene means)

# Simulate gene-specific mean expression levels (e.g., from a Gamma distribution)
gene_means <- rgamma(n_genes, shape = 1/dispersion, scale = mean_expression * dispersion)

# Simulate cell-specific size factors (e.g., to account for sequencing depth differences)
size_factors <- rgamma(n_cells, shape = 1, scale = 1) # Mean 1, but random variation

# Create a matrix to store the simulated data
counts <- matrix(0, nrow = n_genes, ncol = n_cells)

# Simulate counts for each gene and cell
for (i in 1:n_genes) {
  for (j in 1:n_cells) {
    # Simulate Poisson-like expression with overdispersion
    true_mean <- gene_means[i] * size_factors[j]

    # Approximate Poisson with Negative binomial to add dispersion.
    # We approximate a negative binomial with a poisson for simplicity, using rbinom.

    # Simulate counts using rbinom (Bernoulli trials), approximating a poisson.
    # To use rbinom, we need to create a large number of trials, and then count successes.
    # The probability of success will be very low.

    #Using a poisson approximation.
    lambda = true_mean;
    trials = max(1000, round(lambda*10)); #Number of trials. Ensures at least 1000 trials.
    prob = lambda/trials; #Probability of success.

    counts[i, j] <- rbinom(1, trials, prob) #Simulate one cell for one gene.
  }
}

# Optional: Add some zero inflation (dropout)
dropout_rate <- 0.2
dropout_mask <- matrix(runif(n_genes * n_cells) < dropout_rate, nrow = n_genes, ncol = n_cells)
counts[dropout_mask] <- 0

# Name cells and genes
rownames(counts) <- paste0('g',1:n_genes) #Gene names
colnames(counts) <- paste0('c',1:n_cells) #Cell names

# Print the first few rows and columns of the simulated data
print(counts[1:5, 1:5])

hist(counts) # Basic histogram of the raw count data

# Gene-pair interactions
n_interactions = 20 #Number of interactions to simulate
set.seed(0) #Seed for reproducibility
gpairs = matrix(sample(x = 1:n_genes, size = 2*n_interactions), ncol = 2) #Generate random gene pairs
gpairs = as.data.frame(gpairs) #Convert to data frame for easier handling

# Convert counts to a data frame for easier manipulation
scdata <- as.data.frame(t(counts)) # transpose to have cells as rows and genes as columns

for (i in 1:n_interactions){
  g_source = paste0('g',gpairs[i,1]) #Source gene
  g_target = paste0('g',gpairs[i,2]) #Target gene
  if (i%%4 == 0) scdata[g_target] = abs(sin(scdata[g_source]/2)*2+3 + rnorm(n_cells)) #Apply sinusoidal interaction
  if (i%%4 == 1) scdata[g_target] = abs(log1p(scdata[g_source])*5 + rnorm(n_cells)) #Apply logarithmic interaction
  if (i%%4 == 2) scdata[g_target] = abs(cos(scdata[g_source]/2)*2+3 + rnorm(n_cells)) #Apply cosine interaction
  if (i%%4 == 3) scdata[g_target] = abs(scdata[g_source]^2 + scdata[g_source] + rnorm(n_cells)) #Apply quadratic interaction

  #Convert interaction results back to count data (using Poisson approximation)
  lambda = scdata[[g_target]] #Extract the interaction result vector
  trials = pmax(1000, round(lambda*10)) #Number of trials.
  prob = lambda/trials; #Probability of success.

  scdata[g_target] = rbinom(n_cells, trials, prob) #Convert the interaction back to counts.
}

# Transpose back to genes as rows and cells as columns.
counts_with_interactions <- t(as.matrix(scdata))

# Print the first few rows and columns of the simulated data with interactions
print(counts_with_interactions[1:5, 1:5])

write.csv(counts_with_interactions, 'Simulation_scRNAseq/counts.csv') #Save count data to csv
write.csv(gpairs,'Simulation_scRNAseq/gt_GRN.csv') #Save gene interactions to csv

gpairs #Display gene pairs used for interactions

# Plot gene pairs to visualize simulated relationships
plot(counts_with_interactions[14,], counts_with_interactions[96,])
plot(counts_with_interactions[68,], counts_with_interactions[44,])
plot(counts_with_interactions[39,], counts_with_interactions[33,])
plot(counts_with_interactions[1,], counts_with_interactions[35,])

# Calculate and compare linear (Pearson) correlation and non-linear (XICOR) correlation.
cor(counts_with_interactions[14,], counts_with_interactions[96,]) #Pearson correlation
xicor(counts_with_interactions[14,], counts_with_interactions[96,]) #XICOR correlation
xicor(scdata[,51], scdata[,14]) #XICOR Correlation on the transposed data.

cor(counts_with_interactions[68,], counts_with_interactions[44,])
xicor(counts_with_interactions[68,], counts_with_interactions[44,])

cor(counts_with_interactions[39,], counts_with_interactions[33,])
xicor(counts_with_interactions[39,], counts_with_interactions[33,])

cor(counts_with_interactions[1,], counts_with_interactions[35,])
xicor(counts_with_interactions[1,], counts_with_interactions[35,])
