setwd("~/sim-overdispersion/")
source("code/utils.R")
library(castor)
library(ape)
library(umap)
library(dbscan)
library(phytools)
library(colorspace)


tree <- read.tree("trees/squamates_Title_Science2024_ultrametric_constrained.tre")
max_tax <- 350
n_traits <- 10

n_neighbors_values <- c(10, 20, 50)     # Number of neighbors
min_dist_values <- c(0.05, 0.1, 0.2)   # Minimum distance
eps_values <- c(0.5, 1.0, 1.5)         # DBSCAN epsilon
minPts_values <- c(5, 10, 20)          # Minimum points

param_grid <- expand.grid(
  n_neighbors = n_neighbors_values,
  min_dist = min_dist_values,
  eps = eps_values,
  minPts = minPts_values
)

for(i in 1:nrow(param_grid)){
  print(i)
  parallel::mclapply(1:100, function(x) 
    single_run(x, 
      max_taxa = max_tax, 
      n_traits = n_traits, 
      n_neighbors = param_grid$n_neighbors[i], 
      min_dist = param_grid$min_dist[i], 
      eps = param_grid$eps[i], 
      minPts = param_grid$minPts[i]), mc.cores = 10)
}

# Load the file metadata
tables_to_load <- dir("tables/", full.names = TRUE)
file_metadata <- get_file_table(tables_to_load)

# Load all data and merge with file metadata
big_table <- data.frame()
for (file in tables_to_load) {
  data <- read.csv(file)
  fname <- basename(file)
  metadata <- file_metadata[file_metadata$File == fname, ]
  data <- cbind(data, metadata[rep(1, nrow(data)), ])
  big_table <- rbind(big_table, data)
}

# Calculate differences for each trait and simulation
big_table$delta_est <- NA
unique_sim_indices <- unique(big_table$sim_ind)
for (sim in unique_sim_indices) {
  for (trait in 1:n_traits) {
    # Extract estimates for the current simulation and trait
    full_est <- big_table$est[big_table$sim_ind == sim & big_table$type == "full"]
    rndm_est <- big_table$est[big_table$sim_ind == sim & big_table$type == "rndm"]
    clst_est <- big_table$est[big_table$sim_ind == sim & big_table$type == "clst"]
    
    # Compute differences
    big_table$delta_est[big_table$sim_ind == sim & big_table$type == "rndm"] <- rndm_est - full_est
    big_table$delta_est[big_table$sim_ind == sim & big_table$type == "clst"] <- clst_est - full_est
  }
}

# Unique clustering regimes
unique_regimes <- unique(big_table$Regime)

# Plot densities for each regime
par(mfrow = c(ceiling(length(unique_regimes) / 2), 2))  # Arrange plots in grid
for (regime in unique_regimes) {
  regime_data <- big_table[big_table$Regime == regime, ]
  
  rndm_density <- density(regime_data$delta_est[regime_data$type == "rndm"], na.rm = TRUE)
  clst_density <- density(regime_data$delta_est[regime_data$type == "clst"], na.rm = TRUE)
  
  plot(rndm_density, col = "blue", lwd = 2, 
    main = paste("Density for Regime:", regime),
    xlab = "Difference (Subsample - Full)", 
    ylab = "Density",
    xlim = range(c(rndm_density$x, clst_density$x)),
    ylim = range(c(rndm_density$y, clst_density$y)))
  lines(clst_density, col = "red", lwd = 2)
  legend("topright", legend = c("Random", "Cluster-Based"), col = c("blue", "red"), lwd = 2)
}

results <- data.frame()
for (regime in unique_regimes) {
  regime_data <- big_table[big_table$Regime == regime, ]
  
  rndm_diffs <- regime_data$delta_est[regime_data$type == "rndm"]
  clst_diffs <- regime_data$delta_est[regime_data$type == "clst"]
  
  # T-tests
  t_test_rndm <- t.test(rndm_diffs, mu = 0, paired = FALSE)
  t_test_clst <- t.test(clst_diffs, mu = 0, paired = FALSE)
  
  # Wilcoxon tests
  wilcox_rndm <- wilcox.test(rndm_diffs, mu = 0)
  wilcox_clst <- wilcox.test(clst_diffs, mu = 0)
  
  # Additional metrics
  mean_rndm <- mean(rndm_diffs, na.rm = TRUE)
  mean_clst <- mean(clst_diffs, na.rm = TRUE)
  
  sd_rndm <- sd(rndm_diffs, na.rm = TRUE)
  sd_clst <- sd(clst_diffs, na.rm = TRUE)
  
  rmse_rndm <- sqrt(mean(rndm_diffs^2, na.rm = TRUE))
  rmse_clst <- sqrt(mean(clst_diffs^2, na.rm = TRUE))
  
  # Store results
  results <- rbind(results, data.frame(
    Regime = regime,
    Mean_Random = mean_rndm,
    SD_Random = sd_rndm,
    RMSE_Random = rmse_rndm,
    TTest_Random = t_test_rndm$p.value,
    Wilcox_Random = wilcox_rndm$p.value,
    
    Mean_Cluster = mean_clst,
    SD_Cluster = sd_clst,
    RMSE_Cluster = rmse_clst,
    TTest_Cluster = t_test_clst$p.value,
    Wilcox_Cluster = wilcox_clst$p.value
  ))
}

# View results
print(results)
