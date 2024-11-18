setwd("~/sim-overdispersion/")
library(castor)
library(ape)
library(umap)
library(dbscan)
library(phytools)
library(colorspace)

get_weights <- function(tree){
  C = vcv(tree)
  eig <- eigen(C)
  Q <- eig$vectors
  L <- diag((1/sqrt(eig$values)))
  P =  Q %*% L %*% t(Q)
  return(P)
}

mean_lab_color <- function(colors) {
  if (length(colors) == 1) return(colors)  # Single color, no blending
  col_means <- colMeans(cluster_lab[colors, , drop = FALSE])  # Mean in LAB space
  RGB(as(col_means, "RGB"))@coords  # Convert back to RGB
}

paint_tree_branches <- function(tree, tip_colors) {
  # Convert tip colors to LAB color space
  unique_colors <- unique(tip_colors)
  color_lab <- hex2RGB(unique_colors)@coords
  names(color_lab) <- unique_colors
  
  # Initialize branch colors
  branch_colors <- rep(NA, Nedge(tree))  # One color per edge
  edge_names <- tree$edge[, 2]  # Node or tip associated with each edge
  
  # Assign tip colors
  for (tip in 1:Ntip(tree)) {
    tip_label <- tree$tip.label[tip]
    tip_color <- tip_colors[tip_label]
    branch_index <- which(edge_names == tip)
    branch_colors[branch_index] <- tip_color
  }
  
  # Prune tree in postorder traversal
  internal_nodes <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))  # Internal nodes
  for (node in rev(internal_nodes)) {
    descendant_edges <- which(tree$edge[, 1] == node)  # Child edges
    descendant_colors <- branch_colors[descendant_edges]
    
    # Blend colors if multiple descendants have different colors
    if (length(unique(descendant_colors)) == 1) {
      branch_colors[which(edge_names == node)] <- unique(descendant_colors)
    } else {
      valid_colors <- descendant_colors[!is.na(descendant_colors)]  # Remove NAs
      blended_color <- mean_lab_color(valid_colors)
      branch_colors[which(edge_names == node)] <- blended_color
    }
  }
  return(branch_colors)
}


tree <- read.tree("trees/squamates_Title_Science2024_ultrametric_constrained.tre")
max_tax <- 350
n_traits <- 10

single_run <- function(index, max_taxa, n_traits, n_neighbors=15, min_dist = 0.1){
  # simulate
  # D <- get_random_diffusivity_matrix(10, degrees=NULL, V=1)
  print("generating data...")
  D <- diag(rnorm(n_traits)^2, n_traits, n_traits)
  all_states = simulate_bm_model(tree, diffusivity = D)
  tip_states <- all_states$tip_states
  rownames(tip_states) <- tree$tip.label
  
  # tsne
  print("running umap...")
  # tsne_result <- Rtsne(tip_states, dims = 2, perplexity = 50, verbose = FALSE)
  # reduced_data <- tsne_result$Y
  umap_result <- umap(tip_states, n_neighbors = n_neighbors, min_dist = min_dist, n_components = 2)
  reduced_data <- umap_result$layout
  
  
  # hdbscan
  print("clustering data...")
  dbscan_result <- dbscan(reduced_data, eps = .5, minPts = 10)
  clusters <- dbscan_result$cluster
  cluster_sp_list <- split(tree$tip.label, clusters)
  
  # visualize
  print("outputing plot...")
  pdf(file = paste0("plots/cluster_plot_sim_", sprintf("%03d", index), 
    "_traits", n_traits, 
    "_taxa", max_taxa, 
    "_nn", n_neighbors, 
    "_md", min_dist, ".pdf"))
  num_clusters <- length(unique(clusters))
  cluster_colors <- hcl.colors(num_clusters, palette = "Dark 3")  # Choose a good categorical palette
  names(cluster_colors) <- unique(clusters)
  cluster_lab <- hex2RGB(cluster_colors)@coords
  centroids <- aggregate(reduced_data, by = list(cluster = clusters), FUN = mean)
  par(mfrow=c(1,2))
  plot(reduced_data, col = cluster_colors[as.character(clusters)], pch = 16,
    xlab = "UMAP Dimension 1", ylab = "UMAP Dimension 2", main = "UMAP Clustering")
  text(centroids[, 2], centroids[, 3], labels = centroids$cluster, col = "black", font = 2)
  # legend("topright", legend = names(cluster_colors), col = cluster_colors, 
    # pch = 16, ncol = 2, title = "Clusters")
  # legend("topright", legend = unique(clusters), col = unique(clusters + 1), pch = 16, ncol = 2)
  tip_colors <- cluster_colors[as.character(clusters)]
  plot(tree, show.tip.label = FALSE, no.margin = TRUE, direction = "leftwards", type = "fan")
  tiplabels(pch = 16, col = tip_colors, cex = 0.5, offset = 0.5)
  dev.off()
  
  # subsample by cluster
  print("subsample by cluster...")
  n_sample <- round(max_tax/max(clusters))
  subsample_cluster <- lapply(cluster_sp_list, function(x) sample(x, min(c(length(x), n_sample))))
  subtree_cluster <- keep.tip(tree, unlist(subsample_cluster))
  cluster_tip_states <- tip_states[subtree_cluster$tip.label,]
  
  # subsample by random
  print("subsample randomly...")
  subsamples_random <- sample(tree$tip.label, size = Ntip(subtree_cluster))
  subtree_random <- keep.tip(tree, subsamples_random)
  random_tip_states <- tip_states[subtree_random$tip.label,]
  
  print("fitting BM models...")
  fit1=fit2=fit3=list()
  for(i in 1:n_traits){
    cat("\r", i, "out of", n_traits, "...")
    fit1[[i]] = fit_and_compare_bm_models(
      trees1          = tree, 
      tip_states1     = tip_states[,i], 
      trees2          = subtree_cluster,
      tip_states2     = cluster_tip_states[,i],
      Nbootstraps     = 100,
      Nsignificance   = 100)
    fit2[[i]] = fit_and_compare_bm_models(
      trees1          = tree, 
      tip_states1     = tip_states[,i], 
      trees2          = subtree_random,
      tip_states2     = random_tip_states[,i],
      Nbootstraps     = 100,
      Nsignificance   = 100)
    fit3[[i]] = fit_and_compare_bm_models(
      trees1          = subtree_random, 
      tip_states1     = random_tip_states[,i], 
      trees2          = subtree_cluster,
      tip_states2     = cluster_tip_states[,i],
      Nbootstraps     = 100,
      Nsignificance   = 100)
  }
  
  print("saving results.")
  saveRDS(fit1, file = paste0("out/fit1_sim_", sprintf("%03d", index), 
    "_traits", n_traits, 
    "_taxa", max_taxa, 
    "_nn", n_neighbors, 
    "_md", min_dist, ".RDS"))
  saveRDS(fit2, file = paste0("out/fit2_sim_", sprintf("%03d", index), 
    "_traits", n_traits, 
    "_taxa", max_taxa, 
    "_nn", n_neighbors, 
    "_md", min_dist, ".RDS"))
  saveRDS(fit3, file = paste0("out/fit3_sim_", sprintf("%03d", index), 
    "_traits", n_traits, 
    "_taxa", max_taxa, 
    "_nn", n_neighbors, 
    "_md", min_dist, ".RDS"))
  
  print("summarizing results.")
  clst_full_sig <- unlist(lapply(fit1, "[[", "significance"))
  rndm_full_sig <- unlist(lapply(fit2, "[[", "significance"))
  rndm_clst_sig <- unlist(lapply(fit3, "[[", "significance"))
  
  clst_est <- cbind(do.call(rbind, lapply(fit1, function(x) 
    data.frame(type = "clst", sim_ind = index, est=x$fit2$diffusivity, ci_low=x$fit2$CI95lower, ci_upp=x$fit2$CI95upper))), sig = clst_full_sig)
  rndm_est <- cbind(do.call(rbind, lapply(fit2, function(x) 
    data.frame(type = "rndm", sim_ind = index, est=x$fit2$diffusivity, ci_low=x$fit2$CI95lower, ci_upp=x$fit2$CI95upper))), sig = rndm_full_sig)
  full_est <- cbind(do.call(rbind, lapply(fit1, function(x) 
    data.frame(type = "full", sim_ind = index, est=x$fit1$diffusivity, ci_low=x$fit1$CI95lower, ci_upp=x$fit1$CI95upper))), sig = rndm_clst_sig)
  
  sim_est <- rbind(full_est, clst_est, rndm_est)
  write.csv(sim_est, file = paste0("tables/est_table_", "sim_", sprintf("%03d", index), 
    "_traits", n_traits, 
    "_taxa", max_taxa, 
    "_nn", n_neighbors, 
    "_md", min_dist, ".csv"), row.names = FALSE)
  print("Done.")
}

# parallel::mclapply(1:100, function(x) single_run(x, 350, 10, 15, 0.1), mc.cores = 15)

tables_to_load <- dir("tables/", full.names = TRUE)
big_table <- data.frame()
for(i in tables_to_load){
  big_table <- rbind(big_table, read.csv(i))
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

# Compute density for random and cluster-based subsampling
rndm_density <- density(diff_data$delta_est[diff_data$type == "rndm"], na.rm = TRUE)
clst_density <- density(diff_data$delta_est[diff_data$type == "clst"], na.rm = TRUE)

# Plot the densities
plot(rndm_density, col = "blue", lwd = 2, 
  xlab = "Difference (Subsample - Full)", 
  ylab = "Density", 
  main = "Density of Differences in Rate Estimates",
  xlim = range(c(rndm_density$x, clst_density$x)),
  ylim = range(c(rndm_density$y, clst_density$y)))
lines(clst_density, col = "red", lwd = 2)

# Add a legend
legend("topright", legend = c("Random", "Cluster-Based"), 
  col = c("blue", "red"), lwd = 2)


# Merge estimates for plotting
estimates <- merge(
  subset(big_table, type == "full")[, c("sim_ind", "est")],
  subset(big_table, type != "full")[, c("sim_ind", "type", "est")],
  by = c("sim_ind"),
  suffixes = c("_full", "_subsample")
)

# Subset for random and cluster subsampling
rndm_est <- subset(estimates, type == "rndm")
clst_est <- subset(estimates, type == "clst")

# Scatter plot for random subsampling
par(mfrow=c(1,2))
plot(rndm_est$est_full, rndm_est$est_subsample,
  main = "Random Subsampling vs Full Data",
  xlab = "Full Data Estimates",
  ylab = "Random Subsample Estimates",
  pch = 16, col = "blue")
abline(0, 1, lty = 2)  # Reference line (y = x)

# Scatter plot for cluster-based subsampling
plot(clst_est$est_full, clst_est$est_subsample,
  main = "Cluster-Based Subsampling vs Full Data",
  xlab = "Full Data Estimates",
  ylab = "Cluster-Based Subsample Estimates",
  pch = 16, col = "red")
abline(0, 1, lty = 2)  # Reference line (y = x)



# Differences for random and cluster subsampling
rndm_diffs <- diff_data$delta_est[diff_data$type == "rndm"]
clst_diffs <- diff_data$delta_est[diff_data$type == "clst"]

# Paired t-tests
t_test_rndm <- t.test(rndm_diffs, mu = 0, paired = FALSE)
t_test_clst <- t.test(clst_diffs, mu = 0, paired = FALSE)

# Print results
t_test_rndm
t_test_clst

# Wilcoxon signed-rank test
wilcox_rndm <- wilcox.test(rndm_diffs, mu = 0)
wilcox_clst <- wilcox.test(clst_diffs, mu = 0)

# Print results
wilcox_rndm
wilcox_clst

# Mean and SD for random subsampling
mean_rndm <- mean(rndm_diffs, na.rm = TRUE)
sd_rndm <- sd(rndm_diffs, na.rm = TRUE)

# Mean and SD for cluster-based subsampling
mean_clst <- mean(clst_diffs, na.rm = TRUE)
sd_clst <- sd(clst_diffs, na.rm = TRUE)

# Print results
cat("Random Subsampling - Mean:", mean_rndm, "SD:", sd_rndm, "\n")
cat("Cluster Subsampling - Mean:", mean_clst, "SD:", sd_clst, "\n")
