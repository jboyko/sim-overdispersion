setwd("~/sim-overdispersion/")
library(castor)
library(ape)
library(umap)
library(dbscan)

get_weights <- function(tree){
  C = vcv(tree)
  eig <- eigen(C)
  Q <- eig$vectors
  L <- diag((1/sqrt(eig$values)))
  P =  Q %*% L %*% t(Q)
  return(P)
}

tree <- read.tree("trees/squamates_Title_Science2024_ultrametric_constrained.tre", )

# simulate
D <- get_random_diffusivity_matrix(10, degrees=NULL, V=1)
all_states = simulate_bm_model(tree, diffusivity = D)
tip_states <- all_states$tip_states
rownames(tip_states) <- tree$tip.label

# umap
umap_result <- umap(tip_states, n_neighbors = 15, min_dist = 0.1, n_components = 2)
reduced_data <- umap_result$layout

plot(reduced_data)

# dbscan
dbscan_result <- dbscan(reduced_data, eps = 0.4, minPts = 10)
clusters <- dbscan_result$cluster
cluster_sp_list <- split(tree$tip.label, clusters)

# visualize
plot(reduced_data, col = clusters + 1, pch = 16,
  xlab = "UMAP Dimension 1", ylab = "UMAP Dimension 2", main = "UMAP Clustering")
legend("topright", legend = unique(clusters), col = unique(clusters + 1), pch = 16, ncol = 2)

# subsample by cluster
n_sample <- 20
subsample_cluster <- lapply(cluster_sp_list, function(x) sample(x, min(c(length(x), n_sample))))
subtree_cluster <- keep.tip(tree, unlist(subsample_cluster))
cluster_tip_states <- tip_states[subtree_cluster$tip.label,]

# subsample by random
subsamples_random <- sample(tree$tip.label, size = Ntip(subtree_cluster))
subtree_random <- keep.tip(tree, subsamples_random)
random_tip_states <- tip_states[subtree_random$tip.label,]

# subsample by phylorandom
# P <- get_weights(tree)
# X2 <- apply(tip_states_ind, 2, function(x) P %*% x)


# fit and compare BM models between the two data sets
fit = fit_and_compare_bm_models(
  trees1 = subtree_random, 
  tip_states1     = random_tip_states, 
  
  trees2          = subtree_cluster,
  tip_states2     = cluster_tip_states,

  Nbootstraps     = 100,
  Nsignificance   = 100)

# print summary of results
cat(sprintf("Fitted D1 = %g, D2 = %g, significance of log-diff. = %g\n",
  fit$fit1$diffusivity, fit$fit2$diffusivity, fit$significance))


fit$fit1$diffusivity - fit$fit2$diffusivity
