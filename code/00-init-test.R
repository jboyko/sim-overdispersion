setwd("~/sim-overdispersion/")
library(castor)
library(ape)
library(Rtsne)
library(dbscan)
library(phytools)

get_weights <- function(tree){
  C = vcv(tree)
  eig <- eigen(C)
  Q <- eig$vectors
  L <- diag((1/sqrt(eig$values)))
  P =  Q %*% L %*% t(Q)
  return(P)
}

tree <- read.tree("trees/squamates_Title_Science2024_ultrametric_constrained.tre")
max_tax <- 350
n_traits <- 10

single_run <- function(index, max_taxa, n_traits){
  # simulate
  # D <- get_random_diffusivity_matrix(10, degrees=NULL, V=1)
  print("generating data...")
  D <- diag(rnorm(n_traits)^2, n_traits, n_traits)
  all_states = simulate_bm_model(tree, diffusivity = D)
  tip_states <- all_states$tip_states
  rownames(tip_states) <- tree$tip.label
  
  # tsne
  print("running tnse...")
  tsne_result <- Rtsne(tip_states, dims = 2, perplexity = 50, verbose = FALSE)
  reduced_data <- tsne_result$Y
  
  # hdbscan
  print("clustering data...")
  dbscan_result <- dbscan(reduced_data, eps = 1.5, minPts = 10)
  clusters <- dbscan_result$cluster
  cluster_sp_list <- split(tree$tip.label, clusters)
  
  # visualize
  print("outputing plot...")
  pdf(file = paste0("plots/cluster_plot_sim_", index, ".pdf"))
  par(mfrow=c(1,2))
  plot(reduced_data, col = clusters + 1, pch = 16,
    xlab = "TSNE Dimension 1", ylab = "TSNE Dimension 2", main = "TSNE Clustering")
  # legend("topright", legend = unique(clusters), col = unique(clusters + 1), pch = 16, ncol = 2)
  plot(tree, show.tip.label = FALSE, no.margin = TRUE, direction = "leftwards", type = "fan")
  tiplabels(pch = 16, col = clusters + 1, cex = 0.5, offset = 0.5)
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
  saveRDS(fit1, file = paste0("out/fit1_sim", index, ".RDS"))
  saveRDS(fit2, file = paste0("out/fit2_sim", index, ".RDS"))
  saveRDS(fit3, file = paste0("out/fit3_sim", index, ".RDS"))
  
  print("summarizing results.")
  clst_full_sig <- unlist(lapply(fit1, "[[", "significance"))
  rndm_full_sig <- unlist(lapply(fit2, "[[", "significance"))
  rndm_clst_sig <- unlist(lapply(fit3, "[[", "significance"))
  
  clst_est <- do.call(rbind, lapply(fit1, function(x) 
    data.frame(type = "clst", sim_ind = index, est=x$fit2$diffusivity, ci_low=x$fit2$CI95lower, ci_upp=x$fit2$CI95upper)))
  rndm_est <- do.call(rbind, lapply(fit2, function(x) 
    data.frame(type = "rndm", sim_ind = index, est=x$fit2$diffusivity, ci_low=x$fit2$CI95lower, ci_upp=x$fit2$CI95upper)))
  full_est <- do.call(rbind, lapply(fit1, function(x) 
    data.frame(type = "rndm", sim_ind = index, est=x$fit1$diffusivity, ci_low=x$fit1$CI95lower, ci_upp=x$fit1$CI95upper)))
  
  sim_est <- rbind(full_est, clst_est, rndm_est)
  write.csv(sim_est, file = paste0("tables/est_table_", "sim", index, ".csv"), row.names = FALSE)
  print("Done.")
}

parallel::mclapply(1:100, function(x) single_run(x, 350, 10), mc.cores = 10)
