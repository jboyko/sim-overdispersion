get_weights <- function(tree){
  C = vcv(tree)
  eig <- eigen(C)
  Q <- eig$vectors
  L <- diag((1/sqrt(eig$values)))
  P =  Q %*% L %*% t(Q)
  return(P)
}

mean_lab_color <- function(colors) {
  if (length(colors) == 1) return(colors)
  col_means <- colMeans(colors) 
  return(hex(RGB(t(col_means))))
}

paint_tree_branches <- function(tree, tip_colors) {
  # Initialize branch colors
  branch_colors <- rep(NA, Nedge(tree))  # One color per edge
  edge_names <- tree$edge[, 2]  # Node or tip associated with each edge
  
  # Assign tip colors
  for (tip in 1:Ntip(tree)) {
    tip_label <- tree$tip.label[tip]
    tip_color <- tip_colors[tip]
    branch_index <- which(edge_names == tip)
    branch_colors[branch_index] <- tip_color
  }
  
  # Prune tree in postorder traversal
  internal_nodes <- unique(reorder(tree, "postorder")$edge[,1])
  for (node in internal_nodes) {
    descendant_edges <- which(tree$edge[, 1] == node)  # Child edges
    descendant_colors <- branch_colors[descendant_edges]
    
    # Blend colors if multiple descendants have different colors
    if (length(unique(descendant_colors)) == 1) {
      branch_colors[which(edge_names == node)] <- unique(descendant_colors)
    } else {
      valid_colors <- descendant_colors[!is.na(descendant_colors)]  # Remove NAs
      blended_color <- mean_lab_color(hex2RGB(valid_colors)@coords)
      branch_colors[which(edge_names == node)] <- blended_color
    }
  }
  return(branch_colors)
}

get_file_table <- function(tables_to_load) {
  params_table <- do.call(rbind, lapply(basename(tables_to_load), function(fname) {
    matches <- regmatches(fname, regexec(
      "est_table_sim_(\\d+)_traits(\\d+)_taxa(\\d+)_nn(\\d+)_md([0-9.]+)_eps([0-9.]+)_minPts(\\d+)\\.csv", fname
    ))
    if (length(matches[[1]]) > 0) {
      c(
        File = fname,
        Index = as.numeric(matches[[1]][2]),
        Traits = as.numeric(matches[[1]][3]),
        Taxa = as.numeric(matches[[1]][4]),
        Neighbors = as.numeric(matches[[1]][5]),
        MinDist = as.numeric(matches[[1]][6]),
        Eps = as.numeric(matches[[1]][7]),
        MinPts = as.numeric(matches[[1]][8]),
        Regime = paste0("nn", matches[[1]][5], "_md", matches[[1]][6], "_eps", matches[[1]][7], "_minPts", matches[[1]][8])
      )
    } else {
      NULL
    }
  }))
  
  params_df <- as.data.frame(params_table, stringsAsFactors = FALSE)
  params_df$Index <- as.integer(params_df$Index)
  params_df$Traits <- as.integer(params_df$Traits)
  params_df$Taxa <- as.integer(params_df$Taxa)
  params_df$Neighbors <- as.integer(params_df$Neighbors)
  params_df$MinDist <- as.numeric(params_df$MinDist)
  params_df$Eps <- as.numeric(params_df$Eps)
  params_df$MinPts <- as.integer(params_df$MinPts)
  return(params_df)
}

single_run <- function(index, max_taxa, n_traits, n_neighbors=15, min_dist = 0.1, eps = 0.5, minPts = 10){
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
  dbscan_result <- dbscan(reduced_data, eps = eps, minPts = minPts)
  clusters <- dbscan_result$cluster
  cluster_sp_list <- split(tree$tip.label, clusters)
  
  # visualize
  print("outputing plot...")
  pdf(file = paste0("plots/cluster_plot_sim_", sprintf("%03d", index), 
    "_traits", n_traits, 
    "_taxa", max_taxa, 
    "_nn", n_neighbors, 
    "_md", min_dist, 
    "_eps", eps, 
    "_minPts", minPts, ".pdf"))
  num_clusters <- length(unique(clusters))
  cluster_colors <- hcl.colors(num_clusters, palette = "Dark 3")  # Choose a good categorical palette
  names(cluster_colors) <- unique(clusters)
  cluster_lab <- hex2RGB(cluster_colors)@coords
  centroids <- aggregate(reduced_data, by = list(cluster = clusters), FUN = mean)
  par(mfrow=c(1,2))
  plot(reduced_data, col = cluster_colors[as.character(clusters)], pch = 16,
    xlab = "UMAP Dimension 1", ylab = "UMAP Dimension 2", main = "UMAP Clustering")
  text(centroids[, 2], centroids[, 3], labels = centroids$cluster, col = "black", font = 2)
  tip_colors <- cluster_colors[as.character(clusters)]
  branch_colors <- paint_tree_branches(tree, tip_colors)
  plot.phylo(tree, show.tip.label = FALSE, no.margin = TRUE, direction = "leftwards", type = "fan", edge.color = branch_colors, edge.width = 0.25)
  tiplabels(pch = 16, col = tip_colors, cex = 0.25, offset = 0.5)
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
    "_md", min_dist, 
    "_eps", eps, 
    "_minPts", minPts, ".RDS"))
  saveRDS(fit2, file = paste0("out/fit2_sim_", sprintf("%03d", index), 
    "_traits", n_traits, 
    "_taxa", max_taxa, 
    "_nn", n_neighbors, 
    "_md", min_dist, 
    "_eps", eps, 
    "_minPts", minPts, ".RDS"))
  saveRDS(fit3, file = paste0("out/fit3_sim_", sprintf("%03d", index), 
    "_traits", n_traits, 
    "_taxa", max_taxa, 
    "_nn", n_neighbors, 
    "_md", min_dist, 
    "_eps", eps, 
    "_minPts", minPts, ".RDS"))
  
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
  
  sim_est <- cbind(rbind(full_est, clst_est, rndm_est), num_clusters = num_clusters)
  write.csv(sim_est, file = paste0("tables/est_table_", "sim_", sprintf("%03d", index), 
    "_traits", n_traits, 
    "_taxa", max_taxa, 
    "_nn", n_neighbors, 
    "_md", min_dist, 
    "_eps", eps, 
    "_minPts", minPts, ".csv"), row.names = FALSE)
  print("Done.")
}
