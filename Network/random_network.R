library(graphlayouts)
library(SpiecEasi)
library(igraph)
library(dplyr)
library(ggplot2)
library(tidygraph)
library(ggraph)
library(psych)       
library(tidyr)


sample_order <- c(1, 7, 15, 16, 4, 6, 8, 9, 2, 3, 5, 10, 11, 12, 13, 14)
file_names <- sprintf("C:/RStudio/gut/code/data/CD%d_adj.csv", sample_order)
adj_list <- list()
for (i in seq_along(file_names)) {
  file_path <- file.path("C:/RStudio/gut/code/data/adj_matrices", file_names[i])
  adj_list[[i]] <- as.matrix(read.csv(file_path, row.names = 1, check.names = FALSE))
}

align_adjacency <- function(adj1, adj2) {
  species1 <- rownames(adj1)
  species2 <- rownames(adj2)
  common_species <- intersect(species1, species2)
  
  if(length(common_species) < 1) {
    return(list(adj1_common = NULL, adj2_common = NULL))
  }
  
  template <- matrix(0, nrow = length(common_species), 
                     ncol = length(common_species),
                     dimnames = list(common_species, common_species))
  
  adj1_common <- template
  idx1 <- match(common_species, species1)
  adj1_common[, ] <- adj1[idx1, idx1]
  
  adj2_common <- template
  idx2 <- match(common_species, species2)
  adj2_common[, ] <- adj2[idx2, idx2]
  
  return(list(adj1 = adj1_common, adj2 = adj2_common))
}


jaccard_similarity_adv <- function(adj1, adj2) {
  aligned <- align_adjacency(adj1, adj2)
  
  if(is.null(aligned$adj1)) return(NA)  
  
  upper <- upper.tri(aligned$adj1)
  x <- aligned$adj1[upper]
  y <- aligned$adj2[upper]
  
  intersection <- sum(x & y)
  union <- sum(x | y)
  
  return(ifelse(union == 0, 0, intersection/union))
}


generate_random_network <- function(adj_matrix) {
  orig_g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
  orig_nodes <- vcount(orig_g)
  orig_edges <- ecount(orig_g)
  
  random_g <- rewire(orig_g, keeping_degseq(niter = 100 * orig_edges))
  
  if (ecount(random_g) != orig_edges) {

    deg_seq <- degree(orig_g)
    random_g <- sample_degseq(deg_seq, method = "vl")
    
  }
  
  random_adj <- as_adjacency_matrix(random_g, type = "both", sparse = FALSE, 
                                    names = TRUE, attr = NULL)
  
  colnames(random_adj) <- rownames(random_adj) <- rownames(adj_matrix)
  
  return(random_adj)
}

r2a_samples <- c(1, 7, 15, 16)
sample_pairs <- combn(r2a_samples, 2)  
col_names <- apply(sample_pairs, 2, function(x) paste0("CD", x[1], "-CD", x[2]))
results <- matrix(nrow = 1000, ncol = 6)
colnames(results) <- col_names


set.seed(123)  
pb <- txtProgressBar(min = 0, max = 1000, style = 3)

for (sim in 1:1000) {

  rand_nets <- list()
  for (sample_id in r2a_samples) {
    idx <- which(sample_order == sample_id)
    orig_adj <- adj_list[[idx]]
    rand_nets[[as.character(sample_id)]] <- generate_random_network(orig_adj)
  }
  
  sim_vals <- apply(sample_pairs, 2, function(pair) {
    net1 <- rand_nets[[as.character(pair[1])]]
    net2 <- rand_nets[[as.character(pair[2])]]
    jaccard_similarity_adv(net1, net2)
  })
  
  results[sim, ] <- sim_vals
  setTxtProgressBar(pb, sim)
}
close(pb)

write.csv(results, "C:/RStudio/gut/code/data/random_network_pairwise_similarity_R2A.csv", row.names = FALSE)

r_samples <- c(2, 3, 5, 10, 11, 12, 13, 14)
sample_pairs <- combn(r_samples, 2)  
col_names <- apply(sample_pairs, 2, function(x) paste0("CD", x[1], "-CD", x[2]))
results <- matrix(nrow = 1000, ncol = 28)
colnames(results) <- col_names


set.seed(123)  
pb <- txtProgressBar(min = 0, max = 1000, style = 3)

for (sim in 1:1000) {

  rand_nets <- list()
  for (sample_id in r_samples) {
    idx <- which(sample_order == sample_id)
    orig_adj <- adj_list[[idx]]
    rand_nets[[as.character(sample_id)]] <- generate_random_network(orig_adj)
  }
  
  sim_vals <- apply(sample_pairs, 2, function(pair) {
    net1 <- rand_nets[[as.character(pair[1])]]
    net2 <- rand_nets[[as.character(pair[2])]]
    jaccard_similarity_adv(net1, net2)
  })
  
  results[sim, ] <- sim_vals
  setTxtProgressBar(pb, sim)
}
close(pb)


write.csv(results, "C:/RStudio/gut/code/data/random_network_pairwise_similarity_R.csv", row.names = FALSE)

ma_samples <- c(4, 6, 8, 9)
sample_pairs <- combn(ma_samples, 2)  
col_names <- apply(sample_pairs, 2, function(x) paste0("CD", x[1], "-CD", x[2]))
results <- matrix(nrow = 1000, ncol = 6)
colnames(results) <- col_names


set.seed(123)  
pb <- txtProgressBar(min = 0, max = 1000, style = 3)

for (sim in 1:1000) {

  rand_nets <- list()
  for (sample_id in ma_samples) {
    idx <- which(sample_order == sample_id)
    orig_adj <- adj_list[[idx]]
    rand_nets[[as.character(sample_id)]] <- generate_random_network(orig_adj)
  }
  
  sim_vals <- apply(sample_pairs, 2, function(pair) {
    net1 <- rand_nets[[as.character(pair[1])]]
    net2 <- rand_nets[[as.character(pair[2])]]
    jaccard_similarity_adv(net1, net2)
  })
  
  results[sim, ] <- sim_vals
  setTxtProgressBar(pb, sim)
}
close(pb)

write.csv(results, "C:/RStudio/gut/code/data/random_network_pairwise_similarity_mA.csv", row.names = FALSE)