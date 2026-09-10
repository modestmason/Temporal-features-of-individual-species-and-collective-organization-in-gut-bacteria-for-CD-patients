library(graphlayouts)
library(SpiecEasi)
library(igraph)
library(dplyr)
library(ggplot2)
library(tidygraph)
library(ggraph)
library(psych)       
library(tidyr)


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


Sorensen_Dice_adv <- function(adj1, adj2) {

  aligned <- align_adjacency(adj1, adj2)
  
  if(is.null(aligned$adj1)) return(NA)  
  
  upper <- upper.tri(aligned$adj1)
  x <- aligned$adj1[upper]
  y <- aligned$adj2[upper]

  
  return(2 * sum(x & y) / (sum(x) + sum(y)))
}





sample_order <- c(1, 7, 15, 16, 4, 6, 8, 9, 2, 3, 5, 10, 11, 12, 13, 14)
file_names <- sprintf("CD%d_adj.csv", sample_order)
adj_list <- list()
for (i in seq_along(file_names)) {
  file_path <- file.path("C:/RStudio/gut/code/data/adj_matrices", file_names[i])
  adj_list[[i]] <- as.matrix(read.csv(file_path, row.names = 1, check.names = FALSE))
}

n_samples <- length(sample_order)
sim_matrix <- matrix(NA, nrow = n_samples, ncol = n_samples)

for (i in 1:n_samples) {
  for (j in 1:n_samples) {
    sim_matrix[i, j] <- jaccard_similarity_adv(adj_list[[i]], adj_list[[j]])
  }
}

rownames(sim_matrix) <- paste0("CD", sample_order)
colnames(sim_matrix) <- paste0("CD", sample_order)
write.csv(sim_matrix, "jaccard_similarity_matrix.csv", row.names = TRUE)


