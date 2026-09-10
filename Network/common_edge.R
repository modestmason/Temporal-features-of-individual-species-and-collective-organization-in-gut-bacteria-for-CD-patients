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
file_names <- sprintf("CD%d_adj.csv", sample_order)
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

get_common_edges_in_group <- function(adj_list, group_indices, threshold = 0.75) {

  group_adj_list <- adj_list[group_indices]
  all_nodes <- lapply(group_adj_list, rownames)
  common_nodes <- Reduce(intersect, all_nodes)
  
  if (length(common_nodes) < 1) {
    stop("No common node.")
  }

  results <- data.frame(
    Node1 = character(0),
    Node2 = character(0),
    Count = integer(0),
    Proportion = numeric(0),
    Status = character(0)
  )
  
  for (i in 1:(length(common_nodes) - 1)) {
    for (j in (i + 1):length(common_nodes)) {
      node1 <- common_nodes[i]
      node2 <- common_nodes[j]
      edge_count <- 0
      for (adj in group_adj_list) {
        idx1 <- which(rownames(adj) == node1)
        idx2 <- which(colnames(adj) == node2)
        
        if (length(idx1) == 1 && length(idx2) == 1 && adj[idx1, idx2] == 1) {
          edge_count <- edge_count + 1
        }
      }

      proportion <- edge_count / length(group_adj_list)

      status <- ifelse(proportion == 1, "All present",
                       ifelse(proportion >= threshold, "Majority present",
                              "Minority present"))

      results <- rbind(results, data.frame(
        Node1 = node1,
        Node2 = node2,
        Count = edge_count,
        Proportion = proportion,
        Status = status
      ))
    }
  }

  results <- results[order(-results$Count, -results$Proportion), ]

  common_edges <- results[results$Status %in% c("All present", "Majority present"), ]
  
  return(list(
    all_edges = results,
    common_edges = common_edges,
    common_nodes = common_nodes,
    num_networks = length(group_adj_list)
  ))
}

r2a_indices <- c(1, 7, 15, 16)  
r2a_common_results <- get_common_edges_in_group(adj_list, r2a_indices, threshold = 0.75)
write.csv(common_results$r2a_common_edges, "R2A_common_edges.csv", row.names = FALSE)

ma_indices <- c(4, 6, 8, 9)  
ma_common_results <- get_common_edges_in_group(adj_list, ma_indices, threshold = 0.75)
write.csv(common_results$ma_common_edges, "mA_common_edges.csv", row.names = FALSE)

r_indices <- c(2, 3, 5, 10, 11, 12, 13, 14)  
r_common_results <- get_common_edges_in_group(adj_list, r_indices, threshold = 0.75)
write.csv(common_results$r_common_edges, "R_common_edges.csv", row.names = FALSE)