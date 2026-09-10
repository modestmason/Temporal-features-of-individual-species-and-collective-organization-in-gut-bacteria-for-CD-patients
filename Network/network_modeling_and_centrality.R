library(graphlayouts)
library(SpiecEasi)
library(igraph)
library(dplyr)
library(ggplot2)
library(tidygraph)
library(ggraph)
library(psych)       


# Process the data of each patient separately.
abundance_data <- read.csv("C:/RStudio/gut/code/data/CD1.csv", row.names = 1, header = TRUE)
abundance_matrix <- t(as.matrix(abundance_data))
mean_abundance <- colMeans(abundance_matrix)
filtered_species <- names(mean_abundance)[mean_abundance >= 0.0001]
filtered_matrix <- abundance_matrix[, filtered_species]
write.csv(filtered_matrix, "C:/RStudio/gut/code/data/filtered_matrix_CD1.csv", row.names = FALSE)
#Network modeling based on Meinshausen-Bühlman method, supported by SpiecEasi package.
se_mb <- spiec.easi(filtered_matrix, method="mb", lambda.min.ratio=1e-2, nlambda = 50)

adj_matrix_mb <- getRefit(se_mb)
adj_matrix_mb <- as.matrix(adj_matrix_mb)
rownames(adj_matrix_mb) <- colnames(filtered_matrix)
colnames(adj_matrix_mb) <- colnames(filtered_matrix)
sample_id <- tools::file_path_sans_ext(basename("C:/RStudio/gut/code/data/CD1.csv"))
output_file <- file.path("C:/RStudio/gut/code/data/adj_matrices", paste0(sample_id, "_adj.csv"))
write.csv(adj_matrix_mb, file = output_file, row.names = TRUE)


network_unsigned <- graph_from_adjacency_matrix(adj_matrix_mb, 
                                       mode = "undirected",
                                       diag = FALSE)
network_unsigned <- network_unsigned %>%
  as_tbl_graph() %>%
  activate(nodes) %>%         
  filter(degree(., mode = "all") > 0) %>% 
  activate(edges) 


comp <- components(network_unsigned)
largest_comp <- which.max(comp$csize)
network_largest <- induced_subgraph(
  network_unsigned,
  which(comp$membership == largest_comp)
)


closeness_largest <- closeness(network_largest, normalized = TRUE)

closeness_df <- data.frame(
  Species = names(closeness_largest),  
  Closeness = round(closeness_largest, 6),  
  stringsAsFactors = FALSE
)

closeness_df <- closeness_df[order(-closeness_df$Closeness), ]
closeness_df$Rank <- 1:nrow(closeness_df)
write.csv(closeness_df,
          "C:/RStudio/gut/code/Network/closeness_centrality_CD1.csv",
          row.names = FALSE)  
