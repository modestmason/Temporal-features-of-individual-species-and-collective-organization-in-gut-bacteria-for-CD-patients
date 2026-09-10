library(graphlayouts)
library(SpiecEasi)
library(igraph)
library(dplyr)
library(ggplot2)
library(tidygraph)
library(ggraph)

# Network modeling and plot for each patient
abundance_data <- read.csv("C:/RStudio/gut/code/data/CD16.csv", row.names = 1, header = TRUE)

abundance_matrix <- t(as.matrix(abundance_data))
mean_abundance <- colMeans(abundance_matrix)
filtered_species <- names(mean_abundance)[mean_abundance >= 0.0001]
filtered_matrix <- abundance_matrix[, filtered_species]

# Network modeling based on Meinshausen-Bühlman method, supported by SpiecEasi package.
se_mb <- spiec.easi(filtered_matrix, method="mb", lambda.min.ratio=1e-2, nlambda = 50)
adj_matrix_mb <- getRefit(se_mb)
adj_matrix_mb <- as.matrix(adj_matrix_mb)
rownames(adj_matrix_mb) <- colnames(filtered_matrix)
colnames(adj_matrix_mb) <- colnames(filtered_matrix)

network_full <- graph_from_adjacency_matrix(adj_matrix_mb, 
                                            mode = "undirected",
                                            diag = FALSE) %>%
  as_tbl_graph()


network_full <- network_full %>%
  activate(nodes) %>%
  mutate(
    node_color = "#DA6C6C",  
    node_size = 7.5            
  )


set.seed(456)
layout <- create_layout(network_full, 
                        layout = "fr", 
                        weights = NULL)


network_plot <- ggraph(layout) +
  geom_edge_link(edge_color = "grey85", alpha = 1) +
  geom_node_point(aes(color = node_color),  
                  size = 3, 
                  alpha = 1) +
  scale_color_identity() +  
  
  theme_graph(base_family = "sans")


print(network_plot)
ggsave("C:/RStudio/gut/code/Network/CD16_network.tiff", 
       plot = network_plot,
       width = 10, height = 8, dpi = 600)
