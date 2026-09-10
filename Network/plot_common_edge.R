library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Common edges are obtained from common_edge.R.
r2a_edges <- data.frame(
  from = c(
    "s__Butyrivibrio hungatei",
    "s__Bacteroides vulgatus",
    "s__Bacteroides intestinalis",
    "s__Bacteroides intestinalis",
    "s__Bacteroides cellulosilyticus",
    "s__Bacteroides helcogenes",
    "s__Lachnospiraceae bacterium GAM79",
    "s__Veillonella atypica",
    "s__Butyricimonas faecalis",
    "s__Clostridioides difficile",
    "s__Hungatella hathewayi",
    "s__Streptococcus sp. LPB0220"
  ),
  to = c(
    "s__Lachnoclostridium phytofermentans",
    "s__Bacteroides dorei",
    "s__Bacteroides cellulosilyticus",
    "s__Bacteroides helcogenes",
    "s__Bacteroides zoogleoformans",
    "s__Bacteroides heparinolyticus",
    "s__Roseburia hominis",
    "s__Veillonella dispar",
    "s__Paraprevotella xylaniphila",
    "s__Anaerobutyricum hallii",
    "s__Massilistercora timonensis",
    "s__Streptococcus sp. HSISM1"
  )
)

ma_edges <- data.frame(
  from = c(
    "s__Bacteroides helcogenes",
    "s__Veillonella atypica",
    "s__Bacteroides uniformis",
    "s__Flavonifractor plautii",
    "s__Bifidobacterium pseudocatenulatum",
    "s__Bacteroides cellulosilyticus",
    "s__Dysosmobacter welbionis",
    "s__Hungatella hathewayi",
    "s__Lacrimispora saccharolytica"
  ),
  to = c(
    "s__Bacteroides heparinolyticus",
    "s__Veillonella dispar",
    "s__Bacteroides sp. A1C1",
    "s__Clostridiales bacterium CCNA10",
    "s__Bifidobacterium catenulatum",
    "s__Bacteroides heparinolyticus",
    "s__Oscillibacter sp. PEA192",
    "s__Lacrimispora sphenoides",
    "s__Lacrimispora sphenoides"
  )
)

r_edges <- data.frame(
  from = c(
    "s__Bacteroides helcogenes",
    "s__Alistipes megaguti"
  ),
  to = c(
    "s__Bacteroides zoogleoformans",
    "s__Alistipes dispar"
  )
)

all_genera <- c(
  "Bacteroides", "Veillonella", "Flavonifractor", "Bifidobacterium", 
  "Dysosmobacter", "Hungatella", "Lacrimispora", "Clostridiales", 
  "Oscillibacter", "Butyrivibrio", "Lachnoclostridium", "Butyricimonas", 
  "Clostridioides", "Paraprevotella", "Anaerobutyricum", "Massilistercora", 
  "Streptococcus", "Roseburia", "Lachnospiraceae", "Alistipes"
)

# Take the subfigure of common edges in R group for instance.
graph <- graph_from_data_frame(r_edges, directed = FALSE) %>% 
  as_tbl_graph()

genus_names <- V(graph)$name %>% 
  gsub("s__([^ ]+).*", "\\1", .) %>% 
  gsub("\\[|\\]", "", .)  

label_colors <- ifelse(
  genus_names == "Veillonella", "darkred",
  ifelse(genus_names == "Alistipes", "darkblue", "black")
)

graph <- graph %>%
  activate(nodes) %>%
  mutate(
    genus = genus_names,
    genus = factor(genus_names, levels = all_genera),
    label_color = label_colors  
  )

# 创建颜色映射
genus_colors <- c(
  "Bacteroides" = "#8FC1E3",    
  "Veillonella" = "#F7C678",    
  "Flavonifractor" = "#87D68D", 
  "Bifidobacterium" = "#D8A7CA", 
  "Dysosmobacter" = "#7FB3D5",  
  "Hungatella" = "#F0A35C",     
  "Lacrimispora" = "#C6B2D9",   
  "Clostridiales" = "#A2D6C4",  
  "Oscillibacter" = "#68C3D4",  
  "Butyrivibrio" = "#E8C4A1",   
  "Lachnoclostridium" = "#B5D8CC", 
  "Butyricimonas" = "#D9E599",  
  "Clostridioides" = "#A8B8C8", 
  "Paraprevotella" = "#E1C7A5", 
  "Anaerobutyricum" = "#A5D8CC", 
  "Massilistercora" = "#D4C1DA", 
  "Streptococcus" = "#B8C9E6",   
  "Roseburia" = "#F5C7A9",       
  "Lachnospiraceae" = "#AEC9A7", 
  "Alistipes" = "#EDD3A7"        
)

graph <- graph %>%
  activate(nodes) %>%
  mutate(
    fill_color = ifelse(genus %in% names(genus_colors), 
                        genus_colors[genus], 
                        "gray50") 
  )

graph <- graph %>%
  activate(edges) %>%
  mutate(
    same_genus = (genus_names[from] == genus_names[to]),
    weight = ifelse(same_genus, 5, 1) 
  )

layout <- create_layout(
  graph,
  layout = "fr",  
  weights = weight,  
  niter = 10000,      
  start.temp = 0.01, 
  repulserad = vcount(graph)^2.2,
  area = vcount(graph)^1.5 
)

layout$x <- scales::rescale(layout$x, to = c(-1, 1))
layout$y <- scales::rescale(layout$y, to = c(-1, 1))

p <- ggraph(layout) + 
  
  geom_edge_link(
    edge_width = 0.5,   
    edge_color = "gray60",
    alpha = 0.7,
    show.legend = FALSE
  ) + 

  geom_node_point(
    aes(fill = genus),
    size = 6,          
    shape = 21,
    color = "black",
    stroke = 0.5,       
    alpha = 0.9,
    show.legend = TRUE
  ) +

  ggrepel::geom_text_repel(
    aes(x = x, y = y, label = name, color = label_color), 
    size = 2.5,         
    family = "sans",
    point.padding = 0.05,
    box.padding = 0.05,
    segment.color = "grey50",
    segment.size = 0.1,
    min.segment.length = 0,
    max.overlaps = 10,
    bg.color = "white",
    bg.r = 0.1,
    seed = 123,         
    force = 0.1,         
    fontface = "italic"
  ) +

  scale_fill_manual(
    name = "Genus",
    values = genus_colors,
    na.value = "grey50",
    drop = FALSE
  ) +

  scale_color_manual(
    values = c("darkred" = "darkred", "darkblue" = "darkblue", "black" = "black"),
    guide = "none"  
  ) +

  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold", hjust = 0.5),
    legend.text = element_text(face = "italic", size = 10),
    plot.margin = margin(0.5, 0.5, 1.5, 0.5, "cm"),  
    legend.box.margin = margin(0, 0, 0, 0),
    legend.box = "vertical",
    legend.direction = "horizontal"
  ) +
  guides(
    fill = guide_legend(
      nrow = 5,                 
      ncol = 4,                 
      byrow = TRUE,
      override.aes = list(
        size = 4,               
        alpha = 1, 
        stroke = 0.7,
        shape = 21
      ),
      title.position = "top",
      title.hjust = 0.5
    )
  )


p <- p + theme(legend.spacing.y = unit(0, "cm")) 

ggsave("C:/RStudio/gut/code/Figure5B_R.png", p, 
       width = 6, height = 6, dpi = 300, bg = "white")
