library(graphlayouts)
library(SpiecEasi)
library(igraph)
library(dplyr)
library(ggplot2)
library(tidygraph)
library(ggraph)
library(psych)       
library(tidyr)
library(ggrepel)
library(patchwork)

deep_red_species <- c(
  "s__Sutterella faecalis", "s__Veillonella atypica", "s__Clostridium septicum", 
  "s__Enterobacter hormaechei", "s__Shigella dysenteriae", "s__Fusobacterium ulcerans",
  "s__Morganella morganii", "s__Eikenella corrodens", "s__Klebsiella michiganensis",
  "s__Klebsiella oxytoca","s__Klebsiella aerogenes"
)


sample_files <- list.files("C:/RStudio/gut/code/Network/", pattern = "closeness_centrality_", full.names = TRUE)  
all_data <- lapply(sample_files, function(file) {  
  df <- read.csv(file)  
  df$SampleID <- gsub("closeness_centrality_|\\.csv", "", basename(file))  
  return(df)  
})  
all_data <- bind_rows(all_data)

group_mapping <- c(
  "CD2" = "R", "CD3" = "R", "CD5" = "R", "CD10" = "R", 
  "CD11" = "R", "CD12" = "R", "CD13" = "R", "CD14" = "R",
  "CD4" = "mA", "CD6" = "mA", "CD8" = "mA", "CD9" = "mA",
  "CD1" = "R2A", "CD7" = "R2A", "CD15" = "R2A", "CD16" = "R2A"
)

all_data <- all_data %>%
  mutate(
    Group = case_when(
      SampleID %in% names(group_mapping)[group_mapping == "R"] ~ "R",
      SampleID %in% names(group_mapping)[group_mapping == "mA"] ~ "mA",
      SampleID %in% names(group_mapping)[group_mapping == "R2A"] ~ "R2A"
    )
  )


all_species <- unique(all_data$Species)

group_info <- all_data %>%
  select(SampleID, Group) %>%
  distinct(SampleID, .keep_all = TRUE)

complete_data <- all_data %>%
  select(-Rank) %>%
  complete(SampleID, Species = all_species, fill = list(Closeness = 0)) %>%
  select(-any_of("Group")) %>%
  left_join(group_info, by = "SampleID")

library(rstatix)

kruskal_results <- complete_data %>%
  group_by(Species) %>%
  kruskal_test(Closeness ~ Group) %>%
  adjust_pvalue(method = "none") %>%
  add_significance("p.adj")

significant_species <- kruskal_results %>%
  filter(p < 0.05)

dunn_results <- complete_data %>%
  filter(Species %in% significant_species$Species) %>%
  group_by(Species) %>%
  dunn_test(Closeness ~ Group, p.adjust.method = "BH") %>%
  mutate(
    new_group_pair = case_when(
      group1 == "mA" & group2 == "R" ~ "mA_vs_R",
      group1 == "R" & group2 == "mA" ~ "mA_vs_R",
      group1 == "R2A" & group2 == "mA" ~ "R2A_vs_mA",
      group1 == "mA" & group2 == "R2A" ~ "R2A_vs_mA",
      group1 == "R2A" & group2 == "R" ~ "R2A_vs_R",
      group1 == "R" & group2 == "R2A" ~ "R2A_vs_R"
    )
  ) %>%
  separate(new_group_pair, into = c("new_group1", "new_group2"), sep = "_vs_", remove = FALSE) %>%
  mutate(
    group1 = new_group1,
    group2 = new_group2
  ) %>%
  select(-new_group_pair, -new_group1, -new_group2) %>%
  mutate(group_pair = paste0(group1, "_vs_", group2)) %>%
  arrange(factor(group_pair, levels = c("mA_vs_R", "R2A_vs_mA", "R2A_vs_R")))


species_medians <- complete_data %>%
  group_by(Species, Group) %>%
  summarise(median_closeness = median(Closeness, na.rm = TRUE), .groups = "drop")

# 添加log2FC值
dunn_results <- dunn_results %>%
  left_join(species_medians, by = c("Species" = "Species", "group1" = "Group")) %>%
  rename(median_group1 = median_closeness) %>%
  left_join(species_medians, by = c("Species" = "Species", "group2" = "Group")) %>%
  rename(median_group2 = median_closeness) %>%
  mutate(
    log2FC = case_when(
      median_group1 == 0 & median_group2 == 0 ~ 0,
      median_group1 == 0 & median_group2 > 0 ~ log2((median_group2 + 1e-4) / (1e-4)),
      median_group1 > 0 & median_group2 == 0 ~ log2(1e-4 / (median_group1 + 1e-4)),
      TRUE ~ log2((median_group2 + 1e-4) / (median_group1 + 1e-4))
    )
  )

plot_volcano <- function(comparison, data, top_n = 6) {
  # 定义正确的横坐标标签
  comparison_labels <- list(
    "R_vs_mA" = "log₂(R/mA)",
    "mA_vs_R2A" = "log₂(mA/R2A)",
    "R_vs_R2A" = "log₂(R/R2A)"
  )
  
  groups <- strsplit(comparison, "_vs_")[[1]]
  group1 <- groups[1]
  group2 <- groups[2]

  plot_data <- data %>%
    filter(
      Comparison == comparison,
      !is.na(p),
      is.finite(log2FC)
    )
  
  if(nrow(plot_data) == 0) {
    return(ggplot() + 
             geom_blank() + 
             labs(title = paste("Volcano Plot:", gsub("_vs_", " vs ", comparison)),
                  x = comparison_labels[[comparison]],
                  y = "-log₁₀(p-value)") +
             theme_classic())
  }
  

  plot_data <- plot_data %>%
    mutate(
      abs_log2FC = abs(log2FC),
      
      PointType = case_when(
        
        abs_log2FC >= 1 & p.adj <= 0.1 &
          (median_group1 == 0 | median_group2 == 0) ~ 
          case_when(
            median_group1 == 0 & median_group2 == 0 ~ "DeepGray",
            TRUE ~ paste0("Deep", 
                          
                          ifelse(log2FC > 0, group2, group1))
          ),
        
        abs_log2FC >= 1 & p.adj <= 0.1 ~ 
          paste0("Light", 
                 
                 ifelse(log2FC > 0, group2, group1)),
        
        TRUE ~ "Gray"
      )
    )
  

  light_points <- plot_data %>%
    filter(startsWith(PointType, "Light")) %>%
    mutate(Label = Species)
  
  deep_points <- plot_data %>%
    filter(startsWith(PointType, "Deep")) %>%
    arrange(p) %>%  
    slice_head(n = top_n) %>%  
    mutate(Label = Species)
  

  to_label <- bind_rows(light_points, deep_points)
  
  
  to_label <- to_label %>%
    mutate(
      LabelColor = ifelse(Species %in% deep_red_species, "darkred", "black")
    )
  

  plot_data <- plot_data %>%
    left_join(select(to_label, Species, Label, LabelColor), by = "Species") %>%
    mutate(
      Label = coalesce(Label, ""),
      LabelColor = coalesce(LabelColor, "black")
    )
  
  # 定义颜色映射
  color_mapping <- c(
    "DeepR" = "#0065F8",        
    "LightR" = "#7FBCFF",       
    "DeepmA" = "#129990",       
    "LightmA" = "#66C2A5",      
    "DeepR2A" = "#DA6C6C",      
    "LightR2A" = "#FF9D9D",     
    "DeepGray" = "#8B0000",     
    "Gray" = "gray85"           
  )
  

  max_fc <- max(plot_data$abs_log2FC, na.rm = TRUE)
  max_fc <- if(is.finite(max_fc)) max_fc else 1
  size_breaks <- if(max_fc > 0) seq(0, ceiling(max_fc), 1) else c(0, 1)
  

  p <- ggplot(plot_data, aes(x = log2FC, y = -log10(p.adj))) +

    geom_point(
      aes(size = abs_log2FC, color = PointType),
      alpha = 0.8
    ) +

    ggrepel::geom_text_repel(
      data = subset(plot_data, Label != ""),
      aes(
        label = Label,
        color = LabelColor  
      ),
      size = 4.5,
      box.padding = 0.6,
      point.padding = 0.3,
      min.segment.length = 0.1,
      max.overlaps = 100,
      segment.color = "grey40",
      segment.size = 0.2,
      segment.alpha = 0.5,
      show.legend = FALSE,
      max.iter = 50000,
      seed = 123,
      nudge_y = 0.1,
      fontface = "italic" 
    ) +

    scale_color_manual(
      values = c(color_mapping, "darkred" = "darkred", "black" = "black"),
      guide = "none"  
    ) +
    coord_cartesian(
      xlim = c(-12, 12),
      ylim = c(0, 3.5)
    ) +

    scale_size_continuous(
      name = "|log₂FC|",
      range = c(1.5, 6),
      breaks = size_breaks,
      guide = "none"
    ) +
   
    labs(
      x = "log₂(FC)",
      y = "-log₁₀(adjusted p-value)"
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 10)),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      panel.grid.major.y = element_line(color = "grey92"),
      panel.grid.major.x = element_blank()
    )
  
  return(p)
}


dunn_results <- dunn_results %>%
  mutate(
    Comparison = case_when(
      (group1 == "R" & group2 == "mA") ~ "R_vs_mA",
      (group1 == "mA" & group2 == "R") ~ "R_vs_mA",
      (group1 == "mA" & group2 == "R2A") ~ "mA_vs_R2A",
      (group1 == "R2A" & group2 == "mA") ~ "mA_vs_R2A",
      (group1 == "R" & group2 == "R2A") ~ "R_vs_R2A",
      (group1 == "R2A" & group2 == "R") ~ "R_vs_R2A"
    )
  ) %>%
  filter(!is.na(Comparison))


volcano_R_mA <- plot_volcano("R_vs_mA", dunn_results)
volcano_mA_R2A <- plot_volcano("mA_vs_R2A", dunn_results)
volcano_R_R2A <- plot_volcano("R_vs_R2A", dunn_results)


combined_plots <- volcano_R_mA + volcano_mA_R2A + volcano_R_R2A + 
  plot_layout(ncol = 3, guides = "collect") & 
  theme(legend.position = "bottom")


ggsave("C:/RStudio/gut/code/Figure4B.png", 
       combined_plots, width = 20, height = 8, dpi = 300)
