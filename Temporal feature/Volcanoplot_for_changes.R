library(tidyverse)
library(rstatix)
library(ggrepel)
library(patchwork)
library(ggtext)

fluctuation_R <- read_csv("C:/RStudio/gut/code/data/fluctuation_R.csv")  
fluctuation_mA <- read_csv("C:/RStudio/gut/code/data/fluctuation_mA.csv")
fluctuation_R2A <- read_csv("C:/RStudio/gut/code/data/fluctuation_R2A.csv")

fluctuation_R$Group <- "R"
fluctuation_mA$Group <- "mA"
fluctuation_R2A$Group <- "R2A"

combined_data <- bind_rows(fluctuation_R, fluctuation_mA, fluctuation_R2A) %>%
  mutate(Group = factor(Group, levels = c("R", "mA", "R2A"))) %>%  
  pivot_longer(
    cols = -Group,
    names_to = "Species",
    values_to = "Fluctuation"
  )

filtered_species <- combined_data %>%
  group_by(Species) %>%
  summarise(
    Mean_Fluctuation = mean(Fluctuation),
    .groups = "drop"
  ) %>%
  filter(Mean_Fluctuation >= 0) %>%  
  pull(Species)

combined_data_filtered <- combined_data %>%
  filter(Species %in% filtered_species)

kruskal_results_filtered <- combined_data_filtered %>%
  group_by(Species) %>%
  kruskal_test(Fluctuation ~ Group) %>%  
  ungroup() %>%
  adjust_pvalue(method = "none") %>%      
  add_significance("p.adj")



significant_species <- kruskal_results_filtered %>%
  filter(p < 0.05) %>%  
  pull(Species)

calculate_all_pairwise_raw <- function(species_name) {
  data_subset <- combined_data %>%
    filter(Species == species_name, Group %in% c("R", "mA", "R2A"))
  
  if (n_distinct(data_subset$Group) < 2) {
    return(data.frame())
  }
  
  comparisons <- list(
    c("R", "mA"),
    c("R", "R2A"),
    c("mA", "R2A")
  )
  
  map_dfr(comparisons, function(pair) {
    group1 <- pair[1]
    group2 <- pair[2]
    
    pair_data <- data_subset %>% filter(Group %in% c(group1, group2))
    
    dunn_result <- tryCatch({
      pair_data %>% dunn_test(Fluctuation ~ Group, p.adjust.method = "none")
    },
    error = function(e) { return(NULL) }
    )
    
    if (is.null(dunn_result) || nrow(dunn_result) == 0) {
      return(data.frame())
    }
    
    p_raw <- dunn_result %>%
      filter(group1 == !!group1, group2 == !!group2) %>%
      pull(p)
    
    median1 <- median(pair_data$Fluctuation[pair_data$Group == group1], na.rm = TRUE)
    median2 <- median(pair_data$Fluctuation[pair_data$Group == group2], na.rm = TRUE)
    log2fc <- log2((median1 + 1e-10) / (median2 + 1e-10))
    
    data.frame(
      Species = species_name,
      Comparison = paste0(group1, "_vs_", group2),
      log2FC = log2fc,
      p.raw = ifelse(length(p_raw) == 1, p_raw, NA)
    )
  })
}


raw_results <- map_dfr(
  .x = significant_species,
  .f = calculate_all_pairwise_raw
) %>% 
  filter(!is.na(p.raw))

final_results <- raw_results %>%
  mutate(p.adj = p.adjust(p.raw, method = "BH")) %>%
  select(Species, Comparison, log2FC, p.raw, p.adj)



manual_labels <- list(

  "R_vs_mA" = list(
    left = c("s__Sutterella megalosphaeroides", "s__Klebsiella pneumoniae", 
             "s__Veillonella atypica", "s__Veillonella dispar", "s__Clostridium butyricum"),
    right = c("s__Alistipes shahii", "s__Paenibacillus odorifer", "s__Desulfovibrio sp. IOR2")
  ),

  "mA_vs_R2A" = list(
    left = c("s__Acidaminococcus intestini", "s__Cupriavidus pauculus", "s__Suicoccus acidiformans"),
    right = c("s__Bifidobacterium dentium", "s__Cedecea lapagei", "s__Leclercia sp. LSNIH1")
  ),
  
  "R_vs_R2A" = list(
    left = c("s__Dialister hominis", "s__Alistipes shahii", 
             "s__Dehalobacter restrictus", "s__Desulfosporosinus meridiei"),
    right = c("s__Campylobacter concisus", "s__Klebsiella pneumoniae", "s__Veillonella atypica", 
              "s__Clostridium perfringens", "s__Klebsiella oxytoca", "s__Eikenella corrodens")
  )
)


plot_volcano <- function(comparison, data) {

  darkblue_species <- c("s__Alistipes shahii")
  darkred_species <- c(
    "s__Sutterella megalosphaeroides", "s__Klebsiella pneumoniae", 
    "s__Veillonella atypica", "s__Veillonella dispar", "s__Clostridium butyricum",
    "s__Campylobacter concisus", "s__Klebsiella oxytoca", "s__Eikenella corrodens","s__Clostridium perfringens"
  )
  
  current_labels <- manual_labels[[comparison]]
  
  plot_data <- data %>%
    filter(
      Comparison == comparison,
      !is.na(p.adj),
      is.finite(log2FC)
    ) %>%
    mutate(
      abs_log2FC = abs(log2FC),
      ColorGroup = case_when(
        p.adj >= 0.1 ~ "Not Significant",
        comparison == "R_vs_mA" & log2FC > 0 ~ "R",
        comparison == "R_vs_mA" & log2FC < 0 ~ "mA",
        comparison == "mA_vs_R2A" & log2FC > 0 ~ "mA",
        comparison == "mA_vs_R2A" & log2FC < 0 ~ "R2A",
        comparison == "R_vs_R2A" & log2FC > 0 ~ "R",
        comparison == "R_vs_R2A" & log2FC < 0 ~ "R2A",
        TRUE ~ "Not Significant"
      ),

      ManualLabel = case_when(
        Species %in% current_labels$left ~ Species,
        Species %in% current_labels$right ~ Species,
        TRUE ~ ""
      ),

      ItalicLabel = ifelse(ManualLabel != "", paste0("italic('", ManualLabel, "')"), ""),

      LabelColor = case_when(
        Species %in% darkblue_species ~ "darkblue",
        Species %in% darkred_species ~ "darkred",
        TRUE ~ "black"
      )
    )
  
  color_values <- c(
    "R" = "#0065F8", 
    "mA" = "#129990", 
    "R2A" = "#DA6C6C", 
    "Not Significant" = "grey80"
  )
  
  p <- ggplot(plot_data, aes(x = log2FC, y = -log10(p.adj))) +
    geom_point(
      aes(size = abs_log2FC, color = ColorGroup),
      alpha = 0.7
    ) +
    geom_text_repel(
      data = filter(plot_data, ManualLabel != "" & log2FC < 0),
      aes(label = ItalicLabel),
      parse = TRUE,  
      color = filter(plot_data, ManualLabel != "" & log2FC < 0)$LabelColor,  # 应用颜色
      size = 4.5,
      direction = "y",
      nudge_x = -1,
      segment.curvature = -0.1,
      segment.ncp = 3,
      segment.angle = 20,
      segment.size = 0.3,
      point.padding = 0.4,
      box.padding = 0.6,
      force = 5,
      force_pull = 0.5,
      seed = 123,
      min.segment.length = 0,
      max.overlaps = 100,
      show.legend = FALSE
    ) +

    geom_text_repel(
      data = filter(plot_data, ManualLabel != "" & log2FC > 0),
      aes(label = ItalicLabel),
      parse = TRUE,
      color = filter(plot_data, ManualLabel != "" & log2FC > 0)$LabelColor,  # 应用颜色
      size = 4.5,
      direction = "y",
      nudge_x = 1,
      segment.curvature = 0.1,
      segment.ncp = 3,
      segment.angle = 160,
      segment.size = 0.3,
      point.padding = 0.4,
      box.padding = 0.6,
      force = 5,
      force_pull = 0.5,
      seed = 123,
      min.segment.length = 0,
      max.overlaps = 100,
      show.legend = FALSE
    ) +
    scale_color_manual(values = color_values, guide = 'none') +
    scale_size_continuous(name = "|log2FC|", range = c(1.5, 6), guide = "none") +
    coord_cartesian(xlim = c(-15, 15), ylim = c(0, 2)) +
    labs(x = "log₂(FC)", y = "-log₁₀(adjusted p-value)") +
    theme_classic() +
    theme(
      axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 10)),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      text = element_text(family = "")
    )
  
  return(p)
}


volcano_AB <- plot_volcano("R_vs_mA", final_results)
volcano_BC <- plot_volcano("mA_vs_R2A", final_results)
volcano_AC <- plot_volcano("R_vs_R2A", final_results)


combined_plots <- volcano_AB + volcano_BC + volcano_AC + plot_layout(ncol = 3)
ggsave("C:/RStudio/gut/code/Figure2C.png", 
       combined_plots, width = 20, height = 8, dpi = 300)  
