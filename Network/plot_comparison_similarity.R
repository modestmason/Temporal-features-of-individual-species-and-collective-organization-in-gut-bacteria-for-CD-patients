library(rstatix)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(readr)

# The Jaccard indice of real patients calculated from jaccard_similarity_matrix (obtained from network_similarity.R)
df_real <- data.frame(
  Group = rep(c("R", "mA", "R2A"), times = c(28, 6, 6)),
  Value = c(
    0.058823529,0.063660477,0.040178571,0.056701031,0.040462428,0.045731707,0.080275229,
    0.07277628,0.055776892,0.065934066,0.037383178,0.037037037,0.055555556,0.059360731,
    0.044776119,0.049132948,0.077639752,0.059850374,0.025806452,0.028901734,0.018404908,0.053872054,
    0.041958042,0.064285714,0.034351145,0.030508475,0.041297935,0.039215686,
    0.040968343,0.061016949,0.047619048,0.04743083,0.037243948,0.080103359,
    0.073446328,0.060377358,0.087064677,0.037593985,0.052757794,0.069767442
  )
)


read_matrix_to_vector <- function(file_name) {
  as.vector(as.matrix(read_csv(file_name, col_names = FALSE, show_col_types = FALSE)))
}


df_rand <- bind_rows(
  data.frame(
    Group = "R-rand",
    Value = read_matrix_to_vector("C:/RStudio/gut/code/data/random_network_pairwise_similarity_R.csv")
  ),
  data.frame(
    Group = "mA-rand",
    Value = read_matrix_to_vector("C:/RStudio/gut/code/data/random_network_pairwise_similarity_mA.csv")
  ),
  data.frame(
    Group = "R2A-rand",
    Value = read_matrix_to_vector("C:/RStudio/gut/code/data/random_network_pairwise_similarity_R2A.csv")
  )
)

df_combined <- bind_rows(
  df_real %>% mutate(Type = "Actual"),
  df_rand %>% mutate(Type = "Random")
) %>%
  mutate(
    Group = factor(Group, 
                   levels = c("R-rand", "mA-rand","R2A-rand","R", "mA",  "R2A"),
                   labels = c("R-rand", "mA-rand",  "R2A-rand", "R", "mA","R2A")),
    Type = factor(Type, levels = c("Actual", "Random"))
  )


kruskal_results <- df_real %>%
  kruskal_test(Value ~ Group) %>%
  adjust_pvalue(method = "none") %>%
  add_significance("p.adj")


dunn_results <- df_real %>%
  dunn_test(Value ~ Group, p.adjust.method = "BH") 
y_max <- max(df_real$Value) * 1.25
group_positions <- c("R" = 4, "mA" = 5, "R2A" = 6)
dunn_results <- dunn_results %>%
  mutate(

    xmin = group_positions[group1],
    xmax = group_positions[group2],
    temp = ifelse(xmin > xmax, xmin, xmax),
    xmin = ifelse(xmin > xmax, xmax, xmin),
    xmax = temp
  ) %>%
  select(-temp) %>%
  mutate(
    p_label = ifelse(
      p < 0.001, 
      "p < 0.001",
      paste0("p = ", format(round(p, 3), nsmall = 3))
    )
  ) %>%
  arrange(xmax - xmin, xmin) %>%
  mutate(
    y.position = y_max * (1 + 0.1 * row_number()) # 每个比较增加10%的高度
  )

morandi_colors <- c(
  "R-rand" = "#7FBCFF",  
  "mA-rand" = "#66C2A5", 
  "R2A-rand" = "#FF9D9D", 
  "R" = "#0065F8",        
  "mA" = "#129990",       
  "R2A" = "#DA6C6C"       
)


p <- ggplot(df_combined, aes(x = Group, y = Value, color = Group)) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    outlier.size = 2,
    outlier.alpha = 1
  ) +
  geom_jitter(
    data = filter(df_combined, Type == "Actual"),
    width = 0.15,
    size = 3,
    alpha = 1
  ) +
  scale_color_manual(values = morandi_colors) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 15)),  
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 15)),  
    axis.text.x = element_text(size = 16, face = "bold", color = "black"),  
    axis.text.y = element_text(size = 16, face = "bold", color = "black"),  
    plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5)
  ) +
  labs(
    x = "Group",
    y = "Network Similarity",
    title = "Comparison of Network Similarity Across Groups",
    subtitle = paste0("Kruskal-Wallis test: χ²(", kruskal_results$df, ") = ",
                      round(kruskal_results$statistic, 2), 
                      ", p = ", format(kruskal_results$p, scientific = FALSE, digits = 7.5))
  ) +
  scale_y_continuous(
    limits = c(0, max(dunn_results$y.position) * 1.15), # 确保所有p值标签可见
    expand = expansion(mult = c(0, 0.1))
  )

p <- p +
  geom_segment(
    data = dunn_results,
    aes(x = xmin, xend = xmax, y = y.position, yend = y.position),
    color = "black",
    size = 1
  ) +
  geom_text(
    data = dunn_results,
    aes(x = (xmin + xmax)/2, y = y.position * 1.01, label = p_label),
    size = 7.5,
    color = "black",
    vjust = 0
  )


print(p)

ggsave("C:/RStudio/gut/code/Figure5A.png", p, 
       width = 14, height = 6, dpi = 300, bg = "white")







