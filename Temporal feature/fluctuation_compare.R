library(rstatix)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(scales)
# Fluctuations are calcualted by meanfluc.m using matlab.
df <- data.frame(
  Group = rep(c("R", "mA", "R2A"), times = c(8, 4, 4)),
  Value = c(
    9.39029159132e-06,9.67014315852923e-06,8.03402016970377e-06,
    2.08504615732369e-05,6.77914762452736e-06,8.50920333534266e-06,
    1.12392474108170e-05,2.08000027124773e-05,
    1.77581794755876e-05,3.49303752260399e-06,1.20623110307414e-05,9.93139256945585e-06,
    3.18501327822268e-05,3.27558543012141e-05,1.66250242314647e-05,2.63148510397831e-05
  )
)


df$Group <- factor(df$Group, levels = c("R", "mA", "R2A"))

kruskal_results <- df %>%
  kruskal_test(Value ~ Group) %>%
  adjust_pvalue(method = "none") %>%
  add_significance("p.adj")

dunn_results <- df %>%
  dunn_test(Value ~ Group, p.adjust.method = "BH") 

y_max <- max(df$Value) * 1.25

group_positions <- c("R" = 1, "mA" = 2, "R2A" = 3)

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
    y.position = y_max * (1 + 0.1 * row_number()) 
  )


p <- ggplot(df, aes(x = Group, y = Value, color = Group)) +

  geom_boxplot(
    width = 0.6,
    outlier.shape = 20,
    outlier.size = 2,
    outlier.alpha = 1
  ) +

  geom_jitter(
    width = 0.15,
    size = 3,
    alpha = 1
  ) +

  scale_color_manual(values = c("#0065F8", "#129990", "#DA6C6C")) +
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
    y = "Averaged squared first-order difference",
    title = "Comparison of First-order difference",
    subtitle = paste0("Kruskal-Wallis test: χ²(", kruskal_results$df, ") = ",
                      round(kruskal_results$statistic, 2), 
                      ", p = ", format(kruskal_results$p, scientific = FALSE, digits = 7.5))
  ) +
  scale_y_continuous(
    limits = c(0, max(dunn_results$y.position) * 1.15),
    expand = expansion(mult = c(0, 0.1)),
    labels = label_number(accuracy = 0.00001) 
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
    aes(x = (xmin + xmax)/2, y = y.position * 1.02, label = p_label),
    size = 7,
    color = "black",
    vjust = 0
  )

print(p)


ggsave("C:/RStudio/gut/code/Figure2A.png", p, 
       width = 14, height = 6, dpi = 300, bg = "white")