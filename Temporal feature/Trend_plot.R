library(tidyverse)
library(ggtext)
library(ggrepel)
library(scales)

complete_data <- read_csv("C:/RStudio/gut/code/data/Trend_data.csv") %>%
  mutate(
    Group = factor(Group, levels = c("R", "mA", "R2A")),
    log10_p = -log10(p_value)
  )

groups <- list(
  R = c(2, 3, 5, 10, 11, 12, 13, 14),
  mA = c(4, 6, 8, 9),
  R2A = c(1, 7, 15, 16)
)

group_order <- c(
  paste("R", sort(groups$R), sep = "-"),
  paste("mA", sort(groups$mA), sep = "-"),
  paste("R2A", sort(groups$R2A), sep = "-")
)

trend_data <- complete_data %>%
  filter(p_value < 0.05) %>%
  mutate(
    log10_p = -log10(p_value),
    Group_Subject = factor(paste(Group, Patient, sep = "-"))
  )

species_order <- trend_data %>%
  group_by(Species) %>%
  summarise(mean_tau = mean(tau, na.rm = TRUE)) %>%
  arrange(desc(mean_tau)) %>%
  pull(Species)

darkgreen_targets <- c(
  "Alistipes",          
  "Bifidobacterium",    
  "s__Faecalibacterium prausnitzii",  
  "s__Akkermansia muciniphila",
  "s__Roseburia intestinalis"
)

darkred_targets <- c(
  "Klebsiella",        
  "Veillonella",        
  "s__Fusobacterium mortiferum",
  "s__Escherichia coli"  
)

species_labels <- sapply(species_order, function(species) {
  
  genus <- strsplit(sub("^s__", "", species), " ")[[1]][1]
  
  if(species %in% darkgreen_targets || genus %in% darkgreen_targets) {
    
    return(paste0("<span style='color:darkblue'>*", species, "*</span>"))
  } else if(species %in% darkred_targets || genus %in% darkred_targets) {
    
    return(paste0("<span style='color:darkred'>*", species, "*</span>"))
  } else {
   
    return(paste0("*", species, "*"))
  }
})
names(species_labels) <- species_order


plot_data <- complete_data %>%
  filter(Species %in% unique(trend_data$Species)) %>% 
  mutate(
    log10_p = -log10(p_value),
    Group_Subject = paste(Group, Patient, sep = "-"),
    Species = factor(Species, levels = species_order) 
  ) %>%
  mutate(Group_Subject = factor(Group_Subject, levels = group_order))


r_end <- length(groups$R)
mA_end <- r_end + length(groups$mA)
group_breaks <- c(r_end + 0.5, mA_end + 0.5)


bubble_plot <- ggplot(plot_data, aes(x = Group_Subject, y = Species)) +
  
  geom_vline(xintercept = group_breaks, 
             color = "grey50", 
             linewidth = 1.2,
             linetype = "dashed") +
  
  geom_point(
    aes(size = log10_p, fill = tau),
    data = plot_data %>% filter(p_value < 0.05),  
    shape = 21, 
    color = "black",
    stroke = 0.3
  ) +
 
  scale_size_continuous(
    name = "-log10(p-value)",
    range = c(3, 12),
    breaks = c(1.3, 2, 3),
    labels = c("p=0.05", "p=0.01", "p=0.001")
  ) +
  
  scale_fill_gradient2(
    name = "Mann-Kendall tau",
    low = "#1F509A",
    mid = "white",
    high = "#E38E49",
    midpoint = 0,
    limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1)
  ) +
  
  labs(
    title = "Temporal Trends",
    x = "Group - Patient ID",
    y = "Species"
  ) +
  
  scale_y_discrete(
    labels = species_labels[levels(plot_data$Species)]
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 15
    ),
    axis.text.y = ggtext::element_markdown(size = 10),  # Markdown渲染
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
  ) +
  
  scale_x_discrete(limits = group_order) +
  guides(
    size = guide_legend(order = 1),
    fill = guide_colorbar(order = 2)
  )


print(bubble_plot)
ggsave("C:/RStudio/gut/code/Figure1.png", 
       plot = bubble_plot,
       width = 14, 
       height = 14, 
       dpi = 600)
