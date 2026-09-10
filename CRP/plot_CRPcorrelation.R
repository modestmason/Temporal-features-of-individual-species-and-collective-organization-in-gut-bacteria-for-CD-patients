library(tidyverse)
library(ggtext)
library(ggrepel)
library(scales)

data <- read_csv("C:/RStudio/gut/code/data/CRP-spearman.csv") %>%
  mutate(
    Group = factor(Group, levels = c("R", "mA", "R2A")),
    log10_p = -log10(p)
  )

significant_species <- data %>%
  filter(p < 0.05) %>%
  pull(Species) %>%
  unique()

plot_data <- data %>%
  filter(Species %in% significant_species) %>%
  mutate(
    Group_Subject = factor(paste(Group, Subject, sep = "-"),
                           levels = c("R-12", "R-14", "mA-4", "R2A-1", "R2A-7", "R2A-15", "R2A-16"))
  )

darkblue_genera <- c("Bifidobacterium", "Roseburia", "Megamonas", "Blautia")
darkred_species <- "[Ruminococcus] gnavus"

species_order <- plot_data %>%
  group_by(Species) %>%
  summarise(mean_rho = mean(rho, na.rm = TRUE)) %>%
  arrange(desc(mean_rho)) %>%
  pull(Species)

species_labels <- sapply(species_order, function(species) {

  genus <- gsub(" .*", "", species)  
  
  if (any(sapply(darkblue_genera, function(g) grepl(g, genus)))) {
    return(paste0("<span style='color:#0000CC;'>*", species, "*</span>"))
  }

  else if (species == darkred_species) {
    return(paste0("<span style='color:#CC0000;'>*", species, "*</span>"))
  }

  else {
    return(paste0("*", species, "*"))
  }
})
names(species_labels) <- species_order

plot_data$Species <- factor(plot_data$Species, levels = species_order)

group_breaks <- c(
  which(levels(plot_data$Group_Subject) == "R-14") + 0.5,
  which(levels(plot_data$Group_Subject) == "mA-4") + 0.5
)


bubble_plot <- ggplot(plot_data, aes(x = Group_Subject, y = Species)) +

  geom_vline(xintercept = group_breaks, 
             color = "grey50", 
             linewidth = 1.2,
             linetype = "dashed") +

  geom_point(
    aes(size = log10_p, fill = rho),
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
    name = "Spearman's rho",
    low = "#1F509A",   
    mid = "white",     
    high = "#E38E49",  
    midpoint = 0,      
    limits = c(-1, 1), 
    breaks = c(-1, -0.5, 0, 0.5, 1)
  ) +

  labs(
    title = "Microbiota-CRP Associations",
    x = "Group - Patient ID",
    y = "Species"
  ) +

  scale_y_discrete(labels = species_labels) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold", hjust = 0.5, size = 15),
    axis.text.y = ggtext::element_markdown(face = "italic", size = 15),  
    axis.title = element_text(face = "bold", size = 20),
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
    legend.title = element_text(face = "bold", size = 15),
    legend.text = element_text(size = 9),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
  ) +

  guides(
    size = guide_legend(order = 1),
    fill = guide_colorbar(order = 2)
  )


print(bubble_plot)
ggsave("C:/RStudio/gut/code/Figure3.png", 
       plot = bubble_plot,
       width = 14, 
       height = 10, 
       dpi = 600)
