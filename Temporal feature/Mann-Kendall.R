library(Kendall)
library(tidyr)
library(purrr)
library(readr)
library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggrepel)
library(patchwork)
library(dplyr)


setwd("C:/RStudio/gut/code/data")
groups <- list(
  R = c(2, 3, 5, 10, 11, 12, 13, 14),
  mA = c(4, 6, 8, 9),
  R2A = c(1, 7, 15, 16)
)
patients <- data.frame(
  patient_id = unlist(groups),
  group = rep(names(groups), sapply(groups, length))
)
patient_data <- list()
all_species <- c()  
for (i in 1:nrow(patients)) {
  pid <- patients$patient_id[i]
  group <- patients$group[i]
  file_name <- sprintf("CD%d.csv", pid)
  
  if (file.exists(file_name)) {
    df <- read_csv(file_name, show_col_types = FALSE)
    current_species <- df$Species
    all_species <- union(all_species, current_species)
    df_long <- df %>%
      pivot_longer(
        cols = starts_with("Month"),
        names_to = "Month",
        values_to = "Abundance"
      ) %>%
      mutate(
        Month = as.numeric(gsub("Month", "", Month)),
        Patient = as.character(pid),
        Group = group
      ) %>%
      filter(!is.na(Abundance))

    patient_data[[as.character(pid)]] <- df_long
  }
}


all_data <- bind_rows(patient_data)


complete_data <- all_data %>%

  distinct(Group,Patient,Species,.keep_all=FALSE) %>%
  mutate(trend=0,tau=0,p_value=0)


for (i in 1:nrow(complete_data)) {
  group_val <- complete_data$Group[i]
  patient_val <- complete_data$Patient[i]
  species_val <- complete_data$Species[i]
  

  species_data <- all_data %>%
    filter(Group == group_val,
           Patient == patient_val,
           Species == species_val)
  

  mean_abundance <- mean(species_data$Abundance, na.rm = TRUE)
  

  if (mean_abundance > 0.001 ) {
    mk_test <- MannKendall(species_data$Abundance)
    complete_data$p_value[i]<- mk_test$sl
    if (mk_test$sl < 0.05) {

      complete_data$tau[i] <- mk_test$tau
      if (mk_test$tau > 0) {
        complete_data$trend[i] <- 1  
      } else {
        complete_data$trend[i] <- -1 
      }
    } else {
      complete_data$trend[i] <- 0  
    }
  } else {

    complete_data$trend[i] <- 0
    complete_data$tau[i] <- 0
    complete_data$p_value[i] <- 1
  }
}
write_csv(complete_data, "Trend_data.csv")
complete_data <- read_csv("C:/RStudio/gut/code/data/Trend_data.csv")


kruskal_results <- complete_data %>%
  group_by(Species) %>%
  kruskal_test(trend ~ Group) %>%  
  ungroup() %>%
  adjust_pvalue(method = "none") %>%      
  add_significance("p.adj")


significant_species <- kruskal_results %>%
  filter(p < 0.05) %>%  
  pull(Species)

calculate_all_pairwise_raw <- function(species_name) {
  data_subset <- complete_data %>%
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
    

    medians <- pair_data %>%
      group_by(Group) %>%
      summarise(
        median_tau = median(trend, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      )
    
    median1 <- medians %>% filter(Group == group1) %>% pull(median_tau)
    median2 <- medians %>% filter(Group == group2) %>% pull(median_tau)
    

    epsilon <- 1e-8  
    adjusted_median1 <- median1 + sign(median1)*epsilon
    adjusted_median2 <- median2 + sign(median2)*epsilon
    

    log2fc <- log2(abs(adjusted_median1) / abs(adjusted_median2))
    

    if (sign(median1) != sign(median2)) {
      direction <- ifelse(median1 > median2, 1, -1)
      log2fc <- log2fc * direction
    }
    

    dunn_result <- tryCatch({
      pair_data %>% dunn_test(trend ~ Group, p.adjust.method = "none")
    }, error = function(e) { return(NULL) })
    
    if (is.null(dunn_result) || nrow(dunn_result) == 0) {
      return(data.frame())
    }
    
    p_raw <- dunn_result %>%
      filter(group1 == !!group1, group2 == !!group2) %>%
      pull(p)
    

    data.frame(
      Species = species_name,
      Comparison = paste0(group1, "_vs_", group2),
      Group1 = group1,
      Group2 = group2,
      Median1 = median1,
      Median2 = median2,
      log2FC = ifelse(is.finite(log2fc), log2fc, NA),
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
  mutate(
    p.adj = p.adjust(p.raw, method = "none")  
  ) %>%
  select(Species, Comparison, Group1, Group2, Median1, Median2, log2FC, p.raw, p.adj)


final_results <- final_results %>%
  mutate(
    Direction = case_when(
      log2FC > 0 ~ "Increase",
      log2FC < 0 ~ "Decrease",
      TRUE ~ "No Change"
    ),
    Significance = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )




