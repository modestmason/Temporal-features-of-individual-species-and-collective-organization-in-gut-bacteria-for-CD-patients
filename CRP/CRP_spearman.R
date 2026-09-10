
library(tidyverse)
library(coin)
library(doParallel)  

registerDoParallel(cores = 4)
#Process the data of each patient separately, and manually compile the csv file recording the results of Spearman test
full_data <- read_csv("C:/RStudio/gut/code/data/CRP_species_1.csv")

crp <- full_data[[1]]  
species_data <- full_data[, -1]  
mean_abundance <- species_data %>% 
  summarise(across(everything(), ~mean(., na.rm = TRUE))) %>%
  unlist()

valid_species <- species_data[, mean_abundance >= 0.001]
results <- map_df(colnames(valid_species), function(species_name) {

  df <- tibble(
    CRP = crp,
    Abundance = valid_species[[species_name]]
  ) %>%
    filter(complete.cases(.))  

  tryCatch({
    test_res <- spearman_test(
      CRP ~ Abundance, 
      data = df,
      distribution = approximate(nresample = 10000)
    )
    
    rho <- cor(df$CRP, df$Abundance, method = "spearman")
    p_val <- as.numeric(pvalue(test_res))

    tibble(
      Species = species_name,
      Rho = rho,
      P_value = p_val,
      N = nrow(df)
    )
  }, error = function(e) {
    message("Species ", species_name, " Failure: ", e$message)
    NULL
  })
})

results_adj <- results %>%
  filter(!is.na(P_value)) %>%
  mutate(
    Adj_P = p.adjust(P_value, method = "BH"),
    Significance = factor(case_when(
      Adj_P < 0.001 ~ "***",
      Adj_P < 0.01 ~ "**",
      Adj_P < 0.05 ~ "*",
      TRUE ~ ""
    ))
  ) %>%
  arrange(Adj_P)


