
# Install required packages (uncomment if needed)
# install.packages(c("tidyverse", "lme4", "gridExtra", "ggplot2", "ggalluvial", "igraph", "ggraph", "viridis"))

# Load required libraries
library(tidyverse)
library(lme4)
library(gridExtra)
library(ggplot2)

project_dir <- "sim_output/20260914_084346_expl_measles"

## Load Data

config_exp <- 1 ## This is called within run_rStride()

exp_name <- sprintf("exp%04d", config_exp)

output_folder_name <- file.path(project_dir, exp_name)

population_snapshot_filename <- file.path(output_folder_name, dir(output_folder_name, pattern = 'snapshot.csv'))

population_snapshot <- read.table(population_snapshot_filename,header=T,sep=',') #%>%
  # mutate(immunity_status = ifelse(immunity_status == "immune", 1, 0))

###########################

if(1 == 0){
  population_snapshot_output <- analyze_immunity_clustering(population_snapshot)
}

analyze_immunity_clustering <- function(df){
  # Convert immunity_status to binary
  df$immune_binary <- as.numeric(df$immunity_status == "immune")
  
  head(df)
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("CLUSTERING ANALYSIS: Immune/Susceptible Status\n")
  cat(rep("=", 70), "\n\n", sep = "")
  
  cat("Dataset Summary:\n")
  cat("  Total individuals:", nrow(df), "\n")
  cat("  Immune:", sum(df$immune_binary), "\n")
  cat("  Susceptible:", sum(!df$immune_binary), "\n")
  cat("  Households:", n_distinct(df[["household_id"]]), "\n")
  cat("  Communities:", n_distinct(df[["community_weekday_id"]]), "\n\n")
  
  # ============================================================================
  # CHI-SQUARE TEST OF INDEPENDENCE
  # ============================================================================
  
  # Significant p-value suggests immunity status is NOT independent of household/community. (i.e., clustering exists)
  
  # Test for household clustering
  cat("By Household")
  household_table <- table(df$household_id, df$immune_binary)
  chi2_household <- chisq.test(household_table)
  print(chi2_household)
  
  # Test for weekday community clustering
  cat("By Weekday Community")
  community_table <- table(df$community_weekday_id, df$immune_binary)
  chi2_community <- chisq.test(community_table)
  print(chi2_community)
  
  # ============================================================================
  # EFFECT SIZE - CRAMÉR'S V (for chi-square tests)
  # ============================================================================
  
  # (V > 0.1 = small effect, V > 0.3 = medium, V > 0.5 = large)
  
  cramersV <- function(chi2_stat, n, min_dim) {
    sqrt(chi2_stat / (n * (min_dim - 1)))
  }
  
  # For household
  n_households <- nrow(household_table)
  v_household <- cramersV(chi2_household$statistic, nrow(df), 
                          min(ncol(household_table), nrow(household_table)))
  cat("Household Cramér's V:", v_household, "\n")
  
  # For community
  n_communities <- nrow(community_table)
  v_community <- cramersV(chi2_community$statistic, nrow(df), 
                          min(ncol(community_table), nrow(community_table)))
  cat("Community Cramér's V:", v_community, "\n")
  
  # ========================================================================
  # MIXED-EFFECTS MODELS
  # ========================================================================
  
  null_model <- glm(immune_binary ~ 1, family = binomial, data = df)
  
  household_model <- glmer(immune_binary ~ 1 + (1 | household_id), 
                           family = binomial, data = df)
  
  community_model <- glmer(immune_binary ~ 1 + (1 | community_weekday_id), 
                           family = binomial, data = df)
  
  both_model <- glmer(immune_binary ~ 1 + (1 | household_id) + (1 | community_weekday_id), 
                      family = binomial, data = df)
  
  # Extract variance components
  vc_household <- VarCorr(household_model)
  vc_community <- VarCorr(community_model)
  vc_both <- VarCorr(both_model)
  
  var_household_between <- as.numeric(vc_household[[1]][1])
  var_community_between <- as.numeric(vc_community[[1]][1])
  var_household_both <- as.numeric(vc_both[[1]][1])
  var_community_both <- as.numeric(vc_both[[2]][1])
  
  residual_var_logistic <- pi^2 / 3
  
  # Calculate ICC
  icc_household <- var_household_between / (var_household_between + residual_var_logistic)
  icc_community <- var_community_between / (var_community_between + residual_var_logistic)
  icc_household_adj <- var_household_both / (var_household_both + residual_var_logistic)
  icc_community_adj <- var_community_both / (var_community_both + residual_var_logistic)
  
  # ========================================================================
  # SUMMARY STATISTICS BY GROUP
  # ========================================================================
  
  household_summary <- df %>%
    group_by(household_id) %>%
    summarise(
      n = n(),
      pct_immune = mean(immune_binary),
      .groups = 'drop'
    ) %>%
    arrange(desc(pct_immune))
  
  community_summary <- df %>%
    group_by(community_weekday_id) %>%
    summarise(
      n = n(),
      pct_immune = mean(immune_binary),
      .groups = 'drop'
    ) %>%
    arrange(desc(pct_immune))
  
  # ========================================================================
  # CREATE RESULTS LIST
  # ========================================================================
  
  results <- list(
    # Raw data
    data = df,
    
    # Statistical test results
    statistics = list(
      chi2_household = chi2_household,
      chi2_community = chi2_community,
      cramers_v_household = v_household,
      cramers_v_community = v_community,
      icc_household = icc_household,
      icc_community = icc_community,
      icc_household_adjusted = icc_household_adj,
      icc_community_adjusted = icc_community_adj
    ),
    
    # Fitted models
    models = list(
      null_model = null_model,
      household_model = household_model,
      community_model = community_model,
      both_model = both_model
    ),
    
    # Summary tables
    summaries = list(
      household = household_summary,
      community = community_summary
    ),
    
    # AIC comparison
    aic_comparison = tibble(
      Model = c("Null", "Household", "Community", "Both"),
      AIC = c(AIC(null_model), AIC(household_model), AIC(community_model), AIC(both_model)),
      Delta_AIC = c(0, 
                    AIC(household_model) - AIC(null_model),
                    AIC(community_model) - AIC(null_model),
                    AIC(both_model) - AIC(null_model))
    ),
    
    # Column names used (for visualization functions)
    column_names = list(
      outcome_col = outcome_col,
      household_col = household_col,
      community_col = community_col
    )
  )
  
  class(results) <- c("immunity_clustering", "list")
  
  if (verbose) {
    cat("\n✓ Analysis complete!\n")
    cat("Use print(results) for summary, or results$statistics for details\n\n")
  }
  
  return(results)
}

# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

plot_household_distribution <- function(results, bins = 50) {
  
  household_summary <- results$summaries$household
  
  p <- ggplot(household_summary, aes(x = pct_immune)) +
    geom_histogram(bins = bins, fill = "#2b83ba", alpha = 0.7, color = "black") +
    geom_vline(aes(xintercept = mean(pct_immune)), color = "red", linetype = "dashed", size = 1) +
    theme_minimal() +
    labs(
      title = "Distribution of Immunity Rates Across Households",
      subtitle = paste0("Red line = mean (", 
                        round(mean(household_summary$pct_immune) * 100, 1), "%)"),
      x = "Proportion Immune",
      y = "Number of Households"
    ) +
    scale_x_continuous(labels = scales::percent)
  
  return(p)
}

plot_household_size_effect <- function(results) {
  
  df <- results$data
  household_col <- results$column_names$household_col
  
  size_analysis <- df %>%
    group_by(!!sym(household_col)) %>%
    mutate(household_size = n(), immune_pct = mean(immune_binary)) %>%
    group_by(household_size) %>%
    summarise(
      n_households = n_distinct(household_id),
      mean_immune = mean(immune_pct),
      sd_immune = sd(immune_pct),
      .groups = 'drop'
    ) %>%
    filter(household_size <= 15)  # Limit to reasonable sizes for clarity
  
  p <- ggplot(size_analysis, aes(x = household_size, y = mean_immune)) +
    geom_point(aes(size = n_households), alpha = 0.6, color = "#1b9e77") +
    geom_line(color = "#1b9e77", alpha = 0.5) +
    geom_errorbar(aes(ymin = mean_immune - sd_immune, ymax = mean_immune + sd_immune),
                  width = 0.2, alpha = 0.3) +
    theme_minimal() +
    theme(legend.position = "right") +
    labs(
      title = "Immunity Rate by Household Size",
      x = "Household Size (number of people)",
      y = "Mean Proportion Immune",
      size = "# of Households"
    ) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(breaks = 1:15)
  
  return(p)
}

plot_community_distribution <- function(results, bins = 30) {
  
  community_summary <- results$summaries$community
  
  p <- ggplot(community_summary, aes(x = pct_immune)) +
    geom_histogram(bins = bins, fill = "#d95f02", alpha = 0.7, color = "black") +
    geom_vline(aes(xintercept = mean(pct_immune)), color = "red", linetype = "dashed", size = 1) +
    theme_minimal() +
    labs(
      title = "Distribution of Immunity Rates Across Communities",
      subtitle = paste0("Red line = mean (", 
                        round(mean(community_summary$pct_immune) * 100, 1), "%)"),
      x = "Proportion Immune",
      y = "Number of Communities"
    ) +
    scale_x_continuous(labels = scales::percent)
  
  return(p)
}

plot_bubble_chart <- function(results) {
  
  df <- results$data
  household_col <- results$column_names$household_col
  community_col <- results$column_names$community_col
  
  bubble_data <- df %>%
    group_by(household_id, community_weekday_id) %>%
    summarise(
      n = n(),
      pct_immune = mean(immune_binary),
      .groups = 'drop'
    ) %>%
    sample_n(min(nrow(.), 2000))  # Subsample if too many points
  
  p <- ggplot(bubble_data, aes(x = community_weekday_id, y = household_id)) +
    geom_point(aes(size = n, fill = pct_immune), alpha = 0.5, shape = 21, color = "black") +
    scale_fill_gradient2(low = "#d73027", mid = "#fee090", high = "#1a9850",
                        limits = c(0, 1),
                        labels = scales::percent) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = "Household-Community Cross-Clustering",
      subtitle = "Bubble size = number of individuals; Color = % immune",
      x = "Community ID",
      y = "Household (rows)",
      size = "N people",
      fill = "% Immune"
    )
  
  return(p)
}






