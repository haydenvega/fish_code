library(tidyverse)

# Load required libraries
library(tidyverse)
library(broom)
library(ggpubr)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(tidyverse)
library(lme4)
library(broom.mixed)
library(emmeans)


df <- read_csv(here("data", "cleaned_skel.csv"))
df$delta_mass <- df$f.weight-df$i.weight


#look at how delta mass varies with surface area

##############
#plot initial weight agains final weight for unwounded corals 
# Filter data for unwounded corals no fish to think about allometry 
##############


# Calculate change in mass (delta_mass)
df <- df %>%
  mutate(delta_mass = f.weight - i.weight)

# Filter data for unwounded corals and no fish treatment
allometry_df <- df %>%
  filter(wound_raw == "not wounded", treatment == "no fish") %>%
  mutate(
    log_i_weight = log(i.weight),
    log_f_weight = log(f.weight)
  )

# Fit separate regression models to extract species-specific b-values with confidence intervals
b_values_df <- allometry_df %>%
  group_by(species) %>%
  do({
    model <- lm(log_f_weight ~ log_i_weight, data = .)
    tidy_model <- broom::tidy(model, conf.int = TRUE)
    intercept_row <- tidy_model[1, ]
    slope_row <- tidy_model[2, ]
    tibble(
      intercept = intercept_row$estimate,
      intercept_conf.low = intercept_row$conf.low,
      intercept_conf.high = intercept_row$conf.high,
      b_value = slope_row$estimate,
      b_conf.low = slope_row$conf.low,
      b_conf.high = slope_row$conf.high
    )
  })

# Merge regression parameters back to the main data
allometry_df <- allometry_df %>%
  left_join(b_values_df, by = "species")

# Visualization: Initial vs. Final Weight with species-specific regression lines
allometry_df %>%
  ggplot(aes(x = log_i_weight, y = log_f_weight)) +
  geom_point(alpha = 0.7, color = "steelblue") +
  geom_abline(aes(intercept = intercept, slope = b_value), color = "darkred", size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  facet_wrap(~species, scales = "free") +
  labs(
    title = "Species-Specific Allometric Scaling (Unwounded, No Fish Treatment)",
    x = "Log Initial Weight",
    y = "Log Final Weight"
  ) +
  theme_bw(base_size = 14)

# Display species-specific b-values with confidence intervals
print(b_values_df)

# Visualization: delta_mass vs. Surface Area (unwounded, no fish)
ggplot(allometry_df, aes(x = CSA_cm2, y = delta_mass)) +
  geom_point(alpha = 0.7, color = "forestgreen") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  facet_wrap(~species, scales = "free") +
  labs(
    title = "Coral Mass Change vs. Surface Area (Unwounded, No Fish Treatment)",
    x = "Surface Area (cm²)",
    y = "Change in Mass (delta_mass)"
  ) +
  theme_classic(base_size = 14)





##############
#plot initial weight agains final weight for unwounded corals 
# Filter data for unwounded corals no fish to think about allometry 
##############


# Calculate the ratio of final to initial mass as a measure of growth
df <- df %>%
  mutate(growth_ratio = f.weight / i.weight)

# Visualization: Growth Ratio by Treatment and Wound Status
ggplot(df, aes(x = treatment, y = growth_ratio, fill = treatment)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~species + wound_raw, scales = "free") +
  labs(
    title = "Effect of Fish Treatment and Wounding on Coral Growth Ratio",
    x = "Treatment",
    y = "Final/Initial Mass Ratio"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

# Statistical analysis: Linear model to assess effects
model_results <- df %>%
  nest(data = -species) %>%
  mutate(
    model = map(data, ~lm(growth_ratio ~ treatment * wound_raw, data = .x)),
    summary = map(model, tidy)
  ) %>%
  select(species, summary) %>%
  unnest(summary)

# Display model results
print(model_results)



# Calculate the ratio of final to initial mass as a measure of growth
df <- df %>%
  mutate(growth_ratio = f.weight / i.weight)

# Statistical analysis: summarize growth ratio means and confidence intervals
summary_df <- df %>%
  group_by(species, treatment, wound_raw) %>%
  summarize(
    mean_growth = mean(growth_ratio, na.rm = TRUE),
    se_growth = sd(growth_ratio, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    ci_lower = mean_growth - qt(0.975, df = n()-1) * se_growth,
    ci_upper = mean_growth + qt(0.975, df = n()-1) * se_growth
  )

# Visualization: Interaction plot with means and confidence intervals
ggplot(summary_df, aes(x = treatment, y = mean_growth, color = wound_raw, group = wound_raw)) +
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  geom_line(position = position_dodge(width = 0.3), linewidth = 1) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2,
                position = position_dodge(width = 0.3)) +
  facet_wrap(~species, scales = "free") +
  labs(
    title = "Interaction Effect of Fish Treatment and Wounding on Coral Growth",
    x = "Treatment",
    y = "Mean Growth Ratio (Final/Initial Mass)",
    color = "Wound Status"
  ) +
  theme_bw(base_size = 14)





#######
##plot differetn species separately 
######
library(tidyverse)
library(cowplot)

# Define distinct and visually appealing colors for each species
colors_acropora <- c("#0072B2", "#D55E00")  # Blue and Orange for Acropora
colors_porites <- c("#56B4E9", "#009E73")   # Light Blue and Green for Porites
custom_shapes <- c(1, 16)                  # Solid for "wounded", Open for "not wounded"

# Custom theme for consistent style
custom_theme <- theme_classic(base_size = 14) +
  theme(
    legend.position = "right",             # Place legend on the right
    legend.box = "vertical",               # Arrange legend items vertically
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# Plot for Acropora
plot_acropora <- ggplot(df %>% filter(species == "Acropora"), 
                        aes(x = CSA_cm2, y = delta_mass, color = treatment, shape = wound_raw)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, aes(linetype = interaction(treatment, wound_raw)), linewidth = 1) +
  labs(
    title = "Acropora",
    x = "Surface Area (cm²)",
    y = "Change in Mass (g)",
    color = "Treatment",
    shape = "Wounding",
    linetype = "Treatment & Wounding"
  ) +
  scale_color_manual(values = colors_acropora) +
  scale_shape_manual(values = custom_shapes) +
  custom_theme

# Plot for Porites
plot_porites <- ggplot(df %>% filter(species == "Porites"), 
                       aes(x = CSA_cm2, y = delta_mass, color = treatment, shape = wound_raw)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, aes(linetype = interaction(treatment, wound_raw)), linewidth = 1) +
  labs(
    title = "Porites",
    x = "Surface Area (cm²)",
    y = "Change in Mass (g)",
    color = "Treatment",
    shape = "Wounding",
    linetype = "Treatment & Wounding"
  ) +
  scale_color_manual(values = colors_porites) +
  scale_shape_manual(values = custom_shapes) +
  custom_theme

# Combine the two plots into a single figure
final_plot <- plot_grid(
  plot_acropora,
  plot_porites,
  ncol = 2, 
  labels = c("A", "B"), # Panel labels
  label_fontface = "bold",
  align = "v" # Align vertically for better composition
)

# Save the plot as a high-resolution figure
ggsave("publication_quality_plots_legend_right.png", plot = final_plot, width = 16, height = 7, dpi = 300)

# Display the final plot
print(final_plot)




#######################################
##LMER analysis of all species pooled##
#######################################


# Load necessary libraries
library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(MuMIn)

# Fit linear mixed-effects models with decreasing complexity (including random tank effect)
full_model <- lmer(delta_mass ~ CSA_cm2 * treatment * wound_raw * species + (1 | tank), data = df)
model_no_four_way <- lmer(delta_mass ~ CSA_cm2 * treatment * wound_raw + CSA_cm2 * species + treatment * species + wound_raw * species + (1 | tank), data = df)
model_no_three_way <- lmer(delta_mass ~ CSA_cm2 * treatment + CSA_cm2 * wound_raw + CSA_cm2 * species + treatment * wound_raw + treatment * species + wound_raw * species + (1 | tank), data = df)
model_no_interactions <- lmer(delta_mass ~ CSA_cm2 + treatment + wound_raw + species + (1 | tank), data = df)

# Compare models using AICc (corrected for small sample size)
library(MuMIn)
aicc_values <- AICc(full_model, model_no_four_way, model_no_three_way, model_no_interactions)
print(aicc_values)

# Select the best model based on lowest AICc
best_model_name <- rownames(aicc_values)[which.min(aicc_values$AICc)]
best_model <- switch(best_model_name,
                     full_model = full_model,
                     model_no_four_way = model_no_four_way,
                     model_no_three_way = model_no_three_way,
                     model_no_interactions = model_no_interactions)

# Summary of the best model
cat("\n### Summary of the Best Model (Selected by AICc) ###\n")
print(summary(best_model))

# Generate ANOVA table for the best model
cat("\n### ANOVA Table ###\n")
anova_table <- anova(best_model)
print(anova_table)

# Type II ANOVA for marginal p-values
library(car)
cat("\n### Type II ANOVA ###\n")
type_ii_anova <- Anova(best_model, type = "II", test.statistic = "F")
print(type_ii_anova)

# Post-hoc pairwise comparisons if interactions remain
library(emmeans)
if (any(grepl(":", names(fixef(best_model))))) {
  cat("\n### Post-Hoc Pairwise Comparisons ###\n")
  posthoc <- emmeans(best_model, pairwise ~ treatment * wound_raw * species)
  print(posthoc)
}







# Create a new data frame for predictions within observed ranges
prediction_data <- df %>%
  group_by(treatment, wound_raw, species, tank) %>%
  summarize(
    CSA_min = min(CSA_cm2),
    CSA_max = max(CSA_cm2),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  do({
    data.frame(
      CSA_cm2 = seq(.$CSA_min, .$CSA_max, length.out = 100),
      treatment = .$treatment,
      wound_raw = .$wound_raw,
      species = .$species,
      tank = .$tank
    )
  }) %>%
  ungroup()

# Add predictions to the new data frame
prediction_data <- prediction_data %>%
  mutate(
    predicted = predict(best_model, newdata = ., re.form = NULL) # Predicted values
  )

# Plot with prediction lines limited to observed range
csa_lmer<-ggplot(df, aes(x = CSA_cm2, y = delta_mass, color = treatment, shape = wound_raw)) +
  geom_point(size = 3, alpha = 0.8) +  # Raw data points
  geom_line(
    data = prediction_data,
    aes(y = predicted, group = interaction(treatment, wound_raw, tank)),
    linetype = "solid", size = 1
  ) +  # Linear fits limited to observed range
  facet_wrap(~ species, labeller = labeller(species = c("Acropora" = "Acropora", "Porites" = "Porites"))) +
  labs(
    title = "Effect of Surface Area on Coral Growth (Delta Mass)",
    subtitle = "Linear regression fits accounting for tank random effects (limited to observed range)",
    x = "Surface Area (cm²)",
    y = "Change in Mass (g)",
    color = "Treatment",
    shape = "Wounding"
  ) +
  scale_color_manual(values = c("fish" = "#0072B2", "no fish" = "#D55E00")) +  # Distinct colors
  scale_shape_manual(values = c("not wounded" = 16, "wounded" = 1)) +  # Solid/open shapes
  theme_classic(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12)
  )

ggsave("plots/coral_growth_plot.png", plot = csa_lmer, width = 10, height = 7, dpi = 300)


# First, calculate the mean CSA_cm2 from your data for reference
mean_CSA <- mean(df$CSA_cm2, na.rm = TRUE)
# Calculate mean CSA_cm2
mean_CSA <- mean(df$CSA_cm2, na.rm = TRUE)

# Reference dataframe at average CSA
reference_df <- expand.grid(
  CSA_cm2 = mean_CSA,
  treatment = c("no fish", "fish"),
  wound_raw = c("not wounded", "wounded"),
  species = unique(df$species)
)

# Corrected prediction (ignoring random effect explicitly)
reference_df$predicted_growth <- predict(best_model, newdata = reference_df, re.form = NA)

# View the corrected predictions clearly
print(reference_df)



# Calculate mean predicted growth per treatment
treatment_means <- reference_df %>%
  group_by(treatment) %>%
  summarize(mean_growth = mean(predicted_growth))

print(treatment_means)

# Clearly calculate the percentage increase (fish vs no fish)
fish_effect_percentage <- treatment_means %>%
  pivot_wider(names_from = treatment, values_from = mean_growth) %>%
  mutate(
    growth_diff = fish - `no fish`,
    percent_increase = (growth_diff / `no fish`) * 100
  )

# Show explicitly
print(fish_effect_percentage)






################################################
#using LMER to analyze the data separately #####
################################################


# Filter data by species
df_acropora <- df %>% filter(species == "Acropora")
df_porites <- df %>% filter(species == "Porites")

# Define a function to perform model selection and analysis
analyze_species <- function(data, species_name) {
  cat("\n### Analyzing", species_name, "###\n")
  
  # Fit models with different levels of interactions
  full_model <- lmer(delta_mass ~ CSA_cm2 * treatment * wound_raw + (1 | tank), data = data)
  model_no_three_way <- lmer(delta_mass ~ CSA_cm2 * treatment + CSA_cm2 * wound_raw + treatment * wound_raw + (1 | tank), data = data)
  model_no_two_way <- lmer(delta_mass ~ CSA_cm2 + treatment + wound_raw + (1 | tank), data = data)
  
  # Compare models using AIC
  aic_values <- AIC(full_model, model_no_three_way, model_no_two_way)
  print(aic_values)
  
  # Select the best model based on AIC
  best_model_name <- rownames(aic_values)[which.min(aic_values$AIC)]
  best_model <- get(best_model_name)
  
  # Output summary of the best model
  cat("\n### Summary of the Best Model ###\n")
  print(summary(best_model))
  
  # Generate ANOVA table for the best model
  cat("\n### ANOVA Table ###\n")
  anova_table <- anova(best_model)
  print(anova_table)
  
  # Type II ANOVA for marginal p-values
  cat("\n### Type II ANOVA ###\n")
  type_ii_anova <- Anova(best_model, type = "II", test.statistic = "F")
  print(type_ii_anova)
  
  # Perform post-hoc pairwise comparisons if interactions remain in the best model
  if (any(grepl(":", rownames(summary(best_model)$coefficients)))) {
    cat("\n### Post-Hoc Pairwise Comparisons ###\n")
    posthoc <- emmeans(best_model, pairwise ~ treatment * wound_raw)
    print(posthoc)
  }
}

# Analyze Acropora
analyze_species(df_acropora, "Acropora")

# Analyze Porites
analyze_species(df_porites, "Porites")



















#######################################
## Fixed-effects analysis (species pooled)
#######################################

# Load libraries
library(tidyverse)
library(MuMIn)  # For AICc calculation
library(emmeans)

# Fit linear models (no random effects)
full_model <- lm(delta_mass ~ CSA_cm2 * treatment * wound_raw * species, data = df)
model_no_four_way <- lm(delta_mass ~ CSA_cm2 * treatment * wound_raw + CSA_cm2 * species + treatment * species + wound_raw * species, data = df)
model_no_three_way <- lm(delta_mass ~ CSA_cm2 * treatment + CSA_cm2 * wound_raw + CSA_cm2 * species + treatment * wound_raw + treatment * species + wound_raw * species, data = df)
model_no_interactions <- lm(delta_mass ~ CSA_cm2 + treatment + wound_raw + species, data = df)

# Compare models using AICc explicitly
library(MuMIn)
aicc_values <- AICc(full_model, model_no_four_way, model_no_three_way, model_no_interactions)
print(aicc_values)

# Clearly select best model based on lowest AICc
best_model_name <- rownames(aicc_values)[which.min(aicc_values$AICc)]
best_model <- switch(best_model_name,
                     full_model = full_model,
                     model_no_four_way = model_no_four_way,
                     model_no_three_way = model_no_three_way,
                     model_no_interactions = model_no_interactions)

# Summary of the best model
cat("\n### Summary of the Best Model (Selected by AICc) ###\n")
print(summary(best_model))

# ANOVA table
cat("\n### ANOVA Table ###\n")
print(anova(best_model))

# Post-hoc pairwise comparisons if interactions remain
if (any(grepl(":", names(coef(best_model))))) {
  cat("\n### Post-Hoc Pairwise Comparisons ###\n")
  posthoc <- emmeans(best_model, pairwise ~ treatment * wound_raw * species)
  print(posthoc)
}

################################################
# Separate species analyses
################################################

analyze_species <- function(data, species_name) {
  cat("\n### Analyzing", species_name, "###\n")
  
  # Fit models with decreasing complexity
  full_model <- lm(delta_mass ~ CSA_cm2 * treatment * wound_raw, data = data)
  model_no_three_way <- lm(delta_mass ~ CSA_cm2 * treatment + CSA_cm2 * wound_raw + treatment * wound_raw, data = data)
  model_no_two_way <- lm(delta_mass ~ CSA_cm2 + treatment + wound_raw, data = data)
  
  # Compare models using AIC
  aic_values <- AIC(full_model, model_no_three_way, model_no_two_way)
  print(aic_values)
  
  # Select the best model based on lowest AIC
  best_model_name <- rownames(aic_values)[which.min(aic_values$AIC)]
  best_model <- switch(best_model_name,
                       full_model = full_model,
                       model_no_three_way = model_no_three_way,
                       model_no_two_way = model_no_two_way)
  
  # Summary of the best model
  cat("\n### Summary of the Best Model ###\n")
  print(summary(best_model))
  
  # Generate ANOVA table
  cat("\n### ANOVA Table ###\n")
  print(anova(best_model))
  
  # Post-hoc comparisons if interactions remain
  if (any(grepl(":", names(coef(best_model))))) {
    cat("\n### Post-Hoc Pairwise Comparisons ###\n")
    posthoc <- emmeans(best_model, pairwise ~ treatment * wound_raw)
    print(posthoc)
  }
  
  # Plot estimated means and confidence intervals
  emm <- emmeans(best_model, ~ treatment * wound_raw)
  plot_data <- as.data.frame(emm)
  
  ggplot(plot_data, aes(x = treatment, y = emmean, color = wound_raw, group = wound_raw)) +
    geom_point(position = position_dodge(0.3), size = 3) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(0.3)) +
    geom_line(position = position_dodge(0.3)) +
    labs(
      title = paste("Delta Mass by Treatment and Wounding for", species_name),
      x = "Treatment",
      y = "Estimated Mean Delta Mass (g)",
      color = "Wound Status"
    ) +
    theme_bw(base_size = 14) -> p
  
  print(p)
  
  # Display summarized results as table
  cat("\n### Estimated Means and Confidence Intervals ###\n")
  print(plot_data)
}

# Analyze Acropora
analyze_species(df %>% filter(species == "Acropora"), "Acropora")

# Analyze Porites
analyze_species(df %>% filter(species == "Porites"), "Porites")












# Import data
# df <- read_csv("data/cleaned_skel.csv")

# Calculate growth metric (final/initial mass)
df <- df %>%
  mutate(growth_ratio = f.weight / i.weight)

######################################
## Fixed-effects analysis (species pooled)
#######################################

# Fit linear models with decreasing complexity (no random effects)
full_model <- lm(growth_ratio ~ treatment * wound_raw * species, data = df)
model_no_three_way <- lm(growth_ratio ~ treatment * wound_raw + treatment * species + wound_raw * species, data = df)
model_no_interactions <- lm(growth_ratio ~ treatment + wound_raw + species, data = df)

# Compare models using AIC
aic_values <- AIC(full_model, model_no_three_way, model_no_interactions)
print(aic_values)

# Select the best model based on lowest AIC
best_model_name <- rownames(aic_values)[which.min(aic_values$AIC)]
best_model <- switch(best_model_name,
                     full_model = full_model,
                     model_no_three_way = model_no_three_way,
                     model_no_interactions = model_no_interactions)

# Summary of the best model
cat("\n### Summary of Best Model (Species Pooled) ###\n")
print(summary(best_model))

# ANOVA table
cat("\n### ANOVA Table ###\n")
print(anova(best_model))

# Post-hoc tests if interactions exist
if (any(str_detect(names(coef(best_model)), ":"))) {
  cat("\n### Post-Hoc Comparisons ###\n")
  posthoc <- emmeans(best_model, pairwise ~ treatment * wound_raw * species)
  print(posthoc)
}




# Generate estimated marginal means from the best model
emm <- emmeans(best_model, ~ treatment + wound_raw + species)

# Convert emmeans output to a dataframe
emm_df <- as.data.frame(emm)

# Plot raw data points and model predictions
ggplot() +
  # Raw data points
  geom_jitter(data = df, 
              aes(x = treatment, y = growth_ratio, color = wound_raw),
              width = 0.1, height = 0, alpha = 0.5, size = 2) +
  
  # Model predictions (mean ± 95% CI)
  geom_point(data = emm_df, 
             aes(x = treatment, y = emmean, group = wound_raw, color = wound_raw), 
             position = position_dodge(width = 0.5), size = 4, shape = 18) +
  
  geom_errorbar(data = emm_df, 
                aes(x = treatment, ymin = lower.CL, ymax = upper.CL, color = wound_raw), 
                width = 0.2, position = position_dodge(width = 0.5), linewidth = 1) +
  
  facet_wrap(~species,scales="free") +
  labs(
    title = "Model Estimates vs Raw Data (Growth Ratio)",
    subtitle = "Main effects only: Treatment, Wounding, and Species",
    y = "Growth Ratio (Final / Initial Mass)",
    x = "Treatment",
    color = "Wound Status"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )


ggsave("plots/coral_growth_model_vs_data_pub_quality1.png", width = 12, height = 8, dpi = 600)



# Assuming best_model and data (df, emm_df) are already prepared
# Set treatment order explicitly
df$treatment <- factor(df$treatment, levels = c("no fish", "fish"))
emm_df$treatment <- factor(emm_df$treatment, levels = c("no fish", "fish"))

# Create a publication-quality plot
ggplot() +
  # Raw data points
  geom_jitter(data = df, 
              aes(x = treatment, y = growth_ratio, fill = wound_raw),
              width = 0.15, shape = 21, alpha = 0.6, color = "black", size = 2.5) +
  
  # Model estimates (means)
  geom_point(data = emm_df, 
             aes(x = treatment, y = emmean, fill = wound_raw, alpha= 0.6,color = wound_raw, group = wound_raw), 
             position = position_dodge(width = 0.5), size = 4, shape = 23, stroke = 1.2) +
  
  # Error bars (95% CI)
  geom_errorbar(data = emm_df, 
                aes(x = treatment, ymin = lower.CL, ymax = upper.CL, , alpha= 0.6,color = wound_raw, group = wound_raw), 
                width = 0.2, position = position_dodge(width = 0.5), linewidth = 0.8) +
  
  # Vertical facet by species
  facet_wrap(~species, ncol = 1, labeller = labeller(species = label_value),scales="free") +
  
  # Labels
  labs(
    title = "Effect of Treatment and Wounding on Coral Growth",
    subtitle = "Growth ratio (final/initial mass) with model estimates ± 95% CI",
    y = "Growth Ratio (Final / Initial Mass)",
    x = "Treatment",
    fill = "Wound Status",
    color = "Wound Status"
  ) +
  
  # Colorblind-friendly palette for both fill and color
  scale_fill_manual(values = c("not wounded" = "#56B4E9", "wounded" = "#E69F00")) +
  scale_color_manual(values = c("not wounded" = "#56B4E9", "wounded" = "#E69F00")) +
  
  # Classic minimal theme for publication
  theme_classic(base_size = 16) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14)
  )

# Save high-resolution figure
ggsave("plots/coral_growth_model_vs_data_pub_quality.png", width = 8, height = 10, dpi = 600)





library(tidyverse)
library(broom)

# Extract tidy summary with confidence intervals
effects_raw <- tidy(model_no_interactions, conf.int = TRUE)

# Extract the exact terms for fish treatment and wounding
effects_of_interest <- effects_raw %>%
  filter(term %in% c("treatmentfish", "wound_rawwounded")) %>%
  select(term, estimate, conf.low, conf.high, p.value)

# Clearly convert estimates to percentage differences
effects_percentage <- effects_of_interest %>%
  mutate(
    percent_change = estimate * 100,
    percent_low = conf.low * 100,
    percent_high = conf.high * 100
  ) %>%
  select(term, percent_change, percent_low, percent_high, p.value)

# Clear output table for reporting
effects_percentage %>%
  mutate(term = recode(term, 
                       "treatmentfish" = "Fish Treatment (compared to No Fish)",
                       "wound_rawwounded" = "Wounding Effect")) %>%
  rename(
    "Effect" = term,
    "Percent Change (%)" = percent_change,
    "Lower 95% CI (%)" = percent_low,
    "Upper 95% CI (%)" = percent_high,
    "P-value" = p.value
  ) %>%
  print()







##############################
# Separate species analysis
##############################

analyze_species <- function(data, species_name) {
  cat("\n## Analysis for", species_name, "##\n")
  
  # Fit models
  model_full <- lm(growth_ratio ~ treatment * wound_raw, data = data)
  model_simple <- lm(growth_ratio ~ treatment + wound_raw, data = data)
  
  # Compare AIC
  aic_species <- AIC(model_full, model_simple)
  print(aic_species)
  
  # Select best model explicitly
  best_species_model <- if (aic_species$AIC[1] < aic_species$AIC[2]) model_full else model_simple
  
  cat("\n### Best Model for", unique(data$species), "###\n")
  print(summary(best_model_species <- best_species_model))
  
  cat("\n### ANOVA Table ###\n")
  print(anova(best_model_species))
  
  # Post-hoc comparisons if interactions exist
  if (any(str_detect(names(coef(best_model_species)), ":"))) {
    cat("\n### Post-Hoc Comparisons ###\n")
    posthoc <- emmeans(best_model_species, pairwise ~ treatment * wound_raw)
    print(posthoc)
  }
  
  # Plot predicted values with confidence intervals
  emm <- emmeans(best_model_species, ~ treatment * wound_raw)
  plot_data <- as.data.frame(emm)
  
  ggplot(plot_data, aes(x = treatment, y = emmean, color = wound_raw, group = wound_raw)) +
    geom_point(position = position_dodge(width = 0.3), size = 3) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(0.3)) +
    geom_line(position = position_dodge(0.3)) +
    labs(
      title = paste("Growth Ratio by Treatment and Wounding for", species_name),
      x = "Treatment",
      y = "Estimated Mean Growth Ratio",
      color = "Wound Status"
    ) +
    theme_bw(base_size = 14) -> p
  
  print(p)
  
  # Print table
  print(plot_data)
}

# Analyze Acropora
analyze_species(df %>% filter(species == "Acropora"), "Acropora")

# Analyze Porites
analyze_species(df %>% filter(species == "Porites"), "Porites")








####
#some quick stats ignoring tank
####

# Fit a linear model with interaction terms
interaction_model <- lm(delta_mass ~ CSA_cm2 * wound_raw * treatment * species, data = df)

# Summary of the model
summary(interaction_model)

# Generate ANOVA table (Type I)
anova_table <- anova(interaction_model)
print(anova_table)

# Alternatively, generate a Type II or Type III ANOVA table
# Type II
anova_table_type_II <- Anova(interaction_model, type = "II")
print(anova_table_type_II)

# Type III (requires contrast specification, especially for categorical variables)
anova_table_type_III <- Anova(interaction_model, type = "III")
print(anova_table_type_III)
# Add predicted values to the data frame
df$predicted_delta_mass <- predict(interaction_model, newdata = df)

# Plot predictions
ggplot(df, aes(x = CSA_cm2, y = predicted_delta_mass, color = treatment)) +
  geom_line(aes(linetype = wound_raw), size = 1) +
  facet_wrap(~ species, scales = "free") +
  labs(
    title = "Predicted Effects of Wounding and Fish on Coral Growth",
    x = "Surface Area (cm²)",
    y = "Predicted Change in Mass (g)",
    color = "Fish Treatment",
    linetype = "Wounding"
  ) +
  theme_minimal()


#pull out selected columns and rename them 

df1<- df %>% 
  select(unique_id, length_days, treatment, wound_raw, tank, species, growth.g, growth.percm2, data_quality) %>% 
  mutate(tank = factor(tank, levels = c("Fish A", "Fish B", "No Fish A", "No Fish B"))) %>% 
  rename(wound_treatment = wound_raw) %>% 
  rename(fish_treatment = treatment)

#df2 <- df1 %>%
 # select(unique_id, species, wound_treatment, fish_treatment) %>%
  #distinct() %>%
  #group_by(species, fish_treatment,  wound_treatment) %>%
  #summarise(count = n(), .groups = 'drop') 


my_theme2 <- theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = "white"), # Ensure plot background is white
    panel.background = element_rect(fill = "white", color = "white") # Ensure panel background is white
  )

library(ggplot2)
library(RColorBrewer)

  
# Ensure df1 has the right structure and correct factor levels
df1$tank <- factor(df1$tank)
df1$fish_treatment <- factor(df1$fish_treatment)
df1$species <- factor(df1$species)
df1$wound_treatment <- factor(df1$wound_treatment)

#data_frame_for_modeling
#write_csv(df1, "data/clean_modeling_df.csv")

dfhighdq <- df1 %>% 
  filter(is.na(data_quality))

#all_together
 ggplot(df1, aes(x = interaction(species, fish_treatment), y = growth.percm2, group = interaction(species, fish_treatment), color = wound_treatment)) +
  geom_boxplot(outlier.shape = NA, fill = "white") +  # Hides outliers in the boxplots
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +  # Adds jittered points to show individual data points
  facet_wrap(~species, scales = "free_x") +  # Each species gets a separate panel
  xlab("treatment") +
  ylab("growth per cm2") +
  ggtitle("n = 80 skeletal growth g per cm2 after 38 days") +
  scale_color_manual(values = c("wounded" = "#9B241C",  "not wounded" = "#006479")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x labels for readability
  my_theme2  # Assuming my_theme2 is defined correctly

# Print the plo

# scale_color_manual(values = c("Fish A" = "#9B241C", "Fish B" = "#9B740B", "No Fish A" = "#004A8E", "No Fish B" = "#006479")) +

ggsave("plots/all_together_all_data.pdf")


colorthem  <- c("#81CDC1", "#003B32")


dfa<- df1 %>% 
  filter(species == "Acropora")

ggplot(dfa, aes(x = interaction(wound_treatment, fish_treatment), y = growth.percm2, group = interaction(wound_treatment, fish_treatment), color = tank)) +
  geom_boxplot(outlier.shape = NA, fill = "white") +  # Hides outliers in the boxplots
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +  # Adds jittered points to show individual data points
  facet_wrap(~wound_treatment, scales = "free_x") +  # Each species gets a separate panel
  xlab("treatment") +
  ylab("growth percm2") +
  ggtitle("Acropora") +
  scale_color_manual(values = c("#9B241C", "#006479", "#81CDC1", "#003B32")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x labels for readability
  my_theme2 

ggsave("plots/tank_acropora_by_wounded.pdf")


ggplot(dfa, aes(x = interaction(wound_treatment, fish_treatment), y = growth.percm2, group = interaction(wound_treatment, fish_treatment), color = wound_treatment)) +
  geom_boxplot(outlier.shape = NA, fill = "white") +  # Hides outliers in the boxplots
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +  # Adds jittered points to show individual data points
  facet_wrap(~fish_treatment, scales = "free_x") +  # Each species gets a separate panel
  xlab("treatment") +
  ylab("growth per cm2") +
  ggtitle("Acropora") +
  scale_color_manual(values = c("wounded" = "#9B241C",  "not wounded" = "#006479")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x labels for readability
  my_theme2 
ggsave("plots/acropora_by_fish.pdf")

scale_color_manual(values = c("wounded" = "#9B241C",  "not wounded" = "#006479")) +


  dfp<- df1 %>% 
  filter(species == "Porites")
  

ggplot(dfp, aes(x = interaction(wound_treatment, fish_treatment), y = growth.percm2, group = interaction(wound_treatment, fish_treatment), color = wound_treatment)) +
  geom_boxplot(outlier.shape = NA, fill = "white") +  # Hides outliers in the boxplots
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +  # Adds jittered points to show individual data points
  facet_wrap(~wound_treatment, scales = "free_x") +  # Each species gets a separate panel
  xlab("treatment") +
  ylab("growth per cm2") +
  ggtitle("Porites n = 40 ") +
  scale_color_manual(values = c("wounded" = "#9B241C",  "not wounded" = "#006479")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x labels for readability
  my_theme2 

ggsave("plots/n=40_porites_by_wounded.pdf")


ggplot(dfp, aes(x = interaction(wound_treatment, fish_treatment), y = growth.percm2, group = interaction(wound_treatment, fish_treatment), color = tank)) +
  geom_boxplot(outlier.shape = NA, fill = "white") +  # Hides outliers in the boxplots
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +  # Adds jittered points to show individual data points
  facet_wrap(~fish_treatment, scales = "free_x") +  # Each species gets a separate panel
  xlab("treatment") +
  ylab("growth std by initial weight") +
  ggtitle("Porites all") +
  scale_color_manual(values = c("#9B241C", "#006479", "#81CDC1", "#003B32")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x labels for readability
  my_theme2 

ggsave("plots/N=40_porites_by_fish.pdf")





###################### remove low data quality


dfahq<- dfhighdq %>% 
  filter(species == "Acropora")

ggplot(dfahq, aes(x = interaction(wound_treatment, fish_treatment), y = growth.percm2, group = interaction(wound_treatment, fish_treatment), color = tank)) +
  geom_boxplot(outlier.shape = NA, fill = "white") +  # Hides outliers in the boxplots
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +  # Adds jittered points to show individual data points
  facet_wrap(~wound_treatment, scales = "free_x") +  # Each species gets a separate panel
  xlab("treatment") +
  ylab("growth percm2") +
  ggtitle("Acropora highDQ n = 37") +
  scale_color_manual(values = c("#9B241C", "#006479", "#81CDC1", "#003B32")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x labels for readability
  my_theme2 

ggsave("plots/highDQ_TANKacropora_by_wounded.pdf")


ggplot(dfahq, aes(x = interaction(wound_treatment, fish_treatment), y = growth.percm2, group = interaction(wound_treatment, fish_treatment), color = wound_treatment)) +
  geom_boxplot(outlier.shape = NA, fill = "white") +  # Hides outliers in the boxplots
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +  # Adds jittered points to show individual data points
  facet_wrap(~fish_treatment, scales = "free_x") +  # Each species gets a separate panel
  xlab("treatment") +
  ylab("growth per cm2") +
  ggtitle("Acropora n = 37 high DQ") +
  scale_color_manual(values = c("wounded" = "#9B241C",  "not wounded" = "#006479")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x labels for readability
  my_theme2 
ggsave("plots/high_DQ_acropora_by_fish.pdf")

scale_color_manual(values = c("wounded" = "#9B241C",  "not wounded" = "#006479")) +
  
  
  dfphighdq<- dfhighdq %>% 
  filter(species == "Porites")


ggplot(dfphighdq, aes(x = interaction(wound_treatment, fish_treatment), y = growth.percm2, group = interaction(wound_treatment, fish_treatment), color = tank)) +
  geom_boxplot(outlier.shape = NA, fill = "white") +  # Hides outliers in the boxplots
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +  # Adds jittered points to show individual data points
  facet_wrap(~wound_treatment, scales = "free_x") +  # Each species gets a separate panel
  xlab("treatment") +
  ylab("growth per cm2") +
  ggtitle("Porites high DQ n = 33") +
  scale_color_manual(values = c("#9B241C", "#006479", "#81CDC1", "#003B32")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x labels for readability
  my_theme2 

ggsave("plots/TANKhighDQ_porites_by_wounded.pdf")


ggplot(dfphighdq, aes(x = interaction(wound_treatment, fish_treatment), y = growth.percm2, group = interaction(wound_treatment, fish_treatment), color = wound_treatment)) +
  geom_boxplot(outlier.shape = NA, fill = "white") +  # Hides outliers in the boxplots
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +  # Adds jittered points to show individual data points
  facet_wrap(~fish_treatment, scales = "free_x") +  # Each species gets a separate panel
  xlab("treatment") +
  ylab("growth std by initial weight") +
  ggtitle("Porites high DQ n =33") +
  scale_color_manual(values = c("wounded" = "#9B241C",  "not wounded" = "#006479")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x labels for readability
  my_theme2 

ggsave("plots/high_DQ_porites_by_fish.pdf")






