# ===============================
# NOTE TO HAYDEN
# ===============================
# This script performs an allometric analysis of coral growth using buoyant weight measurements.
# It filters and reshapes initial and final weights, merges in tank metadata, 
# tests for initial mass differences across treatments, fits an allometric mixed-effects model 
# (log-log), and evaluates treatment effects (fish and wound) on allometrically scaled growth 
# using likelihood ratio tests.

# ✅ Next Steps:
# 1. Incorporate **time** and **surface area** into the growth metric to reflect standardized
#    growth per day per surface area — consistent with Pete Edmunds’ approach.
#    The formula will be:
#       (final - initial) / (surface area × time)
#    which gives a rate in **mg/cm²/day**.
#
# 2. For surface area:
#    - Use **wax-dipped surface area estimates**.
#    - Apply the **dowel calibration curve**.
#    - Subtract **wound area** from total surface area for a corrected estimate.
#
# 3. Once that’s calculated, rerun the same analysis (Sections 5–7) with this new growth response.
#    That includes:
#    - Visualization of raw and scaled growth.
#    - Allometric adjustment (if needed).
#    - Statistical testing (lmer and likelihood ratio tests).


# ===============================
# 1. Load Required Libraries
# ===============================
library(tidyverse)
library(janitor)
library(here)
library(gt)
library(broom.mixed)
library(ggpubr)
library(lme4)

# ===============================
# 2. Load Buoyant Weight Data
# ===============================
bw_raw <- read_csv(here("data", "fish_regen_buoyantweight (1).csv")) %>%
  clean_names() %>%
  filter(coral_id != 45)%>%
  filter(coral_id != 91)

# Split initial and final weights
bw_initial <- bw_raw %>%
  filter(date == "initial") %>%
  select(coral_id, wound, fish, initial_weight = bouyantweight_g)

bw_final <- bw_raw %>%
  filter(date == "final") %>%
  select(coral_id, final_weight = bouyantweight_g)

# Join initial and final to get one row per coral
bw_wide <- left_join(bw_initial, bw_final, by = "coral_id") %>%
  mutate(
    delta_mass = final_weight - initial_weight,
    coral_id = as.character(coral_id)
  ) 



  #.....................Load Surface Area Data.....................ADRIAN LOOK HERE!
sa <- read_csv(here("data", "fish_regen_mastersheet_wound closure_necrosis_sa - mastersheet.csv")) %>% 
  clean_names() %>% 
  mutate(coral_id = as.character(coral_id)) %>% 
  select(-fish, -date_collected, -date_fragment) %>% 
  mutate(sa_col=ifelse(wound == "Large", sa_cal - 3, #removing 3 cm^2 from large wounds and 1 from small
                       ifelse(wound == "Small", sa_cal - 1, sa_cal))) %>% 
  select(-wound) #removing the column so it doesn't duplicate when we join with bw_wide

bw_sa_wide <- left_join(bw_wide, sa, by = "coral_id") %>% 
  mutate(
    wound = factor(wound, levels = c("No Wound", "Small", "Large")),
    tank = as.factor(tank)
  ) %>% 
  mutate(tank = factor(tank, levels = sort(unique(as.numeric(as.character(tank))))))
  
# ===============================
# 3. Load Tank Metadata
# ===============================
tank_df <- read_csv(here("data", "Copy of fish_regen_mastersheet_wound closure_necrosis - mastersheet.csv")) %>%
  clean_names() %>%
  select(coral_id, tank) %>%
  mutate(coral_id = as.character(coral_id))

# ===============================
# 4. Merge Tank Info into Wide Format
# ===============================
bw_merged <- left_join(bw_wide, tank_df, by = "coral_id") %>%
  mutate(
    wound = factor(wound, levels = c("No Wound", "Small", "Large")),
    tank = as.factor(tank)
  )


# ===============================
# 3. Visualize Raw Growth Data across treatment and tank
# ===============================

#do tanks or treatments have distinguishable sizes initialy? 

bw_plot <- bw_merged %>%
  mutate(tank = factor(tank, levels = sort(unique(as.numeric(as.character(tank))))))

ggplot(bw_plot, aes(x = tank, y = initial_weight, fill = wound)) +
  # Violin plot with dodge
  geom_violin(
    aes(group = interaction(tank, wound)),
    trim = FALSE, color = NA,
    alpha = 0.4,
    position = position_dodge(width = 0.9)
  ) +
  # Boxplot with matching dodge
  geom_boxplot(
    aes(group = interaction(tank, wound)),
    width = 0.15,
    position = position_dodge(width = 0.9),
    outlier.shape = NA,
    color = "black", alpha = 0.7
  ) +
  # Jittered raw points
  geom_jitter(
    aes(color = wound),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
    size = 1.5, alpha = 0.6
  ) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "Initial Coral Mass by Tank and Wound Treatment",
    subtitle = "Violin + boxplot + individual data points",
    x = "Tank",
    y = "Initial Buoyant Weight (g)",
    fill = "Wound",
    color = "Wound"
  ) +
  theme_pubr(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0),
    legend.position = "top"
  )

ggsave(
  here("figures", "growth","initial_mass_by_treatment.png"),
  width = 8,
  height = 6,
  dpi = 300
)
#ask if there are difference statistically THIS IS SMART

# Fit linear model with interaction
model_lm <- lm(initial_weight ~ wound * tank, data = bw_plot)

# ANOVA table
anova(model_lm)

#no so we move forward with our analysis looking at allooometry next 

#...................test for differences in SA...................ADRIAN LOOK HERE!

model_sa <- lm(sa_cal ~ wound * tank, data = bw_sa_wide)

anova(model_sa)

bw_plot <- bw_merged %>%
  mutate(tank = factor(tank, levels = sort(unique(as.numeric(as.character(tank))))))
#No significant differences

# ===============================
# 3. Transform Data for Allometric Modeling
# ===============================
# Create a dataset with log-transformed weights and growth ratio
# Drop NAs to ensure valid model fitting
log_df <- bw_merged %>%
  mutate(
    log_i_weight = log(initial_weight),
    log_f_weight = log(final_weight),
    growth_ratio = final_weight / initial_weight
  ) %>%
  drop_na(log_i_weight, log_f_weight, tank)

# ===============================
# 4. Fit Mixed-Effects Allometric Model
# ===============================
# Fit a linear mixed-effects model:
# Predict final log-weight from initial log-weight
# Include tank as a random intercept to account for tank-level variation

model_lmer <- lmer(log_f_weight ~ log_i_weight + (1 | tank), data = log_df)

# Summary with fixed effect (slope), intercept, and random effects
summary(model_lmer)

# Extract allometric slope (b-value) with 95% confidence interval
b_value <- tidy(model_lmer, effects = "fixed", conf.int = TRUE) %>%
  filter(term == "log_i_weight")
print(b_value)


# Generate predicted log final weights across observed range of initial weights
fit_df <- log_df %>%
  select(log_i_weight) %>%
  distinct() %>%
  arrange(log_i_weight) %>%
  mutate(predicted = predict(model_lmer, newdata = ., re.form = NA))


#  Plot Allometric Fit and Observed Data
ggplot(log_df, aes(x = log_i_weight, y = log_f_weight)) +
  geom_point(size = 2.2, shape = 21, fill = "steelblue", color = "black", alpha = 0.8, stroke = 0.3) +
  geom_line(data = fit_df, aes(x = log_i_weight, y = predicted), color = "darkred", linewidth = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray60", linewidth = 0.8) +
  labs(
    title = "Allometric Scaling of Coral Buoyant Weight",
    subtitle = "Mixed-effects model with tank as a random intercept",
    x = "Log Initial Weight (g)",
    y = "Log Final Weight (g)"
  ) +
  theme_pubr(base_size = 14)


library(gt)

# Extract and prepare the b-value data
b_value_gt <- tidy(model_lmer, effects = "fixed", conf.int = TRUE) %>%
  filter(term == "log_i_weight") %>%
  transmute(
    Term = "Allometric slope (b-value)",
    Estimate = round(estimate, 3),
    `95% CI` = paste0("[", round(conf.low, 3), ", ", round(conf.high, 3), "]")
  )

# Create and format the gt table for the b value 
b_value_gt %>%
  gt() %>%
  tab_header(
    title = "Estimated Allometric Slope",
    subtitle = "Mixed-effects model with tank as a random intercept"
  ) %>%
  cols_label(
    Term = "",
    Estimate = "Estimate",
    `95% CI` = "95% Confidence Interval"
  ) %>%
  tab_options(
    table.font.size = 14,
    heading.title.font.size = 16,
    heading.subtitle.font.size = 13
  )


#so the change in mas  is pretty close to isometric but a little lower

# ===============================
# 5. Visualize Growth (Final - Initial) by Initial Mass and Treatment
# ===============================

# Prep data: drop NAs, ensure proper factor levels
growth_df <- bw_merged %>%
  drop_na(initial_weight, final_weight, wound) %>%
  mutate(
    delta_mass = final_weight - initial_weight,
    wound = factor(wound, levels = c("No Wound", "Small", "Large"))
  )

# Plot: Growth vs. Initial Mass, colored by Wound Treatment
ggplot(growth_df, aes(x = initial_weight, y = delta_mass, color = wound)) +
  geom_point(size = 2.2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.1) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "Change in Coral Mass as a Function of Initial Mass",
    subtitle = "Colored by Wound Treatment",
    x = "Initial Buoyant Weight (g)",
    y = "Change in Mass (g)",
    color = "Wound Treatment"
  ) +
  theme_pubr(base_size = 14)

ggsave(
  here("figures", "growth","growth_vs_initial_mass.png"),
  width = 8,
  height = 6,
  dpi = 300
)

#..................Visualizing SA by Inital Mass and Treatment...................ADRIAN LOOK HERE!

# Prep data: drop NAs, ensure proper factor levels
bw_sa_df <- bw_sa_wide %>%
  drop_na(initial_weight, final_weight, wound) %>%
  filter(coral_id != 103) %>% 
  mutate(
    delta_mass = final_weight - initial_weight,
    wound = factor(wound, levels = c("No Wound", "Small", "Large"))
  )

# Plot: SA vs. Initial Mass, colored by Wound Treatment
ggplot(bw_sa_df, aes(x = initial_weight, y = sa_cal, color = wound)) +
  geom_point(size = 2.2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.1) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "Coral Surface Area as a Function of Initial Mass",
    subtitle = "Colored by Wound Treatment",
    x = "Initial Buoyant Weight (g)",
    y = "Surface Area (cm^2)",
    color = "Wound Treatment"
  ) +
  theme_pubr(base_size = 14)

ggsave(
  here("figures", "growth","sa_vs_initial_mass.png"),
  width = 8,
  height = 6,
  dpi = 300
)

# ===============================
# 6. Statistical Analysis of Allometrically Scaled Growth
# ===============================

# --- Step 1: Prepare Data with Allometric Scaling ---
b_est <- b_value$estimate[1]  # b-value from previous section

growth_stats <- bw_merged %>%
  drop_na(initial_weight, final_weight, fish, wound, tank) %>%
  mutate(
    wound = factor(wound, levels = c("No Wound", "Small", "Large")),
    fish = factor(fish),
    tank = factor(tank),
    growth_scaled = log(final_weight/initial_weight ^ b_est) 
  )

# --- Step 2: Fit Mixed-Effects Models ---

# Full model with interaction
mod_int <- lmer(growth_scaled ~ fish * wound + (1 | tank), data = growth_stats)

# Additive model without interaction
mod_add <- lmer(growth_scaled ~ fish + wound + (1 | tank), data = growth_stats)

# Fish-only model
mod_fish <- lmer(growth_scaled ~ fish + (1 | tank), data = growth_stats)

# Wound-only model
mod_wound <- lmer(growth_scaled ~ wound + (1 | tank), data = growth_stats)

# Null model (intercept + tank)
mod_null <- lmer(growth_scaled ~ 1 + (1 | tank), data = growth_stats)

# --- Step 3: Likelihood Ratio Tests ---
lrt_interaction <- anova(mod_add, mod_int, test = "Chisq")
lrt_fish <- anova(mod_wound, mod_add, test = "Chisq")
lrt_wound <- anova(mod_fish, mod_add, test = "Chisq")

# --- Step 4: Summarize Results in a Table ---


lrt_table <- tibble(
  Test = c("Interaction (Fish × Wound)", "Fish Effect", "Wound Effect"),
  ChiSq = c(
    lrt_interaction$Chisq[2],
    lrt_fish$Chisq[2],
    lrt_wound$Chisq[2]
  ),
  Df = c(
    lrt_interaction$Df[2],
    lrt_fish$Df[2],
    lrt_wound$Df[2]
  ),
  p_value = c(
    lrt_interaction$`Pr(>Chisq)`[2],
    lrt_fish$`Pr(>Chisq)`[2],
    lrt_wound$`Pr(>Chisq)`[2]
  )
) %>%
  mutate(across(where(is.numeric), round, 3))

# Format with gt
lrt_table %>%
  gt() %>%
  tab_header(
    title = "Table 2. Likelihood Ratio Tests for Effects on Coral Growth",
    subtitle = paste("Scaled growth:", expression((final - initial)/initial^b))
  ) %>%
  cols_label(
    Test = "Model Comparison",
    ChiSq = "Chi-squared",
    Df = "Degrees of Freedom",
    p_value = "P-value"
  ) %>%
  fmt_number(columns = vars(ChiSq, p_value), decimals = 3) %>%
  tab_options(
    table.font.size = 14,
    heading.title.font.size = 16,
    heading.subtitle.font.size = 13
  )

#..........Statistical Analysis of Growth/SA/Time Metric.........ADRIAN LOOK HERE!


bw_sa_stats <- bw_sa_df %>%
  drop_na(initial_weight, final_weight, fish, wound, tank) %>%
  mutate(
    wound = factor(wound, levels = c("No Wound", "Small", "Large")),
    fish = factor(fish),
    tank = factor(tank),
    bw_sa_time = (final_weight - initial_weight) / (sa_cal * 21) 
  )

# --- Step 2: Fit Mixed-Effects Models ---

# Full model with interaction
mod_int_sa <- lmer(bw_sa_time ~ fish * wound + (1 | tank), data = bw_sa_stats)

# Additive model without interaction
mod_add_sa <- lmer(bw_sa_time ~ fish + wound + (1 | tank), data = bw_sa_stats)

# Fish-only model
mod_fish_sa <- lmer(bw_sa_time ~ fish + (1 | tank), data = bw_sa_stats)

# Wound-only model
mod_wound_sa <- lmer(bw_sa_time ~ wound + (1 | tank), data = bw_sa_stats)

# Null model (intercept + tank)
mod_null_sa <- lmer(bw_sa_time ~ 1 + (1 | tank), data = bw_sa_stats)

# --- Step 3: Likelihood Ratio Tests ---
lrt_interaction_sa <- anova(mod_add_sa, mod_int_sa, test = "Chisq")
lrt_fish_sa <- anova(mod_wound_sa, mod_add_sa, test = "Chisq")
lrt_wound_sa <- anova(mod_fish_sa, mod_add_sa, test = "Chisq")

# --- Step 4: Summarize Results in a Table ---


lrt_table <- tibble(
  Test = c("Interaction (Fish × Wound)", "Fish Effect", "Wound Effect"),
  ChiSq = c(
    lrt_interaction_sa$Chisq[2],
    lrt_fish_sa$Chisq[2],
    lrt_wound_sa$Chisq[2]
  ),
  Df = c(
    lrt_interaction_sa$Df[2],
    lrt_fish_sa$Df[2],
    lrt_wound_sa$Df[2]
  ),
  p_value = c(
    lrt_interaction_sa$`Pr(>Chisq)`[2],
    lrt_fish_sa$`Pr(>Chisq)`[2],
    lrt_wound_sa$`Pr(>Chisq)`[2]
  )
) %>%
  mutate(across(where(is.numeric), round, 3))

# Format with gt
lrt_table %>%
  gt() %>%
  tab_header(
    title = "Table 2. Likelihood Ratio Tests for Effects on Coral Growth as a function of Surface Area and Time",
    subtitle = paste("Scaled growth:", expression((final - initial)/(surface_area * 21)))
  ) %>%
  cols_label(
    Test = "Model Comparison",
    ChiSq = "Chi-squared",
    Df = "Degrees of Freedom",
    p_value = "P-value"
  ) %>%
  fmt_number(columns = vars(ChiSq, p_value), decimals = 3) %>%
  tab_options(
    table.font.size = 14,
    heading.title.font.size = 16,
    heading.subtitle.font.size = 13
  )

#so with no interaction the final model is going to be just the main effects 


# ===============================
# 7. Visualize Scaled Growth Across Treatments (Faceted by Wound)
# ===============================

# Ensure fish levels are ordered correctly
growth_stats <- growth_stats %>%
  mutate(fish = factor(fish, levels = c("No Fish", "Fish")))

# Plot: Scaled growth by Fish, faceted by Wound
ggplot(growth_stats, aes(x = fish, y = growth_scaled, fill = fish)) +
  geom_boxplot(
    alpha = 0.5,
    width = 0.6,
    outlier.shape = NA,
    color = "black"
  ) +
  geom_jitter(
    aes(color = fish),
    width = 0.15,
    size = 2,
    alpha = 0.6
  ) +
  facet_wrap(~wound) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  labs(
    # title = "Allometrically Scaled Coral Growth by Treatment",
    # subtitle = expression(paste("Growth normalized by ", initial^b, ", grouped by Fish Treatment")),
    x = "Fish Treatment",
    y = paste("Calcification", expression(log(Mf/Mi^b))),
    fill = "Fish",
    color = "Fish"
  ) +
  theme_pubr(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "top"
  )

ggsave(
  here("figures","growth", "scaled_growth_by_treatment.png"),
  width = 8,
  height = 6,
  dpi = 300
)

# Visualizing Growth, SA, Time Metric Across Treatments (facet by wound) ~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ADRIAN LOOK HERE!


# Ensure fish levels are ordered correctly
bw_sa_stats_viz <- bw_sa_stats %>%
  mutate(fish = factor(fish, levels = c("No Fish", "Fish")))

# Plot: Scaled growth by Fish, faceted by Wound
ggplot(bw_sa_stats_viz, aes(x = fish, y = bw_sa_time, fill = fish)) +
  geom_boxplot(
    alpha = 0.5,
    width = 0.6,
    outlier.shape = NA,
    color = "black"
  ) +
  geom_jitter(
    aes(color = fish),
    width = 0.15,
    size = 2,
    alpha = 0.6
  ) +
  facet_wrap(~wound) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  labs(
    # title = "Coral Growth as a Function of SA and Time by Treatment",
    # subtitle = expression(paste("Growth normalized by ", initial^b, ", grouped by Fish Treatment")),
    x = "Fish Treatment",
    y = expression((Final - Initial) / (Surface_Area * Time)),
    fill = "Fish",
    color = "Fish"
  ) +
  theme_pubr(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "top"
  )

ggsave(
  here("figures", "growth","growth_sa_time_by_treatment.png"),
  width = 8,
  height = 6,
  dpi = 300
)
