
#......................Installing packages.......................
library(tidyverse)
library(janitor)
library(here)
library(gt)
library(broom.mixed)
#.............Loading in calibrations and coral data.............
calibration_df <- read_csv(here("data", "wax_dip_calibrations.csv"))
coral_df <- read_csv(here("data", "wax_dip_corals.csv"))

#...Calculating the curve coefficients for slope and intercept...
stnd.curve <- lm(calculated_sa~wax_weight_diff, data=calibration_df)
plot(calculated_sa~wax_weight_diff, data=calibration_df)
stnd.curve$coefficients
summary(stnd.curve)$r.squared

#........Calculate surface area using the standard curve.........
coral_df$SA_cal <- stnd.curve$coefficients[2] * coral_df$wax_diff + stnd.curve$coefficients[1]

#..........................Save as a csv.........................

write.csv(coral_df, "coral_sa.csv", row.names = FALSE)
