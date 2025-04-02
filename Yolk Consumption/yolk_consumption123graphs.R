library(tidyverse)
library(lubridate)
library(ggplot2)
library(plotly)

# Read the CSV
master_growth_df = 'av_slopes_all.csv'
master_slopes = read.csv(master_growth_df)

# Convert relevant columns to factors
master_slopes$brood <- factor(master_slopes$brood)
master_slopes$parent.treatment <- factor(master_slopes$parent.treatment, 
                                         levels = c("control", "low", "medium", "high"))
master_slopes$treatment <- factor(master_slopes$treatment, 
                                  levels = c("control", "treatment"))

# Filter for Brood 1, 2, and 3
brood_slopes <- master_slopes %>% filter(brood %in% c("1", "2", "3"))
view(brood_slopes)
# Create the combined box plot for Brood 1, 2, and 3

ggplot(brood_slopes, 
       aes(x = parent.treatment, 
           y = abs(avg_yolk_slope),  
           fill = treatment, 
           color = brood)) +  
  geom_boxplot(position = position_dodge2(width = 5, preserve = "single"), 
               size = 1,  # Increase outline thickness
               width = 0.9) +  
  labs(title = "Average Growth Rates (Absolute Value)",
       x = "Parental Treatment",
       y = " Average Yolk Consumption (%/day)") +
  scale_fill_manual(name = "Embryo Treatment",
                    values = c("control" = "skyblue3", "treatment" = "tomato"),
                    labels = c("control" = "Control", "treatment" = "Heat Shocked")) +  
  scale_color_manual(name = "Brood",  
                     values = c("1" = "black", "2" = "purple", "3" = "green"),
                     labels = c("1" = "1", "2" = "2", "3" = "3")) +  
  theme_minimal() +
  theme(legend.title = element_text(size = 14),  
        legend.text = element_text(size = 13),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 15),  
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 14))  



library(lme4)
library(lmerTest)  # For p-values




#-------------------------------------------------------------------------------
# by exposures
view(brood_slopes)

exposure_plot <- ggplot(brood_slopes, aes(x = factor(num_hs_exposures), y = abs(avg_yolk_slope), fill = treatment)) +  
  geom_boxplot() +  # Create the box plot
  labs(title = "Average % Heat Shock Exposures",
       x = "Number of Heat Shock Exposures",
       y = "% yolk consumed/ day") +
  scale_fill_manual(values = c("control" = "skyblue3", "treatment" = "tomato")) +  # Colors for control and treatment
  scale_x_discrete(limits = as.character(0:18)) +  # Correct X-axis limits using discrete values
  theme_minimal() +
  theme(legend.title = element_blank(),  # Remove legend title
        plot.title = element_text(hjust = 0.5),  # Center the title
        axis.text.x = element_text(size = 11),  # Increase x-axis text size
        axis.title.x = element_text(size = 12),  # Increase x-axis title size
        axis.title.y = element_text(size = 12))  # Increase y-axis title size



#-------------------------------------------------------------------------------


# Load the individual data
individual_df_read = "yolk_consumption_regression_all.csv"
master_individual = read.csv(individual_df_read)
view(master_individual)
# View the structure of the data to ensure it's loaded correctly

# Assuming your data frame 'master_individual' contains 'slope', 'cage', 'brood', and 'treatment' columns

library(dplyr)

# Modify the yolk_slope column to take the absolute value
master_individual <- master_individual %>%
  mutate(yolk_slope = abs(yolk_slope))

# Now create the plot with the absolute values of yolk_slope
ggplot(master_individual, aes(x = factor(cage), y = yolk_slope, fill = treatment, color = factor(brood))) +
  geom_boxplot(outlier.shape = NA, size = 0.5) +  # Decrease outline thickness (adjusted to 0.5)
  labs(title = "Yolk consumption by Cage, Treatment, and Brood",
       x = "Cage", 
       y = "Yolk %/ day", 
       fill = "Treatment", 
       color = "Brood") +
  scale_fill_manual(values = c("control" = "skyblue3", "treatment" = "tomato")) +  # Fill colors for treatment
  scale_color_manual(values = c("1" = "black", "2" = "purple", "3" = "green")) +  # Outline colors for brood
  scale_x_discrete(limits = as.character(1:48)) +  # Set x-axis limits to range from 1 to 48
  theme_minimal() +
  theme(legend.title = element_blank(),  # Remove legend title
        plot.title = element_text(hjust = 0.5),  # Center the plot title
        axis.text.x = element_text(size = 11),  # Increase x-axis text size
        axis.title.x = element_text(size = 12),  # Increase x-axis title size
        axis.title.y = element_text(size = 12))  # Increase y-axis title size






ggplot(master_individual %>% filter(cage %in% c("4", "8", "7", "19", "10")), 
       aes(x = factor(cage), y = yolk_slope, fill = treatment, color = factor(brood))) +
  geom_boxplot(outlier.shape = NA, size = 0.8) +  # Increase outline thickness
  labs(title = "Yolk Consumption by Cage, Treatment, and Brood (Brood 1)",
       x = "Cage", 
       y = "Yolk consumption (% / day)", 
       fill = "Embryo Treatment", 
       color = "Brood") +  # Keep the brood legend
  scale_fill_manual(values = c("control" = "skyblue3", "treatment" = "tomato"),
                    labels = c("control" = "Control", "treatment" = "Heat Shocked")) +  # Labels updated
  scale_color_manual(values = c("1" = "black", "2" = "purple", "3" = "green"),
                     labels = c("1" = "1", "2" = "2", "3" = "3")) +  # Brood labels
  scale_x_discrete(limits = c("4", "8", "7", "19", "10")) +  # Filter specific cages
  theme_minimal() +
  theme(legend.title = element_text(size = 14),  # Adjust legend title size
        legend.text = element_text(size = 13),  # Adjust legend text size
        plot.title = element_text(hjust = 0.5),  # Center the plot title
        axis.text.x = element_text(size = 15),  # Increase x-axis text size
        axis.title.x = element_text(size = 14),  # Increase x-axis title size
        axis.title.y = element_text(size = 15),  # Increase y-axis title size
        axis.text.y = element_text(size = 14))  # Increase y-axis text size


ggplot(master_individual %>% filter(brood == 1), 
       aes(x = factor(cage), y = yolk_slope, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, size = 0.8) +  # Increase outline thickness
  labs(title = "Yolk Consumption by Cage, Treatment, and Brood 1",
       x = "Cage", 
       y = "Yolk consumed (% / day)", 
       fill = "Embryo Treatment") +  # Legend title updated
  scale_fill_manual(values = c("control" = "skyblue3", "treatment" = "tomato"),
                    labels = c("control" = "Control", "treatment" = "Heat Shocked")) +  # Labels updated
  scale_x_discrete() +  # Use the cage values without needing to specify limits
  theme_minimal() +
  theme(legend.title = element_text(size = 14),  # Adjust legend title size
        legend.text = element_text(size = 13),  # Adjust legend text size
        plot.title = element_text(hjust = 0.5),  # Center the plot title
        axis.text.x = element_text(size = 15, angle = 90, hjust = 1),  # Rotate x-axis text
        axis.title.x = element_text(size = 14),  # Increase x-axis title size
        axis.title.y = element_text(size = 15),  # Increase y-axis title size
        axis.text.y = element_text(size = 14),  # Increase y-axis text size
        legend.key = element_blank(),  # Remove legend key
        strip.text = element_text(size = 14)) +  # Increase size of facet labels (parental groups)
  facet_wrap(~ factor(parent.treatment, levels = c("control", "low", "medium", "high")), 
             scales = "free_x", ncol = 4)  # Set order of parent.treatment and display in 4 columns








#OLD PRELIMINARY STATS
# ANOVA by parent treatment-------------------------------------------------------
# Filter the dataset for low parent treatment
low_parent_data <- master_individual %>%
  filter(parent.treatment == "low")

# Run ANOVA for low parent treatment
anova_low_parent <- aov(yolk_slope ~ brood * treatment, data = low_parent_data)

# Get the summary of the ANOVA for low parent treatment
summary(anova_low_parent)

# Filter the dataset for medium parent treatment
medium_parent_data <- master_individual %>%
  filter(parent.treatment == "medium")

# Run ANOVA for medium parent treatment
anova_medium_parent <- aov(yolk_slope ~ brood * treatment, data = medium_parent_data)

# Get the summary of the ANOVA for medium parent treatment
summary(anova_medium_parent)
#------------------------------------------------------------------------------




