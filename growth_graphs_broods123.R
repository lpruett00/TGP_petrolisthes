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

# Create the combined box plot for Brood 1, 2, and 3
ggplot(brood_slopes, 
       aes(x = parent.treatment, 
           y = avg_growth_slope, 
           fill = treatment, 
           color = brood)) +  
  geom_boxplot(position = position_dodge2(width = 5, preserve = "single"), 
               size = 1, 
               width = .9) +  # Increase box width
  labs(title = "Average growth rates",
       x = "Parental Treatment",
       y = "Average growth rate (mm^2/week)") +
  scale_fill_manual(name = "Embryo Treatment",  # Rename fill legend
                    values = c("control" = "skyblue3", "treatment" = "tomato"),
                    labels = c("control" = "Control", "treatment" = "Heat Shocked")) +  
  scale_color_manual(name = "Brood",  
                     values = c("1" = "black", "2" = "purple", "3" = "green"),
                     labels = c("1" = "Brood 1", "2" = "Brood 2", "3" = "Brood 3")) +  
  theme_minimal() +
  theme(legend.title = element_text(size = 13), 
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 13))



#-------------------------------------------------------------------------------
# by exposures
view(brood_slopes)

# Ensure num_hs_exposures is a factor with levels from 0 to 18
brood_slopes$num_hs_exposures <- factor(brood_slopes$num_hs_exposures, levels = 0:18)

exposure_plot <- ggplot(brood_slopes, aes(x = num_hs_exposures, y = avg_growth_slope, fill = treatment)) +
  geom_boxplot() +
  labs(title = "Average Growth Rates by Number of Heat Wave Exposures",
       x = "Number of Heat Shock Exposures",
       y = "Average Growth Rate (mm²/week)") +
  scale_fill_manual(values = c("control" = "skyblue3", "treatment" = "tomato")) +
  scale_x_discrete(drop = FALSE) +  # Ensures all levels appear even if missing data
  theme_minimal() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

# Display the plot
exposure_plot





# Load the individual data
individual_df_read = "growth_rate_regression.csv"
master_individual = read.csv(individual_df_read)
view(master_individual)



# Plot by cages

ggplot(master_individual, aes(x = factor(cage), y = slope, fill = treatment, color = factor(brood))) +
  geom_boxplot(outlier.shape = NA, size = 0.5) +  # Decrease outline thickness (adjusted to 0.5)
  labs(title = "Growth Rates by Cage, Treatment, and Brood",
       x = "Cage", 
       y = "Growth Rate (mm^2/week)", 
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
       aes(x = factor(cage), y = slope, fill = treatment, color = factor(brood))) +
  geom_boxplot(outlier.shape = NA, size = 1) +  # Increased outline thickness
  labs(title = "Growth Rates by Cage, Treatment, and Brood",
       x = "Cage", 
       y = "Growth Rate (mm²/week)", 
       fill = "Embryo Treatment", 
       color = "Brood") +
  scale_fill_manual(name = "Embryo Treatment",
                    values = c("control" = "skyblue3", "treatment" = "tomato"),
                    labels = c("control" = "Control", "treatment" = "Heat Shocked")) +  
  scale_color_manual(name = "Brood", 
                     values = c("1" = "black", "2" = "purple", "3" = "green")) +  
  scale_x_discrete(limits = c("4", "8", "7", "19", "10")) +  
  theme_minimal() +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 15),  # No angle applied
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 14)) +  
  geom_text(aes(x = factor(cage), y = -Inf, label = parent.treatment), 
            data = master_individual %>% filter(cage %in% c("4", "8", "7", "19", "10")),
            inherit.aes = FALSE, 
            size = 4, angle = 0, hjust = 0.5, vjust = 1.5, color = "black")  

# Perform a three-way ANOVA with embryo nested in cage
anova_model <- aov(slope ~ brood * treatment * parent.treatment + Error(cage/embryo), 
                   data = master_individual)

# Display the ANOVA table
summary(anova_model)




# Calculate the difference between control and treatment slopes (growth rates)
brood_diff <- master_individual %>%
  group_by(brood) %>%
  summarise(control = mean(slope[treatment == "control"], na.rm = TRUE),
            treatment = mean(slope[treatment == "treatment"], na.rm = TRUE)) %>%
  mutate(diff = abs(control - treatment))  # Calculate the absolute difference

# Check the structure of the data to make sure it's correct
str(brood_diff)

# If the diff column is now populated, proceed to convert 'diff' to numeric and 'brood' to factor
brood_diff <- brood_diff %>%
  mutate(diff = as.numeric(diff), 
         brood = as.factor(brood))

# Check the data again
print(brood_diff)

























