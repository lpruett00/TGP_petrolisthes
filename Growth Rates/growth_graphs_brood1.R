# Load necessary libraries
library(tidyverse)
library(lubridate)  # For handling date conversion
library(ggplot2)
library(plotly)

# Read the CSV
master_growth_df = 'master_file_growth_yolk.csv'
master_growth = read.csv(master_growth_df)

# Convert date column to a uniform Date format (MM-DD-YY for display)
master_growth <- master_growth %>%
  mutate(date = mdy(date))  # Converts mixed formats (M/D/YY, MM-DD-YY) to Date
master_growth <- master_growth %>%
  mutate(area_mm = as.numeric(area_mm))  # Keep date as Date object
# Sort the data by embryo and date
master_growth <- master_growth %>%
  arrange(embryo, date)

view(master_growth)

#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------

# Fit a regression model for each embryo
# Filter data to ensure there are no missing values and each embryo has more than one data point



growth_rate_regression <- master_growth %>%
  filter(brood == 1, !is.na(area_mm), !is.na(date)) %>%  # Ensure only brood 1
  group_by(embryo) %>%
  filter(n() > 1) %>%  # Ensure at least two data points per embryo
  do({
    model <- lm(area_mm ~ date, data = .)  # Fit linear regression (size vs. time)
    data.frame(slope = coef(model)[2],  # Extract the slope (growth rate)
               intercept = coef(model)[1],  
               r_squared = summary(model)$r.squared)  
  }) %>%
  ungroup() %>%
  left_join(master_growth %>% 
              filter(brood == 1) %>%
              select(embryo, brood, num_hs_exposures, cage, treatment, parent.treatment), 
            by = "embryo") %>%  # Join only necessary columns
  distinct()  # Remove duplicate rows

view(growth_rate_regression)
library(writexl)

# Save the dataframe as an Excel file
write_xlsx(growth_rate_regression, "growth_rate_regression_embryo.xlsx")
# Average the slopes by cage, embryo treatment, and include parent treatment
average_slopes_by_cage_treatment_parent <- growth_rate_regression %>%
  group_by(cage, treatment, parent.treatment) %>%
  summarize(
    avg_slope = mean(slope, na.rm = TRUE),  # Calculate avg slope per group
    num_hs_exposures = unique(num_hs_exposures),  # Add num_hs_exposures (assuming it’s the same for each embryo)
    .groups = "drop"
  )  # Ensure that the `num_hs_exposures` is added and that the data is grouped correctly



average_slopes_by_cage_treatment_parent <- growth_rate_regression %>%
  # Join with master_growth to get treatment info and parent_treatment
  group_by(cage, treatment, parent.treatment) %>%  # Group by cage, treatment, and parent_treatment
  summarize(
    avg_slope = mean(slope, na.rm = TRUE),  # Calculate the average slope for each group
    sd_slope = sd(slope, na.rm = TRUE),  # Calculate the standard deviation for each group
    num_hs_exposures = unique(num_hs_exposures),  # Add num_hs_exposures (assuming it's the same for each embryo)
    .groups = "drop"
  ) %>%
  ungroup()




# Ensure parental_treatment is a factor with the correct order
average_slopes_by_cage_treatment_parent$parent.treatment <- factor(average_slopes_by_cage_treatment_parent$parent.treatment, 
                                                                   levels = c("control", "low", "medium", "high"))

# Ensure treatment is a factor for proper ordering in the plot
average_slopes_by_cage_treatment_parent$treatment <- factor(average_slopes_by_cage_treatment_parent$treatment, 
                                                            levels = c("control", "treatment"))

# Create the box plot
ggplot(average_slopes_by_cage_treatment_parent, aes(x = parent.treatment, y = avg_slope, fill = treatment)) +
  geom_boxplot() +  # Create the box plot
  labs(title = "Comparison of Average Slopes by Parental Treatment",
       x = "Parental Treatment",
       y = "Average Slope (Growth Rate)") +
  scale_fill_manual(values = c("control" = "skyblue3", "treatment" = "tomato")) +  # Colors for control and treatment
  theme_minimal() +
  theme(legend.title = element_blank(),  # Remove legend title
        plot.title = element_text(hjust = 0.5),  # Center the title
        axis.text.x = element_text(size = 11),  # Increase x-axis text size
        axis.title.x = element_text(size = 12),  # Increase x-axis title size
        axis.title.y = element_text(size = 12))  # Increase y-axis title size

# Fit the mixed-effects model with parental.treatment and treatment as fixed effects, cage as a random effect
model <- lmer(avg_slope ~ treatment * parent.treatment + (1 | cage), data = average_slopes_by_cage_treatment_parent)

# Check the model summary
summary(model)
# Install and load the 'car' package if you don't have it
# install.packages("car")
library(car)

# Perform Type III ANOVA
Anova(model, type = 3)



## Stats

# one way anova: compare growth rates of CONTROL embryos across parental treatment groups
# Subset only control embryos
control_embryos <- subset(average_slopes_by_cage_treatment_parent, treatment == "control")
view(control_embryos)
shapiro.test(control_embryos$avg_slope) #data is noramlly distributed
library(car)
leveneTest(avg_slope ~ parent.treatment, data = control_embryos)
anova_control <- aov(avg_slope ~ parent.treatment, data = control_embryos)
summary(anova_control)


# one way anova: compare growth rates of TREATMENT embryos across parental treatment group

treatment_embryos <- subset(average_slopes_by_cage_treatment_parent, treatment == "treatment")
view(treatment_embryos)
# Check normality assumption
shapiro.test(treatment_embryos$avg_slope)
#Was not normally distributed, so running this test
kruskal.test(avg_slope ~ parent.treatment, data = treatment_embryos)
# no significance 


# two way anova testing interaction effects (parent treatment, embryo treatment, AND interactive effect of parent and embryo treatment)
average_slopes_by_cage_treatment_parent$treatment <- factor(average_slopes_by_cage_treatment_parent$treatment)
average_slopes_by_cage_treatment_parent$parent.treatment <- factor(average_slopes_by_cage_treatment_parent$parent.treatment)

# Two-way ANOVA (parent.treatment * treatment)
shapiro.test(residuals(two_way_anova))

two_way_anova <- aov(avg_slope ~ parent.treatment * treatment, data = average_slopes_by_cage_treatment_parent)

summary(two_way_anova)

tukey_results <- TukeyHSD(two_way_anova)

print(tukey_results)


# Three way anova
kruskal.test(avg_slope ~ parent.treatment, data = average_slopes_by_cage_treatment_parent)
kruskal.test(avg_slope ~ treatment, data = average_slopes_by_cage_treatment_parent)
kruskal.test(avg_slope ~ cage, data = average_slopes_by_cage_treatment_parent)


# Fit a mixed-effects model
library(lme4)
mixed_model <- lmer(avg_slope ~ parent.treatment * treatment + (1 | cage), 
                    data = average_slopes_by_cage_treatment_parent)

# Check summary of the model
summary(mixed_model)







# is there a difference between control and treatment growth rates within parental groups?

control_control <- subset(average_slopes_by_cage_treatment_parent, parent.treatment == "control" & treatment == "control")
treatment_control <- subset(average_slopes_by_cage_treatment_parent, parent.treatment == "control" & treatment == "treatment")

control_low <- subset(average_slopes_by_cage_treatment_parent, parent.treatment == "low" & treatment == "control")
treatment_low <- subset(average_slopes_by_cage_treatment_parent, parent.treatment == "low" & treatment == "treatment")

control_medium <- subset(average_slopes_by_cage_treatment_parent, parent.treatment == "medium" & treatment == "control")
treatment_medium <- subset(average_slopes_by_cage_treatment_parent, parent.treatment == "medium" & treatment == "treatment")

control_high <- subset(average_slopes_by_cage_treatment_parent, parent.treatment == "high" & treatment == "control")
treatment_high <- subset(average_slopes_by_cage_treatment_parent, parent.treatment == "high" & treatment == "treatment")

# Perform t-test for each parental treatment group
t_test_control <- t.test(control_control$avg_slope, treatment_control$avg_slope)
t_test_low <- t.test(control_low$avg_slope, treatment_low$avg_slope)
t_test_medium <- t.test(control_medium$avg_slope, treatment_medium$avg_slope)
t_test_high <- t.test(control_high$avg_slope, treatment_high$avg_slope)

# results
t_test_control # p = .99
t_test_low # p =.27
t_test_medium # p = .97
t_test_high # p = .15




#-------------------------------------------------------------------------------
# BOX PLOT BY EXPOSURES

exposure_plot = ggplot(average_slopes_by_cage_treatment_parent, aes(x = factor(num_hs_exposures, levels = as.character(1:15)), y = avg_slope, fill = treatment)) +
  geom_boxplot() +  
  labs(title = "Comparison of Average Slopes by Number of Heat Shock Exposures (Brood 1)",
       x = "Number of Heat Shock Exposures",
       y = "Average Growth Rate (mm²/week)") +
  scale_fill_manual(name = "Embryo Treatment",
                    values = c("control" = "skyblue3", "treatment" = "tomato"),
                    labels = c("control" = "Control", "treatment" = "Heat Shocked")) +  
  scale_x_discrete(limits = as.character(1:15)) +  # Ensure x-axis goes from 1 to 15
  theme_minimal() +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 14))  

# Display the plot
exposure_plot


#stats
two_way_anova <- aov(avg_slope ~ factor(num_hs_exposures) * treatment, data = average_slopes_by_cage_treatment_parent)
summary(two_way_anova)






#this was made using individual embryo slopes
view(growth_rate_regression)
cage_plot = ggplot(growth_rate_regression, aes(x = factor(cage), y = slope, fill = treatment)) +
  geom_boxplot() +  
  labs(title = "Comparison of Average Slopes by Moms (Brood 1)",
       x = "Cage Number",
       y = "Average growth rate (mm^2/ week)") +
  scale_fill_manual(name = "Embryo Treatment",
                    values = c("control" = "skyblue3", "treatment" = "tomato"),
                    labels = c("control" = "Control", "treatment" = "Heat Shocked")) +  
  theme_minimal() +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 14))  

# Display the plot
cage_plot



one_way_anova_cage <- aov(slope ~ factor(cage), data = growth_rate_regression)
summary(one_way_anova_cage)
TukeyHSD(one_way_anova_cage)
























# **For making one graph specified by brood # and parent treatment that follows embryo growth over time
# not super important to me 
filtered_growth <- master_growth %>%
  filter(parent.treatment == "medium", brood == 1, cage =="10", treatment == "treatment") %>%
  mutate(area_mm = as.numeric(area_mm))  

# Create the plot 
grow_plot = ggplot(filtered_growth, aes(x = date, y = area_mm, color = embryo, group = embryo)) +
  geom_line(size = 1, show.legend = FALSE) +  
  geom_point(size = 2, shape = 21, fill = "white", show.legend = FALSE) +  
  labs(title = "Embryo Growth Over Time: Medium Parent Treatment, Brood 1",
       x = "Date",
       y = "Area (mm²)") +
  scale_y_continuous(limits = c(0, 0.9)) +  
  theme_minimal()

# Display the static ggplot
grow_plot

# Convert the ggplot to an interactive plot
interactive_plot = ggplotly(grow_plot)

# Display the interactive plot
interactive_plot

# -----------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------

#-------------------------------------
#send to stillman
#-------------------------------------
# reading and formatting csv
library(tidyverse)
library(lubridate)  # For handling date conversion
library(ggplot2)
library(plotly)

# Read the CSV
master_growth_df = 'master_file_growth_yolk.csv'
master_growth = read.csv(master_growth_df)

# Convert date column to a uniform Date format (MM-DD-YY for display)
master_growth <- master_growth %>%
  mutate(date = mdy(date))  # Converts mixed formats (M/D/YY, MM-DD-YY) to Date
master_growth <- master_growth %>%
  mutate(area_mm = as.numeric(area_mm))  
master_growth <- master_growth %>%
  arrange(embryo, date)

# Check sorting
view(master_growth)




# Calculate growth rate as change in area per date and filter for brood
master_growth <- master_growth %>%
  filter(brood == 1) %>%
  group_by(embryo) %>%
  arrange(embryo, date) %>%  # Ensure dates are sorted
  mutate(growth_rate = (area_mm - lag(area_mm)) / as.numeric(difftime(date, lag(date), units = "days"))) %>%
  ungroup()
view(master_growth)



# Calculate the average growth rate for each embryo

average_growth_rate_data <- master_growth %>%
  group_by(embryo) %>%
  summarize(
    avg_growth_rate = mean(growth_rate, na.rm = TRUE),
    treatment = first(treatment),  # Retain the first treatment value for each embryo
    brood = first(brood),  # Retain the first brood value for each embryo
    cage = first(cage),  # Retain the first cage value for each embryo
    cage_name = first(cage.name),
    parent_treatment = first(parent.treatment),
    exposures = first(num_hs_exposures)# Retain the first cage.name value for each embryo
  ) %>%
  ungroup()

# View the updated data frame
view(average_growth_rate_data)



# calculate average growth rate of control/ treatment embryos in each cage
cage_avg_growth_rate <- master_growth %>%
  group_by(cage, treatment) %>%
  summarize(
    avg_growth_rate = mean(growth_rate, na.rm = TRUE),
    brood = first(brood),  # Retain the first brood value for each group
    cage_name = first(cage.name),  # Retain the first cage.name value for each group
    parental_treatment = first(parent.treatment),
    exposures = first(num_hs_exposures)# Retain the first parent.treatment value for each group
  ) %>%
  ungroup()

# View the results
view(cage_avg_growth_rate)

# Make a box plot 

cage_avg_growth_rate$parental_treatment <- factor(cage_avg_growth_rate$parental_treatment,
                                                  levels = c("control", "low", "medium", "high"))

# Create the combined box plot with organized parental treatment
average_gr = ggplot(cage_avg_growth_rate, aes(x = parental_treatment, y = avg_growth_rate, fill = treatment)) +
  geom_boxplot() +  # Create the box plot
  labs(title = "Growth Rates (brood 1)",
       x = "Parent Treatment Group",
       y = "Average Growth Rate (mm^2/week)",
       fill = "Embryo Treatment Group") +  # Add legend title
  scale_fill_manual(name = "Embryo Treatment Group", 
                    values = c("control" = "skyblue3", "treatment" = "tomato")) +  # Color for control and treatment
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),  # Center title
        axis.text.x = element_text(size = 11),  # Increase x-axis text size
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))  # Increase axis title size

average_gr


#--------------------------












#----------------------------------------------------------------

















library(plotly)

# Create an interactive plot using Plotly
interactive_plot <- ggplot(growth_rate_regression, aes(x = embryo, y = slope, color = embryo)) +
  geom_point() +  # Plot points for each embryo's growth rate
  labs(title = "Growth Rates (Slopes) for Each Embryo",
       x = NULL,  # Remove x-axis label
       y = "Growth Rate (mm^2/day)") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # Remove x-axis labels
        axis.ticks.x = element_blank())  # Remove x-axis ticks

# Convert ggplot to plotly for interactivity
interactive_plotly <- ggplotly(interactive_plot)

# Show the interactive plot
interactive_plotly






master_growth <- master_growth %>%
  filter(brood == 1) %>%
  group_by(embryo) %>%
  arrange(embryo, date) %>%  # Ensure dates are sorted
  mutate(growth_rate = (area_mm - lag(area_mm)) / as.numeric(difftime(date, lag(date), units = "days"))) %>%
  ungroup()

# Check if growth rates were calculated correctly
#view(master_growth)

# Test one embryo to make sure everything is right
#view(master_growth %>% filter(embryo == "A1.10"))

# Calculate the average growth rate for each embryo

average_growth_rate_data <- master_growth %>%
  group_by(embryo) %>%
  summarize(
    avg_growth_rate = mean(growth_rate, na.rm = TRUE),
    treatment = first(treatment),  # Retain the first treatment value for each embryo
    brood = first(brood),  # Retain the first brood value for each embryo
    cage = first(cage),  # Retain the first cage value for each embryo
    cage_name = first(cage.name),
    parent_treatment = first(parent.treatment),
    exposures = first(num_hs_exposures)# Retain the first cage.name value for each embryo
  ) %>%
  ungroup()

# View the updated data frame
view(average_growth_rate_data)

# Check for a specific embryo (e.g., "H10.19")
view(average_growth_rate_data %>% filter(embryo == "H10.19"))


# Ensure parental_treatment is a factor with the correct ordering
average_growth_rate_data$parental_treatment <- factor(average_growth_rate_data$parent_treatment, 
                                                      levels = c("control", "low", "medium", "high"))

# Create the box plot
ggplot(average_growth_rate_data, aes(x = parental_treatment, y = avg_growth_rate, fill = treatment)) +
  geom_boxplot() +  # Create the box plot
  labs(title = "Average Growth Rate by individual embryos",
       x = "Parental Treatment",
       y = "Average Growth Rate") +
  scale_fill_manual(values = c("control" = "skyblue3", "treatment" = "tomato")) +  # Colors for control and treatment
  theme_minimal() +
  theme(legend.title = element_blank())  # Remove legend title


#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------

# calculate average growth rate in each cage divided by treatment and control
# we are grouping embryos by 12 (control & treatment)

# Calculate average growth rate for control and treatment embryos within each cage, 

cage_avg_growth_rate <- master_growth %>%
  group_by(cage, treatment) %>%
  summarize(
    avg_growth_rate = mean(growth_rate, na.rm = TRUE),
    brood = first(brood),  # Retain the first brood value for each group
    cage_name = first(cage.name),  # Retain the first cage.name value for each group
    parental_treatment = first(parent.treatment),
    exposures = first(num_hs_exposures)# Retain the first parent.treatment value for each group
  ) %>%
  ungroup()

# View the results
view(cage_avg_growth_rate)

# Make a box plot 

# Create the box plot
# Create the combined box plot
# Ensure parental_treatment is a factor with the correct ordering
cage_avg_growth_rate$parental_treatment <- factor(cage_avg_growth_rate$parental_treatment,
                                                  levels = c("control", "low", "medium", "high"))

# Create the combined box plot with organized parental treatment
average_gr = ggplot(cage_avg_growth_rate, aes(x = parental_treatment, y = avg_growth_rate, fill = treatment)) +
  geom_boxplot() +  # Create the box plot
  labs(title = "Growth Rates (brood 1)",
       x = "Parent Treatment Group",
       y = "Average Growth Rate (mm^2/week)",
       fill = "Embryo Treatment Group") +  # Add legend title
  scale_fill_manual(name = "Embryo Treatment Group", 
                    values = c("control" = "skyblue3", "treatment" = "tomato")) +  # Color for control and treatment
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),  # Center title
        axis.text.x = element_text(size = 11),  # Increase x-axis text size
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))  # Increase axis title size

average_gr

#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
# BOX PLOT BY EXPOSURES

# Ensure 'exposures' is a factor (if it's not already)
# Ensure 'exposures' is a factor
cage_avg_growth_rate$exposures <- factor(cage_avg_growth_rate$exposures)
view(cage_avg_growth_rate)


# Create the updated box plot
exposure = ggplot(cage_avg_growth_rate, aes(x = exposures, y = avg_growth_rate, fill = treatment)) +
  geom_boxplot() +  # Create the box plot
  labs(title = "Average Growth Rate by Number of Exposures",
       x = "Number of parental heatshock exposures",
       y = "Average Growth Rate (mm^2/week)") +
  scale_fill_manual(values = c("control" = "skyblue3", "treatment" = "tomato")) +  # Colors for control and treatment
  theme_minimal() +
  theme(
    legend.title = element_blank(),  # Remove legend title
    plot.title = element_text(hjust = 0.5)  # Center title
  )

exposure

#--------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------
# STATS


# one way anova: compare growth rates of CONTROL embryos across parental treatment groups

# Subset only control embryos
control_embryos <- subset(cage_avg_growth_rate, treatment == "control")
view(control_embryos)
# Check normality assumption
shapiro.test(control_embryos$avg_growth_rate)

# Check homogeneity of variance
library(car)
leveneTest(avg_growth_rate ~ parental_treatment, data = control_embryos)

# Run one-way ANOVA
anova_control <- aov(avg_growth_rate ~ parental_treatment, data = control_embryos)
summary(anova_control)

# Post-hoc Tukey test (if ANOVA is significant)
TukeyHSD(anova_control)



# one way anova: compare growth rates of TREATMENT embryos across parental treatment group

treatment_embryos <- subset(cage_avg_growth_rate, treatment == "treatment")
view(treatment_embryos)
# Check normality assumption
shapiro.test(treatment_embryos$avg_growth_rate)

#Was not normally distributed, so running this test
kruskal.test(avg_growth_rate ~ parental_treatment, data = treatment_embryos)
# no significance 




# is there a difference between control and treatment growth rates within parental groups?

# Subset data for each parental treatment group
control_control <- subset(cage_avg_growth_rate, parental_treatment == "control" & treatment == "control")
treatment_control <- subset(cage_avg_growth_rate, parental_treatment == "control" & treatment == "treatment")

control_low <- subset(cage_avg_growth_rate, parental_treatment == "low" & treatment == "control")
treatment_low <- subset(cage_avg_growth_rate, parental_treatment == "low" & treatment == "treatment")

control_medium <- subset(cage_avg_growth_rate, parental_treatment == "medium" & treatment == "control")
treatment_medium <- subset(cage_avg_growth_rate, parental_treatment == "medium" & treatment == "treatment")

control_high <- subset(cage_avg_growth_rate, parental_treatment == "high" & treatment == "control")
treatment_high <- subset(cage_avg_growth_rate, parental_treatment == "high" & treatment == "treatment")

# Perform t-test for each parental treatment group
t_test_control <- t.test(control_control$avg_growth_rate, treatment_control$avg_growth_rate)
t_test_low <- t.test(control_low$avg_growth_rate, treatment_low$avg_growth_rate)
t_test_medium <- t.test(control_medium$avg_growth_rate, treatment_medium$avg_growth_rate)
t_test_high <- t.test(control_high$avg_growth_rate, treatment_high$avg_growth_rate)

# results
t_test_control # p = .97
t_test_low # p =.17
t_test_medium # p = .99
t_test_high # p = .13




