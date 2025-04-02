
library(tidyverse)
library(lubridate)  
library(ggplot2)
library(plotly)


master_growth_df = 'master_file_growth_yolk.csv'
master_growth = read.csv(master_growth_df)
# Convert date column to a uniform Date format 
master_growth <- master_growth %>%
  mutate(date = mdy(date))  # Converts mixed formats (M/D/YY, MM-DD-YY) to Date
master_growth <- master_growth %>%
  mutate(area_mm = as.numeric(area_mm))  
# Sort the data by embryo and date
master_growth <- master_growth %>%
  arrange(embryo, date)

view(master_growth)

#---------------------------------------------------------------------------------------------------------

#find average growth rates of control and treatment embryos by cage

#---------------------------------------------------------------------------------------------------------


# Calculate the average growth rate for each embryo in Brood 1
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

view(growth_rate_regression) #save as xcel and csv


# Check for a specific embryo 
view(growth_rate_regression %>% filter(embryo == "A1.10"))

growth_rate_regression = read.csv('growth_rate_regression_embryo.csv')

average_slopes_by_cage_treatment_parent <- growth_rate_regression %>%
  group_by(cage, treatment, parent.treatment) %>%
  summarize(
    avg_slope = mean(slope, na.rm = TRUE),  # Calculate avg slope per group
    num_hs_exposures = unique(num_hs_exposures),  
    .groups = "drop"
  )  
average_slopes_by_cage_treatment_parent <- growth_rate_regression %>%
  # Join with master_growth to get treatment info and parent_treatment
  group_by(cage, treatment, parent.treatment) %>%  
  summarize(
    avg_slope = mean(slope, na.rm = TRUE),  
    sd_slope = sd(slope, na.rm = TRUE),  
    num_hs_exposures = unique(num_hs_exposures),  
    .groups = "drop"
  ) %>%
  ungroup()


#plot average growth rates

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



#this was made using individual embryo slopes to make growth rates by cage
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


#-----------------------------------------------
#scatter plot
#-----------------------------------------------
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


#-------------------------------------------------------------------------------






#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------
# playing around with stats, not using this data
#-----------------------------------------------

# one way anova: compare growth rates of CONTROL embryos across parental treatment group
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
