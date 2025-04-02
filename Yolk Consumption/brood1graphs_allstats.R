# Load necessary libraries
library(tidyverse)
library(lubridate)  # For handling date conversion
library(ggplot2)
library(plotly)
master_growth_df = 'master_file_growth_yolk.csv'
master_growth = read.csv(master_growth_df)
# Convert date column to a uniform Date format (MM-DD-YY for display)
master_growth <- master_growth %>%
  mutate(date = mdy(date))  # Converts mixed formats (M/D/YY, MM-DD-YY) to Date
master_growth <- master_growth %>%
  arrange(embryo, date)
master_growth <- master_growth %>%
  filter(yolk.cover.percentage != "", !is.na(yolk.cover.percentage), !is.na(date))  # Remove empty strings and NA values
view(master_growth)
#-------------------------------------------------------------------------------

# Convert yolk.cover.percentage to numeric, replacing non-numeric values with NA
master_growth <- master_growth %>%
  mutate(yolk.cover.percentage = as.numeric(yolk.cover.percentage))

# Check for NAs after conversion
sum(is.na(master_growth$yolk.cover.percentage))


# Calculate yolk consumption regression slope for each embryo in brood 1
yolk_consumption_regression <- master_growth %>%
  filter(brood == 1, !is.na(yolk.cover.percentage), 
         !is.nan(yolk.cover.percentage), 
         !is.infinite(yolk.cover.percentage), 
         !is.na(date)) %>%  # Only brood 1, no invalid values
  group_by(embryo) %>%
  filter(n() > 1) %>%  # Ensure at least 2 data points per embryo
  do({
    model <- lm(yolk.cover.percentage ~ date, data = .)  # Fit regression
    # Keep all original columns and add the regression coefficients
    cbind(., 
          slope = coef(model)[2],  # Extract slope
          intercept = coef(model)[1],  
          r_squared = summary(model)$r.squared)  
  }) %>%
  ungroup()

# Remove duplicates based on unique combination of all columns
yolk_consumption_regression <- yolk_consumption_regression %>%
  distinct(embryo, .keep_all = TRUE)

# View the data with duplicates removed
view(yolk_consumption_regression)


#------------------------------------------
#stats for brood 1

library(nlme)

yolk_lmm_lme <- lme(slope ~ treatment * parent.treatment, 
                    random = ~1 | cage/embryo, 
                    data = yolk_consumption_regression)
summary(yolk_lmm_lme)
anova(yolk_lmm_lme)

yolk_lmm_lme <- lme(slope ~ treatment * cage, 
                    random = ~1 | cage/embryo, 
                    data = yolk_consumption_regression)



anova(yolk_lmm_lme)




#-------------------------------------------

# stats across broods 


all_embryo_data = read.csv('yolk_consumption_allembryos.csv')
view(all_embryo_data)
# run stats 

library(lmerTest)

all_embryo_data$brood <- as.factor(all_embryo_data$brood)
anova_model <- lmer(slope ~ brood * treatment * parent.treatment + (1 | cage/embryo), 
                    data = all_embryo_data)
anova(anova_model)

tukey_results <- emmeans(anova_model, ~ brood * treatment * parent.treatment)
tukey_test <- pairs(tukey_results, adjust = "tukey")

# View the Tukey test results
summary(tukey_test)

emm_results <- emmeans(anova_model, ~ brood * treatment * parent.treatment)
pairwise_comparisons <- contrast(emm_results, method = "pairwise")
summary(pairwise_comparisons)


#-------------------------------------------
# graphs       

# Average slopes for yolk consumption by cage, treatment, and parent treatment
average_slopes_by_cage_treatment_parent_yolk <- yolk_consumption_regression %>%
  group_by(cage, treatment, parent.treatment) %>%
  summarize(
    avg_yolk_slope = mean(slope, na.rm = TRUE),  # Calculate average slope for yolk consumption
    num_hs_exposures = unique(num_hs_exposures),  # Assuming num_hs_exposures is the same for each embryo
    .groups = "drop"  # Drop grouping after summarizing
  )

# View the results
view(average_slopes_by_cage_treatment_parent_yolk)



# Ensure parental_treatment is a factor with the correct order
average_slopes_by_cage_treatment_parent_yolk$parent.treatment <- factor(average_slopes_by_cage_treatment_parent_yolk$parent.treatment, 
                                                                   levels = c("control", "low", "medium", "high"))

# Ensure treatment is a factor for proper ordering in the plot
average_slopes_by_cage_treatment_parent_yolk$treatment <- factor(average_slopes_by_cage_treatment_parent_yolk$treatment, 
                                                            levels = c("control", "treatment"))

# Convert avg_slope to positive by taking the absolute value
average_slopes_by_cage_treatment_parent_yolk$avg_yolk_slope <- abs(average_slopes_by_cage_treatment_parent_yolk$avg_yolk_slope)

# Create the box plot
ggplot(average_slopes_by_cage_treatment_parent_yolk, aes(x = parent.treatment, y = avg_yolk_slope, fill = treatment)) +
  geom_boxplot() +  # Create the box plot
  labs(title = "Average yolk consumption rates by Parental Treatment Brood 1",
       x = "Parental Treatment",
       y = "Average Yolk consumption rates (% consumed/ day)") +
  scale_fill_manual(values = c("control" = "skyblue3", "treatment" = "tomato")) +  # Colors for control and treatment
  scale_x_discrete(drop = FALSE) +  # Force display of all factor levels, even empty ones
  theme_minimal() +
  theme(legend.title = element_blank(),        # Remove legend title
        plot.title = element_text(hjust = 0.5),  # Center the title
        axis.text.x = element_text(size = 11),   # Increase x-axis text size
        axis.title.x = element_text(size = 12),    # Increase x-axis title size
        axis.title.y = element_text(size = 12))    # Increase y-axis title size



ggplot(average_slopes_by_cage_treatment_parent_yolk, 
       aes(x = factor(num_hs_exposures, levels = as.character(1:15)), 
           y = avg_yolk_slope, 
           fill = treatment)) +
  geom_boxplot() +  
  labs(title = "Average Yolk Consumption Rates by Number of Heat Shock Exposures (Brood 1)",
       x = "Number of Heat Shock Exposures",
       y = "Average Yolk Consumption Rate (% consumed/day)") +
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




ggplot(average_slopes_by_cage_treatment_parent_yolk, 
       aes(x = factor(cage), 
           y = avg_yolk_slope, 
           fill = treatment)) +
  geom_boxplot() +  
  labs(title = "Average Yolk Consumption Rates by Cage (Brood 1)",
       x = "Cage Number",
       y = "Average Yolk Consumption Rate (% consumed/day)") +
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




# Calculate standard deviations for each parental treatment and embryo treatment
sd_summary <- average_slopes_by_cage_treatment_parent_yolk %>%
  group_by(parent.treatment, treatment) %>%
  summarize(
    sd_yolk_slope = sd(avg_yolk_slope, na.rm = TRUE),
    .groups = "drop"
  )

# View results
print(sd_summary)




#testing varience
library(car)
# Levene's Test for homogeneity of variance
leveneTest(avg_yolk_slope ~ parent.treatment, data = average_slopes_by_cage_treatment_parent_yolk)
# Bartlett's Test for homogeneity of variance (normality assumption)
bartlett.test(avg_yolk_slope ~ parent.treatment, data = average_slopes_by_cage_treatment_parent_yolk)









#OLD PRELIMINARY STATS

# one way anova: compare growth rates of CONTROL embryos across parental treatment groups
# Subset only control embryos
control_embryos <- subset(average_slopes_by_cage_treatment_parent_yolk, treatment == "control")
shapiro.test(control_embryos$avg_yolk_slope) #data is noramlly distributed
library(car)
leveneTest(avg_yolk_slope ~ parent.treatment, data = control_embryos)
anova_control <- aov(avg_slope ~ parent.treatment, data = control_embryos)
summary(anova_control)

treatment_embryos <- subset(average_slopes_by_cage_treatment_parent_yolk, treatment == "treatment")
view(treatment_embryos)
# Check normality assumption
shapiro.test(treatment_embryos$avg_yolk_slope)
leveneTest(avg_yolk_slope ~ parent.treatment, data = treatment_embryos)
anova_treatment <- aov(avg_yolk_slope ~ parent.treatment, data = treatment_embryos)
summary(anova_treatment)



average_slopes_by_cage_treatment_parent_yolk$treatment <- factor(average_slopes_by_cage_treatment_parent_yolk$treatment)
average_slopes_by_cage_treatment_parent_yolk$parent.treatment <- factor(average_slopes_by_cage_treatment_parent_yolk$parent.treatment)

# Two-way ANOVA (parent.treatment * treatment)
two_way_anova <- aov(avg_yolk_slope ~ parent.treatment * treatment, data = average_slopes_by_cage_treatment_parent_yolk)

# View the summary of the ANOVA
summary(two_way_anova)



# is there a difference between control and treatment within each parent group
control_control <- subset(average_slopes_by_cage_treatment_parent_yolk, parent.treatment == "control" & treatment == "control")
treatment_control <- subset(average_slopes_by_cage_treatment_parent_yolk, parent.treatment == "control" & treatment == "treatment")

control_low <- subset(average_slopes_by_cage_treatment_parent_yolk, parent.treatment == "low" & treatment == "control")
treatment_low <- subset(average_slopes_by_cage_treatment_parent_yolk, parent.treatment == "low" & treatment == "treatment")

control_medium <- subset(average_slopes_by_cage_treatment_parent_yolk, parent.treatment == "medium" & treatment == "control")
treatment_medium <- subset(average_slopes_by_cage_treatment_parent_yolk, parent.treatment == "medium" & treatment == "treatment")

control_high <- subset(average_slopes_by_cage_treatment_parent_yolk, parent.treatment == "high" & treatment == "control")
treatment_high <- subset(average_slopes_by_cage_treatment_parent_yolk, parent.treatment == "high" & treatment == "treatment")

# Perform t-test for each parental treatment group
t_test_control <- t.test(control_control$avg_yolk_slope, treatment_control$avg_yolk_slope)
t_test_low <- t.test(control_low$avg_yolk_slope, treatment_low$avg_yolk_slope)
t_test_medium <- t.test(control_medium$avg_yolk_slope, treatment_medium$avg_yolk_slope)
t_test_high <- t.test(control_high$avg_yolk_slope, treatment_high$avg_yolk_slope)

# results
t_test_control # p = .96
t_test_low # p =.26
t_test_medium # p = .11
t_test_high # p = .24











