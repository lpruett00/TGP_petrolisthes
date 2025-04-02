# Load necessary libraries
library(tidyverse)
library(lubridate)  # For handling date conversion
library(ggplot2)
library(plotly)
library(nlme)
library(multcomp)
library(emmeans)

# Read the CSV
master_growth_df = 'growth_rate_regression_embryo.csv'
embryo_growth = read.csv(master_growth_df)
view(embryo_growth)


#--------------------------------------

#------------------------------------------------


# Three way ANOVA (embryo nested in cage) using embryo_growth
model_lme <- lme(slope ~ treatment * parent.treatment, 
                 random = ~1 | cage/embryo, 
                 data = embryo_growth)

# View the model summary
summary(model_lme)
anova(model_lme)

# Tukey test

emm <- emmeans(model_lme, pairwise ~ treatment * parent.treatment, adjust = "tukey")
emm$contrasts


# Perform pairwise comparisons between parental treatments with Tukey adjustment
parent_emm <- emmeans(model_lme, pairwise ~ parent.treatment, adjust = "tukey")
parent_emm$contrasts




# num exposures


# Fit a linear mixed-effects model including num_hs_exposures and treatment as fixed effects
model_lme <- lme(slope ~ num_hs_exposures * treatment, 
                 random = ~1 | cage/embryo, 
                 data = embryo_growth)

# View the model summary
summary(model_lme)

# Run an ANOVA on the model
anova(model_lme)



# by mother

# Fit the linear mixed-effects model with embryos nested within cages
model_lme_mother <- lme(slope ~ 1, 
                        random = ~1 | cage/embryo, 
                        data = embryo_growth)

# View the model summary
summary(model_lme_mother)

# Run an ANOVA on the model to assess the significance of maternal effects
anova(model_lme_mother)




# All broods
all_embryo_regression = read.csv("growth_rate_regression_embryo_allbroods.csv")
view(all_embryo_regression)

# Filter data for specific cages
filtered_data <- all_embryo_regression %>% 
  filter(cage %in% c(4, 8, 7, 19, 10))

str(filtered_data)


# Fit the linear mixed-effects model
model_lmer <- lmer(slope ~ treatment * parent.treatment * brood  + (1 | cage/embryo), 
                   data = filtered_data)

# Display the summary of the model
summary(model_lmer)

# Perform ANOVA on the model
anova(model_lmer)



filtered_data <- all_embryo_regression %>% 
  filter(cage %in% c(4, 8, 7, 19, 10) & brood != 8)

# Check the structure of the filtered data
str(filtered_data)

# Fit the linear mixed-effects model
model_lmer <- lmer(slope ~ treatment * parent.treatment * brood * cage + (1 | cage/embryo), 
                   data = filtered_data)

# Display the summary of the model
summary(model_lmer)

# Perform ANOVA on the model
anova(model_lmer)






