# Load necessary libraries
library(tidyverse)
install.packages("writexl")  
library(writexl)

# Read the data from the CSV file you are converting
file_name_one = "Medium_1.csv"
data_one = read.csv(file_name_one)
data_one <- data_one %>%
  mutate(embryo_area_numeric = as.numeric(as.character(embryo.area..pixels.))) 

# Check for any warnings after conversion
if (any(is.na(data_one$embryo_area_numeric))) {
  warning("Some values in 'embryo.area..pixels.' could not be converted to numeric.")
}

# Create the "area mm" column by converting the area in pixels to mm² based on magnification
data_one <- data_one %>%
  mutate(area_mm = case_when(
    magnification == 13.5 ~ embryo_area_numeric / 580644,
    magnification == 12 ~ embryo_area_numeric / 454276,
    magnification == 10 ~ embryo_area_numeric / 326041,
    magnification == 8 ~ embryo_area_numeric / 208849,
    TRUE ~ NA_real_  # Assign NA if magnification doesn't match any of the conditions
  ))

# Save  data to an Excel file
write_xlsx(data_one, "medium_1_area.xlsx")



