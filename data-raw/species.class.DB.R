library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)

data.file <- "./data/ecological species traits final v5.xlsx"
data <- read_xlsx(data.file)
species.classification <- data %>%
  dplyr::select(acceptedname,temperament,succession_type) %>%
  mutate(species = tolower(acceptedname))

usethis::use_data(species.classification, overwrite = TRUE)
