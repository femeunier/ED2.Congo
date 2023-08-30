pftmapping <- read.csv("data-raw/pftmapping.csv", sep = ";")

usethis::use_data(pftmapping, overwrite = TRUE)
