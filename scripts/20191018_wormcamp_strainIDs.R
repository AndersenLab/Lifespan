#!/usr/bin/env Rscript
#Load necessary packages
library(tidyverse)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load raw strain ids and reorganize
csv <- readr::read_csv("data/20190801_wormcamp_strainIDs.csv") %>%
  dplyr::select(strain, wormcamp_plate_num, well, random_ID = rand_ID, everything(), -`Vol to add for 1000`)
csv[csv == "#N/A"] <- NA
csv <- csv %>%
  dplyr::mutate(strain = ifelse(is.na(strain) & random_ID == "N2", "N2", strain))

# write tsv for data analysis
readr::write_tsv(csv, "data/20190801_wormcamp_strainIDs_processed.tsv")
