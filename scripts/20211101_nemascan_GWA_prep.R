library(tidyverse)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

#======================================================#
#                load data for run on NemaScan
#   commit = dc38bb0e4516ae737c9c4851df6fe75c2917fd1e                      
#======================================================#
ls_df <- data.table::fread("data/20191107_lifespan_GWA_phenotypes.csv") %>%
  dplyr::mutate(group = ifelse(group == "Censor", "censor", group),
                censor = ifelse(group == "censor", TRUE, FALSE))

strain_id_df <- data.table::fread("data/20190801_wormcamp_strainIDs_processed.tsv") %>%
  dplyr::mutate(wellNum = case_when(well == "D06" ~ 1, well == "C06" ~ 2, well == "B06" ~ 3, well == "A06" ~ 4,
                                    well == "D05" ~ 5, well == "C05" ~ 6, well == "B05" ~ 7, well == "A05" ~ 8,
                                    well == "D04" ~ 9, well == "C04" ~ 10,well == "B04" ~ 11,well == "A04" ~ 12,
                                    well == "D03" ~ 13,well == "C03" ~ 14,well == "B03" ~ 15,well == "A03" ~ 16,
                                    well == "D02" ~ 17,well == "C02" ~ 18,well == "B02" ~ 19,well == "A02" ~ 20,
                                    well == "D01" ~ 21,well == "C01" ~ 22,well == "B01" ~ 23,well == "A01" ~ 24),
                plateName = as.character(glue::glue('WI_{wormcamp_plate_num}')))


# Join data to include strain names for censored wells
ls_df_2 <- full_join(ls_df, strain_id_df) %>%
  dplyr::mutate(plateName = factor(plateName, levels = c("WI_1", "WI_2", "WI_3", "WI_4",  "WI_5", 
                                                         "WI_6",  "WI_7",  "WI_8",  "WI_9", "WI_10",
                                                         "WI_11", "WI_12", "WI_13", "WI_14", "WI_15",
                                                         "WI_16", "WI_17", "WI_18", "WI_19", "WI_20")))

# generate trait files for mappings based on sliding censoring threshold. No regression or use of control N2 strains
ls_df_CFY_tf  <- ls_df_2 %>%
  dplyr::filter(censor == FALSE) %>%
  dplyr::mutate(t97_t89_diff = T97 - T89) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(t97_mean = mean(T97),
                t89_mean = mean(T89),
                t97_t89_diff_mean = mean(t97_t89_diff),
                n = n()) %>%
  dplyr::select(strain, t97_mean, t89_mean, t97_t89_diff_mean) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(strain, .keep_all = TRUE) %>%
  dplyr::filter(!(strain %in% c("N2_1", "N2_2", "N2_3", "N2_4", "N2_5")))

# export the lifespan data
rio::export(ls_df_CFY_tf, 'data/traitfile_ls_GWA_CFY.tsv', format = "tsv")


# NemaScan command 20211101
# nextflow run andersenlab/nemascan --vcf 20210121 --traitfile /projects/b1059/projects/Tim/Lifespan/traitfile_ls_GWA_CFY.tsv --sthresh EIGEN -r dc38bb0e4516ae737c9c4851df6fe75c2917fd1e

