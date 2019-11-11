try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
library(tidyverse)
library(survival)
library(flexsurv)
library(ggplot2)
library(ggfortify)
library(cowplot)
library(ggplot2)

# Original processing to reformat data to long form.
# Load csv file

#raw <- as.data.frame(read.csv("20171214_Raw_br1_tr2.csv", header = T))

# reshape to long format

#df <- raw %>%
#dplyr::group_by(strain) %>%
#dplyr::mutate(rep = rep(1:n())) %>%
#tidyr::gather(activity_type, activity, -strain, -life_span, -rep) %>%
#tidyr::separate(col = activity_type, c("activity_type", "time_d"), sep = "_._") %>%
#dplyr::mutate(bio_rep = 1,
#             tech_rep = 2) %>%
#dplyr::select(strain, tech_rep, bio_rep, rep, life_span, activity_type, time_d, activity) %>%
#dplyr::arrange(strain, tech_rep, bio_rep, rep)

#write.csv(df, file = "20171214_Processed_br1_tr2.csv", row.names = F)

# repeat above script for number of files to generate processed files for all biological and technical replicates
# bind individual dataframes together 
df_1 <- as.data.frame(read.csv(file = "20171214_Processed_br1_tr2.csv", header = T))
df_2 <- as.data.frame(read.csv(file = "20171214_Processed_br2_tr2.csv", header = T))
df_3 <- as.data.frame(read.csv(file = "20171214_Processed_br3_tr1.csv", header = T))
df_4 <- as.data.frame(read.csv(file = "20171214_Processed_br4_tr2.csv", header = T))
df_5 <- as.data.frame(read.csv(file = "20171214_Processed_br5_tr1.csv", header = T))
df_6 <- as.data.frame(read.csv(file = "20171214_Processed_br6_tr2.csv", header = T))
df <- bind_rows(df_1, df_2, df_3, df_4, df_5, df_6)
write.csv(df, file = "20170116_FullProcessed.csv", row.names = F)