try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
library(tidyverse)
library(survival)
library(flexsurv)
library(ggplot2)
library(ggfortify)
library(cowplot)



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

# reshape data for survival object, ignore technical replicate and include just biological reps
df_km <- df %>%
  dplyr::distinct(strain, bio_rep, rep, life_span) %>%
  dplyr::mutate(status = ifelse(is.na(life_span), 0, 1)) %>%
  dplyr::mutate(bio_rep = factor(bio_rep))

#Asign functions to calculate 90, 50, 10 percent survival times with Weibull model.
weibull_50 <- function(shape, scale) {
  qweibull(0.5, shape = shape, scale = scale)
}
weibull_10 <- function(shape, scale) {
  qweibull(0.9, shape = shape, scale = scale)
}
weibull_90 <- function(shape, scale) {
  qweibull(0.1, shape = shape, scale = scale)
}

# establish lists to populate with for loop
block_list <- list()
strain_list <- list()

# nested for loop to extract parameters for each strain in each block
for(b in 1:length(unique(df_km$bio_rep))){
  
  block_df <- df_km %>%
    dplyr::filter(bio_rep == (unique(df_km$bio_rep)[b]))

  for(s in 1:length(unique(block_df$strain))){
  
    strain_df <- block_df %>%
      dplyr::filter(strain == (unique(block_df$strain)[s]))
  
    w_fit = flexsurvreg(Surv(life_span, status) ~ 1, dist = "weibull", data = strain_df)
  
    w_90 <- summary(w_fit, fn = weibull_90, t=1)
  
    w_50 <- summary(w_fit, fn = weibull_50, t=1)
  
    w_10 <- summary(w_fit, fn = weibull_10, t=1)
  
    strain_surv <- data.frame(block = unique(block_df$bio_rep),
                  strain = unique(block_df$strain)[s],
                  w_shape_est = w_fit[["res"]][1],
                  w_scale_est = w_fit[["res"]][2],
                  w_surv90_est = w_90[[1]][2],
                  w_surv50_est = w_50[[1]][2],
                  w_surv10_est = w_10[[1]][2])
    strain_surv <- strain_surv %>%
      dplyr::rename('w_surv90_est' = est, 'w_surv50_est' = est.1, 'w_surv10_est' = est.2)

    strain_list[[s]] <- strain_surv
  }
  
  strain_df <- bind_rows(strain_list)
  
  block_list[[b]] <- strain_df
  
}

#combine all data
h_df <- dplyr::bind_rows(block_list) %>%
  dplyr::arrange(strain, block)

#calculate heritability






### Old code not necessary for heritability calculation
# bulid survival object and plot average for all strains
km <- with(df_km, Surv(life_span, status))
km_fit <- survfit(Surv(life_span, status) ~ 1, data=df_km)
summary(km_fit, times = seq(1,34, 1))
autoplot(km_fit)


# replot with individual all individual strains and just CB4856 and N2
km_fit_strains <- survfit(Surv(life_span, status) ~ strain, data=df_km)
summary(km_fit_strains, times = seq(1,34, 1))
km_fit_strains_plot <- autoplot(km_fit_strains)


df_km_cbn2 <- with(df_km %>% dplyr::filter(grepl("N2", strain) | grepl("CB4856", strain)), Surv(life_span, status))
km_fit_cbn2 <- survfit(Surv(life_span, status) ~ strain, data=df_km %>% dplyr::filter(grepl("N2", strain) | grepl("CB4856", strain)))
summary(km_fit_cbn2, times = seq(1,34, 1))
km_fit_cbn2_plot <- autoplot(km_fit_cbn2)

km_summary_plot <- cowplot::plot_grid(km_fit_strains_plot, km_fit_cbn2_plot)

# separate all bio reps and plot N2 and CB
df_km <- df %>%
  dplyr::distinct(strain, tech_rep, bio_rep, rep, .keep_all = T ) %>%
  dplyr::mutate(status = ifelse(is.na(life_span), 0, 1))

df_km_1 <- df_km %>%
  dplyr::filter(grepl("N2", strain) | grepl("CB4856", strain),
                bio_rep == 1)

df_km_2 <- df_km %>%
  dplyr::filter(grepl("N2", strain) | grepl("CB4856", strain),
                bio_rep == 2)

df_km_3 <- df_km %>%
  dplyr::filter(grepl("N2", strain) | grepl("CB4856", strain),
                bio_rep == 3)

df_km_4 <- df_km %>%
  dplyr::filter(grepl("N2", strain) | grepl("CB4856", strain),
                bio_rep == 4)

df_km_5 <- df_km %>%
  dplyr::filter(grepl("N2", strain) | grepl("CB4856", strain),
                bio_rep == 5)

df_km_6 <- df_km %>%
  dplyr::filter(grepl("N2", strain) | grepl("CB4856", strain),
                bio_rep == 6)

km_fit_1 <- survfit(Surv(life_span, status) ~ strain, data = df_km_1)
summary(km_fit_1, times = seq(1,34, 1))
km1 <- autoplot(km_fit_1)

km_fit_2 <- survfit(Surv(life_span, status) ~ strain, data = df_km_2)
summary(km_fit_2, times = seq(1,34, 1))
km2 <- autoplot(km_fit_2)

km_fit_3 <- survfit(Surv(life_span, status) ~ strain, data = df_km_3)
summary(km_fit_3, times = seq(1,34, 1))
km3 <- autoplot(km_fit_3)

km_fit_4 <- survfit(Surv(life_span, status) ~ strain, data = df_km_4)
summary(km_fit_4, times = seq(1,34, 1))
km4 <- autoplot(km_fit_4)

km_fit_5 <- survfit(Surv(life_span, status) ~ strain, data = df_km_5)
summary(km_fit_5, times = seq(1,34, 1))
km5 <- autoplot(km_fit_5)

km_fit_6 <- survfit(Surv(life_span, status) ~ strain, data = df_km_6)
summary(km_fit_6, times = seq(1,34, 1))
km6 <- autoplot(km_fit_6)

km_fullplot <- plot_grid(km1,km2,km3,km4,km5,km6)

# replot with individual strains per bio rep
km_1 <- with(df_km %>% dplyr::filter(bio_rep == 1), Surv(life_span, status))
km_fit_1 <- survfit(Surv(life_span, status) ~ strain, data=df_km %>% dplyr::filter(bio_rep == 1))
summary(km_fit_1, times = seq(1,34, 1))
autoplot(km_fit_1)

# Stefan's way: process data and plot survival curves, assess if there are block effects. Then extract traits from curves.
# this may not be appropriate for individual worms in wells b/c independent reps are not handled correctly.
# furthermore, there are unequal timepoints among experiments/bioreps so some time points have just one value. Weights unequal.

#remove base activity b/c it's not needed here, only sitm activity determines live or dead
df_proc1 <- df %>%
  dplyr::filter(activity_type != "base", !is.na(life_span)) %>%
  dplyr::mutate(alive = ifelse(!is.na(activity), 1, 0)) %>%
  dplyr::group_by(strain, tech_rep, bio_rep, time_d) %>%
  dplyr::summarise(frac_alive = sum(alive)/n(),
                   n = n())

# plot survival curves for all strains
raw_plot_all <- ggplot(df_proc1) +
  aes(x=time_d, y=frac_alive, color=strain) +
  geom_line() +
  labs(x="Day", y="Fraction alive") +
  theme(legend.position = "right") +
  theme_grey() +
  facet_grid(bio_rep ~ tech_rep, labeller = label_both)

# plot survival curves for just N2 and CB4856
raw_plot_subset <- ggplot(df_proc1 %>% dplyr::filter(grepl("N2", strain) | grepl("CB4856", strain))) +
  aes(x=time_d, y=frac_alive, color=strain) +
  geom_line() +
  labs(x="Day", y="Fraction alive") +
  theme(legend.position = "right") +
  theme_grey() + 
  facet_grid(bio_rep ~ tech_rep, labeller = label_both)

# Large block effects are present!
