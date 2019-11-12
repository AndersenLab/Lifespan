#!/usr/bin/env Rscript
#Load necessary packages
library(tidyverse)
library(rio)
library(viridis)
library(boot)
# sommer package is not installing from github -> library(sommer)


# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

########################
### define functions ###
########################
# Heritability
# data is data frame that contains strain and phenotype column
# indicies are used by the boot function to sample from the 'data' data.frame
H2.test.boot <- function(data, indicies){
  
  d <- data[indicies,]
  
  pheno <- as.data.frame(dplyr::select(d, phenotype))[,1]
  strain <- as.factor(d$strain)
  
  reffMod <- lme4::lmer(pheno ~ 1 + (1|strain))
  
  Variances <- as.data.frame(lme4::VarCorr(reffMod, comp = "Variance"))
  
  Vg <- Variances$vcov[1]
  Ve <- Variances$vcov[2]
  H2 <- Vg/(Vg+Ve)
  
  # errors <- sqrt(diag(lme4::VarCorr(reffMod, comp = "Variance")$strain))
  
  return(H2)
}
# data is data frame that contains strain and phenotype column
H2.test <- function(data){
  
  pheno <- as.data.frame(dplyr::select(data, phenotype))[,1]
  strain <- as.factor(data$strain)
  
  reffMod <- lme4::lmer(pheno ~ 1 + (1|strain))
  
  Variances <- as.data.frame(lme4::VarCorr(reffMod, comp = "Variance"))
  
  Vg <- Variances$vcov[1]
  Ve <- Variances$vcov[2]
  H2 <- Vg/(Vg+Ve)
  
  # errors <- sqrt(diag(lme4::VarCorr(reffMod, comp = "Variance")$strain))
  
  return(H2)
}

# df is data frame that contains strain and phenotype column
H2.calc <- function(df, boot = T){
  
  
  if(boot == T){
    # bootstrapping with 1000 replications
    results <- boot(data=df, statistic=H2.test.boot, R=500)
    
    # get 95% confidence interval
    ci <- boot.ci(results, type="bca")
    
    H2_errors <- data.frame(H2 = ci$t0, ci_l = ci$bca[4], ci_r = ci$bca[5])
    
    return(H2_errors)
    
  } else {
    
    H2 <- data.frame(H2 = H2.test(data = df), ci_l = NA, ci_r = NA)
    return(H2)
  }
  
}

# remove outliers function (tukey's fences)
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}
 


#################
### load data ###
#################
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
                
                
# design plot data shaping
ls_df_2_design  <- ls_df_2 %>%
  dplyr::mutate(row = ifelse(wellNum %in% c(1:6), "A",
                             ifelse(wellNum %in% c(7:12), "B",
                                    ifelse(wellNum %in% c(13:18), "C", "D"))),
                column = ifelse(wellNum %in% c(1,7,13,19), 1,
                                ifelse(wellNum %in% c(2,8,14,20), 2,
                                       ifelse(wellNum %in% c(3,9,15,21), 3,
                                              ifelse(wellNum %in% c(4,10,16,22), 4,
                                                     ifelse(wellNum %in% c(5,11,17,23), 5, 6))))),
                column = factor(column, levels = c("6", "5", "4", "3", "2", "1"))) %>%
  dplyr::group_by(censor) %>%
  dplyr::mutate(number_censored_or_not = n())

# plot design with conservative censoring
ww_design_plot <- ggplot(ls_df_2_design) +
  aes(x = column, y = row) +
  geom_tile(aes(fill = factor(censor, levels = c("TRUE","FALSE")))) +
  geom_text(aes(label=strain), size = 1.15) +
  facet_wrap(~plateName) +
  theme_bw() +
  labs(y = "", x = "", fill = "censored", title = "Lifespan GWA (127 uncensored wells)") +
  theme(axis.text.x = element_text(size = 10, face = "plain", color = "black"),
        axis.text.y = element_text(size = 10, face = "plain", color = "black"),
        axis.title.x = element_text(size = 10, face = "bold", color = "black"),
        axis.title.y = element_text(size = 10, face = "bold", color = "black"),
        plot.title = element_text(size = 12, face = "bold", color = "black"),
        strip.text = element_text (size = 10, face = "plain", color = "black"),
        legend.text = element_text(size = 10, face = "plain", color = "black"))

ggsave(ww_design_plot, filename = "plots/20191111_GWA_design_conservative_censoring.png", width = 6.84, height = 5.12, dpi = 300)

# Remove 105 wells did not contain worms. 87 are NA for strain and 18 are labelled as few or none
ls_df_3 <- ls_df_2 %>%
  dplyr::filter(!(notes %in% c("few", "none")) & !is.na(strain))

# Remove 6 wells with total activity less than 0
ls_df_4 <- ls_df_3 %>%
  dplyr::filter(TotAct > 0)

# Remove 20 wells from plate 2 (plate had general FUDR failure, 4 wells already filtered due to negative total activity),
# Remove 7 wells with contamination (growth in wells according to Anthony)
filtering_list <- c("WI_1_11", "WI_1_15", "WI_5_9", "WI_6_9", "WI_7_8", "WI_8_22", "WI_15_8")

ls_df_5 <- ls_df_4 %>%
  dplyr::filter(plateName != "WI_2") %>%
  dplyr::mutate(filter = as.character(glue::glue('{plateName}_{wellNum}'))) %>%
  dplyr::filter(!(filter %in% filtering_list))

# Observe distribution of T97 - T89 for all wells.  those types of wells are almost always empty.
# For example, Censored well WI_1 #16 has LS = HS = 30 days (the limit of the experiment),
# which clearly means there was a persistent/constant noise floor and no actual worms that lived and then died.
ls_df_6 <- ls_df_5 %>%
  dplyr::mutate(t97_t89_diff = T97 - T89)

# remove 9 wells with t97-T89 values that are zero as these can't be informative
ls_df_7 <- ls_df_6 %>%
  dplyr::filter(t97_t89_diff > 0)

# Plot density distribution for T97 - T89 for remaining wells
initial_filtering_t97_t89_density <- ggplot(ls_df_7) +
  aes(t97_t89_diff) + 
  stat_ecdf(geom = "step") +
  geom_vline(xintercept = 2.5, lty = 2, color = "red") +
  theme_bw() +
  xlim(0,30) +
  labs(y = "Cumulative Density", x = "T97 - T89", title = "333 non-filtered wells with worms") +
  theme(axis.text.x = element_text(size = 10, face = "plain", color = "black"),
        axis.text.y = element_text(size = 10, face = "plain", color = "black"),
        axis.title.x = element_text(size = 10, face = "bold", color = "black"),
        axis.title.y = element_text(size = 10, face = "bold", color = "black"),
        strip.text.x = element_text(size = 10, face = "plain", color = "black"),
        strip.text.y = element_text(size = 10, face = "plain", color = "black"),
        plot.title = element_text(size = 12, face = "bold", color = "black"))
  # Getting density values from plot
  initial_filtering_density = ggplot_build(initial_filtering_t97_t89_density)$data[[1]]
  
  
# Plot density distribution for T97 - T89 for wells were strains were not loaed into wells
  ls_df_no_worms_loaded <- ls_df_2 %>%
    dplyr::filter(is.na(strain)) %>%
    dplyr::mutate(t97_t89_diff = T97 - T89) %>%
    dplyr::filter(t97_t89_diff != "NaN")
  
  no_worm_t97_t89_density <- ggplot(ls_df_no_worms_loaded) +
    aes(t97_t89_diff) + 
    stat_ecdf(geom = "step", color = "red") +
    geom_vline(xintercept = 2.5, lty = 2, color = "red") +
    #geom_density() +
    theme_bw() +
    xlim(0,30) +
    labs(y = "Cumulative Density", x = "T97 - T89", title = "87 wells with no worms added") +
    theme(axis.text.x = element_text(size = 10, face = "plain", color = "black"),
           axis.text.y = element_text(size = 10, face = "plain", color = "black"),
           axis.title.x = element_text(size = 10, face = "bold", color = "black"),
           axis.title.y = element_text(size = 10, face = "bold", color = "black"),
           strip.text.x = element_text(size = 10, face = "plain", color = "black"),
           strip.text.y = element_text(size = 10, face = "plain", color = "black"),
           plot.title = element_text(size = 12, face = "bold", color = "black"))

# Getting density values from plot
 no_worm_density = ggplot_build(no_worm_t97_t89_density)$data[[1]]
 
# put above together on one plot
bind_plot <- cowplot::plot_grid(no_worm_t97_t89_density, initial_filtering_t97_t89_density, ncol = 2)
 
ggsave(bind_plot, filename = "plots/20191111_GWA_too_few_worms_threshold.png", width = 6.84, height = 2.56, dpi = 300) 

# noise probabilty table
noise_df <- no_worm_density %>%
  dplyr::filter(is.finite(x)) %>%
  dplyr::mutate(y = 1-y) %>%
  dplyr::rename(`noise probability` = y,
                `T97 - T89 (d)` = x) %>%
  dplyr::select(`T97 - T89 (d)`, `noise probability`) %>%
  dplyr::arrange(`noise probability`)

noise_prob_p <- ggplot(noise_df) +
  aes(x = `T97 - T89 (d)`, y = `noise probability`) +
  geom_line() +
  geom_vline(xintercept = 2.5, lty = 2, color = "red") +
  xlim(0,30) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, face = "plain", color = "black"),
        axis.text.y = element_text(size = 10, face = "plain", color = "black"),
        axis.title.x = element_text(size = 10, face = "bold", color = "black"),
        axis.title.y = element_text(size = 10, face = "bold", color = "black"),
        strip.text.x = element_text(size = 10, face = "plain", color = "black"),
        strip.text.y = element_text(size = 10, face = "plain", color = "black"),
        plot.title = element_text(size = 12, face = "bold", color = "black"))

ggsave(noise_prob_p, filename = "plots/20191111_GWA_noise_prob.png", width = 3.42, height = 2.56, dpi = 300)

# apply 90% threshold censoring (T97-T89 > 3.1) on 333 non-filtered wells with worms
ls_df_90 <- ls_df_7 %>%
  dplyr::mutate(uniq_id = as.character(glue::glue('{plateName}_{wellNum}'))) %>%
  dplyr::filter(t97_t89_diff > 3.1)
 
ls_df_85 <- ls_df_7 %>%
  dplyr::mutate(uniq_id = as.character(glue::glue('{plateName}_{wellNum}'))) %>%
  dplyr::filter(t97_t89_diff > 2.54)

ls_df_80 <- ls_df_7 %>%
  dplyr::mutate(uniq_id = as.character(glue::glue('{plateName}_{wellNum}'))) %>%
  dplyr::filter(t97_t89_diff > 2)

ls_df_75 <- ls_df_7 %>%
  dplyr::mutate(uniq_id = as.character(glue::glue('{plateName}_{wellNum}'))) %>%
  dplyr::filter(t97_t89_diff > 1.33)

# plot box plots for various censoring strategies. ls_df_90 is most conservation
                                                #  ls_df_2 censor filter next
                                                #  ls_df_85
                                                #  ls_df_80
                                                #  ls_df_75
ls_df_90_box <- 

# generate trait files for mappings based on sliding censoring threshold




###############
### tsting ###
###############
# look at distribution of total activity for remaining wells (369)
frequency_plot <- ggplot(ls_df_4) +
  aes(TotAct, stat(density), color = as.character(censor)) +
  geom_freqpoly(bins = 20) +
  xlim(0,1.6) +
  theme_bw() +
  labs(y = "Density", x = "Total Activity", color = "censored", title = "Lifespan GWA censoring (n = 369 wells)") +
  theme(axis.text.x = element_text(size = 14, face = "plain", color = "black"),
        axis.text.y = element_text(size = 14, face = "plain", color = "black"),
        axis.title.x = element_text(size = 16, face = "bold", color = "black"),
        axis.title.y = element_text(size = 16, face = "bold", color = "black"),
        strip.text.x = element_text(size = 16, face = "bold", color = "black"),
        strip.text.y = element_text(size = 16, face = "bold", color = "black"),
        plot.title = element_text(size = 16, face = "bold", color = "black"))

ggsave(frequency_plot, filename = "plots/20191111_GWA_Total_Activity_histogram.png", width = 6.84, height = 5.12, dpi = 300)

# generate two data sets. 1 uses Anthony censoring meathod. 2 uses criteria a(if note says, few, none), b(if strain is NA), c(if strain is not NA, but total activity is ) 


# testing is na = 87
# testing is few, none = 18
test <- ls_df_2 %>% dplyr::filter(notes %in% c("few", "none") & !is.na(strain))
length(test$strain
       )
