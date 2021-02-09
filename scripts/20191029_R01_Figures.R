#!/usr/bin/env Rscript
#Load necessary packages
library(tidyverse)
library(geonames)
library(ggrepel)
library(maps)
library(ggthemes)
library(geosphere)
library(scatterpie)
library(legendMap)
library(rio)
library(viridis)
library(boot)
library(sommer)
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

#######################################################
### Figure 1A: World Map of 12 divergent strain set ###
#######################################################
# set color palette
strain.colours <- c("DL238"="gold2",
                    "N2"="plum4",
                    "JU775"= "darkorange1", 
                    "MY16"="lightskyblue2",
                    "CX11314"="firebrick",
                    "ED3017"= "burlywood3",
                    "LKC34"="gray51", 
                    "JU258"="springgreen4",
                    "JT11398"="lightpink2",
                    "CB4856"="deepskyblue4",
                    "MY23"="mediumpurple4",
                    "EG4725"="chocolate")

# Pull list of stain names for divergent set
strains <- names(strain.colours)

# load in wi strain info for lat long coordinates for map
wi_df <- data.table::fread("data/WI_strain_list.csv")

isolation_info <- wi_df %>%
  dplyr::filter(reference_strain == 1) %>%
  dplyr::select(strain, isotype, long = longitude, lat = latitude, state, country) %>%
  dplyr::filter(lat != "None") %>%
  dplyr::filter(!is.na(lat)) %>%
  dplyr::filter(isotype %in% strains)

isolation_info$lat <- as.numeric(isolation_info$lat)
isolation_info$long <- as.numeric(isolation_info$long)

# # Plot world map collection locations
dev.off()
world <- map_data("world")
world <- world[world$region != "Antarctica",] # intercourse antarctica

world_plot <- ggplot()+ geom_map(data=world, map=world,
                                      aes(x=long, y=lat, map_id=region),
                                      color="black", fill="#E6E6E6", size=0.25)+
  scale_fill_manual(values = strain.colours,name = "Population") +
  theme_map() +
  geom_point(aes(long, lat, fill = isotype, alpha = 0.75),
             data=isolation_info,
             shape = 21, size = 3) +
  geom_label_repel(aes(long, lat, label = isotype, fill = isotype),
                   data=isolation_info,
                   fontface = 'bold', color = 'white',
                   segment.color = 'black', size = 4) +
  theme(legend.position = "none",
        plot.margin = unit(c(0,0.1,0,0), "cm")) +
coord_quickmap(xlim = c(-170,80), ylim = c(80,-75))# coord_quickmap keeps correct aspect ratio
world_plot
###################################################
### Figure 1 - supplemental data Heritabilities ###
###################################################
# load replicated wormwatcher data for divergent set
ww_df <- data.table::fread('data/2019-05-24_AllDataLongFormat.csv', header = T)

# data shaping
ww_df_proc <- ww_df %>%
  dplyr::rename( strain = Strain,
                 Plate_Num = `Plate Num`,
                 ls = `Lifespan est (days)`,
                 hs = `Healthspan est (days)`,
                 t = `Time (days)`,
                 act = `Activity (a.u.)`) %>%
  dplyr::distinct(Plate_Num, Plate_Name, strain, Well_Num, .keep_all = T) %>%
  dplyr::filter(ls != "NaN") %>%
  #dplyr::mutate(indicies = seq(1, length(t), by = 1)) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(reps = n(),
                median_ls = median(ls),
                median_hs = median(hs),
                mean_ls = mean(ls),
                mean_hs = mean(hs),
                stdev_ls = sd(ls),
                sterr_ls = (sd(ls)/sqrt(n())),
                stdev_hs = sd(hs),
                sterr_hs = (sd(hs)/sqrt(n()))) %>%
  dplyr::ungroup() %>%
  tidyr::gather(trait, phenotype, -Plate_Num, -Plate_Name, -Well_Num, -strain, -t, -act, -reps, -median_ls, -median_hs, -mean_ls, -mean_hs, -stdev_ls, -sterr_ls, -stdev_hs, -sterr_hs) %>%
  dplyr::arrange(median_hs)

# remove outliers in new df
ww_df_proc_outliers_removed <- ww_df_proc %>%
  dplyr::group_by(strain, trait) %>%
  dplyr::mutate(phenotype = remove_outliers(phenotype)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(phenotype))

# flag outliers in old data frame
ww_df_proc <- ww_df_proc %>%
  dplyr::group_by(strain, trait) %>%
  dplyr::mutate(outlier = (remove_outliers(phenotype)),
                outlier = ifelse(is.na(outlier), TRUE, FALSE)) %>%
  dplyr::ungroup()

# shaping data for H2 with bootstrapping (broad-sense)
ww_df_H2_ls <- ww_df_proc %>%
  dplyr::filter(trait == "ls") 
ww_df_H2_hs <- ww_df_proc %>%
  dplyr::filter(trait == "hs") 

ww_df_H2_ls_outliers_removed <- ww_df_proc_outliers_removed %>%
  dplyr::filter(trait == "ls") 
ww_df_H2_hs_outliers_removed <- ww_df_proc_outliers_removed %>%
  dplyr::filter(trait == "hs") 


# calculate H2 for both traits
ww_ls_h2 <- H2.calc(ww_df_H2_ls, boot = T)
ww_hs_h2 <- H2.calc(ww_df_H2_hs, boot = T)

ww_ls_h2_outliers_removed <- H2.calc(ww_df_H2_ls_outliers_removed, boot = T)
ww_hs_h2_outliers_removed <- H2.calc(ww_df_H2_hs_outliers_removed, boot = T)

# join H2 plots
ww_ls_h2_full <- bind_rows(ww_ls_h2, ww_ls_h2_outliers_removed) %>%
  dplyr::mutate(trait = as.factor(c("WW_LS","WW_LS no outliers")))

ww_hs_h2_full <- bind_rows(ww_hs_h2, ww_hs_h2_outliers_removed) %>%
  dplyr::mutate(trait = as.factor(c("WW_HS","WW_HS no outliers")))

# ploting data
ww_ls_h2_plot <- ggplot(ww_ls_h2_full) +
  geom_point(aes(x = trait, y = H2), size = 3) +
  geom_errorbar(aes(ymin=ci_l, ymax=ci_r, x = trait), colour="black", width=0.05) +
  labs(x = "", y = "Heritability") +
  ylim(0, 1) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ww_hs_h2_plot <- ggplot(ww_hs_h2_full) +
  geom_point(aes(x = trait, y = H2), size = 3) +
  geom_errorbar(aes(ymin=ci_l, ymax=ci_r, x = trait), colour="black", width=0.05) +
  labs(x = "", y = "Heritability") +
  ylim(0, 1) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# design plot data shaping
ww_design <- ww_df_proc %>%
  dplyr::mutate(row = ifelse(Well_Num %in% c(1:6), "A",
                             ifelse(Well_Num %in% c(7:12), "B",
                                    ifelse(Well_Num %in% c(13:18), "C", "D"))),
                column = ifelse(Well_Num %in% c(1,7,13,19), 1,
                                ifelse(Well_Num %in% c(2,8,14,20), 2,
                                       ifelse(Well_Num %in% c(3,9,15,21), 3,
                                              ifelse(Well_Num %in% c(4,10,16,22), 4,
                                                     ifelse(Well_Num %in% c(5,11,17,23), 5, 6))))))

# plot design
ww_design_plot <- ggplot(ww_design %>% dplyr::filter(trait == "ls")) +
  aes(x = column, y = row) +
  geom_tile(aes(fill = as.numeric(phenotype))) +
  geom_text(aes(label=strain), size = 2) +
  facet_wrap(~Plate_Name) +
  theme_bw() +
  scale_fill_viridis(name="") +
  labs(x="", y="")
ww_design_plot


#############################################################
### Figure 1B-C: Healthspan and lifespan median barcharts ###
#############################################################


# Figure 1B: Healthspan and lifespan barcharts
median_ls_ranks <- ww_df_proc %>%
  dplyr::distinct(strain, .keep_all = T) %>%
  dplyr::arrange(median_ls) %>%
  dplyr::mutate(ls_rank = seq(1:12)) 

median_hs_ranks <- ww_df_proc %>%
  dplyr::distinct(strain, .keep_all = T) %>%
  dplyr::arrange(median_hs) %>%
  dplyr::mutate(hs_rank = seq(1:12))

median_ls_order <- median_ls_ranks$strain

median_hs_order <- median_hs_ranks$strain

# plot median trait values for all strains
median_ww_ls <- ggplot(median_ls_ranks) +
  aes(x=factor(strain, levels = median_ls_order), y=median_ls, fill=strain) +
  scale_fill_manual(values = strain.colours) +
  geom_col(width = 1, color = "black", size = 0.2) +
  annotate(geom = "text", label = expression(H^2 == 0.34), x = 10, y = 20) + # ww_ls_h2_full[1,1]
  labs(x="", y="lifespan (d)") +
  theme_bw() +
  ylim(0, 21) +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black"),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.2,0.2,0,0.2), "cm"))

median_ww_hs <- ggplot(median_hs_ranks) +
  aes(x=factor(strain, levels = median_ls_order), y=median_hs, fill=strain) +
  scale_fill_manual(values = strain.colours) +
  geom_col(width = 1, color = "black", size = 0.2) +
  geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2,
                position=position_dodge(.9))
  annotate(geom = "text", label = expression(H^2 == 0.32), x = 10, y = 15) + # ww_hs_h2_full[1,1]
  ylim(0, 16) +
  labs(x="", y="healthspan (d)") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.text = element_text(colour = "black", size = 8),
        plot.margin = unit(c(0,0.2,0,0.2), "cm"))

median_bars <- cowplot::plot_grid(median_ww_ls, median_ww_hs,
                              rel_heights = c(.85, 1), nrow = 2, labels = c("B", "C"), hjust = 0, vjust = c(1.45, 0.85))


#############################################################
### Figure 1B-C: Healthspan and lifespan mean barcharts   ###
#############################################################


# Figure 1B: Healthspan and lifespan barcharts
mean_ls_ranks <- ww_df_proc %>%
  dplyr::distinct(strain, .keep_all = T) %>%
  dplyr::arrange(mean_ls) %>%
  dplyr::mutate(ls_rank = seq(1:12)) 

mean_hs_ranks <- ww_df_proc %>%
  dplyr::distinct(strain, .keep_all = T) %>%
  dplyr::arrange(mean_hs) %>%
  dplyr::mutate(hs_rank = seq(1:12))

mean_ls_order <- mean_ls_ranks$strain

mean_hs_order <- mean_hs_ranks$strain

# plot median trait values for all strains
mean_ww_ls <- ggplot(mean_ls_ranks) +
  aes(x=factor(strain, levels = mean_ls_order), y=mean_ls, fill=strain) +
  scale_fill_manual(values = strain.colours) +
  geom_col(width = 1, color = "black", size = 0.2) +
  geom_linerange(aes(ymax = mean_ls + stdev_ls, ymin = mean_ls - stdev_ls), size = 0.2) +
  #geom_errorbar(aes(ymax = mean_ls + stdev_ls, ymin = mean_ls), width = .2, position = position_dodge(.9), size = 0.2) +
  annotate(geom = "text", label = expression(H^2 == 0.34), x = 10, y = 20) + # ww_ls_h2_full[1,1]
  labs(x="", y="lifespan (d)") +
  theme_bw() +
  ylim(0, 21) +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black"),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.2,0.2,0,0.2), "cm"))

mean_ww_hs <- ggplot(mean_hs_ranks) +
  aes(x=factor(strain, levels = mean_ls_order), y=mean_hs, fill=strain) +
  scale_fill_manual(values = strain.colours) +
  geom_col(width = 1, color = "black", size = 0.2) +
  geom_linerange(aes(ymax = mean_hs + stdev_hs, ymin = mean_hs - stdev_hs), size = 0.2) +
annotate(geom = "text", label = expression(H^2 == 0.32), x = 10, y = 15) + # ww_hs_h2_full[1,1]
  ylim(0, 16) +
  labs(x="", y="healthspan (d)") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.text = element_text(colour = "black", size = 8),
        plot.margin = unit(c(0,0.2,0,0.2), "cm"))

mean_bars <- cowplot::plot_grid(mean_ww_ls, mean_ww_hs,
                                  rel_heights = c(.85, 1), nrow = 2, labels = c("B", "C"), hjust = 0, vjust = c(1.45, 0.85))
### With standard errors ######
# plot median trait values for all strains
mean_ww_ls <- ggplot(mean_ls_ranks) +
  aes(x=factor(strain, levels = mean_ls_order), y=mean_ls, fill=strain) +
  scale_fill_manual(values = strain.colours) +
  geom_col(width = 1, color = "black", size = 0.2) +
  geom_linerange(aes(ymax = mean_ls + sterr_ls, ymin = mean_ls - sterr_ls), size = 0.2) +
  #geom_errorbar(aes(ymax = mean_ls + stdev_ls, ymin = mean_ls), width = .2, position = position_dodge(.9), size = 0.2) +
  annotate(geom = "text", label = expression(H^2 == 0.34), x = 10, y = 20) + # ww_ls_h2_full[1,1]
  labs(x="", y="lifespan (d)") +
  theme_bw() +
  ylim(0, 21) +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black"),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.2,0.2,0,0.2), "cm"))

mean_ww_hs <- ggplot(mean_hs_ranks) +
  aes(x=factor(strain, levels = mean_ls_order), y=mean_hs, fill=strain) +
  scale_fill_manual(values = strain.colours) +
  geom_col(width = 1, color = "black", size = 0.2) +
  geom_linerange(aes(ymax = mean_hs + sterr_hs, ymin = mean_hs - sterr_hs), size = 0.2) +
  annotate(geom = "text", label = expression(H^2 == 0.32), x = 10, y = 15) + # ww_hs_h2_full[1,1]
  ylim(0, 16) +
  labs(x="", y="healthspan (d)") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.text = element_text(colour = "black", size = 8),
        plot.margin = unit(c(0,0.2,0,0.2), "cm"))

mean_bars_sterr <- cowplot::plot_grid(mean_ww_ls, mean_ww_hs,
                                rel_heights = c(.85, 1), nrow = 2, labels = c("B", "C"), hjust = 0, vjust = c(1.45, 0.85))


# Attempt to put healthspan and lifespan together
# create a dataset
ww_df_proc_distinct <- ww_df_proc %>%
  dplyr::distinct(strain, .keep_all = TRUE) %>%
  dplyr::select(strain, Healthspan = median_hs, Lifespan = median_ls) %>%
  tidyr::gather(trait, value, - strain)

# colors
bar.colors <- c(Healthspan = "#63AABC",
                Lifespan = "#ED3833")
# Grouped
grouped_traits <- ggplot(ww_df_proc_distinct, aes(fill=trait, y=value, x=factor(strain, levels = median_ls_order))) + 
  geom_bar(position="dodge", stat="identity", color = "black", size = 0.3) +
  scale_fill_manual(values = bar.colors) +
  annotate(geom = "text", label = expression('Lifespan '*H^2 == 0.34*''), x = 9, y = 21) + # ww_ls_h2_full[1,1]
  annotate(geom = "text", label = expression('Healthspan '*H^2 == 0.32*''), x = 9.5, y = 18.5) +
  labs(x="", y="Days", fill = "") +
  theme_bw() +
  #ylim(0, 21) +
  theme(legend.position = "bottom",
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.margin = unit(c(0.2,0.2,0,0), "cm"))

# Add all of it together
R01_Figure_v3 <- cowplot::plot_grid(world_plot, mean_bars, ncol = 2, labels = c("A",""), rel_widths = c(1,.4))
#ggsave(R01_Figure_v3, filename = "plots/R01_Figure_v3.pdf", width = 7.5, height = 4, useDingbats = FALSE)
#ggsave(R01_Figure_v3, filename = "plots/R01_Figure_v3.png", width = 7.5, height = 4, dpi = 300)

R01_Figure_v4_sterr <- cowplot::plot_grid(world_plot, mean_bars_sterr, ncol = 2, labels = c("A",""), rel_widths = c(1,.4))
#ggsave(R01_Figure_v4_sterr, filename = "plots/R01_Figure_v4.pdf", width = 7.5, height = 4, useDingbats = FALSE)
#ggsave(R01_Figure_v4_sterr, filename = "plots/R01_Figure_v4.png", width = 7.5, height = 4, dpi = 300)

##############################################################################################
### R01 data analysis: Broad-sense and narrow-sense heritability for Patrick Phillips data ###
##############################################################################################
# Read in data from lifespan machine patrick phillips lab
lm_df1 <- data.table::fread('data/CeNDR_lifespans.csv', header = T) %>%
  dplyr::mutate(strain = ifelse(strain == "EC19191", "ECA191",
                                ifelse(strain == "PL238", "DL238", strain)))
lm_df2 <- data.table::fread('data/CeNDR_071919.csv', header = T) %>%
  dplyr::mutate(`Censored Reason` = as.character(`Censored Reason`),
                `Event Observation Type` = as.character(`Event Observation Type`)) %>%
  dplyr::select(device = Device,
                exp = Experiment,
                plate = `Plate Name`,
                row = `Plate Row`,
                column = `Plate Column`,
                strain = `Strain`,
                temp = `Culturing Temperature`,
                food = `Food Source`,
                event_freq = `Event Frequency`,
                ls = `Age at Death (d)`,
                DNMF = `Duration Not Fast Moving (d)`,
                gap_days = `Longest Gap in Measurement (d)`,
                censor = `Censored`,
                censor_reason = `Censored Reason`,
                event_observation_type = `Event Observation Type`)

# Join lm dataframes and calculate lifespan and healthspan traits per replicate.No regression of device effects or anything
lm_df <- dplyr::full_join(lm_df1, lm_df2) %>%
  dplyr::filter(censor != 1) %>% # remove censored data
  dplyr::mutate(rep = ifelse(strain == "CB4856" & is.na(rep), 5,
                             ifelse(strain == "CX11314" & is.na(rep), 5, rep))) %>% # give CB4856 and CX11314 a fifth replicate. I'm not sure if this is appropriate becuase I don't know what constitutes a replicate
  dplyr::group_by(strain, rep) %>%
  dplyr::mutate(mean_ls = mean(ls),
                sd_ls = sd(ls),
                cv_ls = sd_ls/mean_ls,
                mean_DNMF = mean(DNMF),
                sd_DNMF = sd(DNMF),
                cv_DNMF = sd_DNMF/mean_DNMF) %>% # The DNMF has some impossible numbers, we need to clarify what those are
  dplyr::ungroup()
  

# Build trait file for cegwas2, starting with simple mean lifespan for each strain. This is a mean of means from the 4-5 replicates.
traitfile_small <- lm_df %>%
  dplyr::select(strain, rep, mean_ls, sd_ls, cv_ls) %>%
  dplyr::distinct(strain, rep, .keep_all = TRUE) %>%
  dplyr::select(strain, mean_ls) %>%
  dplyr::arrange(strain) %>%
  dplyr::group_by(strain) %>%
  dplyr::summarise(mean_ls = mean(mean_ls)) %>%
  dplyr::ungroup()

traitfile <- lm_df %>%
  dplyr::select(strain, rep, mean_ls, sd_ls, cv_ls) %>%
  dplyr::distinct(strain, rep, .keep_all = TRUE) %>%
  dplyr::arrange(strain, rep) %>%
  dplyr::select(-rep)

# Write traitfile for cegwas2-nf input
rio::export(traitfile_small, 'data/life_machine_trait_file.tsv')

########################################
### Heritability with sommer package ###
########################################
# load in genotype matrix from cegwas2-nf (WI.20180527.impute.vcf.gz used to generate matrix)
geno_matrix <- data.table::fread("data/life_machine_Genotype_Matrix.tsv")

# df with strain, strain, trait value
df_y <- traitfile_small %>%
  dplyr::rename(value = mean_ls) %>%
  dplyr::mutate(strain1=strain, strain2=strain) %>%
  dplyr::select(-strain)

A <- sommer::A.mat(t(dplyr::select(geno_matrix, -CHROM, -POS, -REF, -ALT)))
E <- sommer::E.mat(t(dplyr::select(geno_matrix, -CHROM, -POS, -REF, -ALT)))

df_H2 <- sommer::mmer(value~1, random=~vs(strain1,Gu=A)+vs(strain2,Gu=E), data=df_y)

(summary(df_H2)$varcomp)

# Broad-sense H2 (additive only). Note we never worry about dominance because we assume homozygous.
pin(df_H2, H2 ~ (V1) / (V1+V2+V3))

# narrow-sense H2 (additive + epistatic variance) / (additive, epistatic, error)
pin(df_H2, H2 ~ (V1+V2) / (V1+V2+V3))
