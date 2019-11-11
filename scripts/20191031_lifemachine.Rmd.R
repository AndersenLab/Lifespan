---
  title: "20190625_lifespan"
author: "Tim C."
date: "6/26/2019"
output:
  html_document:
  toc: true
toc_float:
  collapsed: false
smooth_scroll: false
---
  
  ```{r, echo=F, warning = F, message=F}
library(tidyverse)
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

# Broad-sense H2 (additive only). Not ewe never worry about dominance because we assume homozygous.
pin(df_H2, H2 ~ (V1) / (V1+V2+V3))

# narrow-sense H2 (additive + epistatic variance) / (additive, epistatic, error)
pin(df_H2, H2 ~ (V1+V2) / (V1+V2+V3))
