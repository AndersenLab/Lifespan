#!/usr/bin/env Rscript
#Load necessary packages
library(tidyverse)
library(sommer)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))


########################################
### Heritability with sommer package ###
########################################
# Read traitfile
traitfile_small <- data.table::fread("data/life_machine_trait_file.tsv")

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
sommer::pin(df_H2, H2 ~ (V1) / (V1+V2+V3))

# narrow-sense H2 (additive + epistatic variance) / (additive, epistatic, error)
sommer::pin(df_H2, H2 ~ (V1+V2) / (V1+V2+V3))
