#!/usr/bin/env Rscript
#Load necessary packages
library(tidyverse)
library(cegwas2)
library(rio)
library(viridis)

devtools::install_github("AndersenLab/cegwas2")

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

##########################################
# 1: lload data from burden mapping (VT) #
##########################################
meanT97_VT_dat <- data.table::fread('data/20191112_lifespan_GWA_cegwas2-nf_output/11132019_lifespanEigen_CFY/BURDEN/VT/Data/t97_mean.VariableThresholdPrice.assoc') %>%
  dplyr::mutate(significant = ifelse(PermPvalue < .05/n(), TRUE,FALSE ),
                n = n()) # adding bonferroni significance test to data. Not sure about how mtDNA should be treated here. The plot_burden.R script calculates n() before filtereing mtDNA.

meanT89_VT_dat <- data.table::fread('data/20191112_lifespan_GWA_cegwas2-nf_output/11132019_lifespanEigen_CFY/BURDEN/VT/Data/t89_mean.VariableThresholdPrice.assoc') %>%
  dplyr::mutate(significant = ifelse(PermPvalue < .05/n(), TRUE,FALSE ),
                n = n()) # adding bonferroni significance test to data. Not sure about how mtDNA should be treated here. The plot_burden.R script calculates n() before filtereing mtDNA.

meanT97_T89_VT_dat <- data.table::fread('data/20191112_lifespan_GWA_cegwas2-nf_output/11132019_lifespanEigen_CFY/BURDEN/VT/Data/t97_t89_diff_mean.VariableThresholdPrice.assoc') %>%
  dplyr::mutate(significant = ifelse(PermPvalue < .05/n(), TRUE,FALSE ),
                n = n()) # adding bonferroni significance test to data. Not sure about how mtDNA should be treated here. The plot_burden.R script calculates n() before filtereing mtDNA.

#############################
# 2: load gene descriptions #
#############################
load('/Users/tim/repos/Lifespan/data/gene_descriptions_WS273.Rda')

# join gene descriptions with significant genes
T97_desc <- left_join(meanT97_VT_dat, gene_descriptions %>% dplyr::rename(Gene = wbgene), by = "Gene") %>%
  dplyr::filter(significant == T) %>%
  dplyr::arrange(RANGE)

T89_desc <- left_join(meanT89_VT_dat, gene_descriptions %>% dplyr::rename(Gene = wbgene), by = "Gene") %>%
  dplyr::filter(significant == T) %>%
  dplyr::arrange(RANGE)

T97_T89_desc <- left_join(meanT97_T89_VT_dat, gene_descriptions %>% dplyr::rename(Gene = wbgene), by = "Gene") %>%
  dplyr::filter(significant == T) %>%
  dplyr::arrange(RANGE)

# get list of strains in mapping for pasting into CeNDR and for using in query_vcf
strains<-  tibble(names(data.table::fread('/Users/tim/repos/Lifespan/data/20191112_lifespan_GWA_cegwas2-nf_output/11132019_lifespanEigen_CFY/Genotype_Matrix/Genotype_Matrix.tsv') %>%
                          dplyr::select(-CHROM, -POS, -REF, -ALT))) %>%
  dplyr::rename(names = `names(...)`)

strains_ce <- toString(strains$names)

##############################################
# 3: use query_vcf function to pull variants #
##############################################
# use query_vcf function in cegwas2 to find which strains have variants in these genes. use arguments 
gene_list_T89 <- T89_desc$RANGE

# remove the C45G7.11 b/c it is a ncRNA and causing an error in query_vcf function
gene_list_T89 <- T89_desc %>%
  dplyr::filter(!(RANGE %in% c("IV:2459894-2459964", "X:4244807-4245009"))) %>%
  dplyr::pull(RANGE)

#run query_vcf function
T89_VT_hit_variants <- cegwas2:::query_vcf(gene_list_T89, format = c("GT"), impact = c("LOW", "MODERATE", "HIGH"), samples = strains$names)

T89_VT_hit_variants_count <- T89_VT_hit_variants %>%
  dplyr::filter(genotype == 2) %>%
  dplyr::group_by(query, SAMPLE, impact) %>%
  dplyr::mutate(num_variants_in_sample_gene = n()) %>%
  dplyr::group_by(query, SAMPLE) %>%
  dplyr::mutate(total_num_variants_in_sample_gene = n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(SAMPLE, query, gene_name, impact, num_variants_in_sample_gene, total_num_variants_in_sample_gene) %>%
  dplyr::distinct(SAMPLE, query, impact, .keep_all = T)