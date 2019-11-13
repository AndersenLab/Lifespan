#!/usr/bin/env Rscript
#Load necessary packages
library(tidyverse)
library(rio)
library(viridis)
library(repmis)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# Requires library(repmis) (Need data for chromosome III)!
load('data/haplotype_plot_df.Rda')

# Filter chromosomes appropriately
plot_df_filtered <- plot_df %>%
  dplyr::filter(chromosome == "III") %>%
  dplyr::arrange(chromosome, start, stop) %>%
  # Filter out strains with tiny haplotypes
  dplyr::mutate(filtered_swept_haplotype=
                  (
                    (hap_length > 1E5)
                    &
                      (max_haplotype_shared > 0.03)
                    &
                      (swept_haplotype == TRUE)
                  )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(chromosome, isotype) %>%
  dplyr::mutate(is_swept = (sum(filtered_swept_haplotype) > 0))

# Filter to strains in mapping
strains <- data.table::fread('data/strain_list_CFY.tsv')

plot_df_filtered_mapping <- plot_df_filtered %>%
  dplyr::filter(isotype %in% strains$strain)

###########################################
### build plot that looks like daehan's ###
###########################################
# # X swept (>50%) isotypes
# df_Xhap_swept <- df_Xhap %>%
#   dplyr::filter(filtered_swept_haplotype == "TRUE", max_haplotype_shared >= 0.5) 
# Xswept_strains <- unique(df_Xhap_swept$isotype)
# ## srg-37(del) outcrossed isotypes
# out_strains <- unique(Xswept_strains[!Xswept_strains %in% unique(df_Xhap_NIC2$strain)])  ## over 50% swept, but no srg-37(del)
# 
# df_outcross_haps <- plot_df_filtered %>%
#   dplyr::filter(isotype %in% out_strains)

# outcross_haps <- unique(df_outcross_haps$haplotype)

ranked_by_sharing_oc <- plot_df_filtered_mapping %>%
  dplyr::filter(chromosome== "III") %>%
  dplyr::filter(start<12.8e6 & stop > 12.6e6) %>%
  dplyr::group_by(isotype, haplotype) %>%
  dplyr::mutate(locus_hap_length_1 = stop-start) %>%
  dplyr::mutate(locus_hap_length_2 = stop-12.6e6) %>%
  dplyr::mutate(locus_hap_length_3 = 12.8e6-start) %>%
  dplyr::select(start, stop, isotype, haplotype, locus_hap_length_1, locus_hap_length_2, locus_hap_length_3) %>%
  tidyr::gather(-start, -stop, -isotype, -haplotype, key = "locus_hap_class", value = "length") %>%
  dplyr::group_by(isotype, haplotype) %>%
  dplyr::summarise(locus_haplotype_length = min(length)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(isotype) %>%
  dplyr::mutate(locus_hap_max = ifelse(locus_haplotype_length > 2e5, 2e5, locus_haplotype_length)) %>%
  dplyr::mutate(locus_max = max(locus_hap_max)) %>%
  dplyr::filter(locus_max == locus_hap_max) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(locus_max), haplotype) %>%
  dplyr::mutate(plotpoint=row_number()) %>%
  dplyr::select(isotype, plotpoint) %>%
  dplyr::arrange(plotpoint)

strain_labels_oc <- ranked_by_sharing_oc$isotype

# extract which strains are ref and alt at peak position
refalt_peak <- data.table::fread('data/20191112_lifespan_GWA_cegwas2-nf_output/11132019_lifespanEigen_CFY/Fine_Mappings/Data/t97_mean_snpeff_genes.tsv') %>%
  dplyr::filter(CHROM == "III", POS == 12672440) %>%
  dplyr::select(STRAIN, CHROM, POS, REF, ALT, STRAIN_GENOTYPE) %>%
  dplyr::mutate(genotype_at_peak = ifelse(STRAIN_GENOTYPE == REF, "REF", "ALT"))

ref_strains <- refalt_peak %>%
  dplyr::filter(genotype_at_peak == "REF") %>%
  dplyr::select(strain = STRAIN)

alt_strains <- refalt_peak %>%
  dplyr::filter(genotype_at_peak == "ALT") %>%
  dplyr::select(strain = STRAIN)

plot_oc_hap <- plot_df_filtered_mapping %>%
  dplyr::select(-plotpoint) %>%
  dplyr::filter(chromosome == "III") %>%
  dplyr::left_join(ranked_by_sharing_oc) %>%
  ggplot(.,
         aes(xmin = start/1E6, xmax = stop/1E6,
             ymin = plotpoint - 0.5, ymax = plotpoint + 0.5,
             fill = haplotype)) +
  geom_rect() +
  ggplot2::geom_vline(ggplot2::aes(xintercept = 12575736/1e+06), color = "blue", alpha = 0.7, size = 0.5) + 
  ggplot2::geom_vline(ggplot2::aes(xintercept = 12746972/1e+06), color = "blue", alpha = 0.7, size = 0.5) +
  ggplot2::geom_vline(ggplot2::aes(xintercept = 12672440/1e+06), color = "red", alpha = 0.7, size = 1) +
  #scale_fill_manual(values = mcolor) +
  scale_y_continuous(breaks = 1:length(strain_labels_oc),
                     labels = strain_labels_oc,
                     expand = c(0, 0), position = 'left') +
  xlab("Position (Mb)") +
  labs(title = "T97_mean haplotype structure") +
  coord_cartesian(xlim=c(12, 13.5)) +
  theme_bw() +
  facet_grid(.~chromosome, scales="free", space="free") +
  theme(axis.text.x = element_text(size=12, color = "black"),
        strip.text = element_text(size=12, color = "black"),
        legend.position = 'none',
        axis.title.x = element_text(size=12, color = "black"),
        axis.text.y = element_text(size = 5, color = ifelse(strain_labels_oc %in% ref_strains$strain, "black", "red")),
        plot.margin = margin(b=0.1, l=0.2, r=0.1, t=0.3, unit = "in"))

ggsave(plot_oc_hap, filename = "plots/20191111_GWA_T97_haplotype_structure_CFY.png", width = 6.84, height = 5.12, dpi = 300)
