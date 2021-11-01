library(tidyverse)
library(gggenes)
library(genemodel)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# Filter to strains in mapping
strains <- data.table::fread('data/strain_list_CFY.tsv')

# grab the annotations for taf-7.2 from 20210121 CeNDR release
anno <- data.table::fread("/Users/tim/repos/Lifespan/data/WI.20210121_taf7.2.tsv") %>%
  dplyr::select(chrom = V1, start = V2, ref = V3, alt = V4, consq = V5, wbid = V6, seq_id = V7, bio_class = V8,
                amino_acid = V10, base = V11, variant_strains = V12, BLOSUM = V13, grantham = V14, perc_protein = V15,
                gene_name = V16, impact = V17, divergent_flag = V18)

# grab trait file
trait <- data.table::fread("/Users/tim/repos/Lifespan/data/traitfile_ls_GWA_CFY.tsv")

# expand collapsed strain data
anno2 <- anno %>%
  tidyr::separate_rows(variant_strains)

# find strains in mapping in the anno file only 7 variants total? hmm, there were only two variants in the GWA finemapping plot. 
anno3 <- anno2 %>%
dplyr::filter(variant_strains %in% strains$strain) %>%
dplyr::mutate(id = paste0(start, ref, alt),
              n_variants = length(unique(id)))

#=============================================#
# plot gene and variation
#=============================================#
ggplot(example_genes, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
  geom_gene_arrow() +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3")
view(example_genes)

spl1<-data.frame(
  type=c("5' utr", "coding_region", "intron", "coding_region", "intron", "coding_region","3' utr"), 
  coordinates=c("1-50", "50-100", "100-150", "150-200", "200-250", "250-300","300-350"))

genemodel.plot(model=spl1, start=1, bpstop=350, orientation="reverse", xaxis=T)
mutation.plot(25150214, 25150214, text="P->S", col="black", drop=0, haplotypes=c("red", "blue"))
mutation.plot(75, 75, text="P->S", col="black", drop=0, haplotypes=c("red", "blue"))

# pull taf-7.2 gff model from WS276
gff76_taf <- data.table::fread("data/WS276_taf7.2.gff3")
gff82_taf <- data.table::fread("data/WS282_taf7.2.gff3")

gene_df_276 <- gff76_taf %>%
  dplyr::rename(type = V3, start = V4, stop = V5) %>%
  dplyr::mutate(type2 = case_when(type == "mRNA") )
  dplyr::mutate(coordinates = paste0(start, "-", stop))
  type = V3, 
#>Y111B2A.42
#CCACAATCATTCTTCGACGACACACCAGTAGCATCTTCCGACGATCCACCAGACTTCGAAAGTCACATTGTACTACGTGTACCTGAAGATTGTGTGGGTAGAATCGAGAAAATCATTCAATCGGACGGAAAACACGAGGAATTCTCGTTAAATTTGAATTCAGACGCTCGAAATTCCACAGTCAGAATCGGAAATCAACTGTTAAATGGAAAAATCCTGGATCTTCCCACTATAACAGAAATCCACAAGACATTAGACAACAAAAGCCTGTATAAAGTCGCAGATGTCTCACAGATCCTTGTCTGCACCCATGATTCCATCAATTCAATAGCTTCAAGCTCTGAAGATGCTGCTCAGAAAGCAGCAGCGGCAAAGGCAAAACAATGGCAATACCCGCACGGACTGACGCCTCCCATGAAATCGGCGAGGAAGAAGCGATTTCGAAAGACTAAGAAGAAAAAGTTTATGGATGCTCCAGAGGTTGAGAAGGAGCTCAAGAGACTGCTCCGTGCAGATTTGGAAGCGGATAGTGTGAAATGGGAAATTGTGGAAGGGAATAAGGAAGGTGCGACGGATGAAG

#==============================================#
# Look at mapping
#==============================================#
test <- data.table::fread("data/NemaScan/Analysis_Results-20211101/Mapping/Processed/processed_t97_mean_AGGREGATE_mapping.tsv")
