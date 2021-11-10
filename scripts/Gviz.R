library(tidyverse)
library(biomaRt)
library(Gviz)
library(GenomicRanges)
library(genomation)
library(trackViewer)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

#==================================================#
# Using Gvis to plot gene models as tracks with UCSC
#===================================================#
# make genome track
axTrack <- Gviz::GenomeAxisTrack()

# setup for using UCSC genome browser to get taf-7.2 from ce11
geneTrack <- Gviz::BiomartGeneRegionTrack(genome="ce11",
                                    chromosome="chrIII", start=12670689, end=12673423,
                                    transcriptAnnotation="symbol", name= "Genes")

# get Ideogram: a schematic representation of chromosome, showing their relative size and their cytostain banding pattern.
itrack <- Gviz::IdeogramTrack(genome = "ce11", chromosome = "chrIII")

# plot taf-7.2
Gviz::plotTracks(list(itrack, axTrack, geneTrack), chromosome="chrIII", from=12670689, to=12673423)

#=================================================================================#
# make a GeneRegionTrack from local gff file
#=================================================================================#
# make genome track
axTrack <- Gviz::GenomeAxisTrack()

# make GRanges object from Gff file
gene_GRs <- genomation::gffToGRanges(gff.file = "data/WS282_taf7.2.gff3")
# try to make TxDbF object from Gff file
gene_TxDb <- GenomicFeatures::makeTxDbFromGFF(file = "data/WS282_taf7.2.gff3", format = "gff3")

# make GeneRegionTrack object for Gvis
options(ucscChromosomeNames=FALSE)
gene_GRT <- as(gene_GRs, "GeneRegionTrack")
gene_GRT_TxDb <- Gviz::GeneRegionTrack(range = gene_TxDb, name = "taf-7.2", gene = "taf-7.2")

# plot taf-7.2
Gviz::plotTracks(list(axTrack, gene_GRT), from=12670689, to=12673423)
Gviz::plotTracks(list(axTrack, gene_GRT_TxDb))

#========================================================================================#
# get release gff from WB, grep to string, convert to GRanges object, make GeneRegionTrack
#=========================================================================================#
# make GeneRegionTrack from data.frame
taf72_gff <- data.table::fread("data/WS282_taf7.2.gff3") %>%
  dplyr::mutate()



#A data.frame object: the data.frame needs to contain at least the two mandatory columns start and end with the range coordinates.
#It may also contain a chromosome and a strand column with the chromosome and strand information for each range.
#If missing, this information will be drawn from the constructor's chromosome or strand arguments.
#In addition, the feature, exon, transcript, gene and symbol data can be provided as columns in the data.frame.
#The above comments about potential default values also apply here.


#========================================================================================#
# get release gff from WB, grep to string, convert to GRanges object, make GeneRegionTrack
#=========================================================================================#



#==========================================#
# old stuff
#=========================================#

#===================================================#
# get gff from wormbase or QUEST? or Local
#===================================================#
wormbase <- useMart(biomart = "parasite_mart", 
                    host = "https://parasite.wormbase.org", 
                    port = 443)
listDatasets(wormbase)
wormbase2 <- useDataset(mart = wormbase, dataset = "wbps_gene")

listAttributes(wormbase2)
listFilters(wormbase2)

test <- getBM(attributes = c("chromosome_name", 
  "display_name_1010",
                             "assembly_accession_1010",
                             "exon_chrom_start",
                             "wbps_gene_id",
                             "transcript_db_name",
                             "transcript_count",
                             "transcript_5_utr_start",
                             "transcript_5_utr_end",
                             "transcript_3_utr_start",
                             "transcript_3_utr_end",
                             "external_gene_id",
                             "wbps_transcript_id",
                             "transcript_biotype"), 
              filters = "wbps_transcript_id", 
              values = c("Y111B2A.16.1"), 
              mart = wormbase2)
test


#=======================================================#
# Using Gvis to plot gene models as tracks with gff3 files
#=======================================================#
tets <- data.table::fread("data/WS282_taf7.2.gff3")

# try to make if from TxDbF
gene2 <- GenomicFeatures::makeTxDbFromGFF(file = "data/WS282_taf7.2.gff3", format = "gff3")

# try to make it from GenomicRanges
gene3 <- GenomicRanges::GRanges()

gffRangedData<-import.gff("myFile.gff")
myGranges<-as(gffRangedData, "GRanges")

## Load some sample data
data(cyp2b10)

## Directly from the data.frame
gene2 <- gff82_taf %>%
  dplyr::mutate(width = NA,
                strand = NA,
                gene = NA,
                exon = NA,
                transcript = NA,
                symbol = NA) %>%
  dplyr::rename(chromosome = V1, start = V4, end = V5,feature = V3) %>%
  dplyr::filter(feature %in% c("five_prime_UTR", "CDS", "three_prime_UTR")) %>%
  dplyr::mutate(feature = case_when(feature == "five_prime_UTR" ~ "utr5",
                                    feature == "CDS" ~ "coding_region",
                                    feature == "three_prime_UTR" ~ "utr3",
                                    TRUE ~ NA_character_)) %>%
  dplyr::filter(grepl("Transcript:Y111B2A.16.1", V9)) %>%
  dplyr::mutate(width = end - start,
                strand = V7)

glimpse(cyp2b10)
glimpse(gene2)
grTrack <- GeneRegionTrack(cyp2b10)

gene2track <- GeneRegionTrack(range = gene2)

#========================================
gtfFile <- system.file("extdata","GTF_files","Aedes_aegypti.partial.gtf",
                       package="GenomicFeatures")
chrominfo <- data.frame(chrom = c('supercont1.1','supercont1.2'),
                        length=c(5220442, 5300000),
                        is_circular=c(FALSE, FALSE))
metadata <- data.frame(name="Resource URL",
                       value=paste0("ftp://ftp.ensemblgenomes.org/pub/metazoa/",
                                    "release-13/gtf/aedes_aegypti/"))
txdb2 <- GenomicFeatures::makeTxDbFromGFF(file=gtfFile,
                         chrominfo=chrominfo,
                         dataSource="ensemblgenomes",
                         organism="Aedes aegypti",
                         metadata=metadata)


#-------------------------------------------------------------------

## The empty object
GeneRegionTrack()

## Load some sample data
data(cyp2b10)

## Construct the object
grTrack <- GeneRegionTrack(start=26682683, end=26711643,
                           rstart=cyp2b10$start, rends=cyp2b10$end, chromosome=7, genome="mm9",
                           transcript=cyp2b10$transcript, gene=cyp2b10$gene, symbol=cyp2b10$symbol,
                           feature=cyp2b10$feature, exon=cyp2b10$exon,
                           name="Cyp2b10", strand=cyp2b10$strand)

## Directly from the data.frame
grTrack <- GeneRegionTrack(cyp2b10)
