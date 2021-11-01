library(tidyverse)
# format cM file for Tim ce popgen PCA analysis
# need: (1) chrom (2) start pos (3) end pos (4) cM

# download gff from wormbase
# wget ftp://wormbase.org/pub/wormbase/releases/WS276/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS76.annotations.gff3.gz

# pull out physical map
# zcat c_elegans.PRJNA13758.WS276.annotations.gff3.gz | grep pmap > c_elegans.PRJNA13758.WS276.annotations.cM.gff3

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# re-format in R
rm_row <- function(df) {
    
    df %<>% dplyr::arrange(pos)
    
    keep_fw <- vector()
    keep_rv <- vector()
    
    if (df[1, "cM"] > df[2, "cM"] | df[nrow(df), "cM"] < df[nrow(df)-1, "cM"] ) {
        #        print("cM in first/last position needs investigation")
        df = df[-nrow(df), ]
    } 
    
    keep_fw[1] <- TRUE
    keep_fw[nrow(df)] <- TRUE
    
    keep_rv[1] <- TRUE
    keep_rv[nrow(df)] <- TRUE
    
    prev_pos <- df[1, "cM"] # for irritarate downward
    last_pos <- df[nrow(df), "cM"] # for irritarate backward
    
    for (n in 2:(nrow(df)-1)) {
        keep_fw[n] <- (prev_pos < df[n, "cM"]) & (df[n, "cM"] < df[n+1, "cM"])
        if (keep_fw[n] == TRUE) { 
            prev_pos <- df[n, "cM"] 
        }
    }
    
    for (n in (nrow(df)-1):2) {
        keep_rv[n] <- (last_pos > df[n, "cM"]) & (df[n, "cM"] > df[n-1, "cM"])
        if (keep_rv[n] == TRUE) { 
            last_pos <- df[n, "cM"] 
        }
    }
    
    if (sum(keep_fw) >= sum(keep_rv)) {
        keep <- keep_fw
    } else { keep <- keep_rv }
    
    
    df2 <- df[keep,] 
    out <- df2 %>% 
        dplyr::mutate(V3=".") %>% 
        dplyr::mutate(cM = cM - min(df2$cM, na.rm = T)) %>% 
        dplyr::arrange(pos) %>% 
        dplyr::select(V1, V3, cM, pos)
    
    if (out$cM[1] == 0) {
        chr <- out$V1[1]
        # write.table(out, paste0("chr", chr, ".map"), quote=F, sep=" ", row.names = F, col.names = F)
    }
    
    return(out)
}

gff <- read.delim("data/c_elegans.PRJNA13758.WS276.annotations.cM.gff3", header=FALSE, stringsAsFactors=FALSE)
gff2 <- gff %>% 
    dplyr::filter(str_detect(V2, "abs")) %>% 
    tidyr::separate(V9, c("others","pmap"), sep = "Note=") %>% 
    tidyr::separate(pmap, c("center","span"), sep=" cM \\(\\+\\/\\- ") %>% 
    dplyr::mutate(span=str_replace(span, " cM\\)", "")) %>% 
    dplyr::select(V1, V2, V4, V5, center, span) %>%  
    dplyr::mutate(center=as.numeric(center), span=as.numeric(span)) %>%  
    dplyr::mutate(left=center-span, right=center+span)

start_pos <- gff2 %>%
    dplyr::select(V1, V2, V4, left) %>% 
    dplyr::rename(pos=V4, cM=left) %>% 
    dplyr::arrange(V1, pos)

# output for beagle
df_new <- start_pos %>% 
    dplyr::filter(cM != 0) %>% 
    dplyr::distinct(cM, .keep_all = T) %>% 
    dplyr::distinct(pos, .keep_all = T) %>% 
    dplyr::group_by(V1) #%>% 
    #dplyr::do( rm_row(.) ) 

# output for cepopgen
df_cm <- df_new %>%
    dplyr::mutate(startpos = dplyr::lag(pos)) %>%
    dplyr::mutate(startpos = ifelse(is.na(startpos), 1, startpos)) %>%
    dplyr::select(chrom = V1, startpos, endpos = pos, cM)

# find our region
lifespan_qtl <- df_cm %>%
    dplyr::filter(chrom == "III" & (startpos > 11000000 & startpos < 12500000))
#readr::write_tsv(df_cm, "Ce_Genetic_Map_WS276.bed", col_names = FALSE)    
    
# 12186411 - 11746766 physical
d1 <- 12.1556 - 15.9044 
d2 <- 14.8405 - 15.9044
d3 <- 14.6970 - 15.9044


# 11505488 - 12052420 physical = ORDER
test <- 10.5024 - 14.1221
