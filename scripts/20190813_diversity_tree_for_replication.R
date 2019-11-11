library(tidyverse)
library(ggtree)
library(ape)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load tree
tree <- ape::read.tree(glue::glue("Data/Trees/330_genome.raxml.bestTree"))

# CeNDR sets
set1 <- c("BRC20067", "ECA246", "ECA251", "CX11271", "CX11276", "CX11285", "DL200", "DL226", "ECA36",
          "ED3040", "ED3048", "ED3049", "ED3052", "JU1088", "JU1172", "JU1213","JU1242", "JU1409",
          "JU1440", "JU1568", "JU1581", "JU2007", "JU2466", "JU2519", "JU311", "JU323", "JU346", "JU360",
          "JU440", "JU774", "JU830", "KR314", "MY1", "NIC1", "NIC199", "NIC2", "NIC252", "NIC256", "NIC260",
          "NIC261", "NIC274", "NIC275", "NIC277", "NIC3", "PB303", "QG2075", "QX1794", "WN2002")

set2 <- c("AB1",
          "ECA243", "CB4854", "ECA248", "CB4932", "ED3005", "ED3046", "ED3073", "ED3077", "EG4349",
          "JU1200", "JU1212", "JU1246", "JU1400", "JU1530", "JU1580", "JU1586", "JU2001", "JU2464",
          "JU2513", "JU2522", "JU310", "JU367", "JU393", "JU394", "JU406", "JU778", "JU792", "JU847",
          "MY10", "MY18", "NIC166", "NIC195", "NIC207", "NIC236", "NIC242", "NIC251", "NIC255",
          "NIC265", "NIC266", "NIC271", "NIC272", "NIC276", "ECA259", "PX179", "QG557", "QW947","QX1212")

set3 <- c("CB4852", "ECA250", "ECA252", "CX11254", "CX11262", "CX11264", "CX11292",
          "CX11307", "CX11315", "ED3011", "ED3012", "EG4347", "EG4724", "EG4946", "GXW1",
          "JU1395", "JU1491", "JU1652", "JU1896", "JU2316", "JU2526", "JU397", "JU561", "JU642", "JU751",
          "JU782", "LSJ1", "NIC231", "NIC258", "NIC259", "NIC262", "NIC267", "NIC268", "NIC269", "PS2025",
          "QG536", "QG556", "QX1211", "QX1233", "QX1791", "QX1792", "QX1793", "RC301", "WN2001", "NIC526", 
          "NIC527", "NIC528", "NIC529")

lifespan_GWA_replication_set <- c("CX11315","EG4724","JU2526", "JU397","NIC268","QX1211","QX1791")
# define colors
highlight_color <- "#BE0032"
background_color <- "#848482"

# highlight branches for strains of interest
branch_strains <- list(CONNECT = lifespan_GWA_replication_set)

tree_pt_h <- ggtree::groupOTU(tree, branch_strains)

plot <- ggtree(tree_pt_h,
       branch.length="rate", 
       aes(color=group)) + 
  geom_tiplab(align = F) +
  scale_color_manual(values=c(background_color, highlight_color), 
                     name = "Presence of TALT", 
                     breaks=c("0", "TALT"),
                     labels=c("FALSE", "TRUE")) + 
  theme(legend.position="right")+
  theme_tree2() 

ggsave('plots/lifespan_GWA_replication_set.pdf', height = 30, width = 7.5)
