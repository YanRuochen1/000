rm(list=ls())   
options(stringsAsFactors = F)
library(tidyverse)
library(multiMiR)
gene_diff <- c("CPT2")
gene2mir <- get_multimir(org     = 'hsa',
                         target  = gene_diff,
                         table   = 'validated',
                         summary = TRUE,
                         predicted.cutoff.type = 'n',
                         predicted.cutoff      = 500000)
ez = gene2mir@data[gene2mir@data$database=="tarbase",]
ez1 = gene2mir@data[gene2mir@data$database=="mirtarbase",]
ez <- rbind(ez,ez1)
dim(ez)
miRNAs = unique(ez$mature_mirna_id)  
mRNAs = unique(ez$target_symbol)
node1 <- ez[,c("mature_mirna_id","target_symbol")]
colnames(node1) <- c("gene1", "gene2")
head(node1)
write.csv(node1,"Table S2 node_tarbase&mitarbase_m_mi.csv")
starbase <- read.csv("starBaseV3_lncRNA.csv")
dim(starbase)
miRNAname1 <- c("hsa-miR-30b-5p","hsa-miR-106a-5p","hsa-miR-30a-5p")
lnc_mi = starbase[starbase$miRNAname %in% miRNAname1,]
length(unique(lnc_mi$miRNAname))
idx = lnc_mi$pancancerNum >0 & lnc_mi$clipExpNum>1
table(idx)
lncRNAs = unique(lnc_mi$geneName[idx])
node2 <- lnc_mi[lnc_mi$geneName %in% lncRNAs,]
node2 <- node2[,c("miRNAname","geneName")] %>% 
  as.data.frame() %>% 
  distinct()
colnames(node2) <- c("gene1", "gene2")

