rm(list = ls())
library(tidyverse)
library(DESeq2)
library(limma)
library(factoextra)
library(FactoMineR)
library(RColorBrewer)
res1 <- read.csv("DEG.csv") %>% 
  column_to_rownames("Symbol")
res_deseq_select1 <- res1 %>% 
  arrange(padj) 
log2FoldChange_cutoff <- 1 
type1 = (res_deseq_select1$padj)<0.05 & (res_deseq_select1$log2FoldChange < -log2FoldChange_cutoff)
type2 = (res_deseq_select1$padj)<0.05 & (res_deseq_select1$log2FoldChange > log2FoldChange_cutoff)  
res_deseq_select1$change = ifelse(type1, "down",ifelse(type2,"up","not"))
table(res_deseq_select1$change)
library(ggplot2)
pdf("Figure 2A_mrna_Volcano.pdf")
vocano_result <- res1 %>%
  mutate(Type = if_else(padj > 0.05 | abs(log2FoldChange) < 1, 
                        'No_Significant', 
                        if_else(log2FoldChange >= 1, 
                                'Up_Regulated', 
                                'Down_Regulated'))) 
ggplot(vocano_result,
       aes(log2FoldChange,
           -log10(padj))) +
  geom_point(size = 5, 
             alpha = 0.65,
             aes(color = Type),
             show.legend = T) +
  scale_color_manual(values = c("steelblue", "gray70", "darkred")) +
  ylim(0, 250) +
  xlim(-10, 12) +
  labs(x = "Log2(fold change)", y = "-log10(p-value)") +
  geom_hline(yintercept = -log10(0.05),  
             linetype = 'dotdash', 
             color = 'gray70') + 
  geom_vline(xintercept = c(-1, 1), 
             linetype = 'dotdash', 
             color = 'gray70')
dev.off()
down <- res_deseq_select1 %>% 
  dplyr::filter(change == "down") %>% 
  rownames_to_column("Symbol")
up <- res_deseq_select1 %>%  
  dplyr::filter(change == "up") %>% 
  rownames_to_column("Symbol")
LMRGS_gene <- read.csv("LMRGS.csv",header = F) %>% 
  t() %>% 
  unique() %>% 
  as.data.frame()
rownames(LMRGS_gene) <- NULL
colnames(LMRGS_gene) <- "Symbol"
down_gene <- down %>% 
  filter(Symbol %in% intersect(down$Symbol,LMRGS_gene$Symbol))
up_gene <- up %>% 
  filter(Symbol %in% intersect(up$Symbol,LMRGS_gene$Symbol))
deg_LMRGS <- rbind(down_gene,up_gene)
write.csv(deg_LMRGS,file = "deg_LMRGS_gene_1.csv") 
library(ggvenn)
x <- list(Down_Regulated = down$Symbol,LMRGs = LMRGS_gene$Symbol)   
pdf("Figure 2B_degs_down vs LMRGS.pdf")
ggvenn(
  x,
  columns = NULL,
  show_elements = FALSE,
  show_percentage = F,
  fill_color = c("steelblue", "darkred"),
  fill_alpha = 0.4,
  stroke_color = "black",
  stroke_alpha = 0.5,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = "black",
  set_name_size = 4,
  text_color = "black",
  text_size = 4,
  label_sep = ","
)
dev.off()
x <- list(Up_Regulated = up$Symbol,LMRGs = LMRGS_gene$Symbol) 
pdf("Figure 2B_degs_up vs LMRGS.pdf")
ggvenn(
  x,
  columns = NULL,
  show_elements = FALSE,
  show_percentage = F,
  fill_color = c("steelblue", "darkred"),
  fill_alpha = 0.4,
  stroke_color = "black",
  stroke_alpha = 0.5,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = "black",
  set_name_size = 4,
  text_color = "black",
  text_size = 4,
  label_sep = ","
)
dev.off()
library(R.utils)
R.utils::setOption("clusterProfiler.download.method",'auto')
gene <- read.csv("deg_LMRGS_gene_1.csv")
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
genelist <- bitr(gene$Symbol,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = "org.Hs.eg.db")
gene_f <- gene %>% 
  left_join(genelist, by = c("Symbol"="SYMBOL" )) 
kk <- enrichKEGG(gene         =  gene_f$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.1,
                 qvalueCutoff =0.1)
save(kk, file = "0301KEGG_TCGA_COAD.Rdata")
pdf("Figure 2Ckegg_barplot.pdf")
barplot(kk, showCategory = 10,color = "pvalue")
dev.off()
gene <- read.csv("deg_LMRGS_gene_1.csv")
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
genelist <- bitr(gene$Symbol,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = "org.Hs.eg.db")
gene_f <- gene %>% 
  left_join(genelist, by = c("Symbol"="SYMBOL" )) 
ego <- enrichGO(gene = gene_f$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)
save(ego, file = "0301GO_TCGA_COAD.Rdata")
pdf("Figure 2D_go_barplot.pdf")
barplot(ego, drop = TRUE, showCategory =5,split="ONTOLOGY" ) + 
  facet_grid(ONTOLOGY~., scale='free')
dev.off()


