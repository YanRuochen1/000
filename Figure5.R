rm(list = ls())
library(tidyverse)
library(IOBR)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(pheatmap)
library(cowplot)
blue <- "#4682B4"
red <- "#8B0000"
load("Figure5data.Rda")
cpt2_cibersort_s <- cpt2_cibersort[,c(1,7)]
cpt2_epic_s <- cpt2_epic[,c(1,5,6,10)]
cpt2_mcp_s <- cpt2_mcp[,c(1,8,13)]
cpt2_timer_s <- cpt2_timer[,c(1,4,6,9)]
cpt2_xcell_s <- cpt2_xcell[,c(1,6,9,19,28,29,31,35,37,43,44,47,51,52,57,62)]
tme_combine <- cpt2_cibersort_s %>% 
  inner_join(cpt2_mcp_s, "ID") %>% 
  inner_join(cpt2_xcell_s, "ID") %>%
  inner_join(cpt2_epic_s, "ID") %>% 
  inner_join(cpt2_timer_s, "ID") 
head(tme_combine)
tme_combine1 <-  column_to_rownames(tme_combine, "ID")
tme_combine2 <- scale(tme_combine1,
                      center = TRUE, scale = TRUE)
tme_combine2[tme_combine2 > 2] = 2
tme_combine2[tme_combine2 < (-2)] = -2
Type <- rbind(data.frame(ID = colnames(cpt2_cibersort_s)[2:ncol(cpt2_cibersort_s)], method = "CIBERSORT"),
              data.frame(ID = colnames(cpt2_mcp_s)[2:ncol(cpt2_mcp_s)], method = "MCP_counter"),
              data.frame(ID = colnames(cpt2_xcell_s)[2:ncol(cpt2_xcell_s)], method = "xCell"),
              data.frame(ID = colnames(cpt2_epic_s)[2:ncol(cpt2_epic_s)], method = "EPIC"),
              data.frame(ID = colnames(cpt2_timer_s)[2:ncol(cpt2_timer_s)], method = "TIMER"))
Type <- column_to_rownames(Type, "ID")
annCol <- target.gene.expr1 %>% column_to_rownames("ID")
colnames(annCol) <- c("CPT2","Group")
annRow <- Type %>%
  arrange(method)
library(RColorBrewer)
display.brewer.all()
brewer.pal(7, "Set1")
annColors <- list(method = c("CIBERSORT" = "#E41A1C", 
                             "EPIC" = "#377EB8",
                             "MCP_counter" = "#984EA3",
                             "TIMER" = "#FFFF33",
                             "xCell" = "#A65628"),
                  "CPT2" = colorRampPalette(c(blue, "white",red))(100), 
                  "Group" = c("high" = red, "low" = blue))
Sam_order <- rownames(annCol[order(annCol$CPT2),])
plot_data <- t(tme_combine2)
plot_data <- plot_data[rownames(annRow), Sam_order]
annCol <- annCol[Sam_order, ,drop = F]
pdf("Figure 5A_immune.pdf", 10, 12)
pheatmap(plot_data,
         color = colorRampPalette(c(blue, "white",red))(100), 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE, 
         show_colnames = FALSE, 
         fontsize_row = 5,
         annotation_col = annCol,
         annotation_row = annRow,
         annotation_colors = annColors, 
         gaps_col = table(annCol$Group)[2], 
         gaps_row = cumsum(table(annRow$method)) 
)
dev.off()
box_inhi1$Group <- factor(box_inhi1$Group,levels = c("low","high"))
ggplot(na.omit(box_inhi1), 
       aes(variable, value, fill = Group)) +
  geom_boxplot(alpha = 0.85) +
  scale_y_continuous(name = "Relative Expression (Log2+1)") +
  scale_fill_manual(values = c(blue, red))+
  scale_x_discrete(name = "") +
  theme_bw() +
  theme(text = element_text(size = 10),
        legend.position = "top",
        axis.text.x  = element_text(angle = 60, vjust = 1, hjust = 1)) + 
  stat_compare_means(aes(label = ..p.signif..))
ggsave("Figure 5B_Immunoinhibitor.pdf", width = 10, height = 6)
dev.off()
library(tidyverse)
library(UCSCXenaTools)
library(maftools)
library(IOBR)
library(cowplot)
High <- target.gene.expr[which(target.gene.expr$Group == "high"), ]
High_maf <- maf_TCGA1 %>%
  dplyr::filter(Tumor_Sample_Barcode %in% High$ID)
Low <- target.gene.expr[which(target.gene.expr$Group == "low"), ]
Low_maf <- maf_TCGA1 %>%
  dplyr::filter(Tumor_Sample_Barcode %in% Low$ID)
High <- read.maf(maf = High_maf, isTCGA = TRUE)
pdf("Figure 5C_SNP_High_cpt2.pdf")
oncoplot(maf = High,
         top = 20, 
         fontSize = 0.8,
         legendFontSize = 1,
         annotationFontSize = 1)
dev.off()
Low <- read.maf(maf = Low_maf, isTCGA = TRUE)
pdf("Figure 5D_SNP_Low_cpt2.pdf")
oncoplot(maf = Low,
         top = 20, 
         fontSize = 0.8,
         legendFontSize = 1,
         annotationFontSize = 1)
dev.off()
