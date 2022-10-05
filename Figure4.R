rm(list=ls())   
options(stringsAsFactors = F)
library(tidyverse)
library(beeswarm)
library(ggplot2) 
library(reshape2)
library(ggpubr)
load("figure4data.RDA")
expr$Group <- factor(expr$Group,levels = c("Normal","Tumor"))
exp$Group <- factor(exp$Group,
                    levels = c("Normal","Tumor"))
compaired <- list(c("Normal","Tumor")) 
library(ggpubr)
pdf("Figure 4B_TvsN_CPT2.pdf")
ggboxplot(exp,
          x = "Group", y = "CPT2",
          fill = "Group", palette = c("steelblue","darkred"),alpha = 0.6)+
  stat_compare_means(comparisons = compaired,
                     method = "wilcox.test",   
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  labs(x = "", y = "CPT2 Expression Levels (log2(TPM+1))") +
  theme_bw()
dev.off()
pdf("Figure 4C_Pair_CPT2.pdf")
compaired <- list(c("Normal","Tumor"))
ggplot(expr, aes(x = Group, y = CPT2, group = names)) +
  geom_line(color="black", size=1,alpha =0.4) +
  geom_point(aes(fill=Group), shape = 21, size = 3,alpha = 0.6)+
  labs(x = "", y = "CPT2 Expression Levels (log2(TPM+1))") +
  theme_bw()+
  scale_fill_manual(values = c("steelblue","darkred"))+
  geom_signif(comparisons = list(c("Normal","Tumor")),
              map_signif_level = T)
dev.off()
library(RColorBrewer)
library(ggpubr)
library(tidyverse)
library(dplyr)
table(exp$Pathologic_N)
pdf("Figure 4D_1_Pathologic_N_CPT2.pdf")
exp1$Pathologic_N <- ifelse(exp1$Pathologic_N == "N0","N0","N1 & N2")
exp1$Pathologic_N <- factor(exp1$Pathologic_N,levels = c("N0","N1 & N2"))
compaired <- list(c("N0","N1 & N2")) 
ggboxplot(exp1,
          x = "Pathologic_N", y = "CPT2",
          fill = "Pathologic_N", palette = c("steelblue" ,"darkred"),alpha=0.6)+
  theme(legend.position = "right")+
  stat_compare_means(comparisons = compaired,
                     method = "wilcox.test",   
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))  +
  theme_bw()+
  labs(x = "", y = "CPT2 Expression Levels (log2(TPM+1))") 
dev.off()
pdf("Figure 4D_2_Pathologic_M_CPT2.pdf")
exp1$Pathologic_M <- factor(exp1$Pathologic_M, levels = c('M0','M1'))
compaired <- list(c('M0','M1')) 
ggboxplot(exp1,
          x = "Pathologic_M", y = "CPT2",
          fill = "Pathologic_M", palette = c("steelblue" ,"darkred"),alpha=0.6)+
  theme(legend.position = "right")+
  stat_compare_means(comparisons = compaired,
                     method = "wilcox.test",   
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))  +
  theme_bw()+
  labs(x = "", y = "CPT2 Expression Levels (log2(TPM+1))") 
dev.off()
table(exp1$Tumor_Stage)
pdf("Figure 4D_3_Tumor_Stage_CPT2.pdf")
exp1$Tumor_Stage <- ifelse(exp1$Tumor_Stage %in% c("III","IV"),"III & IV","I & II")
exp1$Tumor_Stage <- factor(exp1$Tumor_Stage, levels = c("I & II",'III & IV'))
compaired <- list(c("I & II",'III & IV'))
ggboxplot(exp1,
          x = "Tumor_Stage", y = "CPT2",
          fill = "Tumor_Stage", palette = c("steelblue" ,"darkred"),alpha=0.6)+
  theme(legend.position = "right")+
  stat_compare_means(comparisons = compaired,
                     method = "wilcox.test",   
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))  +
  theme_bw()+
  labs(x = "", y = "CPT2 Expression Levels (log2(TPM+1))") 
dev.off()
pdf("Figure 4D_4__Lymphatic_Invasion_CPT2.pdf")
exp1$Lymphatic_Invasion <- factor(exp1$Lymphatic_Invasion,
                                 levels = c('NO','YES'))
compaired <- list(c('NO','YES')) 
ggboxplot(exp1,
          x = "Lymphatic_Invasion", y = "CPT2",
          fill = "Lymphatic_Invasion", palette = c("steelblue" ,"darkred"),alpha=0.6)+
  theme(legend.position = "right")+
  stat_compare_means(comparisons = compaired,
                     method = "wilcox.test", 
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))  +
  theme_bw()+
  labs(x = "", y = "CPT2 Expression Levels (log2(TPM+1))") 
dev.off()

rm(list = ls())
library(ggpubr)
library(RColorBrewer)
library(tidyverse)
library(survival)
library(survminer)
library(dplyr)
load("figure4data.RDA")
exp1$OS.Time <- exp1$OS.Time/365
exp1$CPT2_Group <- ifelse( exp1$CPT2 >median( exp1$CPT2),"High","Low")
exp1$CPT2_Group <- factor(exp1$CPT2_Group,levels = c("Low","High"))
point <- surv_cutpoint(exp1,"OS.Time","Event",variables = "CPT2")
summary(point)
res.cat <- surv_categorize(point) 
colnames(res.cat)[3] <- "CPT2"
res.cat$CPT2 <- ifelse(res.cat$CPT2 == "low","Low","High")
exp1$Pathologic_N <- ifelse(exp1$Pathologic_N == "N0","N0","N1 & N2")
exp1$Pathologic_N <- factor(exp1$Pathologic_N,levels = c("N0","N1 & N2"))
res.cat_Pathologic_N <- cbind(res.cat,exp1["Pathologic_N"])
fit4 <- survfit(Surv(OS.Time, Event) ~CPT2+Pathologic_N, data = res.cat_Pathologic_N)
pdf("Figure 4E_1_CPT2+Pathologic_N.pdf",width = 13.5,height = 8)
ggsurvplot(fit4,
           data = res.cat_Pathologic_N,
           risk.table = TRUE,
           palette = c("darkred","FireBrick","CornflowerBlue", "steelblue"),
           pval = T,
           legend=c(0.85,0.85))
dev.off()
exp1$Pathologic_M <- factor(exp1$Pathologic_M, levels = c('M0','M1'))
res.cat_Pathologic_M <- cbind(res.cat,exp1["Pathologic_M"])
fit3 <- survfit(Surv(OS.Time, Event) ~CPT2+Pathologic_M, data = res.cat_Pathologic_M)
pdf("Figure 4E_2_CPT2+Pathologic_M.pdf",width = 13.5,height = 8)
ggsurvplot(fit3,
           data = res.cat_Pathologic_M,
           risk.table = TRUE,
           palette = c("darkred","FireBrick","CornflowerBlue", "steelblue"),
           pval = T,
           legend=c(0.85,0.85))
dev.off()
exp1$Tumor_Stage <- ifelse(exp1$Tumor_Stage %in% c("III","IV"),"III & IV","I & II")
exp1$Tumor_Stage <- factor(exp1$Tumor_Stage, levels = c("I & II",'III & IV'))
res.cat_Tumor_Stage <- cbind(res.cat,exp1["Tumor_Stage"])
pdf("Figure 4E_3_CPT2+Tumor_Stage.pdf",width = 13.5,height = 8)
fit2 <- survfit(Surv(OS.Time, Event) ~CPT2+Tumor_Stage, data = res.cat_Tumor_Stage)
ggsurvplot(fit2,
           data = res.cat_Tumor_Stage,
           risk.table = TRUE,
           palette = c("darkred","FireBrick","CornflowerBlue", "steelblue"),
           pval = T,
           legend=c(0.85,0.85))
dev.off()
pdf("Figure 4E_4_CPT2+Lymphatic_Invasion.pdf",width = 13.5,height = 8)
res.cat_lym <- cbind(res.cat,exp1[10])
fit1 <- survfit(Surv(OS.Time, Event) ~CPT2+Lymphatic_Invasion, data = res.cat_lym)
ggsurvplot(fit1,
           data = res.cat_lym,
           risk.table = TRUE,
           palette = c("darkred","FireBrick","CornflowerBlue", "steelblue"),
           pval = T,
           legend=c(0.85,0.85))
dev.off()






