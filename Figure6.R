rm(list = ls())
library(tidyverse)
library(ggpubr)
library(dplyr)
load("Figure6data.Rda")
 methy450_s_tumor1$CPT2_Group <- ifelse( methy450_s_tumor1$CPT2 >median( methy450_s_tumor1$CPT2),"High","Low")
 methy450_s_tumor1$CPT2_Group <- factor(methy450_s_tumor1$CPT2_Group,levels = c("Low","High"))
 compaired <- list(c("Low","High")) 
 ggboxplot(methy450_s_tumor1,
           x = "CPT2_Group", y = "cg01810926",
           fill = "CPT2_Group", palette = c("steelblue" ,"darkred"),alpha=0.6)+
   theme(legend.position = "right")+
   stat_compare_means(comparisons = compaired,
                      method = "wilcox.test",   
                      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                       symbols = c("***", "**", "*", "ns")))  +
   theme_bw()+
   labs(x = "", y = "cg01810926") 
 ggsave("Figure 6B_cg01810926 vs CPT2.pdf")
 dev.off()
 ggplot(data = methy450_s_tumor1, aes(x = CPT2, y = cg01810926)) + 
    geom_point(shape = 21, fill = "steelblue",
               colour = "black", size = 4) +
    geom_smooth(method = "lm", fill = "gray70",
                color = "red") +
    ggtitle("Correlation between CPT2 and cg01810926") +
    theme(plot.title = element_text(hjust = 0.50,
                                    color="black",
                                    size = 15)) +
    ggpubr::stat_regline_equation(label.x = 3.5, label.y = 1,
                                  size = 4.5, fontface = "bold") +
    ggpubr::stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
                     label.x = 3.5, label.y = 0.9,
                     size = 4.5, fontface = "bold")
 ggsave("Figure 6C_correlation_cg01810926 vs CPT2.pdf",width = 8., height =8)
 dev.off()
 library(survival)
 library(survminer)
 pdf("Figure 6D_cg01810926_CPT2.pdf")
 point <- surv_cutpoint(methy450_s_tumor,"OS.Time","event",variables = "cg01810926")
 summary(point)
 res.cat <- surv_categorize(point) 
 colnames(res.cat)[3] <- "cg01810926"
 fit1 <- survfit(Surv(OS.Time, event) ~cg01810926, data = res.cat)
 ggsurvplot(fit1,
            data = res.cat,
            risk.table = TRUE,
            palette = c("darkred", "steelblue"),
            pval = T)
 dev.off()
 