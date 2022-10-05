rm(list=ls())
library(VennDiagram)  
options(stringsAsFactors = FALSE) 
library(tidyverse)
MCC <- read.csv("MCC.csv") [2] %>% 
  setNames("MCC") %>% 
  as.matrix()
EPC <- read.csv("EPC.csv")[2] %>% 
  setNames("EPC") %>% 
  as.matrix()
DMNC <- read.csv("DMNC.csv")[2] %>% 
  setNames("DMNC") %>% 
  as.matrix()
DEGREE <- read.csv("DEGREE.csv")[2] %>% 
  setNames("DEGREE") %>% 
  as.matrix()
data <- list(MCC,EPC,DMNC,DEGREE)
names(data) <- c("MCC","EPC","DMNC","DEGREE")
library(ggvenn)
pdf("Figure 3B_veenplot.pdf")
ggvenn(
  data,
  columns = NULL,
  show_elements = FALSE,
  show_percentage = FALSE,
  fill_color = c("blue", "yellow", "green", "red"),
  fill_alpha = 0.5,
  stroke_color = "black",
  stroke_alpha = 1,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = "black",
  set_name_size = 6,
  text_color = "black",
  text_size = 4,
  label_sep = ","
)
dev.off()
library(survival)
library(survminer)
library(dplyr)
cli_deg1 <- read.csv("COAD_tumor_clinical1.csv")
outTab=data.frame()
for(i in colnames(cli_deg1[,3:ncol(cli_deg1)])){
  cox <- coxph(Surv(OS.Time, event) ~ cli_deg1[,i], data = cli_deg1)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.csv(outTab,"Table S1_all_uncox.csv")
outTab <- outTab %>%
  mutate(res= ifelse(outTab$pvalue<0.05,"pass","no"))
outTab_pass <- outTab %>%
  filter(res=="pass")
cli_deg3 <- cli_deg1 %>%
  dplyr::select(c("OS.Time","event",outTab_pass$id))
multiCox=coxph(Surv(OS.Time, event) ~ ., data = cli_deg3)
multiCoxSum=summary(multiCox)
outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"]
)
library(forestplot)
pdf("Figure 3C_multi_forest.pdf")
ggforest(multiCox,
         main = "HR",
         cpositions = c(0.02,0.22, 0.4), 
         fontsize = 0.7, 
         refLabel = "reference", 
         noDigits = 4)
dev.off()
outTab=cbind(id=row.names(outTab),outTab) %>% 
  as.data.frame() 
outTab_s <- outTab %>%
  dplyr::filter(pvalue<0.05)
genes <- intersect(MCC,intersect(DMNC,intersect(EPC,DEGREE)))
genes <- as.data.frame(genes)
degs <- intersect(outTab_s$id,genes$genes)
library(ggvenn)
x <- list(Independent_Prognostic_Factor =outTab_s$id ,Hub_genes = genes$genes)  
pdf("Figure 3D_Independent_Prognostic_Factor_cytoscape.pdf")
ggvenn(
  x,
  columns = NULL,
  show_elements = FALSE,
  show_percentage = F,
  fill_color = c("steelblue", "darkred"),
  fill_alpha = 0.4,
  stroke_color = "black",
  stroke_alpha = 0.1,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = "black",
  set_name_size = 4,
  text_color = "black",
  text_size = 4,
  label_sep = ","
)
dev.off()
