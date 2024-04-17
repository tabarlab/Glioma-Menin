### Required Libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(UCell)
set.seed(1234)
library(future)
plan(multicore, workers = 10)
options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM

######## Read Seurat Object
RNA_all <- readRDS('./scRNA_YY_combined.rds')

######## Progenitor Gene Module from Ruiz-Moreno et al. 2022 (https://doi.org/10.1101/2022.08.27.505439)
progenitor_genes <- c('CDK4', 'TUBA1A', 'TUBB2B', 'SOX4', 'MARCKSL1', 
                        'BCAN', 'HES6', 'STMN1', 'CKB', 'RP11-620J15.3', 
                        'GPM6A', 'MTRNR2L8', 'NNAT', 'RBP1', 'SCRG1', 
                        'PTPRZ1', 'NOVA1', 'FXYD6', 'TUBB', 'PFN2', 
                        'MLLT11', 'OLIG1', 'STMN2', 'PTN', 'UCHL1', 
                        'DLL3', 'MAP2', 'IGFBP2', 'SOX11', 'H4C3', 'OS9', 
                        'BEX1', 'TSFM', 'NFIB', 'MIR9-1HG', 'MEG3', 
                        'OLIG2', 'TCF4', 'TUBB2A', 'PCSK1N', 'MAP1B', 
                        'MARCHF9', 'TUBA1B', 'TSPAN31', 'SCD5')

######## Calculate Module Score 
RNA_all <- AddModuleScore_UCell(RNA_subset, features = list(progenitor_genes), name = "progenitor")

######## Figure 6B
RNA_all[[]] %>% 
  ggplot(aes(x=cell.line, y=signature_1progenitor, fill=MEN1.sensitivity)) +
  geom_boxplot()+
  scale_x_discrete(limits=c("MSK-20","MSK-19", "MSK-28","MSK-23", 'MSK-18', "MSK-04")) + 
  scale_fill_manual(values = c('#EB466F', "#26547C")) +
  theme_bw() +
  ylim(0.05,0.6)+
  xlab("Cell Line") + 
  ylab("Progenitor Score") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size=20)) -> p1

ggsave(filename = 'boxplot_cell_line_progenitorScore.pdf', plot = p1, height = 10, width = 10, dpi = 300)

######## Extended Figure 6C
RNA_all[[]] %>% 
  ggplot(aes(x=cell.line, y=signature_1progenitor, fill=MEN1.sensitivity)) +
  #geom_hline(yintercept = 0, color='red', size=.5) +
  geom_boxplot()+
  scale_x_discrete(limits=c("MSK-19", "MSK-19KO")) +
  scale_fill_manual(values = c('#EB466F', "#FEA82F")) + 
  theme_bw() +
  ylim(0.05,0.6)+
  xlab("Cell Line") + 
  ylab("Progenitor Score") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size=20)) -> p2

ggsave(filename = 'boxplot_WT_vs_KO_progenitorScore.pdf', plot = p2, height = 7, width = 6, dpi = 300)