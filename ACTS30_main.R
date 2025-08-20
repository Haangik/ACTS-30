### ACTS-30 single cell RNA-seq analysis
rm(list=ls())

# Library preparation
require(Seurat)
require(dplyr)
require(scRepertoire)
require(ggplot2)
require(ggpubr)
require(pbapply)
require(parallel)
require(tidyverse)
require(RColorBrewer)
require(spacexr)
require(ggrepel)
require(SeuratWrappers)
require(monocle3)
require(tradeSeq)
require(slingshot)
require(viridis)
require(cowplot)
require(gridExtra)
require(scales)
require(PseudotimeDE)
require(genekitr)
require(org.Hs.eg.db)
require(clusterProfiler)
require(CellChat)

# Loading required data
### Treatment outcomes of patients
patients_pCR<-c('ACTSCRHBO00100T', 'ACTSCRHBO00500T', 'ACTSCRHBO01300T',
                'ACTSCRHBO01800T', 'ACTSCRHBO02000T', 'ACTSCRHBO02200T',
                'ACTSCRHBO02400T', 'ACTSCRHBO02500T')
patients_nonpCR<-c('ACTSCRHBO00600T', 'ACTSCRHBO00700T', 'ACTSCRHBO00900T',
                   'ACTSCRHBO01000T', 'ACTSCRHBO01400T', 'ACTSCRHBO01500T',
                   'ACTSCRHBO01700T', 'ACTSCRHBO02300T', 'ACTSCRHBO02700T',
                   'ACTSCRHBO02800T', 'ACTSCRHBO02900T')

### Data was pre-processed to .RDS files
Total <- readRDS('Data/ACTS30_Total.RDS') # Total scRNA-seq data object
CD4T<-readRDS("Data/CD4T.RDS") # Re-clustered CD4+ T cell subset
CD8T<-readRDS("Data/CD8T.RDS") # Re-clustered CD8+ T cell subset
myeloid<-readRDS("Data/ACTS_Myeloid.RDS")
CAF<-readRDS("Data/ACTS_CAF.RDS")

### Generation of merged T cell metadata
Total_T <- merge(
  CD4T,
  y = list(CD8T),
  add.cell.ids = NULL
)
Total_T@meta.data<-droplevels(Total_T@meta.data)
Total_T@meta.data$sample_id_TCR<-sapply(Total_T@meta.data$sample_id, FUN=function(x){
  sub('RHB', 'JHB', x)
})

### Pre-processed T cell receptor contig list
contig_list<-readRDS("Data/TCR_Contig.RDS")
tcell_barcode<-intersect(Total_T@meta.data$cell, unlist(lapply(contig_list, FUN=function(x){x$barcode})))

### metadata summarized by patient level
sm <- Total_T@meta.data %>%
  mutate(response_2 = ifelse(main_group1=="pCR", "pCR", "Non-pCR")) %>%
  mutate(response_3 = main_group1)
sm <- sm[,c("sample_id", "response_2", "response_3")] %>% distinct

### Replacement of sample number in TCR sequencing data
barcode_sample_mapping<-Total_T@meta.data %>% dplyr::select(sample_id, sample_num) %>% unique

### Spatial transcriptomic data
slide.seq_PCR <- readRDS('Data/slide.seq_PCR.RDS')
slide.seq_NonPCR <- readRDS('Data/slide.seq_NonPCR.RDS')
RCTD_PCR <- readRDS('Data/RCTD_main_PCR.RDS')
RCTD_NonPCR <- readRDS('Data/RCTD_main_NonPCR.RDS')
RCTD_PCR_CD8 <- readRDS('Data/RCTD_CD8T_PCR.RDS')
RCTD_NonPCR_CD8 <- readRDS('Data/RCTD_CD8T_NonPCR.RDS')
RCTD_PCR_CD4 <- readRDS('Data/RCTD_CD4T_PCR.RDS')
RCTD_NonPCR_CD4 <- readRDS('Data/RCTD_CD4T_NonPCR.RDS')
RCTD_PCR_merge <- readRDS('RCTD_Merge_PCR.RDS')
RCTD_NonPCR_merge <- readRDS('RCTD_Merge_NonPCR.RDS')

# Article figures
source("src/Figure1.R") # Figure 1
source("src/Figure2.R") # Figure 2
source("src/Figure3.R") # Figure 3
source("src/Figure4.R") # Figure 4
source("src/Figure5.R") # Figure 5
source("src/Figure6.R") # Figure 6
