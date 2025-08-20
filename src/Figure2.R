set.seed(501)

# Figure 2A. UMAP visualization of intra-tumoral CD8+ T cells
use_color <- c('CD8 Tn'= '#F8766D',
               'CD8 Tcm'= '#C49A00',
               'CD8 Tem' = '#53B400',
               'CD8 Trm'='#00BD8F', 
               'CD8 Temra'= '#00B6EB', 
               'CD8 Tpex'= '#9F82FF',
               'CD8 Tex' = '#FB58D4')
Figure2A<-DimPlot(CD8T, cols=use_color, group.by = 'T_Celltype1s', label = T, repel = T)+
  ggtitle("Intra-tumoral CD8+ T cells")+
  theme(
    title = element_text(color = "black", size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5)
  ) 
Figure2A
dev.copy(png, filename = "Figure2A.png",
         width = 2000, height = 2000, res = 300)
dev.off()

# Figure 2B. CD8+ T cell abundance
dat_Tcell_abundance <- table(Total_T$sample_id, Total_T$T_Celltype1s)
dat_Tcell_abundance <- prop.table(table(Total_T$sample_id, Total_T$T_Celltype1s),1) %>% as.data.frame
colnames(dat_Tcell_abundance) <- c("Batch", "Celltype", "Freq")

main_group=distinct(Total_T@meta.data[,c('sample_id','main_group1')])
dat_Tcell_abundance$main_group<-sapply(dat_Tcell_abundance$Batch, FUN=function(x){
  res<-as.character(main_group$main_group1[match(x, main_group$sample_id)])
  if(res != 'pCR'){
    res<-"Non-pCR"
  }
  return(res)
})
dat_Tcell_abundance$Celltype <- factor(dat_Tcell_abundance$Celltype)
dat_Tcell_abundance$main_group <- factor(dat_Tcell_abundance$main_group)

pval_df <- dat_Tcell_abundance %>%
  group_by(Celltype) %>%
  summarise(
    p_val = wilcox.test(Freq ~ main_group)$p.value/2,
    .groups = "drop"
  ) %>%
  filter(p_val < 0.1) %>%
  mutate(
    label = sprintf("P=%.3f", p_val),
    y.position = tapply(dat_Tcell_abundance$Freq, dat_Tcell_abundance$Celltype, max, na.rm = TRUE)[Celltype] + 0.05
  )

Figure2B<-ggplot(dat_Tcell_abundance %>% 
                   dplyr::filter(str_detect(as.character(Celltype), "CD8")), 
                 aes(x = Celltype, y = Freq, fill = main_group)) +
  geom_boxplot(width = 0.5, color = "black", outlier.shape = NA,
               position = position_dodge(width = 0.75)) +
  geom_text(data = pval_df %>% 
              dplyr::filter(str_detect(as.character(Celltype), "CD8")),
            aes(x = Celltype, y = y.position, label = label),
            inherit.aes = FALSE, vjust = 0) +
  ylab("Relative frequency") +
  xlab("") +
  theme_classic()+
  ggtitle("Relative abundance of intra-tumoral CD8+ T-cells")+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size = 9, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black", size = 11, face = "bold"),
    legend.title = element_text(color = "black", size = 12, face = "bold"),
    legend.text = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    title = element_text(color = "black", size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5)
  )+
  scale_fill_manual(values = c("Non-pCR" = "#1B63A5", 
                               "pCR"    = "#D90016"))+
  labs(fill="")

Figure2B
dev.copy(png, filename = "Figure2B.png",
         width = 2000, height = 1400, res = 300)
dev.off()

# Figure 2C-D. Differentially expressed genes
# P-value plot was obtained from ClueGO enrichment analysis
# TEM Differential Gene Expression
CD8Tem <- subset(Total_T, T_Celltype1s == 'CD8 Tem')
CD8Tem@meta.data <- droplevels(CD8Tem@meta.data)

Idents(CD8Tem) <- 'pcR'
markers.CD8Tem <- FindMarkers(CD8Tem,
                                ident.1 = "Yes",
                                ident.2 = "No",
                                min.pct = 0.25,
                                logfc.threshold = 0.25, verbose = FALSE)
#write.csv(markers.CD8Tem, file = 'Output/CD8Tem_pCR_VS_NonpCR.csv')

# TPEX Differential Gene Expression
CD8Tpex <- subset(Total_T, T_Celltype1s == 'CD8 Tpex')
CD8Tpex@meta.data <- droplevels(CD8Tpex@meta.data)
Idents(CD8Tpex) <- 'pcR'

markers.CD8Tpex <- FindMarkers(CD8Tpex, 
                               ident.1 = "Yes",
                               ident.2 = "No",
                               min.pct = 0.25,
                               logfc.threshold = 0.25, verbose = FALSE)
#write.csv(markers.CD8Tem, file = 'Output/CD8Tem_pCR_VS_NonpCR.csv')

# Figure 2E, CD8+ T cells
target_celltype <- as.character(unique(Total_T$T_Celltype1s))
target_celltype <- target_celltype[grep("CD8", target_celltype)]
combined.TCR<-combineTCR(contig_list, samples=names(contig_list),removeNA = FALSE,
                         removeMulti = FALSE, 
                         filterMulti = T)


mtdata_origin= Total_T@meta.data
for(x in 1:length(combined.TCR)){
  combined.TCR[[x]]$cell<-sapply(combined.TCR[[x]]$barcode, FUN=function(barcode){
    unlist(strsplit(barcode, "_"))[2]
  })
  
  combined.TCR[[x]]$celltype<-sapply(combined.TCR[[x]]$cell, FUN=function(barcode){
    Total_T@meta.data$T_Celltype1s[match(barcode, rownames(Total_T@meta.data))]
  })
  
  combined.TCR[[x]]<-combined.TCR[[x]] %>% 
    dplyr::filter(cell %in% tcell_barcode) %>% 
    dplyr::filter(celltype %in% target_celltype)
  
  combined.TCR[[x]]$celltype<-factor(as.character(combined.TCR[[x]]$celltype),
                                     levels = target_celltype)
  names(combined.TCR)[x]<-as.character(mtdata_origin$sample_id[match(combined.TCR[[x]]$cell, rownames(mtdata_origin))]) %>% unique
}

# Adding clinical outcome variables
## response_3: Non-MPR, MPR, pCR
response_3<-sapply(names(combined.TCR), FUN=function(x){
  id_tcr=Total_T@meta.data %>% dplyr::select(sample_id, sample_id_TCR, main_group1) %>% unique
  return(id_tcr$main_group1[match(x, id_tcr$sample_id)])
})# %>% as.character

## response_2: Non-pCR, pCR
response_2<-sapply(response_3, FUN=function(x){
  if(x=="pCR"){
    return("pCR")
  }else{
    return("Non-pCR")
  }
})

combined.TCR<-addVariable(combined.TCR,
                          variable.name = "response_2",
                          variables=response_2)

cloneCall="gene";chain="both"

TCR_diversity<-clonalDiversity(combined.TCR, 
                               cloneCall = cloneCall,chain = chain, 
                               n.boots=5000, exportTable = T)
TCR_diversity$Response_2<-sm$response_2

Figure2Ea<-ggplot(TCR_diversity, aes(x = Response_2, y = chao1
                                     ,fill=Response_2)) +
  geom_boxplot() +
  geom_point(alpha = 0.6, color = "black") +
  labs(x = "", y = "Diversity")+
  stat_compare_means(
    method = "wilcox.test",
    method.args = list(alternative = "greater"),
    hide.ns = TRUE,   
    label = "p.format",
    label.x = 1.5,
    label.y = 280
  )+
  theme_classic()+ labs(fill = "")+
  scale_fill_manual(values = c("Non-pCR" = "#1B63A5", 
                               "pCR"    = "#D90016"))

Figure2Eb<-ggplot(TCR_diversity, aes(x = Response_2, y = 1-norm.entropy, 
                                     fill=Response_2)) +
  geom_boxplot() +
  geom_point(alpha = 0.6, color = "black") +
  labs(x = "", y = "Clonality")+
  stat_compare_means(
    method = "wilcox.test",
    method.args = list(alternative = "less"),
    hide.ns = TRUE,   
    label = "p.format",
    label.x = 1.5,
    label.y = 0.5
  )+
  theme_classic()+
  labs(fill = "")+
  scale_fill_manual(values = c("Non-pCR" = "#1B63A5", 
                               "pCR"    = "#D90016"))

Figure2Ea+Figure2Eb

dev.copy(png, filename = "Figure2Ea.png",
         width = 2300, height = 1200, res = 300)
dev.off()

TCR_diversity_cell<-do.call(rbind, lapply(target_celltype, FUN=function(c){
  combined.TCR.filter<-lapply(combined.TCR, FUN=function(x){
    x %>% dplyr::filter(celltype == c)
  })
  
  remove_ind<-which(sapply(combined.TCR.filter, nrow)<3)
  if(length(remove_ind)>0){
    TCR_diversity<-clonalDiversity(combined.TCR.filter[-remove_ind], 
                                   cloneCall = "gene",chain = chain, 
                                   n.boots=5000, exportTable = T)
    TCR_diversity$Response_2<-sm$response_2[-remove_ind]  
  }else{
    TCR_diversity<-clonalDiversity(combined.TCR.filter, 
                                   cloneCall = "gene",chain = chain,
                                   n.boots=5000, exportTable = T)
    TCR_diversity$Response_2<-sm$response_2
  }
  
  TCR_diversity$Celltype=c
  return(TCR_diversity)
}))

Figure2Ec<-ggplot(TCR_diversity_cell, 
                  aes(x = Celltype, y = 1-norm.entropy, fill = Response_2)) +
  geom_boxplot() +
  labs(x = "", y = "Clonality") +
  stat_compare_means(
    method = "wilcox.test",
    method.args = list(alternative = "less"),
    hide.ns = TRUE,   
    label = "p.format",
    label.y = 0.52
  ) +
  scale_fill_manual(values = c("Non-pCR" = "#1B63A5", 
                               "pCR"    = "#D90016")) +
  theme_classic()+
  labs(fill = "")
Figure2Ec

dev.copy(png, filename = "Figure2Eb.png",
         width = 2400, height = 1200, res = 300)
dev.off()
