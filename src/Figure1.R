set.seed(501)

# Figure 1B. UMAP visualization of the distributions of 8 major cell types from all single-cell transcriptomic profiles. 
Figure1B<-DimPlot(Total, group.by = 'Main_celltype1s', label = T, repel = T)+
  ggtitle("Intra-tumoral cell types")+
  theme(
    title = element_text(color = "black", size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5)
  ) 
Figure1B

dev.copy(png, filename = "Figure1B.png",
         width = 2000, height = 2000, res = 300)
dev.off()

# Figure 1C. Boxplots showing the relative proportions of major cell types in tumors
dat <- prop.table(table(Total$sample_id, Total$Main_celltype1s),1)
dat <- as.data.frame(dat)
colnames(dat) <- c("sample_id", "Celltype", "Freq")

dat$main_group <- "none"
dat$main_group[dat$sample_id %in% patients_pCR] <- 'pCR'
dat$main_group[dat$sample_id %in% patients_nonpCR] <- 'Non-pCR'
dat$main_group <- factor(dat$main_group, levels = c('Non-pCR', 'pCR'))

Figure1C<-ggplot(dat, aes(x = Celltype, y = Freq, fill = main_group)) + 
  geom_boxplot(width = 0.5, color = "black", outlier.shape = NA, position = position_dodge(width = 0.75)) + 
  ylab("Relative frequency") + 
  xlab("") +
  theme_classic() +
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
  ) +
  labs(fill = "")+
  stat_compare_means(
    method = "wilcox.test",
    hide.ns = TRUE,
    label = "p.format",
    label.y = 1
  ) +
  scale_fill_manual(values = c("Non-pCR" = "#1B63A5", 
                                  "pCR"    = "#D90016"))+
  ggtitle("Relative abundance of intra-tumoral cell types")

Figure1C
dev.copy(png, filename = "Figure1C.png",
         width = 2400, height = 1400, res = 300)
dev.off()

# Figure 1D. TCR repertoire analysis of the entire T cells
target_celltype <- as.character(unique(Total_T$T_Celltype1s))
combined.TCR<-combineTCR(contig_list, samples=names(contig_list),removeNA = FALSE,
                         removeMulti = FALSE, 
                         filterMulti = T)

mtdata_origin= Total_T@meta.data
for(x in 1:length(combined.TCR)){
  combined.TCR[[x]]$cell<-sapply(combined.TCR[[x]]$barcode, FUN=function(barcode){
    unlist(strsplit(barcode, "_"))[2]
  })
  
  combined.TCR[[x]]$celltype<-sapply(combined.TCR[[x]]$cell, FUN=function(barcode){
    mtdata_origin$T_Celltype1s[match(barcode, rownames(Total_T@meta.data))]
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
  id_tcr=mtdata_origin %>% dplyr::select(sample_id, sample_id_TCR, main_group1) %>% unique
  return(id_tcr$main_group1[match(x, id_tcr$sample_id)])
})

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

set.seed(43)
TCR_diversity<-clonalDiversity(combined.TCR, 
                               cloneCall = cloneCall,chain = chain,
                               exportTable = T)
TCR_diversity$Response_2<-sm$response_2

Figure1Da<-ggplot(TCR_diversity, aes(x = Response_2, y = 1-norm.entropy,
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
    label.y = 0.45
  )+
  theme_classic()+
  ggtitle("TCR clonality")+
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(color = "black", size = 14, face = "bold"),)+
  scale_fill_manual(values = c("Non-pCR" = "#1B63A5", 
                               "pCR"    = "#D90016"))+
  labs(fill="")


Figure1Db<-clonalHomeostasis(combined.TCR, 
                             cloneCall = cloneCall, chain=chain, group.by = "response_2")+
  ggtitle("TCR repertoires")+
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(color = "black", size = 14, face = "bold"),)


Figure1D<-Figure1Da+Figure1Db
Figure1D

dev.copy(png, filename = "Figure1D.png",
         width = 2300, height = 1300, res = 300)
dev.off()

# Figure 1E. Spatial transcriptomic maps in representative tumor sections.
use_color <- c('T_NK cells'= '#FF7F00',
               'Epithelial_Aneuploid'= 'red',
               'Epithelial_Diploid' = 'deepskyblue2',
               'Fibroblast cells'='#4DAF4A', 
               'Myeloid cells'= '#984EA3', 
               'Endothelial cells'= '#AEDFAD',
               'Mast cells'='#FFFF33',
               'B_Plasma cells' = '#F781BF')

slide.seq_PCR <- AddMetaData(slide.seq_PCR, metadata = RCTD_PCR@results$results_df)
slide.seq_NonPCR <- AddMetaData(slide.seq_NonPCR, metadata = RCTD_NonPCR@results$results_df)
PCR <- SpatialDimPlot(slide.seq_PCR, group.by = "second_type",
                      images = NULL, image.alpha = 0, cols = use_color, pt.size.factor = 3000)
ggsave("spatial_plot_PCR_umap.png", plot = PCR, width = 1.8, height = 1.8, dpi = 300)
img <- image_read("spatial_plot_PCR_umap.png")
Figure1Ea <- image_rotate(img, -90)  

Non_PCR <- SpatialDimPlot(slide.seq_NonPCR, group.by = "second_type",  images = NULL,  image.alpha = 0, cols = use_color, pt.size.factor = 3000)
ggsave("spatial_plot_NonPCR_umap.png", plot = Non_PCR, width = 1.8, height = 1.8, dpi = 300)
img <- image_read("spatial_plot_NonPCR_umap.png")
Figure1Eb <- image_rotate(img, -90)  
