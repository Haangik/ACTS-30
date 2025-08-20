set.seed(501)

# Figure 4A. UMAP visualization of intra-tumoral CD4+ T cells
Figure4A<-DimPlot(CD4T, group.by = 'New', label = T, repel = T)+
  ggtitle("Intra-tumoral CD4+ T cells")+
  theme(
    title = element_text(color = "black", size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5)
  ) 

# Figure 4B. TCR analysis of CD4+ T cells
target_celltype <- as.character(unique(Total_T$T_Celltype1s))
target_celltype <- target_celltype[grep("CD4", target_celltype)]
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
                               cloneCall = cloneCall,
                               chain = chain,
                               n.boots=5000,
                               exportTable = T)
TCR_diversity$Response_2<-sm$response_2

Figure4Ba<-ggplot(TCR_diversity, aes(x = Response_2, y = 1-norm.entropy,
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
    label.y = 0.15
  )+
  theme_classic()+
  ggtitle("Per pathological response")+
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(color = "black", size = 12, face = "bold"))+
  scale_fill_manual(values = c("Non-pCR" = "#1B63A5", 
                               "pCR"    = "#D90016"))+
  labs(fill="")

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

Figure4Bb<-ggplot(TCR_diversity_cell, aes(x = Celltype, y = 1-norm.entropy,
                                          fill=Response_2)) +
  geom_boxplot() +
  labs(x = "", y = "Clonality") +
  stat_compare_means(method="wilcox.test", 
                     method.args = list(alternative = "less"),
                     hide.ns = TRUE,
                     label = "p.format",
                     label.y = 0.2)+
  theme_classic()+
  scale_fill_manual(values = c("Non-pCR" = "#1B63A5", 
                               "pCR"    = "#D90016"))+
  ggtitle("Per CD4 T-cell phenotype and stratified for response")+
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(color = "black", size = 12, face = "bold"))+
  labs(fill="")

Figure4Ba
dev.copy(png, filename = "Figure4Ba.png",
         width = 1200, height = 1200, res = 300)
dev.off()

Figure4Bb
dev.copy(png, filename = "Figure4Bb.png",
         width = 2400, height = 1200, res = 300)
dev.off()

# Figure 4C. Differentially expressed genes
# P-value plot was obtained from ClueGO enrichment analysis
CD4Tem <- subset(CD4T, New == 'CD4_Tem')
CD4Tem@meta.data <- droplevels(CD4Tem@meta.data)
Idents(CD4Tem) <- 'pcR'
markers.CD4Tem <- FindMarkers(CD4Tem, 
                                ident.1 = "Yes",
                                ident.2 = "No", 
                                min.pct = 0.25, 
                                logfc.threshold = 0.25, verbose = FALSE)
#write.csv(markers.CD4Tem, file = 'Output/CD4Tem_pCR_VS_NonpCR.csv')

# Figure 4D. Pseudotime plot 
cds <- as.cell_data_set(CD4T)
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, num_dim = 100)
cds <- align_cds(cds, alignment_group = "sample_id")
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)
cds@clusters$UMAP$partitions <- reacreate.partition
list_cluster <- CD4T@active.ident
cds@clusters$UMAP$clusters <- list_cluster
cds@int_colData@listData$reducedDims$UMAP <- CD4T@reductions$umap@cell.embeddings
cds <- learn_graph(cds, use_partition = F)
cds <- order_cells(cds)

Figure4D<-plot_cells(cds,
           color_cells_by = 'pseudotime',
           cell_size=0.15,
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = T,
           label_leaves = FALSE,
           group_label_size = 0)

# Figure 4E. Density plot 
new_dats <- data.frame(
  X = reducedDims(cds)$UMAP [, 1],  
  Y = reducedDims(cds)$UMAP [, 2],  
  Response= colData(cds)$pcR
)

Figure4Ea<-ggplot(new_dats, aes(x = X, y = Y)) + 
  geom_point(color = "grey", size = 0.5) + theme_bw(base_size = 14) + 
  geom_pointdensity(data = new_dats[new_dats$Response == "Yes",], adjust = 0.7) + 
  scale_color_viridis() + theme(panel.border = element_rect(color = "black", linewidth = 1),
                                panel.grid.major = element_blank(),  
                                panel.grid.minor = element_blank(),
                                axis.text = element_text(size = 15,face = "bold", color = "black"))

Figure4Eb<-ggplot(new_dats, aes(x = X, y = Y)) + 
  geom_point(color = "grey", size = 0.5) + theme_bw(base_size = 14) + 
  geom_pointdensity(data = new_dats[new_dats$Response == "No",], adjust = 0.8) + 
  scale_color_viridis() + theme(panel.border = element_rect(color = "black", linewidth = 1),
                                panel.grid.major = element_blank(),  
                                panel.grid.minor = element_blank(),
                                axis.text = element_text(size = 15,face = "bold", color = "black"))

Figure4Ea+Figure4Eb

# Figure 4F. CD4+ T cell abundance
# T cell abundance dataset was generated from Figure 2B
Figure4F<-ggplot(dat_Tcell_abundance %>% 
                   filter(str_detect(as.character(Celltype), "CD4")), 
                 aes(x = Celltype, y = Freq, fill = main_group)) +
  geom_boxplot(width = 0.5, color = "black", outlier.shape = NA,
               position = position_dodge(width = 0.75)) +
  geom_text(data = pval_df%>% 
              filter(str_detect(as.character(Celltype), "CD4")), aes(x = Celltype, y = y.position, label = label),
            inherit.aes = FALSE, vjust = 0) +
  ylab("Relative frequency") +
  xlab("") +
  theme_classic()+
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
  ggtitle("Relative abundance of intra-tumoral CD4+ T-cells")+
  scale_fill_manual(values = c("Non-pCR" = "#1B63A5", 
                               "pCR"    = "#D90016"))

Figure4F
dev.copy(png, filename = "Figure4F.png",
         width = 2000, height = 1400, res = 300)
dev.off()

# Figure 4G. Spatial transcriptomic maps of CD4+ T cells
use_color <- c(
  'CD4 Treg' = "red",
  'CD4 CXCL13' = "#FF7F00",
  'CD4 Tn' = "#C4C8F7",
  'CD4 Tcm' = "#39B600",
  'CD4 IL17' = "#F781BF",
  'CD4 Tem' = "royalblue")

slide.seq_PCR <- AddMetaData(slide.seq_PCR, metadata = RCTD_PCR_CD4@results$results_df)
slide.seq_NonPCR <- AddMetaData(slide.seq_NonPCR, metadata = RCTD_NonPCR_CD4@results$results_df)

PCR <- SpatialDimPlot(slide.seq_PCR, group.by = "first_type",  images = NULL, image.alpha = 0, cols = use_color, pt.size.factor = 3000)
ggsave("spatial_plot_PCR_CD4T_umap.png", plot = PCR, width = 2.8, height = 2.8, dpi = 300)
img <- image_read("spatial_plot_PCR_CD4T_umap.png")
Figure4Ga <- image_rotate(img, -90) 

Non_PCR <- SpatialDimPlot(slide.seq_NonPCR, group.by = "first_type",  images = NULL,  image.alpha = 0, cols = use_color, pt.size.factor = 3000)
ggsave("spatial_plot_NonPCR_CD4T_umap.png", plot = Non_PCR, width = 2.8, height = 2.8, dpi = 300)
img <- image_read("spatial_plot_NonPCR_CD4T_umap.png")
Figure4Gb <- image_rotate(img, -90)  

# Figure 4H
Treg <- subset(Total_T, T_Celltype1s == 'CD4_Treg')
exhaustion_genes <- c("LAG3", "TIGIT", "PDCD1", "HAVCR2",
                      "CTLA4", "LAYN", "ENTPD1", "TNFRSF1B", "TNFRSF8",
                      "IL10", "TNFRSF18", "TNFRSF9", "TNFRSF4")
Treg <- AddModuleScore(object = Treg, features = list(exhaustion_genes), name = "ExhaustionScore")
Idents(Treg) <- 'pcR1'
Figure4H <- VlnPlot(ACTS_Treg, features = "ExhaustionScore1", group.by = "pcR1", pt.size = 0) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.format",
                     label.x = 1.5) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.6, colour = "black") +
  theme_minimal() +
  scale_fill_manual(values = c("Non-pCR" = "#1B63A5", "pCR" = "#D90016")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(face = "bold", colour = "black", size = 10), 
        axis.text.y = element_text(face = "bold", colour = "black", size = 10),
        title = element_text(color = "black", size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5)
  )+
  ggtitle("Exhausted signature (Treg)")
Figure4H

dev.copy(png, filename = "Figure4H.png",
         width = 1500, height = 1300, res = 300)
dev.off()
