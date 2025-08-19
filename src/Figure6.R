# Figure 6A-B. CellChat signaling network analysis
# Generating CellChat objects
## Subsetting target cell types
CD4_subset <- subset(CD4T, subset = New %in% c("CD4_Tem", "CD4_Treg"))
CD8_subset <- subset(CD8T, subset = T_Celltype1s %in% c("CD8 Tem", "CD8 Tpex", "CD8 Temra"))
myeloid_subset <- subset(myeloid, subset = M_Celltype1 %in% 
                           c("FABP4+ Macrophage", "FOLR2+ Macrophage", "CD16+ Monocyte", "OLR1+ Monocyte"))
apCAF<-subset(CAF, subset = simple %in% c("Antigen-presenting CAF"))
B_Plasma <- subset(Total, subset = Main_celltype1 == "B_Plasma cells")

B_Plasma$cellchat_group<-as.character(B_Plasma$Main_celltype1)
apCAF$cellchat_group<-as.character(apCAF$simple)
CD4_subset$cellchat_group<-as.character(CD4_subset$New)
CD8_subset$cellchat_group<-as.character(CD8_subset$T_Celltype1s)
myeloid_subset$cellchat_group<-as.character(myeloid_subset$M_Celltype1)

B_Plasma <- RenameCells(B_Plasma, add.cell.id = "BPlasma")
apCAF <- RenameCells(apCAF, add.cell.id = "apCAF")
CD4_subset <- RenameCells(CD4_subset, add.cell.id = "CD4")
CD8_subset <- RenameCells(CD8_subset, add.cell.id = "CD8")
myeloid_subset <- RenameCells(myeloid_subset, add.cell.id = "Myeloid")

obj <- merge(
  B_Plasma,
  y = list(apCAF, CD4_subset, CD8_subset, myeloid_subset),
  add.cell.ids = NULL
)
obj_pcr<-subset(obj, pcR == "Yes")
obj_nonpcr<-subset(obj, pcR != "Yes")
obj_all<-obj
rm(obj)

# Designating cellchat DB
CellChatDB.use<-CellChatDB.human # useCellChatDB.mous eif running on mouse data

## pCR cellchat
cellchat_pcr <- createCellChat(object = obj_pcr, group.by = "cellchat_group", assay = "RNA")
cellchat_pcr<-setIdent(cellchat_pcr,ident.use= "cellchat_group") # set "labels" as default cell identity
cellchat_pcr@DB<-CellChatDB.use
cellchat_pcr<-subsetData(cellchat_pcr) # This step is necessary even if using the whole database
cellchat_pcr<-identifyOverExpressedGenes(cellchat_pcr)
cellchat_pcr<-identifyOverExpressedInteractions(cellchat_pcr)
cellchat_pcr<-smoothData(cellchat_pcr, adj = PPI.human)
cellchat_pcr<-computeCommunProb(cellchat_pcr)
cellchat_pcr<-filterCommunication(cellchat_pcr, min.cells= 10)
cellchat_pcr<-computeCommunProbPathway(cellchat_pcr)
cellchat_pcr<-aggregateNet(cellchat_pcr)

## Non-pCR
cellchat_nonpcr <- createCellChat(object = obj_nonpcr, group.by = "cellchat_group", assay = "RNA")
cellchat_nonpcr<-setIdent(cellchat_nonpcr,ident.use= "cellchat_group") # set "labels" as default cell identity
cellchat_nonpcr@DB<-CellChatDB.use
cellchat_nonpcr<-subsetData(cellchat_nonpcr) # This step is necessary even if using the whole database
cellchat_nonpcr<-identifyOverExpressedGenes(cellchat_nonpcr)
cellchat_nonpcr<-identifyOverExpressedInteractions(cellchat_nonpcr)
cellchat_nonpcr<-smoothData(cellchat_nonpcr, adj = PPI.human)
cellchat_nonpcr<-computeCommunProb(cellchat_nonpcr)
cellchat_nonpcr<-filterCommunication(cellchat_nonpcr, min.cells= 10)
cellchat_nonpcr<-computeCommunProbPathway(cellchat_nonpcr)
cellchat_nonpcr<-aggregateNet(cellchat_nonpcr)

# Cell-cell signaling network analysis
cellchat_pcr <- netAnalysis_computeCentrality(cellchat_pcr)
cellchat_nonpcr <- netAnalysis_computeCentrality(cellchat_nonpcr)
target_pathway=c("MHC-I", "MHC-II", "MIF", "FN1", "GALECTIN", "COLLAGEN", "APP", "CD99", 
                 "ICAM", "CypA", "CD45", "ADGRE", "ANNEXIN", "ApoE", "TGFb", "SPP1", "LCK",
                 "CXCL", "PECAM1", "Prostaglandin", "CCL", "SELPLG")
p1<-netAnalysis_signalingRole_heatmap(cellchat_pcr, signaling = target_pathway, pattern = "outgoing",
                                      title="pCR", width=5)
p2<-netAnalysis_signalingRole_heatmap(cellchat_nonpcr,signaling = target_pathway, pattern = "outgoing",
                                      title="non-pCR", width=5)
Figure6A<-p1+p2

p3<-netAnalysis_signalingRole_heatmap(cellchat_pcr, signaling = target_pathway, pattern = "incoming",
                                      title="pCR", width=5)
p4<-netAnalysis_signalingRole_heatmap(cellchat_nonpcr,signaling = target_pathway, pattern = "incoming",
                                      title="non-pCR", width=5)
Figure6B<-p3+p4

### Figure 6C. CCL bubble plot
pairLR.use <- extractEnrichedLR(cellchat_nonpcr, signaling = c("CCL"))
source=sort(unique(obj_pcr$cellchat_group))[-c(2,4,5, 6, 7, 8)]
target=sort(unique(obj_pcr$cellchat_group))[-c(1,3,4,9, 10, 11)]
p5=netVisual_bubble(cellchat_pcr, sources.use=source, targets.use = target,
                    pairLR.use = pairLR.use, remove.isolate = TRUE, title.name = "pCR")
p6=netVisual_bubble(cellchat_nonpcr, sources.use=source, targets.use = target,
                    pairLR.use = pairLR.use, remove.isolate = TRUE, title.name= "non-pCR")

Figure6C<-p5+p6

##Figure6D. Information flow summary chart
df_pCR<-subsetCommunication(cellchat_pcr, 
                            sources.use=source,
                            targets.use=target,
                            signaling="CCL") %>% mutate(group="pCR")
df_nonpCR<-subsetCommunication(cellchat_nonpcr,
                               sources.use=source,
                               targets.use=target,
                               signaling="CCL") %>% mutate(group="NonpCR")
df_combined<-rbind(df_pCR, df_nonpCR)

info_flow_summary <- df_combined %>%
  group_by(interaction_name, group) %>%
  summarise(sum_prob = sum(prob), .groups = "drop")

info_flow_wide <- info_flow_summary %>%
  tidyr::pivot_wider(names_from = group, values_from = sum_prob, values_fill = 0)

info_flow_wide <- info_flow_wide %>%
  mutate(
    total = pCR + NonpCR,
    pCR_ratio = ifelse(total > 0, pCR / total, 0),
    NonpCR_ratio = ifelse(total > 0, NonpCR / total, 0)
  )

info_flow_plot <- info_flow_wide %>%
  select(interaction_name, pCR_ratio, NonpCR_ratio) %>%
  tidyr::pivot_longer(cols = c(pCR_ratio, NonpCR_ratio),
                      names_to = "group", values_to = "ratio")

info_flow_plot$group <- factor(info_flow_plot$group,
                               levels = c("NonpCR_ratio", "pCR_ratio"),
                               labels = c("Non_pCR", "pCR"))

pval_table <- df_combined %>%
  group_by(interaction_name) %>%
  summarise(
    p_value = if (length(unique(group)) == 2 &&
                  sum(group == "pCR") > 0 &&
                  sum(group == "NonpCR") > 0) {
      wilcox.test(prob[group == 'pCR'], prob[group == 'NonpCR'])$p.value
    } else {
      0
    },
    .groups = "drop"
  ) %>% arrange(p_value)

ordered_levels <- rev(as.character(pval_table$interaction_name))
info_flow_plot$interaction_name <- factor(info_flow_plot$interaction_name, levels = ordered_levels)

pval_labels <- pval_table %>%
  mutate(y_pos = 1.03)

Figure6D<-ggplot(info_flow_plot, aes(x = interaction_name, y = ratio, fill = group)) +
  geom_bar(stat = "identity") +
  geom_text(
    data = pval_labels,
    aes(x = interaction_name, y = y_pos, label = sprintf("%.3g", p_value)),
    inherit.aes = FALSE,
    size = 3,
    hjust = 0
  ) +
  annotate("text", 
           x = length(unique(info_flow_plot$interaction_name)) + 0.5,
           y = 1.01, 
           label = "p-value", 
           hjust = 0,
           vjust = -0.5,
           size = 4,
           fontface = "bold") +
  coord_flip(clip = "off") +
  labs(title = paste0("CCL", " signaling - Normalized info flow"),
       y = "Relative contribution (0–1 scale)",
       x = "Ligand -> Receptor") +
  scale_fill_manual(values = c("Non_pCR" = "#2b3f8c", "pCR" = "#e31a1c")) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.2)),
    labels = function(x) ifelse(x > 1, "", x)
  ) +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.margin = margin(5.5, 60, 5.5, 5.5, "pt")
  )

# Figure6E. Dotplot of selected genes and cell types
obj_all$cellchat_group[which(obj_all$cellchat_group=="Antigen-presenting CAF")]<-"apCAF"
obj_all$pcR<-ifelse(obj_all$pcR=="Yes", "pCR", "Non-pCR")
subsetDotPlot<-function(obj, idents_name, idents, features, title=NULL){
  Idents(obj)<-idents_name
  obj_subset=subset(obj, idents=idents)
  p<-DotPlot(obj_subset, group.by = idents_name,  features = features)+
    labs(title = title,
         color = "Average expression",
         size = "Percent expressed")+
    xlab("")+ylab("")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5))+ 
    guides(
      color = guide_colorbar(order = 1),
      size = guide_legend(order = 2)
    )
  return(p)
}

p9=subsetDotPlot(obj=obj_all, idents_name="cellchat_group",
                 idents=c("OLR1+ Monocyte", "FOLR2+ Macrophage", "FABP4+ Macrophage", "CD16+ Monocyte","apCAF"),
                 features=c("CCL20", "CCL5"))
p10=subsetDotPlot(obj=obj_all, idents_name="cellchat_group",
                  idents=c("CD8 Tpex", "CD8 Temra", "CD8 Tem", "CD4_Treg", "CD4_Tem", "B_Plasma cells"),
                  features=c("CCR4", "CCR5", "CCR6"))

Figure6E<-p9+p10

# Figure 6F. Spatial Cell-cell interaction plot
use_color <- c("CD8 Tem"= "#AEDFAD",
               "CD8 Temra"= "#C4C8F7" , 
               "CD8 Tpex" = "red",
               "CD4 Tem"="#DCAE5A", 
               "CD4 Treg"= "#F9F961", 
               'FABP4+ Macrophage' = "#984EA3",
               'FOLR2+ Macrophage' ="deepskyblue2",
               'OLR1+ Monocyte' = "#FF7F00",
               'CD16+ Monocyte' = "#39B600",
               'Antigen-presenting CAF' =  "#F781BF")

slide.seq_PCR <- AddMetaData(slide.seq_PCR, metadata = RCTD_PCR_merge@results$results_df)
slide.seq_NonPCR <- AddMetaData(slide.seq_NonPCR, metadata = RCTD_NonPCR_merge@results$results_df)

PCR <- SpatialDimPlot(slide.seq_PCR, group.by = "second_type",  images = NULL, image.alpha = 0, cols = use_color, pt.size.factor = 3000)
ggsave("spatial_plot_PCR_Merge_umap.png", plot = PCR, width = 2.8, height = 2.8, dpi = 300)
img <- image_read("spatial_plot_PCR_Merge_umap.png")
Figure6Fa <- image_rotate(img, -90) 

Non_PCR <- SpatialDimPlot(slide.seq_NonPCR, group.by = "second_type",  images = NULL,  image.alpha = 0, cols = use_color, pt.size.factor = 3000)
ggsave("spatial_plot_NonPCR_Merge_umap.png", plot = Non_PCR, width = 2.8, height = 2.8, dpi = 300)
img <- image_read("spatial_plot_NonPCR_Merge_umap.png")
Figure6Fb <- image_rotate(img, -90)  

# Spatial cellchat interaction plot 
Merge_pCR_Spatial <- readRDS('Data/Merge_pCR_SpatialCellchat.rds')
Merge_NonpCR_Spatial <- readRDS('Data/Merge_NonpCR_SpatialCellchat.rds')

cell_colors1 <- c('CD8 Tem' = "#AEDFAD",
                  'CD8 Temra'= "#C4C8F7" , 
                  'CD8 Tpex' = "red",
                  'CD4 Tem' = "#DCAE5A", 
                  'CD4 Treg' = "#F9F961", 
                  'FABP4+ Macrophage' = "#984EA3",
                  'FOLR2+ Macrophage' = "deepskyblue2",
                  'OLR1+ Monocyte' = "#FF7F00",
                  'CD16+ Monocyte' = "#39B600",
                  'Antigen-presenting CAF' = "#F781BF")

cells_PCR1 <- rownames(Merge_pCR_Spatial@net$count)
color_use_PCR1 <- cell_colors1[match(cells_PCR1, names(cell_colors1))]

cells_NonPCR1 <- rownames(Merge_NonpCR_Spatial@net$count)
color_use_NonPCR1 <- cell_colors1[match(cells_NonPCR1, names(cell_colors1))]

Merge_pCR_Spatial1 <- netAnalysis_computeCentrality(Merge_pCR_Spatial, slot.name = "netP") 
Merge_NonpCR_Spatial1 <- netAnalysis_computeCentrality(Merge_NonpCR_Spatial, slot.name = "netP")

pathways.show <- c('CCL') 

png("spatial_plot_Merge_pCR_CCL.png", width = 800, height = 800, res = 90)
netVisual_aggregate(Merge_pCR_Spatial1, signaling = pathways.show, color.use = color_use_PCR1, layout = "spatial", edge.width.max = 1, alpha.image = 0.05, vertex.weight = "incoming", vertex.size.max = 4, vertex.label.cex = 3)
dev.off()
img_pCR <- image_read("spatial_plot_Merge_pCR_CCL.png")
Figure6Fc <- image_rotate(img_pCR, 270) 

png("spatial_plot_Merge_NonpCR_CCL.png", width = 800, height = 800, res = 90)
netVisual_aggregate(Merge_NonpCR_Spatial1, signaling = pathways.show, color.use = color_use_NonPCR1, layout = "spatial", edge.width.max = 1, alpha.image = 0.05, vertex.weight = "incoming", vertex.size.max = 4, vertex.label.cex = 3)
dev.off()
img_NonpCR <- image_read("spatial_plot_Merge_NonpCR_CCL.png")
Figure6Fd <- image_rotate(img_NonpCR, 270) 
