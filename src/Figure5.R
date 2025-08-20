# Figure 5A. UMAP visualization of intra-tumoral myeloid cells
Figure5A<-DimPlot(myeloid, group.by = 'M_Celltype1', label = T, repel = T)
Figure5A


# Figure 5B. Boxplots showing the relative proportions of myeloid cell subtypes in tumors
dat <- prop.table(table(myeloid$sample_id, myeloid$M_Celltype1),1)
dat <- as.data.frame(dat)
colnames(dat) <- c("sample_id", "Celltype", "Freq")

dat$main_group <- "none"
dat$main_group[dat$sample_id %in% patients_pCR] <- 'pCR'
dat$main_group[dat$sample_id %in% patients_nonpCR] <- 'Non-pCR'
dat$main_group <- factor(dat$main_group, levels = c('Non-pCR', 'pCR'))

Figure5B<-ggplot(dat, aes(x = Celltype, y = Freq, fill = main_group)) + 
  geom_boxplot(width = 0.5, color = "black", outlier.shape = NA, position = position_dodge(width = 0.75)) + 
  ylab("Relative frequency") + 
  xlab("") + 
  geom_signif(test = "wilcox.test", comparisons = list(c('Non-pCR', 'pCR')),  map_signif_level = TRUE, 
              step_increase = 0.1,  y_position = c(0,0.7),   color = "black",  aes(group = Celltype)) +
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
    title = element_text(color = "black", size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(fill = "")+
  scale_fill_manual(values = c("Non-pCR" = "#1B63A5", 
                               "pCR"    = "#D90016"))+
  ggtitle("Relative abundance of intra-tumoral myeloid cells")+
  stat_compare_means(
    method = "wilcox.test",
    hide.ns = TRUE,
    label = "p.format",
    label.y = 0.77
  )
Figure5B
dev.copy(png, filename = "Figure5B.png",
         width = 2700, height = 1400, res = 300)
dev.off()

# Figure5C. Myeloid Total DEG 
Idents(myeloid) <- 'pcR'
markers.myeloid<- FindMarkers(myeloid, ident.1 = "Yes",
                                ident.2 = "No",
                                min.pct = 0.25,
                                logfc.threshold = 0.25, verbose = FALSE)

write.csv(markers.myeloid, file = 'output/Myeloid_Total_pCR_VS_NonpCR.csv')

# Figure 5D. Myeloid subtype Defense_response_score, External_stimulus_score plot 
Defense_response_score <- list(c('AKNA', 'APOL1', 'BIRC3', 'CASP8', 'CD48', 'CD74', 'CFD', 'CIITA', 'CST3', 'CST7', 
                                 'CXCL16', 'EIF2AK4', 'FCGR1A', 'FCGR1B', 'FEM1C', 'FES', 'FFAR4', 'FPR2',
                                 'FYN', 'GBP1', 'GBP2', 'GBP4', 'GBP5', 'HDAC9', 'HLA-DPA1', 'HLA-F', 'IFI27', 
                                 'IFIT2', 'IFITM1', 'IFITM2', 'IL15', 'IL16', 'INHBA', 'IRF7', 'IRF8', 'ISG20',
                                 'ITGAL', 'JAK2', 'LGALS3BP', 'LY75', 'LYST', 'MALT1', 'MAPK7', 'MEF2C', 'MPEG1',
                                 'NCF1', 'NLRC5', 'NLRP1', 'NOP53', 'OASL', 'PTK2', 'PTPN22', 'RPL39', 'SERPING1',
                                 'SIGLEC10', 'SLAMF7', 'SP140', 'SPN', 'STAT1', 'TRIM22', 'TRIM38'))

External_stimulus_score <- list(c('APP', 'CALR', 'CCL3', 'CCL4', 'CCR1', 'CD1D', 'CSF1R', 'CXCL8',
                                  'EREG', 'IL1B', 'LGALS1', 'LGMN', 'MAPK13', 'NFKBIA', 'NFKBIZ', 
                                  'NINJ1', 'NLRP3', 'OSM', 'P2RX4', 'PLA2G7', 'PLAU', 'PRKCA',
                                  'PTGS2', 'RNASE1', 'S100A9', 'SERPINE1', 'TNF', 'TNIP1', 'TREM2', 
                                  'VEGFA'))

myeloid$ind <- "FALSE"
myeloid$ind[myeloid$M_Celltype1 %in% c('FABP4+ Macrophage', 'OLR1+ Monocyte', 'FOLR2+ Macrophage', 'CD16+ Monocyte')] <- TRUE
myeloid_subset <- subset(myeloid, subset = ind == TRUE)
myeloid_subset@meta.data <- droplevels(myeloid_subset@meta.data)

myeloid_subset$M_Celltype1 <- factor(myeloid_subset$M_Celltype1, 
                                     levels = c('CD16+ Monocyte', 'FABP4+ Macrophage', 
                                                'OLR1+ Monocyte', 'FOLR2+ Macrophage'))

myeloid_subset <- AddModuleScore(myeloid_subset, features = Defense_response_score, name = 'Defense_response_score')
myeloid_subset <- AddModuleScore(myeloid_subset, features = External_stimulus_score, name = 'External_stimulus_score')

Figure5D<-DotPlot(myeloid_subset,
        features = c('Defense_response_score1', 'External_stimulus_score1'),
        group.by = "M_Celltype1",
        dot.scale = 8, scale.max = 120, scale.min = 25) +
  scale_color_gradient2(low = "gray90", mid = "blue", high = "darkblue", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
Figure5D
dev.copy(png, filename = "Figure5B.png",
         width = 1400, height = 1400, res = 300)
dev.off()

# Figure 5E
use_color <- c(
  'FABP4+ Macrophage' = "royalblue",
  'FOLR2+ Macrophage' ="#E76BF3",
  'SPP1+ Macrophage' = "#984EA3",
  'cDC1' = "#39B600",
  'cDC2' = "#94DFFF",
  'CD14+ Monocyte' = "#F87D74",
  'CD16+ Monocyte' = "#C99100",
  'OLR1+ Monocyte' = "#8ECC7C")

slide.seq_PCR <- AddMetaData(slide.seq_PCR, metadata = RCTD_PCR_myeloid@results$results_df)
slide.seq_NonPCR <- AddMetaData(slide.seq_NonPCR, metadata = RCTD_NonPCR_myeloid@results$results_df)

PCR <- SpatialDimPlot(slide.seq_PCR, group.by = "first_type",  images = NULL, image.alpha = 0, cols = use_color, pt.size.factor = 3000)
ggsave("spatial_plot_PCR_CD8T_umap.png", plot = PCR, width = 2.8, height = 2.8, dpi = 300)
img <- image_read("spatial_plot_PCR_CD8T_umap.png")
Figure5Ea <- image_rotate(img, -90) 

Non_PCR <- SpatialDimPlot(slide.seq_NonPCR, group.by = "first_type",  images = NULL,  image.alpha = 0, cols = use_color, pt.size.factor = 3000)
ggsave("spatial_plot_NonPCR_CD8T_umap.png", plot = Non_PCR, width = 2.8, height = 2.8, dpi = 300)
img <- image_read("spatial_plot_NonPCR_CD8T_umap.png")
Figure5Eb <- image_rotate(img, -90)  

# Figure 5F
CD4T$cellchat_group<-as.character(CD4T$T_Celltype1s)
myeloid$cellchat_group<-as.character(myeloid$M_Celltype1)
CD4_Myeloid<- merge(
  CD4T,
  y = list(myeloid),
  add.cell.ids = NULL
)

obj_NonPCR <- subset(CD4_Myeloid, pcR == 'No')
obj_PCR <- subset(CD4_Myeloid, pcR == 'Yes')

obj_NonPCR@meta.data <- droplevels(obj_NonPCR@meta.data)
obj_PCR@meta.data <- droplevels(obj_PCR@meta.data)


cellchat <- createCellChat(object = obj_NonPCR, group.by = "cellchat_group", assay = "RNA")

cellchat<-setIdent(cellchat,ident.use= "cellchat_group") 
levels(cellchat@idents)
groupSize<-as.numeric(table(cellchat@idents)) 
CellChatDB<-CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB$interaction$annotation
CellChatDB.use=CellChatDB
cellchat@DB<-CellChatDB.use
cellchat<-subsetData(cellchat) 
cellchat<-identifyOverExpressedGenes(cellchat)
cellchat<-identifyOverExpressedInteractions(cellchat)
cellchat<-smoothData(cellchat, adj = PPI.human)
cellchat<-computeCommunProb(cellchat) 
cellchat<-filterCommunication(cellchat,min.cells= 10)
cellchat<-computeCommunProbPathway(cellchat)
cellchat<-aggregateNet(cellchat)
groupSize<-as.numeric(table(cellchat@idents))

cellchat_NonPCR <- cellchat
rm(cellchat)

cellchat <- createCellChat(object = obj_PCR, group.by = "cellchat_group", assay = "RNA")

cellchat<-setIdent(cellchat,ident.use= "cellchat_group") 
levels(cellchat@idents)
groupSize<-as.numeric(table(cellchat@idents)) 
CellChatDB<-CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB$interaction$annotation
CellChatDB.use=CellChatDB
cellchat@DB<-CellChatDB.use
cellchat<-subsetData(cellchat) 
cellchat<-identifyOverExpressedGenes(cellchat)
cellchat<-identifyOverExpressedInteractions(cellchat)
cellchat<-smoothData(cellchat, adj = PPI.human)
cellchat<-computeCommunProb(cellchat) 
cellchat<-filterCommunication(cellchat,min.cells= 10)
cellchat<-computeCommunProbPathway(cellchat)
cellchat<-aggregateNet(cellchat)
groupSize<-as.numeric(table(cellchat@idents))

cellchat_PCR <- cellchat
rm(cellchat)

# Figure 5F. Circos plot 
sourcedata <- c('cDC1', 'cDC2', 'FABP4+ Macrophage', 'FOLR2+ Macrophage', 'CD14+ Monocyte', 'CD16+ Monocyte', 
                'Proliferating Macrophage', 'SPP1+ Macrophage', 'OLR1+ Monocyte')

targetdata <- c('CD4_TN', 'CD4_Tcm', 'CD4_Tem', 'CD4_CXCL13', 'CD4_Treg', 'CD4_IL17')

pathways.show <- c("TNF")
Figure5Fa<-netVisual_aggregate(cellchat_NonPCR, signaling = pathways.show, layout = "circle", edge.width.max = 7,  
                    vertex.label.cex = 0.8, arrow.size = 0.7,, weight.scale = TRUE, 
                    sources.use = sourcedata,  targets.use = targetdata,  vertex.receiver = vertex.receiver)

Figure5Fb<-netVisual_aggregate(cellchat_PCR, signaling = pathways.show, layout = "circle", edge.width.max = 7,
                    vertex.label.cex = 0.8, arrow.size = 0.7, weight.scale = TRUE, sources.use = sourcedata,  targets.use = targetdata,
                    vertex.receiver = vertex.receiver)
Figure5Fa+Figure5Fb

# Figure 5G. Ligand-receptor pairs Score plot 
pathway_name <- "TNF"

TNF_df_PCR <- subsetCommunication(cellchat_PCR, signaling = pathway_name)
TNF_df_NonPCR <- subsetCommunication(cellchat_NonPCR, signaling = pathway_name)

TNF_df_PCR$group <- "pCR"
TNF_df_NonPCR$group <- "non-pCR"

TNF_df_combined <- rbind(TNF_df_PCR, TNF_df_NonPCR)

info_flow_summary <- TNF_df_combined %>%
  group_by(interaction_name, group) %>%
  summarise(sum_prob = sum(prob), .groups = "drop")

info_flow_wide <- info_flow_summary %>%
  tidyr::pivot_wider(names_from = group, values_from = sum_prob, values_fill = 0)

info_flow_wide <- info_flow_wide %>%
  mutate(
    total = `pCR` + `non-pCR`,
    pCR_ratio = ifelse(total > 0, `pCR` / total, 0),
    nonpCR_ratio = ifelse(total > 0, `non-pCR` / total, 0)
  )

info_flow_plot <- info_flow_wide %>%
  select(interaction_name, pCR_ratio, nonpCR_ratio) %>%
  tidyr::pivot_longer(cols = c(pCR_ratio, nonpCR_ratio),
                      names_to = "group", values_to = "ratio")

info_flow_plot$group <- factor(info_flow_plot$group,
                               levels = c("nonpCR_ratio", "pCR_ratio"),
                               labels = c("non-pCR", "pCR"))

Figure5G<-ggplot(info_flow_plot, aes(x = reorder(interaction_name, -ratio), y = ratio, fill = group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste0(pathway_name, " ligand-receptor pairs "),
       y = "Relative contribution (0–1 scale)",
       x = "Ligand -> Receptor") +
  scale_fill_manual(values = c("non-pCR" = "#2b3f8c", "pCR" = "#e31a1c")) +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()
  )

pval_table <- TNF_df_combined %>%
  group_by(interaction_name) %>%
  summarise(
    p_value = wilcox.test(prob[group == "pCR"], prob[group == "non-pCR"])$p.value
  ) %>%
  arrange(p_value)
pval_table


# Figure 5H. Ligand-receptor pairs DotPlot 
CD4_Myeloid$ind <- "FALSE"
CD4_Myeloid$ind[CD4_Myeloid$CellType1 %in%c('cDC1', 'cDC2', 'FABP4+ Macrophage', 'FOLR2+ Macrophage', 'CD14+ Monocyte',
                                                        'CD16+ Monocyte','Proliferating Macrophage', 'SPP1+ Macrophage','OLR1+ Monocyte')] <- TRUE
ACTS30_CD4Myeloid_L <- subset(CD4_Myeloid, subset = ind == TRUE)
ACTS30_CD4Myeloid_L@meta.data <- droplevels(ACTS30_CD4Myeloid_L@meta.data)

CD4_Myeloid$ind <- "FALSE"
CD4_Myeloid$ind[CD4_Myeloid$CellType1 %in%c('CD4_TN', 'CD4_Tcm', 'CD4_Tem', 'CD4_CXCL13', 'CD4_Treg', 'CD4_IL17')] <- TRUE
ACTS30_CD4Myeloid_R <- subset(CD4_Myeloid, subset = ind == TRUE)
ACTS30_CD4Myeloid_R@meta.data <- droplevels(ACTS30_CD4Myeloid_R@meta.data)

Idents(ACTS30_CD4Myeloid_L) <- 'CellType1'
Idents(ACTS30_CD4Myeloid_R) <- 'CellType1'

Figure5Ha<-DotPlot(ACTS30_CD4Myeloid_L, features = 'TNF')
Figure5Hb<-DotPlot(ACTS30_CD4Myeloid_R, features = 'TNFRSF1B')

#Figure 5I. 
# Spatial Cellchat interaction plot 
MyeloidCD4_pCR_Spatial <- readRDS('MyeloidCD4_pCR_SpatialCellchat.rds')
MyeloidCD4_NonpCR_Spatial <- readRDS('MyeloidCD4_NonpCR_SpatialCellchat.rds')

cell_colors1 <-  c(
  'FABP4+ Macrophage' = "#F87D74",
  'FOLR2+ Macrophage' ="#F88008",
  'SPP1+ Macrophage' = "#D7D32D",
  'cDC1' = "#39B600",
  'cDC2' = "#C4C8F7",
  'CD14+ Monocyte' = "#A7FFC1",
  'CD16+ Monocyte' = "#984EA3",
  'OLR1+ Monocyte' = "#1E76A2",
  'CD4_Tcm' = "#97D1DD",
  'CD4_IL17' = "#E76BF3",
  'CD4_CXCL13' = "#FB3B9B",
  'CD4_Tem' = "royalblue",
  'CD4_Treg' = "red",
  'CD4_TN' = "#F9F961")

cells_PCR1 <- rownames(MyeloidCD4_pCR_Spatial@net$count)
color_use_PCR1 <- cell_colors1[match(cells_PCR1, names(cell_colors1))]

cells_NonPCR1 <- rownames(MyeloidCD4_NonpCR_Spatial@net$count)
color_use_NonPCR1 <- cell_colors1[match(cells_NonPCR1, names(cell_colors1))]

groupSize_PCR <- as.numeric(table(MyeloidCD4_pCR_Spatial@idents))
groupSize_NonPCR <- as.numeric(table(MyeloidCD4_NonpCR_Spatial@idents))

MyeloidCD4_pCR_Spatial1 <- netAnalysis_computeCentrality(MyeloidCD4_pCR_Spatial, slot.name = "netP") 
MyeloidCD4_NonpCR_Spatial1 <- netAnalysis_computeCentrality(MyeloidCD4_NonpCR_Spatial, slot.name = "netP")

pathways.show <- c('TNF') 

png("spatial_plot_pCR_umap.png", width = 800, height = 800, res = 90)
netVisual_aggregate(MyeloidCD4_pCR_Spatial1, signaling = pathways.show, color.use = color_use_PCR1, layout = "spatial", edge.width.max = 1, alpha.image = 0.05, vertex.weight = "incoming", vertex.size.max = 4, vertex.label.cex = 3)
dev.off()
img_pCR <- image_read("spatial_plot_pCR_umap.png")
img_rotated_pCR <- image_rotate(img_pCR, 270) 
img_rotated_pCR

png("spatial_plot_NonpCR_umap.png", width = 800, height = 800, res = 90)
netVisual_aggregate(MyeloidCD4_NonpCR_Spatial1, signaling = pathways.show, color.use = color_use_NonPCR1, layout = "spatial", edge.width.max = 1, alpha.image = 0.05, vertex.weight = "incoming", vertex.size.max = 4, vertex.label.cex = 3)
dev.off()
img_NonpCR <- image_read("spatial_plot_NonpCR_umap.png")
img_rotated_NonpCR <- image_rotate(img_NonpCR, 270) 
img_rotated_NonpCR
