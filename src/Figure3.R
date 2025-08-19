# Figure 3A. Pseudotime plot 
cds <- as.cell_data_set(CD8T)
cds <- estimate_size_factors(cds)
cds <- pCRprocess_cds(cds, num_dim = 100)
cds <- align_cds(cds, alignment_group = "sample_id")
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)
cds@clusters$UMAP$partitions <- reacreate.partition
list_cluster <- CD8T@active.ident
cds@clusters$UMAP$clusters <- list_cluster
cds@int_colData@listData$reducedDims$UMAP <- CD8T@reductions$umap@cell.embeddings
cds <- learn_graph(cds, use_partition = F)
cds <- order_cells(cds)

Figure3A<-plot_cells(cds,
           color_cells_by = 'pseudotime',
           cell_size=0.18,
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = T,
           label_leaves = FALSE,
           group_label_size = 0)

# Figure 3B. Density plot 
new_dats <- data.frame(
  X = reducedDims(cds)$UMAP[, 1],  
  Y = reducedDims(cds)$UMAP[, 2],  
  pcR = colData(cds)$pcR
)


Figure3Ba<-ggplot(new_dats, aes(x = X, y = Y)) + 
  geom_point(color = "grey", size = 0.5) + theme_bw(base_size = 14) + 
  geom_pointdensity(data = new_dats[new_dats$pcR == "Yes",], adjust = 0.7) + 
  scale_color_viridis() + theme(panel.border = element_rect(color = "black", linewidth = 1),
                                panel.grid.major = element_blank(),  
                                panel.grid.minor = element_blank(),
                                axis.text = element_text(size = 15,face = "bold", color = "black"))

Figure3Bb<-ggplot(new_dats, aes(x = X, y = Y)) + 
  geom_point(color = "grey", size = 0.5) + theme_bw(base_size = 14) + 
  geom_pointdensity(data = new_dats[new_dats$pcR == "No",], adjust = 1) + 
  scale_color_viridis() + theme(panel.border = element_rect(color = "black", linewidth = 1),
                                panel.grid.major = element_blank(),  
                                panel.grid.minor = element_blank(),
                                axis.text = element_text(size = 15,face = "bold", color = "black"))

Figure3Ba+Figure3Bb

# Figure 3C. Lineage plot 
sc <- readRDS('Sling_sc.rds')
sds <- readRDS('Sling_sds.rds') 

FeaturePlot(sc, 
            feature = 'slingshot_pseudotime_curve1', 
            cols = c('red', 'yellow')) + 
  ggtitle('CD8 Tex trajectory') +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank()) +
  labs(colour = 'Pseudotime')

FeaturePlot(sc, 
            feature = 'slingshot_pseudotime_curve19', 
            cols = c('red', 'yellow')) + 
  ggtitle('CD8 Temra trajectory') +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank()) +
  labs(colour = 'Pseudotime')

# Antigen presentation score plot 

antigen_gene <- list(c("HLA-DMB", "HLA-DMA", "CD74", "HLA-DQB1", "HLA-DQA2", "HLA-DRA", "HLA-DRB5", "HLA-DPA1", "HLA-DRB1", "MFSD6", "HLA-DQA1", "HLA-DPB1"))
sc <- AddModuleScore(sc, features = antigen_gene, 
                     name = "Antigen_Score")

df <- sc@meta.data
psts <- slingPseudotime(sds) %>%
  as.data.frame() %>%
  mutate(cells = rownames(.),
         conditions = df$pcR,
         antigen_score = df$Antigen_Score1)  %>%
  pivot_longer(starts_with("Lineage"), values_to = "pseudotime", names_to = "lineages")

psts1 <- subset(psts, lineages == 'Lineage1')
psts19 <- subset(psts, lineages == 'Lineage19')

psts1_pCR <- psts1[,c("cells", "pseudotime", "antigen_score", "conditions")]
psts1_pCR <- subset(psts1_pCR, conditions %in% "Yes")
psts1_NonpCR <- psts1[,c("cells","pseudotime", "antigen_score", "conditions")]
psts1_NonpCR<- subset(psts1_NonpCR, conditions %in% "No")

colnames(psts1_pCR) <- c("cells","x", "y", "group")
colnames(psts1_NonpCR) <- c("cells","x", "y", "group")
rownames(psts1_pCR) <- psts1_pCR$cells
rownames(psts1_NonpCR) <- psts1_NonpCR$cells
psts1_pCR <- na.omit(psts1_pCR)
psts1_NonpCR <- na.omit(psts1_NonpCR)

summarize_data <- function(data, bin_width = 0.4) {
  data %>%
    mutate(x_bin = floor(x / bin_width) * bin_width) %>%
    group_by(x_bin) %>%
    summarise(
      x = mean(x, na.rm = TRUE),
      weighted_mean = sum(y) / n(),
      .groups = 'drop'
    ) %>%
    na.omit()
}

binned_pCR <- summarize_data(psts1_pCR, bin_width = 0.4)
binned_NonpCR <- summarize_data(psts1_NonpCR, bin_width = 0.4)

plot_tcr_line <- function(data1, data2, interval = 1, x_limits = c(0, 10)) {
  x_values <- seq(x_limits[1], x_limits[2], by = interval)
  
  points_data1 <- data1 %>%
    mutate(nearest_x = sapply(x, function(val) x_values[which.min(abs(x_values - val))])) %>%
    filter(abs(x - nearest_x) < interval)
  
  points_data2 <- data2 %>%
    mutate(nearest_x = sapply(x, function(val) x_values[which.min(abs(x_values - val))])) %>%
    filter(abs(x - nearest_x) < interval)
  
  p <- ggplot() +
    geom_line(data = data1, aes(x = x, y = weighted_mean), color = "red", size = 1) +
    geom_line(data = data2, aes(x = x, y = weighted_mean), color = "blue", size = 1) +
    xlab("Pseudotime") + ylab("Antigen Psentation Score (Mean)") +
    scale_x_continuous(limits = x_limits) +  
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
  
  p <- p + 
    geom_point(data = points_data1, aes(x = x, y = weighted_mean), color = "red", size = 3) +
    geom_point(data = points_data2, aes(x = x, y = weighted_mean), color = "blue", size = 3) +
    theme(legend.position = "bottom",
          panel.grid.major = element_blank(),
          axis.text.x = element_text(face = "bold", size = 13, color = "black"),
          axis.text.y = element_text(face = "bold", size = 13, color = "black"),
          axis.line = element_line(size = 1.2, color = "black") 
    ) + xlim(0,10) + ylim(-0.5, 1)
  
  return(p)
}

p_line <- plot_tcr_line(binned_pCR, binned_NonpCR, interval = 1, x_limits = c(0, 10))
print(p_line)

psts19_pCR <- psts19[,c("cells", "pseudotime", "antigen_score", "conditions")]
psts19_pCR <- subset(psts19_pCR, conditions %in% "Yes")
psts19_NonpCR <- psts19[,c("cells","pseudotime", "antigen_score", "conditions")]
psts19_NonpCR<- subset(psts19_NonpCR, conditions %in% "No")

colnames(psts19_pCR) <- c("cells","x", "y", "group")
colnames(psts19_NonpCR) <- c("cells","x", "y", "group")
rownames(psts19_pCR) <- psts19_pCR$cells
rownames(psts19_NonpCR) <- psts19_NonpCR$cells
psts19_pCR <- na.omit(psts19_pCR)
psts19_NonpCR <- na.omit(psts19_NonpCR)

summarize_data <- function(data, bin_width = 0.4) {
  data %>%
    mutate(x_bin = floor(x / bin_width) * bin_width) %>%
    group_by(x_bin) %>%
    summarise(
      x = mean(x, na.rm = TRUE),
      weighted_mean = sum(y) / n(),
      .groups = 'drop'
    ) %>%
    na.omit()
}

binned_pCR <- summarize_data(psts19_pCR, bin_width = 0.4)
binned_NonpCR <- summarize_data(psts19_NonpCR, bin_width = 0.4)

plot_tcr_line <- function(data1, data2, interval = 1, x_limits = c(0, 10.5)) {
  x_values <- seq(x_limits[1], x_limits[2], by = interval)
  
  points_data1 <- data1 %>%
    mutate(nearest_x = sapply(x, function(val) x_values[which.min(abs(x_values - val))])) %>%
    filter(abs(x - nearest_x) < interval)
  
  points_data2 <- data2 %>%
    mutate(nearest_x = sapply(x, function(val) x_values[which.min(abs(x_values - val))])) %>%
    filter(abs(x - nearest_x) < interval)
  
  p <- ggplot() +
    geom_line(data = data1, aes(x = x, y = weighted_mean), color = "red", size = 1) +
    geom_line(data = data2, aes(x = x, y = weighted_mean), color = "blue", size = 1) +
    xlab("Pseudotime") + ylab("Antigen Psentation Score (Mean)") +
    scale_x_continuous(limits = x_limits) +  # x축 범위 고정
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
  
  p <- p + 
    geom_point(data = points_data1, aes(x = x, y = weighted_mean), color = "red", size = 3) +
    geom_point(data = points_data2, aes(x = x, y = weighted_mean), color = "blue", size = 3) +
    theme(legend.position = "bottom",
          panel.grid.major = element_blank(),
          axis.text.x = element_text(face = "bold", size = 13, color = "black"),
          axis.text.y = element_text(face = "bold", size = 13, color = "black"),
          axis.line = element_line(size = 1.2, color = "black") 
    ) + xlim(0,10.5) + ylim(-0.5, 1)
  
  return(p)
}

Figure3C <- plot_tcr_line(binned_pCR, binned_NonpCR, interval = 1, x_limits = c(0, 10.5))
Figure3C

# Figure 3D. Trade-seq DEG pathway
sce1 <- readRDS('Trade_DEG.rds')

deg <- diffEndTest(sce1)

deg <- deg[!is.na(deg$pvalue) & deg$pvalue < 0.05, ]

up_genes <- deg[deg$logFC1_2 > 0, ]
down_genes <- deg[deg$logFC1_2 < 0, ]

up_genes <- rownames(up_genes)
down_genes <- rownames(down_genes)

up_entrez <- bitr(up_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
down_entrez <- bitr(down_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

ego_up <- enrichGO(gene = up_entrez$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pvalueCutoff = 0.05,
                   readable = TRUE)

ego_down <- enrichGO(gene = down_entrez$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     pvalueCutoff = 0.05,
                     readable = TRUE)

selected_pathways_up <- c("lymphocyte differentiation",
                          "T cell differentiation",
                          "regulation of leukocyte proliferation",
                          "myeloid cell differentiation",
                          "B cell activation")

selected_pathways_down <- c("immune response-inhibiting signal transduction",
                            "cell chemotaxis",
                            "regulation of endopeptidase activity",
                            "response to lipopolysaccharide",
                            "regulation of peptidase activity"
)

up_df <- ego_up@result %>%
  filter(Description %in% selected_pathways_up) %>%
  mutate(EnrichmentScore = -log10(pvalue),
         Group = "pCR_lineage")

down_df <- ego_down@result %>%
  filter(Description %in% selected_pathways_down) %>%
  mutate(EnrichmentScore = log10(pvalue),  
         Group = "NonpCR_lineage")

plot_df <- rbind(up_df, down_df)

Figure3D<-ggplot(plot_df, aes(x = reorder(Description, EnrichmentScore), 
                    y = EnrichmentScore, fill = Group)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(limits = c(-5, 5)) + 
  scale_fill_manual(values = c("NonpCR_lineage" = "#1B63A5", "pCR_lineage" = "#D90016")) +
  labs(x = NULL, y = "Enrichment Score (log10 P-value)", 
       title = "Lineage dependent pathway") +
  theme_minimal(base_size = 12) +
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(size = 13, face = "bold", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )

# Figure 3E. TCR similarity between CD8+ Tpex and Tex cells
target_celltype <- c("CD8 Tpex", "CD8 Tex")
combined.TCR<-combineTCR(contig_list, samples=names(contig_list),removeNA = FALSE,
                         removeMulti = FALSE, 
                         filterMulti = T)
combined.TCR<-addVariable(combined.TCR,
                          variable.name = "response_2",
                          variables=response_2)

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

clu<-makeCluster(50)
patientinfo_list<-lapply(1:length(combined.TCR), FUN=function(k){
  patient_TCR<-combined.TCR[[k]]
  patient_TCR$celltype<-sapply(patient_TCR$cell, FUN=function(barcode){
    return(as.character(Total_T@meta.data$T_Celltype1s[match(barcode, rownames(Total_T@meta.data))]))
  })
  return(patient_TCR)
})

TCR_similarity_singlecell<-list()
for(k in 1:length(combined.TCR)){
  cat(k, "\n")
  patient_TCR<-combined.TCR[[k]]
  patient_TCR$celltype<-sapply(patient_TCR$cell, FUN=function(barcode){
    return(as.character(Total_T@meta.data$T_Celltype1s[match(barcode, rownames(Total_T@meta.data))]))
  })
  patientinfo_list[[k]]<-patient_TCR
  clusterExport(clu, c("patient_TCR"))
  TCR_similarity_singlecell[[k]]<-pblapply(1:nrow(patient_TCR), FUN=function(i){
    seq1<-patient_TCR$cdr3_aa2[i]
    if(is.na(seq1)){
      score<-rep(0, nrow(patient_TCR))
      score[i]<-1
      return(score)
    }
    score<-sapply(1:nrow(patient_TCR), FUN=function(j){
      seq2<-patient_TCR$cdr3_aa2[j]
      if(is.na(seq2)){
        return(0)
      }
      p12<-pwalign::pairwiseAlignment(seq1, seq2, substitutionMatrix="BLOSUM62",
                                      gapOpening=5, gapExtension=2,
                                      scoreOnly=TRUE, type="local")
      
      p1<-pwalign::pairwiseAlignment(seq1, seq1, substitutionMatrix="BLOSUM62",
                                     gapOpening=5, gapExtension=2,
                                     scoreOnly=TRUE, type="local")
      p2<-pwalign::pairwiseAlignment(seq2, seq2, substitutionMatrix="BLOSUM62",
                                     gapOpening=5, gapExtension=2,
                                     scoreOnly=TRUE, type="local")
      
      return(p12/sqrt(p1*p2))
    })
    return(score)
  }, cl=clu)
}


tcr_similarity<-lapply(list(c("CD8 Tpex", "CD8 Tex")), FUN=function(x){
  cell1=x[1]
  cell2=x[2]
  
  similarity.df<-lapply(1:19, FUN=function(k){
    cell1.ind<-which(patientinfo_list[[k]]$celltype==cell1)
    cell2.ind<-which(patientinfo_list[[k]]$celltype==cell2)
    cell.TRA=sapply(strsplit(combined.TCR[[k]]$TCR1, split=".", fixed=T), FUN=function(g){
      g[1]
    })
    cell.TRB=sapply(strsplit(combined.TCR[[k]]$TCR2, split=".", fixed=T), FUN=function(g){
      g[1]
    })
    
    normalized_similarity<-TCR_similarity_singlecell[[k]]
    
    similarity.df<-do.call(rbind, lapply(c(cell1.ind, cell2.ind), FUN=function(i){
      if(patientinfo_list[[k]]$celltype[i]==cell1){
        score=normalized_similarity[[i]][cell2.ind]
        data.frame(patient.ind=k, similarity=score, 
                   cell1.TRA=cell.TRA[i], cell2.TRA=cell.TRA[cell2.ind],
                   cell1.TRB=cell.TRB[i], cell2.TRB=cell.TRB[cell2.ind])
      }else{
        score=normalized_similarity[[i]][cell1.ind]
        data.frame(patient.ind=k, similarity=score, 
                   cell1.TRA=cell.TRA[i], cell2.TRA=cell.TRA[cell1.ind],
                   cell1.TRB=cell.TRB[i], cell2.TRB=cell.TRB[cell1.ind])
      }
    }))
    return(similarity.df)
  })
  pcr_ind<-which(sm$response_2=="pCR")
  nonpcr_ind<-which(sm$response_2!="pCR")
  
  for(i in pcr_ind){
    similarity.df[[i]]$label="pCR"
  }
  
  for(i in nonpcr_ind){
    similarity.df[[i]]$label="Non-pCR"
  }
  
  return(do.call(rbind, similarity.df))
})


tcr_comparison_list<-list(
  c("CD8 Tpex", "CD8 Tex")
)


df_plot<-do.call(rbind, lapply(tcr_comparison_list, FUN=function(x){
  cell1=x[1]
  cell2=x[2]
  
  similarity<-lapply(1:19, FUN=function(k){
    cell1.ind<-which(patientinfo_list[[k]]$celltype==cell1)
    cell2.ind<-which(patientinfo_list[[k]]$celltype==cell2)
    normalized_similarity<-TCR_similarity_singlecell[[k]]
    
    
    similarity<-sapply(c(cell1.ind, cell2.ind), FUN=function(i){
      if(patientinfo_list[[k]]$celltype[i]==cell1){
        mean(normalized_similarity[[i]][cell2.ind])
      }else{
        mean(normalized_similarity[[i]][cell1.ind])
      }
    })
    return(similarity)
  })
  pcr_ind<-which(sm$response_2=="pCR")
  nonpcr_ind<-which(sm$response_2!="pCR")
  
  df<-rbind(data.frame(AA_similarity=unlist(similarity[pcr_ind]),
                       label="pCR"),
            data.frame(AA_similarity=unlist(similarity[nonpcr_ind]),
                       label="Non-pCR")
  )
  df$Comparison=paste0(x, collapse = " - ")
  
  return(df)
}))


df_plot$Comparison=factor(df_plot$Comparison, 
                          levels = c("CD8 Tpex - CD8 Tex"))

Figure3Ea<-ggplot(df_plot, aes(x = label, y = AA_similarity, fill = label)) +
  geom_boxplot(width = 0.8, notch = TRUE) +
  theme_classic() +
  stat_compare_means(
    method = "wilcox.test",
    method.args = list(alternative = "less"),
    hide.ns = TRUE,
    label = "p.format",
    label.x = 1.5,
    label.y = max(df_plot$AA_similarity) * 1.05
  ) +
  xlab("") + 
  ylab("Similarity")+
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(color = "black", size = 12, face = "bold"))+
  ggtitle(
    expression(bold("Clones of " * T[PEX] * " - " * T[EX]))
  ) +
  scale_fill_manual(values = c("Non-pCR" = "#1B63A5", 
                               "pCR"    = "#D90016"))+
  labs(fill="") 

Figure3Eb<-clonalHomeostasis(combined.TCR, 
                             cloneCall = cloneCall, chain=chain, group.by = "response_2")+
  ggtitle("Clonal homeostasis of TCR repertoires")+
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(color = "black", size = 12, face = "bold"))

Figure3Ea+Figure3Eb

dev.copy(png, filename = "Figure3E.png",
         width = 2400, height = 1200, res = 300)
dev.off()

# Figure 3F. Spatial transcriptomic maps of CD8+ T cells
use_color <- c(
  'CD8 Tn' = "#C4C8F7",
  'CD8 Tcm' = "#FF7F00",
  'CD8 Tem' = "#39B600",
  'CD8 Trm' = "#E76BF3",
  'CD8 Temra' = "#A7FFC1",
  'CD8 Tpex' = "royalblue",
  'CD8 Tex' = "red")

slide.seq_PCR <- AddMetaData(slide.seq_PCR, metadata = RCTD_PCR_CD8@results$results_df)
slide.seq_NonPCR <- AddMetaData(slide.seq_NonPCR, metadata = RCTD_NonPCR_CD8@results$results_df)

PCR <- SpatialDimPlot(slide.seq_PCR, group.by = "first_type",  images = NULL, image.alpha = 0, cols = use_color, pt.size.factor = 3000)
ggsave("spatial_plot_PCR_CD8T_umap.png", plot = PCR, width = 2.8, height = 2.8, dpi = 300)
img <- image_read("spatial_plot_PCR_CD8T_umap.png")
Figure3Fa <- image_rotate(img, -90) 

Non_PCR <- SpatialDimPlot(slide.seq_NonPCR, group.by = "first_type",  images = NULL,  image.alpha = 0, cols = use_color, pt.size.factor = 3000)
ggsave("spatial_plot_NonPCR_CD8T_umap.png", plot = Non_PCR, width = 2.8, height = 2.8, dpi = 300)
img <- image_read("spatial_plot_NonPCR_CD8T_umap.png")
Figure3Fb <- image_rotate(img, -90)  
