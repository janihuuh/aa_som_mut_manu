

fhrb1680 <- readRDS("results/scrnaseq/FHRB1680/FHRB1680_2.rds")
Idents(fhrb1680) <- Idents(fhrb1680) %>% extractClusterNumber() %>% as.numeric %>% as.factor() %>% getClusterPhenotypesFHRB1680()
fhrb1680$cluster <- Idents(fhrb1680)
fhrb1680 <- fhrb1680 %>% fixSeurat
fhrb1680_seurat <- fhrb1680

dir.create("results/manuscript/", showWarnings = F)
dir.create("results/manuscript/som_mut/", showWarnings = F)
dir.create("results/manuscript/som_mut/fhrb1680/", showWarnings = F)


## Plot QC
fhrb1680_seurat@meta.data %>% ggplot(aes(cluster, nCount_RNA, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette5(22)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/manuscript//som_mut/fhrb1680//vln_nCount_rna.pdf", width = 7, height = 4)

fhrb1680_seurat@meta.data %>% ggplot(aes(cluster, nFeature_RNA, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette5(22)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/manuscript//som_mut/fhrb1680//vln_nFeature_RNA.pdf", width = 7, height = 4)

fhrb1680_seurat@meta.data %>% ggplot(aes(cluster, percent.mt, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette5(22)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/manuscript//som_mut/fhrb1680//vln_percent.mt.pdf", width = 7, height = 4)

fhrb1680_seurat@meta.data %>% ggplot(aes(cluster, percent.ribo, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette5(22)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/manuscript//som_mut/fhrb1680//vln_percent.ribo.pdf", width = 7, height = 4)

fhrb1680_seurat@meta.data %>% ggplot(aes(cluster, percent.cycle, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette5(22)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/manuscript//som_mut/fhrb1680//vln_percent.cycle.pdf", width = 7, height = 4)

## Plot most notable markers
p <- DotPlot(fhrb1680_seurat, features = rev(unique(big_markers2)), cols = "RdBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "") #+ theme(legend.position = "top")
ggsave(plot = p, "results/manuscript//som_mut/fhrb1680//dotplot_big_markers2.png", width = 9, height = 6)

fhrb1680_seurat <- fhrb1680_seurat %>% fixSeurat()
p <- DimPlot(fhrb1680_seurat, group.by = "cluster", reduction = "latent_umap", cols = getPalette5(22), label = T, repel = T) + theme_void(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP1", y = "UMAP2")
ggsave(plot = p, "results/manuscript//som_mut/fhrb1680//latent_umap.png", width = 6, height = 5)

p <- DotPlot(fhrb1680_seurat, features = rev(unique(big_markers)), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "")
ggsave(plot = p, "results/manuscript//som_mut/fhrb1680//dotplot_big_markers.png", width = 14, height = 6)

p <- DotPlot(fhrb1680_seurat, features = rev(unique(c("CD3E", "NCAM1", "GZMK", "SELL", "XCL1", "XCL2", "KLRC1", "IL7R", "LTB", "FCGR3A", "GZMA", "GZMB", "GZMH", "GZMM", "KLRC2", "ZEB2", "KLF2", "PRDM1", "GZMH", "LAG3"))), cols = "RdYlBu") +
  ggpubr::rotate_x_text(angle = 90) + labs(x = "", y  = "")
ggsave(plot = p, "results/manuscript//som_mut/fhrb1680//dotplot_dufva_markers.png", width = 9, height = 6)

guo_genes <- guo_markers %>% do.call(what = "c") %>% as.character() %>% rev %>% unique()
p <- DotPlot(fhrb1680_seurat, features = guo_genes, cols = "RdYlBu") +
  ggpubr::rotate_x_text(angle = 90) + labs(x = "", y  = "")
ggsave(plot = p, "results/manuscript//som_mut/fhrb1680//dotplot_guo_markers.png", width = 12, height = 6)

p <- DotPlot(fhrb1680_seurat, features = rev(zhang_cd8_markers_genes), cols = "RdYlBu") +
  ggpubr::rotate_x_text(angle = 90) + labs(x = "", y  = "")
ggsave(plot = p, "results/manuscript//som_mut/fhrb1680//dotplot_zhang_cd8_markers.png", width = 12, height = 6)




##  Per time point
fhrb1680_seurat@meta.data %>% group_by(timepoint, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  ggplot(aes(timepoint, prop, fill = cluster)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette3(22)) + theme_classic(base_size = 12) + ggpubr::rotate_x_text(45)
ggsave("results/manuscript//som_mut/fhrb1680//bar_cluster.pdf", width = 7, height = 4)

fhrb1680_seurat@meta.data %>% group_by(timepoint, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  ggplot(aes(timepoint, prop, color = cluster, group = cluster)) + geom_point() + geom_path() + scale_color_manual(values = getPalette3(22)) + theme_classic(base_size = 12) + ggpubr::rotate_x_text(45)
ggsave("results/manuscript//som_mut/fhrb1680//line_cluster.pdf", width = 7, height = 4)

fhrb1680_seurat@reductions$latent_umap@cell.embeddings %>% as.data.frame %>% bind_cols(fhrb1680_seurat@meta.data) %>%
  ggplot(aes(latentumap_1, latentumap_2)) + stat_density_2d(aes(fill = ..level..), geom = "polygon") + facet_wrap(~timepoint) + scale_fill_distiller(direction=1) + theme(legend.position = "none") + theme_bw(base_size = 12) + labs(x = "UMAP 1", y = "UMAP 2") + facets_nice
ggsave("results/manuscript//som_mut/fhrb1680///umap_dens.png", width = 10, height = 4)









## Clonalities
fhrb1680_seurat$new_clonotypes_id[fhrb1680_seurat$new_clonotypes_id == ""] <- NA
clusters <- fhrb1680_seurat@meta.data %>% mutate(clusters = Idents(fhrb1680_seurat)) %>% group_by(clusters, !is.na(new_clonotypes_id)) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% filter(`!is.na(new_clonotypes_id)` == T) %>% filter(prop > 0.4) %>% pull(clusters) %>% as.character() %>% unique()

fhrb1680_seurat@meta.data %>%
  mutate(cluster = Idents(fhrb1680_seurat)) %>%
  filter(cluster %in% clusters) %>%
  filter(!is.na(new_clonotypes_id)) %>%
  group_by(cluster, new_clonotypes_id) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  summarise(diversity = vegan::diversity(prop), gini = ineq::Gini(prop)) %>%

  ggplot(aes(gini,diversity,label=cluster)) + geom_point(shape = 21, fill = "lightgrey", size = 3) + ggrepel::geom_text_repel() + labs(x = "Gini index", y = "Shannon diveristy")
ggsave("results/manuscript//som_mut/fhrb1680//scatter_gini_clonality.pdf", width = 5, height = 4)

## Prop of TCRs
fhrb1680_seurat@meta.data %>% mutate(clusters = Idents(fhrb1680_seurat)) %>%
  group_by(clusters, !is.na(new_clonotypes_id)) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  filter(`!is.na(new_clonotypes_id)` == T) %>%

  ggplot(aes(reorder(clusters, prop), prop)) + geom_bar(stat = "identity") + ggpubr::rotate_x_text(45) +
  geom_hline(yintercept = 0.5, linetype = "dotted") + labs(x = "", y = "proportion of cells with TCR") + ylim(c(0,1))
ggsave("results/manuscript//som_mut/fhrb1680//bar_tcrab_cluster.pdf", width = 5, height = 4)

fhrb1680_seurat@meta.data %>% mutate(clusters = Idents(fhrb1680_seurat)) %>%
  group_by(clusters) %>% summarise(n = n()) %>%

  ggplot(aes(reorder(clusters, n), n, fill = clusters, label = n)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values = getPalette(24)) + theme(legend.position = "none") + labs(x = "", y = "nCells") +
  geom_text()
ggsave("results/manuscript//som_mut/fhrb1680//bar_n_cluster.pdf", width = 5, height = 4)




## Find DE-genes between all the clusters
all_markers <- FindAllMarkers(fhrb1680_seurat, test.use = "t", verbose = T, max.cells.per.ident = 1e3)
fwrite(all_markers, "results/manuscript//som_mut/fhrb1680//all_markers_1e3.txt", sep = "\t", quote = F, row.names = F)
all_markers <- fread("results/manuscript//som_mut/fhrb1680//all_markers_1e3.txt") %>% mutate(cluster = cluster %>% extractClusterNumber() %>% as.numeric %>% as.factor() %>% getClusterPhenotypesAA())


set.seed(122)
top10  <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

cells <- fhrb1680_seurat@meta.data %>% mutate(cluster = Idents(fhrb1680_seurat)) %>% group_by(cluster) %>% sample_n(50) %>% pull(barcode)
p <- DoHeatmap(fhrb1680_seurat, cells = cells, features = top10$gene, angle = 90, group.colors = getPalette(24)) + theme(legend.position = "none")
ggsave(plot = p, "results/manuscript//som_mut/fhrb1680//heatmap_top10_markers.png", width = 16, height = 20)



#############################################

dir.create("results/scrnaseq/FHRB1680/manuscript", showWarnings = F)

## Figure A: the UMAP
DimPlot(fhrb1680, reduction = "latent_umap", label = T, repel = T, cols = getPalette(22), label.size = 3) + theme_void() + theme(legend.position = "none")#  + theme(axis.text=element_text(size=12))
ggsave("results/scrnaseq/FHRB1680/manuscript/umap.pdf", width = 5, height = 4)


## Figure B; the volcano plot
cd8     <- subset(fhrb1680, idents = cd8_idents)
cd8     <- getLatentUMAP(cd8)
cd8     <- fixSeurat(cd8)
nColors <- length(cd8_idents)
Idents(cd8)     <- Idents(cd8) %>% extractClusterNumber() %>% getClusterPhenotypesFHRB1680()

p <- DimPlot(cd8, reduction = "latent_umap", label = T, cols = getPalette(nColors)) + theme_void() + add_guide
ggsave(plot = p, "results/scrnaseq/FHRB1680/cluster_markers/latent_umap_cd8.png", width = 12, height = 8)

p <- DotPlot(cd8, features = zhang_cd8_markers_genes,  cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/scrnaseq/FHRB1680/cluster_markers/dotplot_zhang_cd8.pdf", width = 18, height = 10)

p <- DotPlot(cd8, features = guo_markers_genes,  cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/scrnaseq/FHRB1680/cluster_markers/dotplot_guo_cd8.pdf", width = 18, height = 10)

cd8_markers     <- FindAllMarkers(cd8, test.use = "t", verbose = T)
cd8_top_markers <- cd8_markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
fwrite(cd8_markers, "results/scrnaseq/FHRB1680/cluster_markers/cd8_markers.txt", sep = "\t", quote = F, row.names = F)
fwrite(cd8_top_markers, "results/scrnaseq/FHRB1680/cluster_markers/cd8_top_markers.txt", sep = "\t", quote = F, row.names = F)




geens          <- fread("/Users/janihuuh/Dropbox/aplastic_anemia_sc/results/scrnaseq/FHRB1680/scvi/genes_to_use.csv")
de_genes       <- fread("results/scrnaseq/FHRB1680/cluster_markers/cd8_markers.txt") %>% mutate(dir = ifelse(avg_logFC > 0, "over", "under"))
stat3_genes    <- fread("results/scrnaseq/FHRB1680/stat3mt/stat3_genes.txt")[-1,] ## downloaded from MSigDB
stat3_de_genes <- de_genes %>% filter(cluster == "9 CD8+ EMRA")

ggplot(data = stat3_de_genes, aes(avg_logFC, -log10(p_val), color = dir)) +
  # geom_point(shape = 21, size = 2) +
  geom_point() +
  ggrepel::geom_text_repel(data = subset(stat3_de_genes, log10(avg_logFC) > 0.1 & avg_logFC > 5), aes(avg_logFC, -log10(p_val), label = gene), color = "black", fontface = 3) +
  geom_vline(xintercept = -5, linetype = "dotted") +
  geom_vline(xintercept = 5, linetype = "dotted") + theme_bw() + xlim(c(-160,160)) + scale_color_manual(values = c("salmon", "dodgerblue")) + theme(legend.position = "none") + labs(x = "average logFC", y = "-log10(pval)")
ggsave("results/scrnaseq/FHRB1680/manuscript/volcano_total_cd8_clusters_all.pdf", width = 6, height = 5)

ggplot(data = stat3_de_genes, aes(avg_logFC, -log10(p_val), color = dir)) +
  geom_point() +
  # ggrepel::geom_text_repel(data = subset(stat3_de_genes, gene %in% stat3_genes$STAT3_02), aes(avg_logFC, -log10(p_val), label = gene), color = "black", fontface = 3) +
  ggrepel::geom_text_repel(data = subset(stat3_de_genes, log10(avg_logFC) > 0.1), aes(avg_logFC, -log10(p_val), label = gene), color = "black", fontface = 3) +
  xlim(c(-5,5)) +
  theme_bw(base_size = 17) + scale_color_manual(values = c("salmon", "dodgerblue")) + theme(legend.position = "none") + labs(x = "average logFC", y = "-log10(pval)")
ggsave("results/scrnaseq/FHRB1680/manuscript/volcano_total_cd8_clusters_lim5.pdf", width = 6, height = 5)





## Figure C; STAT3 expression violin plot
fhrb1680 <- AddModuleScore(fhrb1680, features = list("STAT3"), name = "STAT3", nbin = 20)

df <- fhrb1680@meta.data %>% filter(celltype == "CD8+") %>% filter(timepoint == "1") %>% group_by(cluster) %>% summarise(median = median(STAT31))

fhrb1680@meta.data %>% filter(celltype == "CD8+") %>% filter(timepoint == "1") %>%
  left_join(df) %>%

  ggplot(aes(reorder(cluster, median), STAT31, fill = cluster)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75)) +
  theme_classic(base_size = 12) +
  ggpubr::stat_compare_means(label = "p.format") +
  scale_fill_manual(values = getPalette(7)) + theme(legend.position = "none") +
  # labs(x = "", y = expression(paste("scaled ", italic("STAT3"), "\nexpression", sep = " "))) +
  labs(x = "", y = "scaled STAT3 \nexpression", sep = " ") +

  ggpubr::rotate_x_text(angle = 45)
ggsave("results/scrnaseq/FHRB1680/manuscript/violin_stat3_expression.pdf", width = 5, height = 4)





## Figure D; STAT3 target expression violin plot
stat3_genes <- rownames(fhrb1680)[rownames(fhrb1680) %in% as.vector(stat3_genes$STAT3_02)]
fhrb1680 <- AddModuleScore(fhrb1680, features = list(stat3_genes), name = "STAT3_targets", nbin = 20)

df <- fhrb1680@meta.data %>% filter(celltype == "CD8+") %>% filter(timepoint == "1") %>% group_by(cluster) %>% summarise(median = median(STAT3_targets1))

fhrb1680@meta.data %>% filter(celltype == "CD8+") %>% filter(timepoint == "1") %>%
  left_join(df) %>%

  ggplot(aes(reorder(cluster, median), STAT3_targets1, fill = cluster)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75)) +
  #  geom_boxplot(width = 0.4, outlier.shape = NA) +
  ggsignif::geom_signif(comparisons = list(
    #   c("9 CD8+ effector exhausted", "3 CD8+ EMRA"),
    #   c("9 CD8+ effector exhausted", "4 CD8+ EMRA clonal"),
    #   c("9 CD8+ effector exhausted", "7 CD8+ activated"),
    #   c("9 CD8+ effector exhausted", "8 CD8+ naive CM"),
    #   c("9 CD8+ effector exhausted", "10 CD8+ EM effector"),
    c("9 CD8+ activated", "7 CD8+ activated")), step_increase = 0.1, map_signif_level = T) +
  theme_classic(base_size = 12) +
  ggpubr::stat_compare_means(label = "p.format") +
  scale_fill_manual(values = getPalette(7)) + theme(legend.position = "none") +
  # labs(x = "", y = expression(paste("scaled ", italic("STAT3"), "\ntarget expression score", sep = " "))) +
  labs(x = "", y = "scaled STAT3 \n target score", sep = " ") +

  ggpubr::rotate_x_text(angle = 45)
ggsave("results/scrnaseq/FHRB1680/manuscript/violin_stat3_target_expression.pdf", width = 5, height = 4)










## Figure E; the evolution of CD8+ clusters
df <- fhrb1680@meta.data %>% filter(celltype == "CD8+") %>%
  mutate(timepoint = plyr::revalue(as.factor(timepoint), replace = c("1" = "dg", "2" = "remission 1", "3" = "remission 2"))) %>%
  # mutate(timepoint = as.character(timepoint)) %>%
  group_by(timepoint, cluster) %>% summarise(n = n()) %>% mutate(prop = n /sum(n))

ggplot(df, aes(timepoint, prop, fill = cluster, group = cluster)) +
  geom_path() +
  geom_point(shape = 21, size = 5) +
  ggrepel::geom_text_repel(data = subset(df, timepoint == "dg"), aes(timepoint,prop,label=cluster), hjust = 1,  direction = "y",  nudge_x = -5, segment.size = 0.5) +
  # ggrepel::geom_label_repel(data = subset(df, timepoint == "dg"), aes(timepoint,prop,label=cluster), hjust = 1, direction = "y", segment.size = 1) +
  # ggrepel::geom_text_repel(data = subset(viz_df, timepoint == "3"), aes(timepoint,prop,label=cluster), hjust = 0, nudge_x = 0.1, direction = "y") +
  scale_fill_manual(values = getPalette(7)) + theme_classic(base_size = 12) + labs(y = "freq of CD8+ cells") + add_guide + theme(legend.position = "none") + ggpubr::rotate_x_text(angle = 45) + labs(x = "")
ggsave("results/scrnaseq/FHRB1680/manuscript/lineplot_cd8_clusters.pdf", width = 5, height = 4)



## Figure E; DEG in evolution
deg_time <- fread("results/scrnaseq/FHRB1680/effectOnClusters/deg_cluster.txt") %>% filter(cluster == "9 CD8+ EMRA")  %>% filter(direction == "up")


genes_to_plot = c("PRF1", "GZMB", "GZMH", "GNLY", "KLRG1", "STAT3")
# genes_to_plot = c("PRF1", "GZMB", "GZMH", "GNLY")
fhrb1680$timepoint2 <- plyr::revalue(as.factor(fhrb1680$timepoint), replace = c("1" = "dg", "2" = "remission 1", "3" = "remission 2"))

VlnPlot(fhrb1680, idents = c("9 CD8+ IFNg exhausted"), features = genes_to_plot, group.by = "timepoint2", cols = getPalette(4), pt.size = 0.3, combine = T) + labs(x = "time point")
ggsave("results/scrnaseq/FHRB1680/manuscript/vlnplot_cd8_clusters_cytotoxic_genes.pdf", width = 8, height = 6)

assay_data <- fhrb1680@assays$RNA@data[rownames(fhrb1680@assays$RNA@data) %in% genes_to_plot, ]
df <- data.frame(fhrb1680@meta.data, assay_data %>% t())

df %>% filter(cluster == "9 CD8+ IFNg exhausted") %>% dplyr::select(timepoint2, genes_to_plot) %>% melt(id = "timepoint2") %>%
  ggplot(aes(timepoint2, value, fill = timepoint2)) + geom_violin(adjust = 2, draw_quantiles = c(0.25,0.5,0.75)) + facet_wrap(~variable, scales = "free_y") + scale_fill_manual(values = getPalette(4)) + theme_classic(base_size = 12) + ggpubr::rotate_x_text(angle = 45) + facets_nice + labs(x = "", y = "expression") +
  ggpubr::stat_compare_means(label = "p.format") + theme(legend.position = "none") + geom_jitter(size = 0.3)
ggsave("results/scrnaseq/FHRB1680/manuscript/vlnplot_cd8_clusters_cytotoxic_genes.pdf", width = 6, height = 4)

df %>% filter(cluster == "9 CD8+ IFNg exhausted") %>% dplyr::select(timepoint2,  STAT3_targets1) %>% melt(id = "timepoint2") %>%
  ggplot(aes(timepoint2, value, fill = timepoint2)) + geom_violin(adjust = 2, draw_quantiles = c(0.25,0.5,0.75)) + facet_wrap(~variable, scales = "free_y") + scale_fill_manual(values = getPalette(4)) + theme_classic(base_size = 12) + ggpubr::rotate_x_text(angle = 45) + facets_nice + labs(x = "", y = "expression") +
  ggpubr::stat_compare_means(label = "p.format") + theme(legend.position = "none") + geom_jitter(size = 0.3)
ggsave("results/scrnaseq/FHRB1680/manuscript/vlnplot_cd8_clusters_STAT3_targets1.pdf", width = 4, height = 3)




## Suppl: STAT3 target in time points
fhrb1680@meta.data %>% filter(celltype == "CD8+") %>% filter(cluster == "9 CD8+ effector exhausted") %>%

  ggplot(aes(as.factor(timepoint), testing1, fill = cluster)) +
  geom_violin() +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  facet_wrap(~cluster, scales = "free") +

  ggsignif::geom_signif(comparisons = list(c("1", "2"), c("2", "3")), map_signif_level = T, step_increase = 0.1) +
  theme_classic(base_size = 17) +
  scale_fill_manual(values = getPalette(7)) + theme(legend.position = "none") + labs(x = "", y = expression(paste("scaled ", italic("STAT3"), " target expression", sep = " "))) + ggpubr::rotate_x_text(angle = 45)



  ## Individual umaps on the time points
  runSeurat <- function(seurat_object, cells.to.use){

    clonotype_seurat_meta           <- filter(seurat_object@meta.data, barcode %in% cells.to.use)
    rownames(clonotype_seurat_meta) <- clonotype_seurat_meta$barcode

    clonotype_seurat <- subset(seurat_object, cells = clonotype_seurat_meta$barcode)
    clonotype_seurat <- clonotype_seurat %>% NormalizeData(verbose = F)

    clonotype_seurat <- FindVariableFeatures(clonotype_seurat, selection.method = "vst", nfeatures = 1e3)
    seurat_hvg       <- VariableFeatures(clonotype_seurat)

    ## Remove clonality and unwanted genes
    clonality_genes <- getClonalityGenes(clonotype_seurat)
    unwanted_genes  <- getUnwantedGenes(clonotype_seurat)

    seurat_hvg      <- seurat_hvg[!seurat_hvg %in% clonality_genes]
    seurat_hvg      <- seurat_hvg[!seurat_hvg %in% unwanted_genes]

    clonotype_seurat <- clonotype_seurat %>% ScaleData(verbose = F, features = seurat_hvg)
    clonotype_seurat <- RunPCA(clonotype_seurat, features = seurat_hvg, npcs = 50)
    nPCs <- sum(clonotype_seurat[["pca"]]@stdev > 2)

    clonotype_seurat <- RunUMAP(clonotype_seurat, dims = 1:nPCs)

    return(clonotype_seurat)

  }

  fhrb1680_1_cells <- fhrb1680@meta.data %>% filter(timepoint == "1") %>% pull(barcode)
  fhrb1680_2_cells <- fhrb1680@meta.data %>% filter(timepoint == "2") %>% pull(barcode)
  fhrb1680_3_cells <- fhrb1680@meta.data %>% filter(timepoint == "3") %>% pull(barcode)

  fhrb1680_1 <- runSeurat(fhrb1680, cells = fhrb1680_1_cells)
  fhrb1680_2 <- runSeurat(fhrb1680, cells = fhrb1680_2_cells)
  fhrb1680_3 <- runSeurat(fhrb1680, cells = fhrb1680_3_cells)


  ## Visualize
  DimPlot(fhrb1680_1, reduction = "umap", label = T, repel = T, cols = getPalette(nClusters), label.size = 7) + theme_void() + theme(legend.position = "none")
  ggsave("results/scrnaseq/FHRB1680/manuscript/suppl_umap_1.png", width = 12, height = 10)

  DimPlot(fhrb1680_2, reduction = "umap", label = T, repel = T, cols = getPalette(nClusters), label.size = 7) + theme_void() + theme(legend.position = "none")
  ggsave("results/scrnaseq/FHRB1680/manuscript/suppl_umap_2.png", width = 12, height = 10)

  DimPlot(fhrb1680_3, reduction = "umap", label = T, repel = T, cols = getPalette(nClusters), label.size = 7) + theme_void() + theme(legend.position = "none")
  ggsave("results/scrnaseq/FHRB1680/manuscript/suppl_umap_3.png", width = 12, height = 10)
