
dir.create("results/scrnaseq/fhrb1641/manuscript", showWarnings = F)

fhrb1641 <- readRDS("results/scrnaseq/FHRB1641/FHRB1641_seurat.rds")
fhrb1641_seurat <- fhrb1641

dir.create("results/manuscript/", showWarnings = F)
dir.create("results/manuscript/som_mut/", showWarnings = F)
dir.create("results/manuscript/som_mut/fhrb1641/", showWarnings = F)

Idents(fhrb1641_seurat) <- Idents(fhrb1641_seurat) %>% extractClusterNumber() %>% as.numeric %>% as.factor() %>% getClusterPhenotypesFHRB1641()
fhrb1641_seurat$cluster <- Idents(fhrb1641_seurat)
fhrb1641_seurat <- fhrb1641_seurat %>% fixSeurat()

## Plot QC
fhrb1641_seurat@meta.data %>% ggplot(aes(cluster, nCount_RNA, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette5(22)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/manuscript/som_mut/fhrb1641/vln_nCount_rna.pdf", width = 7, height = 4)

fhrb1641_seurat@meta.data %>% ggplot(aes(cluster, nFeature_RNA, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette5(22)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/manuscript/som_mut/fhrb1641/vln_nFeature_RNA.pdf", width = 7, height = 4)

fhrb1641_seurat@meta.data %>% ggplot(aes(cluster, percent.mt, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette5(22)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/manuscript/som_mut/fhrb1641/vln_percent.mt.pdf", width = 7, height = 4)

fhrb1641_seurat@meta.data %>% ggplot(aes(cluster, percent.ribo, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette5(22)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/manuscript/som_mut/fhrb1641/vln_percent.ribo.pdf", width = 7, height = 4)

fhrb1641_seurat@meta.data %>% ggplot(aes(cluster, percent.cycle, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette5(22)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/manuscript/som_mut/fhrb1641/vln_percent.cycle.pdf", width = 7, height = 4)



## Plot most notable markers
DimPlot(fhrb1641_seurat, group.by = "cluster", reduction = "latent_umap", cols = getPalette4(22), label = T, repel = T) + theme_void(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP1", y = "UMAP2")
ggsave("results/manuscript/som_mut/fhrb1641/latent_umap.png", width = 6, height = 5)

p <- DotPlot(fhrb1641_seurat, features = rev(unique(big_markers)), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "")
ggsave(plot = p, "results/manuscript/som_mut/fhrb1641/dotplot_big_markers.png", width = 14, height = 6)

p <- DotPlot(fhrb1641_seurat, features = rev(unique(big_markers2)), cols = "RdBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "") #+ theme(legend.position = "top")
ggsave(plot = p, "results/manuscript//som_mut/fhrb1641//dotplot_big_markers2.png", width = 9, height = 6)

p <- DotPlot(fhrb1641_seurat, features = rev(unique(c("CD3E", "NCAM1", "GZMK", "SELL", "XCL1", "XCL2", "KLRC1", "IL7R", "LTB", "FCGR3A", "GZMA", "GZMB", "GZMH", "GZMM", "KLRC2", "ZEB2", "KLF2", "PRDM1", "GZMH", "LAG3"))), cols = "RdYlBu") +
  ggpubr::rotate_x_text(angle = 90) + labs(x = "", y  = "")
ggsave(plot = p, "results/manuscript/som_mut/fhrb1641/dotplot_dufva_markers.png", width = 9, height = 6)

guo_genes <- guo_markers %>% do.call(what = "c") %>% as.character() %>% rev %>% unique()
p <- DotPlot(fhrb1641_seurat, features = guo_genes, cols = "RdYlBu") +
  ggpubr::rotate_x_text(angle = 90) + labs(x = "", y  = "")
ggsave(plot = p, "results/manuscript/som_mut/fhrb1641/dotplot_guo_markers.png", width = 12, height = 6)

p <- DotPlot(fhrb1641_seurat, features = (zhang_cd8_markers_genes), cols = "RdYlBu") +
  ggpubr::rotate_x_text(angle = 90) + labs(x = "", y  = "")
ggsave(plot = p, "results/manuscript/som_mut/fhrb1641/dotplot_zhang_cd8_markers.png", width = 14, height = 6)

##  Per time point
fhrb1641_seurat@meta.data %>% group_by(timepoint, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  ggplot(aes(timepoint, prop, fill = cluster)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette3(22)) + theme_classic(base_size = 12) + ggpubr::rotate_x_text(45)
ggsave("results/manuscript/som_mut/fhrb1641/bar_cluster.pdf", width = 7, height = 4)

fhrb1641_seurat@meta.data %>% group_by(timepoint, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  ggplot(aes(timepoint, prop, color = cluster, group = cluster)) + geom_point() + geom_path() + scale_color_manual(values = getPalette3(22)) + theme_classic(base_size = 12) + ggpubr::rotate_x_text(45)
ggsave("results/manuscript/som_mut/fhrb1641/line_cluster.pdf", width = 7, height = 4)

fhrb1641_seurat@reductions$latent_umap@cell.embeddings %>% as.data.frame %>% bind_cols(fhrb1641_seurat@meta.data) %>%
  ggplot(aes(latentumap_1, latentumap_2)) + stat_density_2d(aes(fill = ..level..), geom = "polygon") + facet_wrap(~timepoint) + scale_fill_distiller(direction=1) + theme(legend.position = "none") + theme_bw(base_size = 12) + labs(x = "UMAP 1", y = "UMAP 2") + facets_nice
ggsave("results/manuscript/som_mut/fhrb1641//umap_dens.png", width = 10, height = 4)

## Prop of TCRab
fhrb1641_seurat@meta.data %>% mutate(clusters = Idents(fhrb1641_seurat)) %>%
  group_by(clusters, !is.na(new_clonotypes_id)) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  filter(`!is.na(new_clonotypes_id)` == T) %>%

  ggplot(aes(reorder(clusters, prop), prop)) + geom_bar(stat = "identity") + ggpubr::rotate_x_text(45) +
  geom_hline(yintercept = 0.5, linetype = "dotted") + labs(x = "", y = "proportion of cells with TCR") + ylim(c(0,1))
ggsave("results/manuscript/som_mut/fhrb1641/bar_tcrab_cluster.pdf", width = 5, height = 4)

## Clonalities
fhrb1641_seurat$new_clonotypes_id[fhrb1641_seurat$new_clonotypes_id == ""] <- NA
clusters <- fhrb1641_seurat@meta.data %>% mutate(clusters = Idents(fhrb1641_seurat)) %>% group_by(clusters, !is.na(new_clonotypes_id)) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% filter(`!is.na(new_clonotypes_id)` == T) %>% filter(prop > 0.4) %>% pull(clusters) %>% as.character() %>% unique()

fhrb1641_seurat@meta.data %>%
  mutate(cluster = Idents(fhrb1641_seurat)) %>%
  filter(cluster %in% clusters) %>%
  filter(!is.na(new_clonotypes_id)) %>%
  group_by(cluster, new_clonotypes_id) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  summarise(diversity = vegan::diversity(prop), gini = ineq::Gini(prop)) %>%

  ggplot(aes(gini,diversity,label=cluster)) + geom_point(shape = 21, fill = "lightgrey", size = 3) + ggrepel::geom_text_repel() + labs(x = "Gini index", y = "Shannon diveristy")
ggsave("results/manuscript/som_mut/fhrb1641/scatter_gini_clonality.pdf", width = 5, height = 4)

## Find DE-genes between all the clusters
all_markers <- FindAllMarkers(fhrb1641_seurat, test.use = "t", verbose = T, max.cells.per.ident = 1e3)
fwrite(all_markers, "results/manuscript/som_mut/fhrb1641/all_markers_1e3.txt", sep = "\t", quote = F, row.names = F)
all_markers <- fread("results/manuscript/som_mut/fhrb1641/all_markers_1e3.txt") %>% mutate(cluster = cluster %>% extractClusterNumber() %>% as.numeric %>% as.factor() %>% getClusterPhenotypesAA())

set.seed(122)
top10  <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC
cells <- fhrb1641_seurat@meta.data %>% mutate(cluster = Idents(fhrb1641_seurat)) %>% group_by(cluster) %>% sample_n(50) %>% pull(barcode)
p <- DoHeatmap(fhrb1641_seurat, cells = cells, features = top10$gene, angle = 90, group.colors = getPalette(24)) + theme(legend.position = "none")
ggsave(plot = p, "results/manuscript/som_mut/fhrb1641/heatmap_top10_markers.png", width = 16, height = 20)





## Focus on CD8+ cells
fhrb1641_seurat$cluster  <- Idents(fhrb1641_seurat)
fhrb1641_seurat$celltype <- fhrb1641_seurat$cluster %>% extractCoarsePhenotype() %>% do.call(what = "c")
cells.to.keep            <- fhrb1641_seurat@meta.data %>% filter(celltype == "CD8+") %>% pull(barcode)

fhrb1641_cd8_seurat <- subset(fhrb1641_seurat, cells = cells.to.keep)
fhrb1641_cd8_seurat <- fhrb1641_cd8_seurat %>% getLatentUMAP() %>% fixSeurat()

DimPlot(fhrb1641_cd8_seurat, reduction = "latent_umap", label = T, repel = T, cols = getPalette3(7)) + theme_void() + theme(legend.position = "none")
ggsave("results/manuscript/som_mut/fhrb1641/umap_cd8.png", width = 5, height = 4)

## Find DEGs
all_cd8_markers <- FindAllMarkers(fhrb1641_cd8_seurat, test.use = "t")
cd8_markers <- all_cd8_markers %>% filter(avg_logFC > 0 & p_val_adj < 0.05)
fwrite(cd8_markers, "results/manuscript/som_mut/fhrb1641/cd8_markers.txt", sep = "\t", quote = F, row.names = F)








## Figure B; the CD8+ UMAP
cd8_tcrab_df         <- fhrb1641_cd8_seurat@meta.data %>% filter(celltype == "CD8+") %>% filter(!is.na(new_clonotypes_id)) %>% mutate(is_kras = ifelse(trb_cdr3s_aa == "CASSPHRNTEAFF", "KRAS_mt", "KRAS_wt"))
kras_tcrab_df        <- cd8_tcrab_df %>% filter(is_kras == "KRAS_mt")

fhrb1641_cd8_seurat$cluster <- Idents(fhrb1641_cd8_seurat)
fhrb1641_cd8_seurat$is_kras <- ifelse(fhrb1641_cd8_seurat$new_clonotypes_id %in% kras_tcrab_df$new_clonotypes_id, "kras_mt", "kras_wt")

cd8_viz_df <- cbind(fhrb1641_cd8_seurat@meta.data, fhrb1641_cd8_seurat@reductions$latent_umap@cell.embeddings) %>% dplyr::rename(umap_1 = latentumap_1, umap_2 = latentumap_2)

umap_mean <- data.frame(aggregate(umap_1 ~ cluster, cd8_viz_df, median), umap_2 = aggregate(umap_2 ~ cluster, cd8_viz_df, median)[,2])
nClusters <- levels(cd8_viz_df$cluster) %>% length

cd8_viz_df$timepoint2 <- plyr::revalue(as.factor(cd8_viz_df$timepoint), c("1" = "diagnosis", "2" = "active disease", "3" = "relapse"))

## Plot UMAPs with KRAS highlighted
ggplot() +
  geom_point(data = cd8_viz_df, aes(x = umap_1, y = umap_2), color = "gray90", size = 0.8) +
  stat_ellipse(data = cd8_viz_df, aes(x = umap_1, y = umap_2, group = cluster), color = "black", size = 0.5, linetype = "dotted") +
  geom_point(data = subset(cd8_viz_df, is_kras == "kras_mt"), aes(x = umap_1, y = umap_2, fill = as.factor(timepoint2)), size = 3, shape = 21) +

  ggrepel::geom_label_repel(data = umap_mean, aes(x = umap_1, y = umap_2, color = cluster, label = cluster), size = 3, color = "black") +
  theme_void(base_size = 12) + #theme(legend.position = "none") +
  scale_color_manual(values = getPalette(nClusters)) +
  scale_fill_manual(values = getPalette3(4)[c(1,3,2)]) + labs(fill = "Vbeta 7.2 clonotype \ntime point") + add_guide
ggsave("results/scrnaseq/FHRB1641/manuscript/latent_umap_vbeta7.2_evolution.png", width = 7, height = 4)






## Figure C; the KRASmt bar plot
kras_tcrab_df$celltype <- kras_tcrab_df$cluster %>% extractCoarsePhenotype()
kras_tcrab_df %>%
  filter(celltype == "CD8+") %>%
  mutate(timepoint2 = plyr::revalue(as.factor(timepoint), c("1" = "diagnosis", "2" = "active disease", "3" = "relapse"))) %>%

  tidyr::complete(timepoint2, cluster, fill = list(z = 0)) %>%
  group_by(timepoint2, cluster) %>% summarise(n = n()) %>% mutate(prop = n/sum(n)) %>%
  ggplot(aes(timepoint2, prop, fill = cluster, group=cluster)) + geom_area() + theme_bw() + scale_fill_manual(values = getPalette3(9)) +
  labs(x = "", y = "freq of Vbeta 7.2 \nclonotype cells") + theme_classic(base_size = 17) + ggpubr::rotate_x_text(45)
ggsave("results/scrnaseq/FHRB1641/manuscript/bar_evolution.pdf", width = 7, height = 5)





## Figure D; the KRASmt vs KRASwt at relapse
cd8_3      <- subset(fhrb1641_cd8_seurat, subset = timepoint == 3)
cd8_3_emra <- subset(cd8_3, idents = "1 CD8+ EMRA")

Idents(cd8_3_emra) <- cd8_3_emra$is_kras %>% as.factor()
kras_markers3 <- FindMarkers(object = cd8_3_emra, ident.1 = "kras_mt", ident.3 = "kras_wt", only.pos = F, min.pct = 1e-5, logfc.threshold = 1e-5, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))
kras_markers3 <- kras_markers3 %>% mutate(direction = ifelse(avg_logFC > 0, "up_in_kras_mt", "down_in_kras_mt")) %>% mutate(direction = ifelse(abs(avg_logFC) > 0.35 & p_val_adj < 0.05, direction, "unsigf"))
fwrite(kras_markers3, "results/scrnaseq/FHRB1641/kras_clone/volcano_kras_markers_3.txt", sep = "\t", quote = F, row.names = F)

kras_markers3 <- fread("results/scrnaseq/FHRB1641/kras_clone/volcano_kras_markers_3.txt")

ggplot(kras_markers3, aes(avg_logFC, -log10(p_val), color = direction)) + geom_point() + xlim(values = c(-5, 5)) + theme_bw() + scale_color_manual(values = c("dodgerblue", "gray90", "salmon")) +
  labs(x = "average logFC", y = "-log10(pval)") +
  # geom_point(data = subset(kras_markers3, gene == "KRAS"), aes(avg_logFC, -log10(p_val)), color = "black") +
  ggrepel::geom_text_repel(data = subset(kras_markers3, gene == "KRAS"), aes(label = gene), fontface = "italic", size = 5, color = "black") +
  ggrepel::geom_text_repel(data = subset(kras_markers3, direction == "up_in_kras_mt" & -log10(p_val) > 2.5), aes(label = gene), fontface = "italic", segment.size = 0.5) +
  theme_bw(base_size = 12) + theme(legend.position = "none")
ggsave("results/scrnaseq/FHRB1641/manuscript/volcano_kras_markers_3.pdf", width = 5, height = 4)
