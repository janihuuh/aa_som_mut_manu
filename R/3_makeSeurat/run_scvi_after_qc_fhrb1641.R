

fhrb1641 <- readRDS("results/scrnaseq/scrnaseq/FHRB1641/raw_object/FHRB1641.rds")


##### Write count files for scVI
## By time points, without HVGs (but remove genes with 0 expression across all cells)
fhrb1641_1   <- fhrb1641@assays$RNA@counts[, fhrb1641$timepoint == "1"]
fhrb1641_1   <- fhrb1641_1[Matrix::rowSums(fhrb1641_1) != 0, ]  %>% as.data.frame

fhrb1641_2   <- fhrb1641@assays$RNA@counts[, fhrb1641$timepoint == "2"]
fhrb1641_2   <- fhrb1641_2[Matrix::rowSums(fhrb1641_2) != 0, ] %>% as.data.frame

fhrb1641_3   <- fhrb1641@assays$RNA@counts[, fhrb1641$timepoint == "3"]
fhrb1641_3   <- fhrb1641_3[Matrix::rowSums(fhrb1641_3) != 0, ] %>% as.data.frame

genes.to.use <- intersect(intersect(rownames(fhrb1641_1), rownames(fhrb1641_2)), rownames(fhrb1641_3))

clonality_genes <- c(grep("^TRAV", genes.to.use, value = T), grep("^TRBV", genes.to.use, value = T),
                     grep("^TRGV", genes.to.use, value = T), grep("^TRDV", genes.to.use, value = T),
                     grep("^IGLV", genes.to.use, value = T), grep("^IGLC", genes.to.use, value = T),
                     grep("^IGLL", genes.to.use, value = T), grep("^IGKV", genes.to.use, value = T),
                     grep("^IGHV", genes.to.use, value = T), grep("^IGKC", genes.to.use, value = T),
                     grep("^IGH",  genes.to.use, value = T),  grep("^IGK", genes.to.use, value = T))

genes.to.use <- genes.to.use[!genes.to.use %in% clonality_genes]

fhrb1641_1   <- fhrb1641_1[genes.to.use, ]
fhrb1641_2   <- fhrb1641_2[genes.to.use, ]
fhrb1641_3   <- fhrb1641_3[genes.to.use, ]

fwrite(fhrb1641_1, file = "results/scrnaseq/FHRB1641/scvi/fhrb1641_1_qc.csv", sep = ",", quote = F, row.names = T)
fwrite(fhrb1641_2, file = "results/scrnaseq/FHRB1641/scvi/fhrb1641_2_qc.csv", sep = ",", quote = F, row.names = T)
fwrite(fhrb1641_3, file = "results/scrnaseq/FHRB1641/scvi/fhrb1641_3_qc.csv", sep = ",", quote = F, row.names = T)





## Load the scvi and put into seurat object if the cells match
set.seed(123)
fhrb1641_qc_os_latent <- fread("results/scrnaseq/FHRB1641/scvi/oneshot_qc/fhrb1641qc_oneshot_latent.csv")
fhrb1641_qc_os_batch  <- fread("results/scrnaseq/FHRB1641/scvi/oneshot_qc/fhrb1641qc_oneshot_indices.csv")

fhrb1641_qc_os_latent <- fread("results/scrnaseq/FHRB1641/scvi/oneshot_qc/fhrb1641qc_oneshot_latent_wo_tcr.csv")
fhrb1641_qc_os_batch  <- fread("results/scrnaseq/FHRB1641/scvi/oneshot_qc/fhrb1641qc_oneshot_indices_wo_tcr.csv")

fhrb1641_qc_os_latent_umap <- umapr::umap(fhrb1641_qc_os_latent)


## Put embeddings into Seurat object
fhrb1641_latent      <- as.matrix(fhrb1641_qc_os_latent)
fhrb1641_latent_umap <- as.matrix(select(fhrb1641_qc_os_latent_umap, UMAP1:UMAP2))

rownames(fhrb1641_latent)      <- colnames(fhrb1641)
rownames(fhrb1641_latent_umap) <- colnames(fhrb1641)
rownames(fhrb1641_latent) == rownames(fhrb1641_latent_umap)

latent_dim_red            <- CreateDimReducObject(key = "latent", embeddings = as.matrix(x = fhrb1641_latent))
latent_umap_dim_red       <- CreateDimReducObject(key = "latent_umap", embeddings = as.matrix(x = fhrb1641_latent_umap))

fhrb1641[['latent']]      <- latent_dim_red
fhrb1641[['latent_umap']] <- latent_umap_dim_red


## Clustering on the latent space
res                <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
fhrb1641      <- FindNeighbors(fhrb1641, reduction = "latent", dims = 1:30)
fhrb1641      <- FindClusters(object = fhrb1641, resolution = res)
clustering_columns <- grep("res", colnames(fhrb1641@meta.data), value = T)

## Name clusters
fhrb1641$RNA_snn_res.1 <- plyr::revalue(fhrb1641$RNA_snn_res.1, c("0" = "0 NK CD56dim1",
                                                                            "1" = "1 CD8+ EMRA",
                                                                            "2" = "2 CD8+ EM",
                                                                            "3" = "3 CD8+ cytotoxic",
                                                                            "4" = "4 CD4+ naive/CM",
                                                                            "5" = "5 CD8+ EM",

                                                                            "6" = "6 CD4+ TRM",
                                                                            "7" = "7 myeloid monocytes",
                                                                            "8" = "8 NK CD56dim2",
                                                                            "9" = "9 NK NKT",
                                                                           "10" = "10 NK CD56dim_adaptive",

                                                                           "11" = "11 B-cells",
                                                                           "12" = "12 CD8+ naive",
                                                                           "13" = "13 unspecific 1",
                                                                           "14" = "14 MAIT 1",
                                                                           "15" = "15 unspecific 2",

                                                                           "16" = "16 B-cells proBcell_cycling",
                                                                           "17" = "17 myeloid pDC",
                                                                           "18" = "18 junk",
                                                                           "19" = "19 junk erythrocytes",
                                                                           "20" = "20 B-cell plasma_cells"))
Idents(fhrb1641) <- fhrb1641$RNA_snn_res.1



saveRDS(fhrb1641, "results/scrnaseq/FHRB1641/FHRB1641_scvi.rds")
