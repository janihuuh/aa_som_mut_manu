
fhrb1680 <- readRDS("results/scrnaseq/FHRB1680/FHRB1680.rds")

## === Make scVI count files

## By time points, without HVGs (but remove genes with 0 expression across all cells)
fhrb1680_1   <- fhrb1680@assays$RNA@counts[, fhrb1680$timepoint == "1"]
fhrb1680_1   <- fhrb1680_1[Matrix::rowSums(fhrb1680_1) != 0, ]  %>% as.data.frame

fhrb1680_2   <- fhrb1680@assays$RNA@counts[, fhrb1680$timepoint == "2"]
fhrb1680_2   <- fhrb1680_2[Matrix::rowSums(fhrb1680_2) != 0, ] %>% as.data.frame

fhrb1680_3   <- fhrb1680@assays$RNA@counts[, fhrb1680$timepoint == "3"]
fhrb1680_3   <- fhrb1680_3[Matrix::rowSums(fhrb1680_3) != 0, ] %>% as.data.frame


genes.to.use <- intersect(intersect(rownames(fhrb1680_1), rownames(fhrb1680_2)), rownames(fhrb1680_3))

clonality_genes <- c(grep("^TRAV", genes.to.use, value = T), grep("^TRBV", genes.to.use, value = T),
                     grep("^TRGV", genes.to.use, value = T), grep("^TRDV", genes.to.use, value = T),
                     grep("^IGLV", genes.to.use, value = T), grep("^IGLC", genes.to.use, value = T),
                     grep("^IGLL", genes.to.use, value = T), grep("^IGKV", genes.to.use, value = T),
                     grep("^IGHV", genes.to.use, value = T), grep("^IGKC", genes.to.use, value = T),
                     grep("^IGH",  genes.to.use, value = T),  grep("^IGK", genes.to.use, value = T))

genes.to.use <- genes.to.use[!genes.to.use %in% clonality_genes]
fwrite(data.frame(genes.to.use), "results/scrnaseq/fhrb1680/scvi/genes_to_use.csv", sep = ",", quote = F, row.names = T)

fhrb1680_1   <- fhrb1680_1[genes.to.use, ]
fhrb1680_2   <- fhrb1680_2[genes.to.use, ]
fhrb1680_3   <- fhrb1680_3[genes.to.use, ]

fwrite(fhrb1680_1, file = "results/scrnaseq/fhrb1680/scvi/fhrb1680_1_rm_wo_tcr_qc.csv", sep = ",", quote = F, row.names = T)
fwrite(fhrb1680_2, file = "results/scrnaseq/fhrb1680/scvi/fhrb1680_2_rm_wo_tcr_qc.csv", sep = ",", quote = F, row.names = T)
fwrite(fhrb1680_3, file = "results/scrnaseq/fhrb1680/scvi/fhrb1680_3_rm_wo_tcr_qc.csv", sep = ",", quote = F, row.names = T)




## Get latent representation from scVI; UMAP it
set.seed(123)

## One shot
fhrb1680_latent      <- fread("results/scrnaseq/fhrb1680/scvi/oneshot_rm_wo_tcr_qc/fhrb1680rm_oneshot_latent_wo_tcr_qc.csv", header = F)
fhrb1680_batch       <- fread("results/scrnaseq/fhrb1680/scvi/oneshot_rm_wo_tcr_qc/fhrb1680rm_oneshot_indices_wo_tcr_qc.csv", header = F)
fhrb1680_latent_umap <- umapr::umap(fhrb1680_latent)

fwrite(fhrb1680_latent_umap, "results/scrnaseq/fhrb1680/scvi/oneshot_rm_wo_tcr_qc/fhrb1680rm_oneshot_latent_umap_wo_tcr_qc.csv", sep = ",", quote = F, row.names = T)


## Put embeddings into Seurat object
fhrb1680_latent      <- as.matrix(fhrb1680_latent)
fhrb1680_latent_umap <- as.matrix(dplyr::select(fhrb1680_latent_umap, UMAP1:UMAP2))

rownames(fhrb1680_latent)      <- colnames(fhrb1680)
rownames(fhrb1680_latent_umap) <- colnames(fhrb1680)

rownames(fhrb1680_latent) == rownames(fhrb1680_latent_umap)

colnames(fhrb1680_latent_umap) <- c("latent_umap_1", "latent_umap_2")
latent_dim_red            <- CreateDimReducObject(key = "latent", embeddings = as.matrix(x = fhrb1680_latent))
latent_umap_dim_red       <- CreateDimReducObject(key = "latent_umap", embeddings = as.matrix(x = fhrb1680_latent_umap))

fhrb1680[['latent']]      <- latent_dim_red
fhrb1680[['latent_umap']] <- latent_umap_dim_red


## Clustering
res           <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
fhrb1680      <- FindNeighbors(fhrb1680, reduction = "latent", dims = 1:30)
fhrb1680      <- FindClusters(object = fhrb1680, resolution = res)
clustering_columns <- grep("res", colnames(fhrb1680@meta.data), value = T)




#### Decide clustering on res.07 as it agrees best with the latent umap
Idents(fhrb1680) <- fhrb1680$RNA_snn_res.0.7

seurat_file_path = "results/scrnaseq/FHRB1680/FHRB1680_stringent_qc2.rds"
saveRDS(fhrb1680, "results/scrnaseq/FHRB1680/FHRB1680_scVI.rds")
