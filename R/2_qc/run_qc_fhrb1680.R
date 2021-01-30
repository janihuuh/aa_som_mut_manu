
## FHRB1680 was sampled from peripheral blood
dir.create("results/scrnaseq/FHRB1680/qc/")

fhrb1680 <- readRDS("results/scrnaseq/FHRB1680/raw_object/FHRB1680.rds")

min_mito     <- 0
max_mito     <- 10

min_ribo     <- 10
max_ribo     <- 50

min_features <- 300
max_features <- 4500

# min_counts   <- 1e3
min_counts   <- 2e3
max_counts   <- 30e3

## In total, we remove with the following conditions:
percent_mito_outlier <- viz_df %>% dplyr::filter(percent.mt > max_mito | percent.mt < min_mito)                                   %>% pull(barcode) %>% as.character()
percent_ribo_outlier <- viz_df %>% dplyr::filter(percent.ribo > max_ribo | percent.ribo < min_ribo)                               %>% pull(barcode) %>% as.character()
features_outlier     <- viz_df %>% dplyr::filter(total_features < min_features | total_features > max_features)                   %>% pull(barcode) %>% as.character()
umis_outlier         <- viz_df %>% dplyr::filter(total_counts > max_counts | total_counts < min_counts)                           %>% pull(barcode) %>% as.character()


## Remove the cells from Seurat-object and save a new seurat-object
cells.to.use <- colnames(fhrb1680)[!colnames(fhrb1680) %in% outlier_df$barcode]
fhrb1680     <- subset(fhrb1680, cells = cells.to.use)

saveRDS(fhrb1680, "results/scrnaseq/FHRB1680/FHRB1680.rds")
