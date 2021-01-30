
fhrb1641 <- readRDS("results/scrnaseq/scrnaseq/FHRB1641/raw_object/FHRB1641.rds")

###################

## In total, we remove with the following conditions:

min_mito     <- 0
max_mito     <- 10

min_ribo     <- 10
max_ribo     <- 50

min_features <- 500
max_features <- 4500

min_counts   <- 1e3
max_counts   <- 30e3


##############################################################################

percent_mito_outlier <- viz_df %>% dplyr::filter(percent.mt > max_mito | percent.mt < min_mito)                                   %>% pull(barcode) %>% as.character()
percent_ribo_outlier <- viz_df %>% dplyr::filter(percent.ribo > max_ribo | percent.ribo < min_ribo)                               %>% pull(barcode) %>% as.character()
features_outlier     <- viz_df %>% dplyr::filter(total_features < min_features | total_features > max_features)                   %>% pull(barcode) %>% as.character()
umis_outlier         <- viz_df %>% dplyr::filter(total_counts > max_counts | total_counts < min_counts)                           %>% pull(barcode) %>% as.character()

##############################################################################

## Remove the cells from Seurat-object and save a new seurat-object
cells.to.use <- colnames(fhrb1641)[!colnames(fhrb1641) %in% outlier_df$barcode]
fhrb1641     <- subset(fhrb1641, cells = cells.to.use)

saveRDS(fhrb1641, "results/scrnaseq/scrnaseq/FHRB1641/raw_object/FHRB1641.rds")
