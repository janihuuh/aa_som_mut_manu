
cat("Create individual seurat-objects ...")

FHRB1641_1 <- Read10X(data.dir = "data/scRNAseq/FHRB1641_BM_jul12/")  %>% CreateSeuratObject(project = "FHRB1641", min.cells = 0, min.features = 0)
FHRB1641_2 <- Read10X(data.dir = "data/scRNAseq/FHRB1641_BM_sep13/")  %>% CreateSeuratObject(project = "FHRB1641", min.cells = 0, min.features = 0)
FHRB1641_3 <- Read10X(data.dir = "data/scRNAseq/FHRB1641_BM_oct15/")  %>% CreateSeuratObject(project = "FHRB1641", min.cells = 0, min.features = 0)

FHRB1680_1 <- Read10X(data.dir = "data/scRNAseq/FHRB1680_PB_sep12/")  %>% CreateSeuratObject(project = "FHRB1680", min.cells = 0, min.features = 0)
FHRB1680_2 <- Read10X(data.dir = "data/scRNAseq/FHRB1680_PB_feb13/")  %>% CreateSeuratObject(project = "FHRB1680", min.cells = 0, min.features = 0)
FHRB1680_3 <- Read10X(data.dir = "data/scRNAseq/FHRB1680_PB_oct15/")  %>% CreateSeuratObject(project = "FHRB1680", min.cells = 0, min.features = 0)


## Merge time points from individual patients
FHRB1641   <- merge(FHRB1641_1, list(FHRB1641_2, FHRB1641_3), add.cell.ids = c("FHRB1641_1", "FHRB1641_2", "FHRB1641_3"), min.cells = 0, min.genes = 0, do.normalize = F, do.scale = F, do.center = F)
FHRB1680   <- merge(FHRB1680_1, list(FHRB1680_2, FHRB1680_3), add.cell.ids = c("FHRB1680_1", "FHRB1680_2", "FHRB1680_3"), min.cells = 0, min.genes = 0, do.normalize = F, do.scale = F, do.center = F)

## Basic QC
cycle.genes  <- c("ANLN", "ASPM","BIRC5","CCNA2","CCNB1","CCNB2","CCND1","CD63","CDC20","CDCA8","CDKN3","CENPE","CENPF",
                  "CEP55","CKAP2L","DLGAP5","FOXM1","GTSE1","H2AFZ","HIST1H1B", "HIST1H1C", "HIST1H1D", "HIST1H1E", "HIST1H2AJ",
                  "HIST1H4C", "HJURP", "HMGB1", "HMGB2", "HMMR", "KIF11", "KIF14", "KIF15", "KIF2C", "LMNA",
                  "MCM3", "MKI67", "NCAPG", "NUSAP1", "PCNA", "PLK1", "PRC1", "RRM2", "SMC4", "STMN1", "TK1", "TOP2A", "TPX2", "TUBA1B",
                  "TUBB", "TYMS", "UBE2C")

FHRB1641  <- PercentageFeatureSet(FHRB1641, pattern = "^MT-", col.name = "percent.mt")
FHRB1641  <- PercentageFeatureSet(FHRB1641, pattern = "^RP", col.name = "percent.ribo")
FHRB1641  <- PercentageFeatureSet(FHRB1641, features = cycle.genes, col.name = "percent.cycle")
FHRB1641@meta.data$barcode   <- colnames(FHRB1641)
FHRB1641@meta.data$timepoint <- substr(colnames(FHRB1641), 10, 10)

FHRB1680  <- PercentageFeatureSet(FHRB1680, pattern = "^MT-", col.name = "percent.mt")
FHRB1680  <- PercentageFeatureSet(FHRB1680, pattern = "^RP", col.name = "percent.ribo")
FHRB1680  <- PercentageFeatureSet(FHRB1680, features = cycle.genes, col.name = "percent.cycle")
FHRB1680@meta.data$barcode   <- colnames(FHRB1680)
FHRB1680@meta.data$timepoint <- substr(colnames(FHRB1680), 10, 10)

# cat("Save individual seurat-objects ...")
# saveRDS(FHRB1641, "results/scrnaseq/FHRB1641/raw_object/FHRB1641.rds")
# saveRDS(FHRB1680, "results/scrnaseq/FHRB1680/raw_object/FHRB1680.rds")
