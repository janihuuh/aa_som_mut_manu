
cat("Preprocess TCRab...")
## Read in the sequencing results; the most important files are all_contig and clonotypes -files

fhrb1641_1_contig     <- read.delim("data/scRNAseq+TCRseq/FHRB1641_BM_jul12_TCR/all_contig_annotations.csv", sep = ",")
fhrb1641_2_contig     <- read.delim("data/scRNAseq+TCRseq/FHRB1641_BM_sep13_TCR/all_contig_annotations.csv", sep = ",")
fhrb1641_3_contig     <- read.delim("data/scRNAseq+TCRseq/FHRB1641_BM_oct15_TCR/all_contig_annotations.csv", sep = ",")

fhrb1680_1_contig     <- read.delim("data/scRNAseq+TCRseq/FHRB1680_PB_sep12_TCR/all_contig_annotations.csv", sep = ",")
fhrb1680_2_contig     <- read.delim("data/scRNAseq+TCRseq/FHRB1680_PB_feb13_TCR/all_contig_annotations.csv", sep = ",")
fhrb1680_3_contig     <- read.delim("data/scRNAseq+TCRseq/FHRB1680_PB_oct15_TCR/all_contig_annotations.csv", sep = ",")

fhrb1641_1_clonotype  <- read.delim("data/scRNAseq+TCRseq/FHRB1641_BM_jul12_TCR/clonotypes.csv", sep = ",")
fhrb1641_2_clonotype  <- read.delim("data/scRNAseq+TCRseq/FHRB1641_BM_sep13_TCR/clonotypes.csv", sep = ",")
fhrb1641_3_clonotype  <- read.delim("data/scRNAseq+TCRseq/FHRB1641_BM_oct15_TCR/clonotypes.csv", sep = ",")

fhrb1680_1_clonotype  <- read.delim("data/scRNAseq+TCRseq/FHRB1680_PB_sep12_TCR/clonotypes.csv", sep = ",")
fhrb1680_2_clonotype  <- read.delim("data/scRNAseq+TCRseq/FHRB1680_PB_feb13_TCR/clonotypes.csv", sep = ",")
fhrb1680_3_clonotype  <- read.delim("data/scRNAseq+TCRseq/FHRB1680_PB_oct15_TCR/clonotypes.csv", sep = ",")



## Pre-process the data. The function produces two files, barcoded and clonotyped
preprocess_10X_TCR(contig_file = fhrb1641_1_contig, clonotype_file = fhrb1641_1_clonotype, prefix = "data/scRNAseq+TCRseq/preprocessed/fhrb1641_1")
preprocess_10X_TCR(contig_file = fhrb1641_2_contig, clonotype_file = fhrb1641_2_clonotype, prefix = "data/scRNAseq+TCRseq/preprocessed/fhrb1641_2")
preprocess_10X_TCR(contig_file = fhrb1641_3_contig, clonotype_file = fhrb1641_3_clonotype, prefix = "data/scRNAseq+TCRseq/preprocessed/fhrb1641_3")

preprocess_10X_TCR(contig_file = fhrb1680_1_contig, clonotype_file = fhrb1680_1_clonotype, prefix = "data/scRNAseq+TCRseq/preprocessed/fhrb1680_1")
preprocess_10X_TCR(contig_file = fhrb1680_2_contig, clonotype_file = fhrb1680_2_clonotype, prefix = "data/scRNAseq+TCRseq/preprocessed/fhrb1680_2")
preprocess_10X_TCR(contig_file = fhrb1680_3_contig, clonotype_file = fhrb1680_3_clonotype, prefix = "data/scRNAseq+TCRseq/preprocessed/fhrb1680_3")



## Read in the just processed files
fhrb1641_1_clonotyped <- read.delim("data/scRNAseq+TCRseq/preprocessed/fhrb1641_1_clonotyped.txt")
fhrb1641_2_clonotyped <- read.delim("data/scRNAseq+TCRseq/preprocessed/fhrb1641_2_clonotyped.txt")
fhrb1641_3_clonotyped <- read.delim("data/scRNAseq+TCRseq/preprocessed/fhrb1641_3_clonotyped.txt")

fhrb1680_1_clonotyped <- read.delim("data/scRNAseq+TCRseq/preprocessed/fhrb1680_1_clonotyped.txt")
fhrb1680_2_clonotyped <- read.delim("data/scRNAseq+TCRseq/preprocessed/fhrb1680_2_clonotyped.txt")
fhrb1680_3_clonotyped <- read.delim("data/scRNAseq+TCRseq/preprocessed/fhrb1680_3_clonotyped.txt")

fhrb1641_1_barcode    <- read.delim("data/scRNAseq+TCRseq/preprocessed/fhrb1641_1_barcoded.txt") %>% mutate(barcode_uniq = paste0("FHRB1641_1_", barcode)) %>% mutate(barcode_uniq = substr(barcode_uniq, 1, nchar(barcode_uniq) - 2))
fhrb1641_2_barcode    <- read.delim("data/scRNAseq+TCRseq/preprocessed/fhrb1641_2_barcoded.txt") %>% mutate(barcode_uniq = paste0("FHRB1641_2_", barcode)) %>% mutate(barcode_uniq = substr(barcode_uniq, 1, nchar(barcode_uniq) - 2))
fhrb1641_3_barcode    <- read.delim("data/scRNAseq+TCRseq/preprocessed/fhrb1641_3_barcoded.txt") %>% mutate(barcode_uniq = paste0("FHRB1641_3_", barcode)) %>% mutate(barcode_uniq = substr(barcode_uniq, 1, nchar(barcode_uniq) - 2))

fhrb1680_1_barcode    <- read.delim("data/scRNAseq+TCRseq/preprocessed/fhrb1680_1_barcoded.txt") %>% mutate(barcode_uniq = paste0("FHRB1680_1_", barcode)) %>% mutate(barcode_uniq = substr(barcode_uniq, 1, nchar(barcode_uniq) - 2))
fhrb1680_2_barcode    <- read.delim("data/scRNAseq+TCRseq/preprocessed/fhrb1680_2_barcoded.txt") %>% mutate(barcode_uniq = paste0("FHRB1680_2_", barcode)) %>% mutate(barcode_uniq = substr(barcode_uniq, 1, nchar(barcode_uniq) - 2))
fhrb1680_3_barcode    <- read.delim("data/scRNAseq+TCRseq/preprocessed/fhrb1680_3_barcoded.txt") %>% mutate(barcode_uniq = paste0("FHRB1680_3_", barcode)) %>% mutate(barcode_uniq = substr(barcode_uniq, 1, nchar(barcode_uniq) - 2))

## Make new clonotype id:s for each of the patient profiled
fhrb1641_tot_barcode <- rbind(fhrb1641_1_barcode, fhrb1641_2_barcode, fhrb1641_3_barcode) %>% newClonotype_df() %>% mutate(new_clonotypes_id = paste0("FHRB1641_", new_clonotypes_id))
fhrb1680_tot_barcode <- rbind(fhrb1680_1_barcode, fhrb1680_2_barcode, fhrb1680_3_barcode) %>% newClonotype_df() %>% mutate(new_clonotypes_id = paste0("FHRB1680_", new_clonotypes_id))


## Write down results
write.table(fhrb1641_tot_barcode, "data/scRNAseq+TCRseq/preprocessed/fhrb1641_tcrab.txt",  sep = "\t", row.names = F, quote = F)
fhrb1641_tot_barcode[!duplicated(fhrb1641_tot_barcode$new_clonotypes_id), ] %>% write.table("data/scRNAseq+TCRseq/preprocessed/fhrb1641_uniq_tcrab.txt", sep = "\t", row.names = F, quote = F)

write.table(fhrb1680_tot_barcode, "data/scRNAseq+TCRseq/preprocessed/fhrb1680_tcrab.txt",  sep = "\t", row.names = F, quote = F)
fhrb1680_tot_barcode[!duplicated(fhrb1680_tot_barcode$new_clonotypes_id), ] %>% write.table("data/scRNAseq+TCRseq/preprocessed/fhrb1680_uniq_tcrab.txt", sep = "\t", row.names = F, quote = F)
