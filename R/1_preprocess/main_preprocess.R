
## Run the preprocessing pipeline
setwd("/Users/hru/Dropbox/aplastic_anemia")

source("src/R/scRNAseq/main.R")
source("src/R/scRNAseq/preprocess/run_preprocessRNA.R")
source("src/R/scRNAseq/preprocess/run_preprocessTCRab.R")
source("src/R/scRNAseq/preprocess/run_preprocessRNATCRab.R")

cat("Fin")
