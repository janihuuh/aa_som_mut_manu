
library(dplyr)
library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(gridExtra)
library(data.table)

me=system("whoami", intern = TRUE)
setwd(paste0("/Users/", me, "/Dropbox/aplastic_anemia_sc/"))


## Set global ggplot2-themes and coloring schemes
theme_set(theme_classic())

add_guide   <- guides(colour = guide_legend(override.aes = list(size=5)))
getPalette  <- colorRampPalette(brewer.pal(9, "Set1"))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))
getPalette4 <- colorRampPalette(brewer.pal(9, "Pastel2"))
getPalette5 <- colorRampPalette(brewer.pal(8, "Spectral"))
getPalette5 <- colorRampPalette(brewer.pal(8, "Pastel2"))
facets_nice <- theme(strip.background = element_rect(fill="grey96"), strip.text = element_text(colour = 'black'))


## Run all fun_* codes
for(code in list.files("src/R/scRNAseq/", "fun", full.names = T, recursive = T)){

  print(code)
  source(code)

}
