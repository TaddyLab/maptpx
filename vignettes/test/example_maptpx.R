

###############  maptpx solutions ########################


library(devtools)
library(CountClust)
library(singleCellRNASeqMouseDeng2014)
library(GTExV6Brain)
gtex.counts <- Biobase::exprs(GTExV6Brain)
gtex.meta_data <- Biobase::pData(GTExV6Brain)
gtex.gene_names <- rownames(gtex.counts)

library(maptpx)

out <- topics(t(gtex.counts),  K=4, tol=0.1)
