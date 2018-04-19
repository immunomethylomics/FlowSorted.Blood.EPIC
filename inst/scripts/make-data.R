### =========================================================================
### FlowSorted.Blood.EPIC RGChannelSet 
### -------------------------------------------------------------------------
###
##Local processing
setwd("D:/Dropbox (Christensen Lab)/Deconvolution libraries/80 Dartmouth Deconvolution Libraries/80 Dartmouth Deconvolution Libraries")
library(minfi)

sheet<-read.metharray.sheet(getwd())
FlowSorted.Blood.EPIC <- read.metharray.exp(targets = sheet,extended = FALSE)
save(FlowSorted.Blood.EPIC, file = "FlowSorted.Blood.EPIC.rda", version = 3, compress = "xz")

## For future reference as an alternative to build this package database, you can download the files from GEO
## Files are available as GSE110554
#library(minfi)
#library(GEOquery)
#getGEOSuppFiles("GSE110554")
#untar("GSE110554/GSE110554_RAW.tar", exdir = "GSE110554/idat")
#idatFiles <- list.files("GSE110554/idat", pattern = "idat.gz$", full = TRUE)
#sapply(idatFiles, gunzip, overwrite = TRUE)
#FlowSorted.Blood.EPIC <- read.metharray.exp("GSE110554/idat")
#geoMat <- getGEO("GSE110554")
#pD <- pData(geoMat[[1]])
#sampleNames(FlowSorted.Blood.EPIC) <- sub(".*_5", "5", sampleNames(FlowSorted.Blood.EPIC))
#rownames(pD) <- pD$title
#pD <- pD[sampleNames(FlowSorted.Blood.EPIC),]
#pData(FlowSorted.Blood.EPIC) <- pD
#FlowSorted.Blood.EPIC
#save(FlowSorted.Blood.EPIC, file = "FlowSorted.Blood.EPIC.rda", version = 3, compress = "xz")