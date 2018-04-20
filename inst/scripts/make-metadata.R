### =========================================================================
### FlowSorted.Blood.EPIC metadata 
### -------------------------------------------------------------------------
###

meta <- data.frame(
    Title = c("FlowSorted.Blood.EPIC"),
    Description = c(paste0("The FlowSorted.Blood.EPIC package contains ",
"Illumina HumanMethylationEPIC (“EPIC”)) DNA methylation microarray data from ",
"the immunomethylomics group (manuscript submitted), consisting of 37 magnetic ",
"sorted blood cell references and 12 samples, formatted as an RGChannelSet ",
"object (minfi) for integration and normalization using most of the existing ",
"Bioconductor packages.","RGChannelSet R data representation derived from ",
                         "GEO accession GSE110554.")),
    BiocVersion = c("3.7"),
    Genome = rep("hg19", 1), 
    SourceType = rep("tar.gz", 1), 
    SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110554",
    SourceVersion = "Apr 28 2018",
    Species = "Homo sapiens",
    TaxonomyId = 9606,
    Coordinate_1_based = TRUE,
    DataProvider = "GEO",
    Maintainer = "Lucas A Salas <lucas.a.salas.diaz@dartmouth.edu>",
    RDataClass = c("RGChannelSet") ,
    DispatchClass = c(rep("Rda",1)),
RDataPath = c(paste0("FlowSorted.Blood.EPIC/FlowSorted.Blood.EPIC.rda")),
Tags = "",
Notes = c("")
)

write.csv(meta, file="metadata.csv", row.names=FALSE)
