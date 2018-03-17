#' FlowSorted.Blood.EPIC
#'
#' Illumina Human Methylation data from EPIC on immunomagnetic sorted adult 
#' blood cell populations. The FlowSorted.Blood.EPIC package contains Illumina
#' HumanMethylationEPIC (\dQuote{EPIC})) DNA methylation microarray data
#' from the immunomethylomics group (manuscript submitted), consisting of 37 
#' magnetic sorted blood cell references and 12 samples, formatted as an
#' RGChannelSet object for  integration and normalization using
#' most of the existing Bioconductor packages.
#'
#' This package contains data similar to the FlowSorted.Blood.450k
#' package consisting of data from peripheral blood samples generated from
#' adult men and women. However, when using the newer EPIC microarray minfi estimates
#' of cell type composition using the FlowSorted.Blood.450k package
#' are less precise compared to actual cell counts. Hence, this package
#' consists of appropriate data for deconvolution of adult blood samples
#' used in for example EWAS relying in the newer EPIC technology.
#'
#' Researchers may find this package useful as these samples represent
#' different cellular populations ( T lymphocytes (CD4+ and CD8+), B cells
#' (CD19+), monocytes (CD14+), NK cells (CD56+) and Neutrophils of cell
#' sorted blood generated with high purity estimates. As a test of accuracy
#' 12 experimental mixtures were reconstructed using fixed amounts of DNA from
#' purified cells. These data can be integrated with
#' the minfi Bioconductor package to estimate cellular composition in
#' users' whole blood Illumina EPIC samples using a modified version of
#' the algorithm constrained projection/quadratic programming described
#' in Houseman et al. 2012. For a slightly more accurate estimations we also
#' offered an IDOL optimized CpG selection for cell deconvolution as the object 
#' IDOLOptimizedCpGs.
#'
#'
#' @format A class: RGChannelSet, dimensions: 1051815 49 
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110554}
#' The FlowSorted.Blood.EPIC objects is based in samples assayed
#' by Brock Christensen and colleagues; manuscript submmited.
#' @seealso
#' References \enumerate{
#' \item DC Koestler et al. (2016). \emph{Improving cell mixture deconvolution by
#' identifying optimal DNA methylation libraries (IDOL)}. BMC bioinformatics.
#' 17, 120.

#' \item EA Houseman et al. (2012) \emph{DNA methylation arrays as surrogate
#' measures of cell mixture distribution}. BMC Bioinformatics 13, 86.
#' \item \pkg{minfi} package for tools for estimating cell type
#' composition in blood using these data

#' }
#' 
#' # Do not run
#' # data("FlowSorted.Blood.EPIC")
#' 
"FlowSorted.Blood.EPIC"

