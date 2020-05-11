#' IDOL Optimized CpGs matrix for adult blood DNA methylation deconvolution EPIC
#' 
#' @description 
#'     This object is a matrix of dimensions 450 x 6 consisting of the average 
#'     DNA methylation values fo the probes included in the IDOL optimized 
#'     CpGs per each of the six cell types available.  These CpGs are used as 
#'     the backbone for deconvolution and were selected because their 
#'     methylation signature differs across the six normal leukocyte subtypes.
#' 
#' 
#' @format An object of class "matrix" of dimensions 450 x 6.
#' 
#'         The format is:
#'         num [1:450, 1:6] 0.197  0.105  0.135  0.654  0.246 ...
#' 
#' @references LA Salas et al. (2018). \emph{An optimized library for 
#' reference-based deconvolution of whole-blood biospecimens assayed using the 
#' Illumina HumanMethylationEPIC BeadArray}. Genome Biology 19, 64. doi:
#' \href{https://dx.doi.org/10.1186/s13059-018-1448-7}{10.1186/s13059-018-1448-7}
#' @references DC Koestler et al. (2016). \emph{Improving cell mixture 
#' deconvolution by identifying optimal DNA methylation libraries (IDOL)}. 
#' BMC bioinformatics. 17, 120. doi: 
#' \href{https://dx.doi.org/10.1186/s13059-018-1448-7}{10.1186/s13059-018-1448-7}.
#' 
#' @examples
#' # Do not run
#' # data("IDOLOptimizedCpGs.compTable")
#' # head(IDOLOptimizedCpGs.compTable)
"IDOLOptimizedCpGs.compTable"
