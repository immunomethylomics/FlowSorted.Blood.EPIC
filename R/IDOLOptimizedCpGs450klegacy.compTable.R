#' IDOL Optimized CpGs matrix for adult blood DNA methylation deconvolution EPIC
#' 
#' @description 
#'     This object is a matrix of dimensions 350 x 6 consisting of the average 
#'     DNA methylation values fo the probes included in the IDOL optimized 
#'     CpGs per each of the six cell types available to use in the older 450K 
#'     platform.  These CpGs are used as the backbone for deconvolution and 
#'     were selected because their methylation signature differs across the six 
#'     normal leukocyte subtypes. 
#' 
#' @format An object of class "matrix" of dimensions 450 x 6.
#' 
#'         The format is:
#'         num [1:350, 1:6] 0.904  0.122  0.633  0.841  0.135 ...
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
#' # data("IDOLOptimizedCpGs450klegacy.compTable")
#' # head(IDOLOptimizedCpGs450klegacy.compTable)
"IDOLOptimizedCpGs450klegacy.compTable"