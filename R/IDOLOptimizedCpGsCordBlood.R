#'     IDOL Optimized CpGs for umbilical cord blood DNA methylation 
#'     deconvolution 
#' 
#' @description 
#'     This object is a vector of length 517 consisting of the names of the IDOL
#'     optimized CpGs.  These CpGs are used as the backbone for deconvolution
#'     and were selected because their methylation signature differs across the 
#'     six normal leukocyte subtypes and the nucleated red blood cells.
#' 
#' 
#' @format An object of class "character" of length 517.
#' 
#'         The format is:
#'         chr [1:517] "cg12603453" "cg24765783" "cg06975018" "cg19708055" ...
#' 
#' @references K Gervin, LA Salas et al. (2019) \emph{Systematic evaluation and 
#' validation of references and library selection methods for deconvolution of 
#' cord blood DNA methylation data}. Clin Epigenetics 11,125. doi:
#' 10.1186/s13148-019-0717-y
#' @references LA Salas et al. (2018). \emph{An optimized library for 
#' reference-based deconvolution of whole-blood biospecimens assayed using the 
#' Illumina HumanMethylationEPIC BeadArray}. Genome Biology 19, 64. doi:
#' 10.1186/s13059-018-1448-7.
#' @references DC Koestler et al. (2016). \emph{Improving cell mixture 
#' deconvolution by identifying optimal DNA methylation libraries (IDOL)}. 
#' BMC bioinformatics. 17, 120. doi: 10.1186/s12859-016-0943-7.
#' 
#' @examples
#' #data ("IDOLOptimizedCpGsCordBlood")
#' #head(IDOLOptimizedCpGsCordBlood)
#' #See ?estimateCellCounts2 for deconvolution examples
#' @usage 
#' 
#' #data ("IDOLOptimizedCpGsCordBlood")
#' #head(IDOLOptimizedCpGsCordBlood)
#' #See ?estimateCellCounts2 for deconvolution examples
"IDOLOptimizedCpGsCordBlood"
