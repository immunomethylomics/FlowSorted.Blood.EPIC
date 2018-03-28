#'     IDOL Optimized CpGs for adult blood DNA methylation deconvolution EPIC
#' 
#' 
#'     This object is nothing more than a vector of length 450 consisting of the names of 
#'     the IDOL optimized CpGs.  These CpGs are used as the backbone for deconvolution
#'     and were selected because their methylation signature differs across the six
#'     normal leukocyte subtypes.
#' 
#' 
#' @format An object of class "character" of length 450.
#' 
#'         The format is:
#'         chr [1:450] "cg08769189" "cg07661835" "cg00219921" "cg13468685" "cg04329870" "cg14085952" "cg09318840" "cg02133939" "cg07215281" "cg07879474" "cg00482026" "cg20538211" "cg00948513" ...
#' 
#' @references D Koestler et al. (2016). \emph{Improving cell mixture 
#' deconvolution by identifying optimal DNA methylation libraries (IDOL)}. 
#' BMC bioinformatics. 17, 120.
#' 
#' @examples
#' # Do not run
#' # data("IDOLOptimizedCpGs")
"IDOLOptimizedCpGs"
