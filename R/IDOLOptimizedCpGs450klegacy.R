#'     IDOL Optimized CpGs for adult blood DNA methylation deconvolution EPIC
#' 
#' 
#'     This object is nothing more than a vector of length 350 consisting of the
#'     names of the IDOL optimized CpGs to use in the older 450K platform.  
#'     These CpGs are used as the backbone for deconvolution and were selected 
#'     because their methylation signature differs across the six
#'     normal leukocyte subtypes.
#' 
#' 
#' @format An object of class "character" of length 350.
#' 
#'         The format is:
#'         chr [1:350] "cg14232368" "cg15087459" "cg20538211" "cg11944101" "cg00219921" "cg23752651" "cg10978799" "cg22442430" "cg12294633" "cg05492387" "cg24336674" "cg25939861" ...
#' 
#' @references LA Salas et al. (Under review). \emph{An optimized library for 
#' reference-based deconvolution of whole-blood biospecimens assayed using the 
#' Illumina HumanMethylationEPIC BeadArray}. 
#' @references D Koestler et al. (2016). \emph{Improving cell mixture 
#' deconvolution by identifying optimal DNA methylation libraries (IDOL)}. 
#' BMC bioinformatics. 17, 120. doi: 10.1186/s12859-016-0943-7.
#' 
#' @examples
#' # Do not run
#' # data("IDOLOptimizedCpGs450klegacy")
"IDOLOptimizedCpGs450klegacy"