#' estimateCellCounts2 
#' @description 
#' estimateCellCounts2 function allows the use of customized reference   
#' datasets and IDOL probes L-DMR lists
#' @import minfi
#' @import SummarizedExperiment
#' @import S4Vectors
#' @import IlluminaHumanMethylationEPICanno.ilm10b4.hg19
#' @import ExperimentHub
#' @importFrom utils data
#' @importFrom utils read.csv
#' @importFrom utils memory.limit
#' @importFrom  graphics legend
#' @importFrom  graphics plot
#' @importFrom  stats as.formula
#' @importFrom  stats model.matrix
#' @importFrom  stats  lm
#' @importFrom  stats  pf
#' @importFrom  stats  vcov
#' @importFrom genefilter rowFtests
#' @importFrom genefilter rowttests
#' @importFrom quadprog solve.QP
#' @importFrom nlme lme
#' @importFrom nlme getVarCov
#' @examples
#' # Step 1: Load the reference library to extract the artificial mixtures
#' 
#' library(ExperimentHub)
#' hub <- ExperimentHub()
#' query(hub, "FlowSorted.Blood.EPIC")
#' FlowSorted.Blood.EPIC <- hub[["EH1136"]]
#' FlowSorted.Blood.EPIC
#' 
#' # Step 2 separate the reference from the testing dataset if you want to run 
#' # examples for estimations for this function example
#' 
#' RGsetTargets <- FlowSorted.Blood.EPIC[,
#'              FlowSorted.Blood.EPIC$CellType == "MIX"]
#'              
#' sampleNames(RGsetTargets) <- paste(RGsetTargets$CellType,
#'                             seq_len(dim(RGsetTargets)[2]), sep = "_")
#' RGsetTargets
#' 
#' # Step 3: use your favorite package for deconvolution.
#' # Deconvolute a target data set consisting of EPIC DNA methylation 
#' # data profiled in blood, using your prefered method.
#' 
#' # You can use our IDOL optimized DMR library for the EPIC array.  This object
#' # contains a vector of length 450 consisting of the IDs of the IDOL optimized
#' # CpG probes.  These CpGs are used as the backbone for deconvolution and were
#' # selected because their methylation signature differs across the six normal 
#' # leukocyte subtypes.
#' 
#' data (IDOLOptimizedCpGs)
#' head (IDOLOptimizedCpGs)
#' # If you need to deconvolute a 450k legacy dataset use 
#' # IDOLOptimizedCpGs450klegacy instead
#' 
#' # We recommend using Noob processMethod = "preprocessNoob" in minfi for the 
#' # target and reference datasets. 
#' # Cell types included are "CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"
#' 
#' # To use the IDOL optimized list of CpGs (IDOLOptimizedCpGs) use 
#' # estimateCellCounts2 an adaptation of the popular estimateCellCounts in 
#' # minfi. This function also allows including customized reference arrays. 
#' 
#' # Do not run with limited RAM the normalization step requires a big amount 
#' # of memory resources
#' 
#' if (memory.limit()>8000){
#'  countsEPIC<-estimateCellCounts2(RGsetTargets, compositeCellType = "Blood", 
#'                                 processMethod = "preprocessNoob",
#'                                 probeSelect = "IDOL", 
#'                                 cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
#'                                 "Mono", "Neu"), 
#'                                 referencePlatform = 
#'                                 "IlluminaHumanMethylationEPIC",
#'                                 referenceset = NULL,
#'                                 IDOLOptimizedCpGs =IDOLOptimizedCpGs, 
#'                                 returnAll = FALSE)
#'                                 
#' head(countsEPIC$counts)
#' }
#' 
#' # If you prefer CIBERSORT or RPC deconvolution use EpiDISH or similar
#' 
#' # Example not to run
#' 
#' # countsEPIC<-estimateCellCounts2(RGsetTargets, compositeCellType = "Blood", 
#' #                                processMethod = "preprocessNoob",
#' #                                probeSelect = "IDOL", 
#' #                                cellTypes = c("CD8T", "CD4T", "NK", 
#' #                                "Bcell", "Mono", "Neu"), 
#' #                                referencePlatform = 
#' #                                "IlluminaHumanMethylationEPIC",
#' #                                referenceset = NULL,
#' #                                IDOLOptimizedCpGs =IDOLOptimizedCpGs, 
#' #                                returnAll = TRUE)
#' 
#' # library(EpiDISH)
#' # RPC <- epidish(getBeta(countsEPIC2$normalizedData), 
#' # as.matrix(countsEPIC2$compTable[IDOLOptimizedCpGs, 3:8]), method = "RPC")
#' # RPC$estF#RPC count estimates
#' 
#' # CBS <- epidish(getBeta(countsEPIC2$normalizedData), 
#' # as.matrix(countsEPIC2$compTable[IDOLOptimizedCpGs, 3:8]), method = "CBS")
#' # CBS$estF#CBS count estimates
#' 
#' #For umbilical cord blood please use the following
#' 
#' # Do not run
#' # library (FlowSorted.Blood.EPIC)
#' # data("IDOLOptimizedCpGsCordBlood")
#' # Step 1: Load the reference library to extract the umbilical cord samples
#' # library(ExperimentHub)
#' # hub <- ExperimentHub()
#' # myfiles <- query(hub, "FlowSorted.CordBloodCombined.450k")
#' # FlowSorted.CordBloodCombined.450k <- myfiles[[1]]
#' # FlowSorted.CordBloodCombined.450k
#' 
#' # Step 2 separate the reference from the testing dataset if you want to run 
#' # examples for estimations for this function example
#' 
#' # RGsetTargets <- FlowSorted.CordBloodCombined.450k[,
#' # FlowSorted.CordBloodCombined.450k$CellType == "WBC"]
#' # sampleNames(RGsetTargets) <- paste(RGsetTargets$CellType,
#' #                               seq_len(dim(RGsetTargets)[2]), sep = "_")
#' # RGsetTargets
#' 
#' # Step 3: use your favorite package for deconvolution.
#' # Deconvolute a target data set consisting of 450K DNA methylation 
#' # data profiled in blood, using your prefered method.
#' # You can use our IDOL optimized DMR library for the Cord Blood,  This object
#' # contains a vector of length 517 consisting of the IDs of the IDOL optimized
#' # CpG probes.  These CpGs are used as the backbone for deconvolution and were
#' # selected because their methylation signature differs across the six normal 
#' # leukocyte subtypes plus the nucleated red blood cells.
#' 
#' # data (IDOLOptimizedCpGsCordBlood)
#' # head (IDOLOptimizedCpGsCordBlood)
#' # We recommend using Noob processMethod = "preprocessNoob" in minfi for the 
#' # target and reference datasets. 
#' # Cell types included are "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", 
#' # "nRBC"
#' # To use the IDOL optimized list of CpGs (IDOLOptimizedCpGsCordBlood) use 
#' # estimateCellCounts2 from FlowSorted.Blood.EPIC. 
#' # Do not run with limited RAM the normalization step requires a big amount 
#' # of memory resources. Use the parameters as specified below for 
#' # reproducibility.
#' # 
#' # if (memory.limit()>8000){
#' #     countsUCB<-estimateCellCounts2(RGsetTargets, 
#' #                                     compositeCellType = 
#' #                                                "CordBloodCombined", 
#' #                                     processMethod = "preprocessNoob",
#' #                                     probeSelect = "IDOL", 
#' #                                     cellTypes = c("CD8T", "CD4T", "NK",  
#' #                                     "Bcell", "Mono", "Gran", "nRBC"), 
#' #                                     referencePlatform = 
#' #                                         "IlluminaHumanMethylation450k",
#' #                                     referenceset = 
#' #                                      "FlowSorted.CordBloodCombined.450k",
#' #                                     IDOLOptimizedCpGs =
#' #                                       IDOLOptimizedCpGsCordBlood, 
#' #                                     returnAll = FALSE)
#' #     
#' #     head(countsUCB$counts)
#' # }
#' 
#' @references LA Salas et al. (2018). \emph{An optimized library for 
#' reference-based deconvolution of whole-blood biospecimens assayed using the 
#' Illumina HumanMethylationEPIC BeadArray}. Genome Biology 19, 64. doi:
#' 10.1186/s13059-018-1448-7.
#' @references DC Koestler et al. (2016). \emph{Improving cell mixture 
#' deconvolution by identifying optimal DNA methylation libraries (IDOL)}. 
#' BMC bioinformatics. 17, 120. doi: 10.1186/s12859-016-0943-7.
#' @references EA Houseman, et al.(2012)  \emph{DNA methylation arrays as 
#' surrogate measures of cell mixture distribution}. BMC bioinformatics  13:86. 
#' doi:10.1186/1471-2105-13-86.
#' @references AE Jaffe and RA Irizarry.(2014) \emph{Accounting for cellular 
#' heterogeneity is critical in epigenome-wide association studies}. 
#' Genome Biology  15:R31. doi:10.1186/gb-2014-15-2-r31.
#' @references K Gervin, LA Salas et al. (2019) \emph{Systematic evaluation and 
#' validation of references and library selection methods for deconvolution of 
#' cord blood DNA methylation data}. Clin Epigenetics 11,125. doi:
#' 10.1186/s13148-019-0717-y
#' @references KM Bakulski, et al. (2016) \emph{DNA methylation of cord blood 
#' cell types: Applications for mixed cell birth studies}. Epigenetics 11:5. 
#' doi:10.1080/15592294.2016.1161875.
#' @references AJ Titus, et al. (2017). \emph{Cell-type deconvolution from DNA 
#' methylation: a review of recent applications}. Hum Mol Genet 26: R216-R224.
#' @references AE Teschendorff, et al. (2017). \emph{A comparison of 
#' reference-based algorithms for correcting cell-type heterogeneity in 
#' Epigenome-Wide Association Studies}. BMC Bioinformatics 18: 105.
#' @param 
#' rgSet           The input RGChannelSet or raw MethylSet for the procedure.
#' @param
#' compositeCellType   Which composite cell type is being deconvoluted. 
#'                      Should be one of "Blood", "CordBloodCombined",
#'                      "CordBlood", "CordBloodNorway", "CordTissueAndBlood",
#'                      or "DLPFC".
#'                      See details for preferred approaches.
#' @param
#' processMethod Joint normalization/background correction for user and 
#'                reference data 
#'                
#'                Default input, "preprocessNoob" in minfi, you can use "auto" 
#'                for preprocessQuantile for Blood and DLPFC in 450K datasets 
#'                and preprocessNoob otherwise, according to existing 
#'                literature. 
#'                
#'                For MethylSet objects only "preprocessQuantile"  
#'                is available. Set it to any minfi preprocessing function as a 
#'                character if you want to override it, like 
#'                "preprocessFunnorm"
#' @param 
#' probeSelect    How should probes be selected to distinguish cell types? 
#' 
#'                Options include:
#'                1)  "IDOL", for using a customized set of probes obtained from
#'                 IDOL optimization, available for Blood and Umbilical Cord 
#'                 Blood 
#'                2) "both", which selects an equal number (50) of probes (with 
#'                F-stat p-value < 1E-8) with the greatest magnitude of effect 
#'                from the hyper- and hypo-methylated sides, and 
#'                3) "any", which selects the 100 probes (with F-stat p-value 
#'                < 1E-8) with the greatest magnitude of difference regardless 
#'                of direction of effect.  
#'                
#'                
#'                This according to minfi algorithm. Default input "auto" in  
#'                minfi will use "any" for cord blood and "both" otherwise.  
#'                Please see references for more details.
#' @param 
#' cellTypes     A vector of length K that contains the cell type names.  
#'                Default: c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu").
#'                Please notice that this library use Neutrophils instead 
#'                of Granulocytes. See details for your library.
#' @param                 
#' referencePlatform The platform for the reference dataset; if the input  
#'                    rgSet belongs to another platform, it will be converted  
#'                    using minfi function convertArray.
#' @param
#' referenceset It is NULL by default.  
#'             
#'             A custom reference RGChannelSet object (in quotes) if it is not 
#'             a package installed. This option also allows the user to perform
#'             the deconvolution in closed computing clusters without internet
#'             access to ExperimentHub. For that download and save the 
#'             reference and input the resulting object here. If using an 
#'             installed reference package set to NULL.
#' @param
#' IDOLOptimizedCpGs a vector of probe names for cell deconvolution. For 
#'                 EPIC datasets it should be equal to IDOLOptimizedCpGs (no 
#'                 quotes). For legacy 450K datasets you can use 
#'                 IDOLOptimizedCpGs450klegacy (no quotes) if you want to use 
#'                 FlowSorted.EPIC.Blood for the deconvolution. For umbilical
#'                 cord blood use IDOLOptimizedCpGsCordBlood in the 
#'                 FlowSorted.CordBloodCombined.450k package for both 450k and
#'                 EPIC arrays. See details in help for each IDOL list.
#' @param 
#' returnAll    Should the composition table and the normalized user supplied 
#'              data be return? Default is False.
#' @param              
#' verbose Should the function be verbose?
#' @param
#' meanPlot    Whether to plots the average DNA methylation across the  
#'             cell-type discriminating probes within the mixed and sorted 
#'             samples.
#' @param 
#' lessThanOne Should the predictions be constrained to exactly one, 
#'             in minfi default is FALSE, now you can select the option
#' @param 
#' ...  Other arguments for preprocessquantile or other normalizations
#'@return
#' This function will return a list containing matrix of cell counts (counts), 
#' if returnAll=FALSE, or a list containing the counts, mean methylation per 
#' cellType, and the normalized betas (if returnAll is set to TRUE). These 
#' objects are important if you decide to use a different deconvolution 
#' algorithm such as CIBERSORT or robust partial correlation (RPC).
#' @export
estimateCellCounts2 <- function(rgSet, compositeCellType = "Blood", 
                                processMethod = "preprocessNoob", 
                                probeSelect = c("auto", "any", "IDOL"), 
                                cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                                "Mono", "Neu"), 
                        referencePlatform = c("IlluminaHumanMethylation450k", 
                                                "IlluminaHumanMethylationEPIC", 
                                                "IlluminaHumanMethylation27k"), 
                                referenceset = NULL, IDOLOptimizedCpGs = NULL, 
                                returnAll = FALSE, meanPlot = FALSE, 
                                verbose = TRUE, lessThanOne = FALSE,
                                ...) {
    if ((!is(rgSet, "RGChannelSet")) && (!is(rgSet, "MethylSet")))  
        stop(strwrap(sprintf("object is of class '%s', but needs to be of 
                                class 'RGChannelSet' 'RGChannelSetExtended' or 
                                'MethylSet' to use this function", 
                            class(rgSet)), width = 80, prefix = " ", 
                    initial = ""))
    if (!is(rgSet, "RGChannelSet") && 
        (processMethod[1] != "preprocessQuantile")) 
        stop(strwrap(sprintf("object is of class '%s', but needs to be of 
                                class 'RGChannelSet' or 'RGChannelSetExtended' 
                                to use other methods different to 
                                'preprocessQuantile'", class(rgSet)),
                    width = 80, prefix = " ", initial = ""))
    if (is(rgSet, "MethylSet") && 
        (processMethod[1] == "preprocessQuantile")) 
        message(strwrap("[estimateCellCounts2] The function will assume that
                            no preprocessing has been performed. Using 
                            'preprocessQuantile' in prenormalized data is 
                            experimental and it should only be run under the 
                            user responsibility",width = 80, prefix = " ", 
                        initial = ""))
    if (is(rgSet, "RGChannelSetExtended"))
        rgSet <- as(rgSet, "RGChannelSet")
    referencePlatform <- match.arg(referencePlatform)
    rgPlatform <- sub("IlluminaHumanMethylation", "",
                annotation(rgSet)[which(names(annotation(rgSet)) == "array")])
    platform <- sub("IlluminaHumanMethylation", "", referencePlatform)
    if ((compositeCellType == "CordBlood" | 
        compositeCellType == "CordBloodCombined") && (!"nRBC" %in% cellTypes)) 
        message(strwrap("[estimateCellCounts2] Consider including 'nRBC' in 
                        argument 'cellTypes' for cord blood estimation.\n",
                        width = 80, prefix = " ", initial = ""))
    if ((compositeCellType == "Blood") && 
        (referencePlatform == "IlluminaHumanMethylationEPIC") && 
        ("Gran" %in% cellTypes)) 
        message(strwrap("[estimateCellCounts2] Replace 'Gran' for 'Neu' in 
                        argument 'cellTypes' for EPIC blood estimation.\n",
                        width = 80, prefix = " ", initial = ""))
    if ((compositeCellType != "Blood") && 
        ("Neu" %in% cellTypes)) 
        message(strwrap("[estimateCellCounts2] Check whether 'Gran' or 'Neu' is 
                        present in your reference and adjust argument 
                        'cellTypes' for your estimation.\n",
                        width = 80, prefix = " ", initial = ""))
    if (compositeCellType == "CordBloodCombined")
        platform<="450k"
    referencePkg <- sprintf("FlowSorted.%s.%s", compositeCellType, platform)
    subverbose <- max(as.integer(verbose) - 1L, 0L)
    if (!is.null(referenceset)){
        referenceRGset<- get(referenceset)
        if (!is(rgSet, "RGChannelSet"))
            referenceRGset<-preprocessRaw(referenceRGset)
    } else{ 
        if (!require(referencePkg, character.only = TRUE)) 
            stop(strwrap(sprintf("Could not find reference data package for 
                                compositeCellType '%s' and referencePlatform 
                                '%s' (inferred package name is '%s')", 
                                compositeCellType, platform, referencePkg),
                        width = 80, prefix = " ", initial = ""))
        if((referencePkg!="FlowSorted.Blood.EPIC")  && 
            (referencePkg!="FlowSorted.CordBloodCombined.450k")){
            referenceRGset <- get(referencePkg)
        } else if (referencePkg=="FlowSorted.Blood.EPIC"){
            hub <- ExperimentHub()
            referenceRGset <-hub[["EH1136"]]
        } else if (referencePkg=="FlowSorted.CordBloodCombined.450k"){
            hub <- ExperimentHub()
            referenceRGset <-hub[["EH2256"]]            
        }
        if (!is(rgSet, "RGChannelSet"))
            referenceRGset<-preprocessRaw(referenceRGset) 
    }   
    if (rgPlatform != platform) {
        rgSet <- convertArray(rgSet, outType = referencePlatform, 
                                verbose = TRUE)
    }
    if (!"CellType" %in% names(colData(referenceRGset))) 
        stop(strwrap(sprintf("the reference sorted dataset (in this case '%s') 
                            needs to have a phenoData column called 
                            'CellType'"), names(referencePkg),
                    width = 80, prefix = " ", initial = ""))
    if (sum(colnames(rgSet) %in% colnames(referenceRGset)) > 0) 
        stop(strwrap("the sample/column names in the user set must not be in 
                    the reference data ", width = 80, prefix = " ", 
                    initial = ""))
    if (!all(cellTypes %in% referenceRGset$CellType)) 
        stop(strwrap(sprintf("all elements of argument 'cellTypes' needs to be 
                            part of the reference phenoData columns 'CellType' 
                            (containg the following elements: '%s')", 
                            paste(unique(referenceRGset$cellType), 
                                    collapse = "', '")), width = 80, 
                    prefix = " ", initial = ""))
    if (length(unique(cellTypes)) < 2) 
        stop("At least 2 cell types must be provided.")
    if ((processMethod == "auto") && 
        (compositeCellType %in% c("Blood", "DLPFC"))) 
        processMethod <- "preprocessQuantile"
    if ((processMethod == "auto") && 
        (!compositeCellType %in% c("Blood","DLPFC")) && 
        (is(rgSet, "RGChannelSet"))) 
        processMethod <- "preprocessNoob"
    processMethod <- get(processMethod)
    if ((probeSelect == "auto") && 
        (compositeCellType %in% c("CordBloodCombined", "CordBlood", 
                                    "CordBloodNorway", "CordTissueAndBlood"))) {
        probeSelect <- "any"
    }
    if ((probeSelect == "auto") && 
        (!compositeCellType %in% c("CordBloodCombined", "CordBlood", 
                                    "CordBloodNorway", "CordTissueAndBlood"))) {
        probeSelect <- "both"
    }
    if (verbose) 
        message(strwrap("[estimateCellCounts2] Combining user data with 
                        reference (flow sorted) data.\n", width = 80, 
                        prefix = " ", initial = ""))
    newpd <- DataFrame(sampleNames = c(colnames(rgSet), 
                                        colnames(referenceRGset)), 
                        studyIndex = rep(c("user", "reference"), 
                                        times = c(ncol(rgSet), 
                                                    ncol(referenceRGset))))
    if(is.null(rgSet$CellType))
        rgSet$CellType<-rep("NA", dim(rgSet)[2])
    if(is.null(rgSet$Age))
        rgSet$Age<-rep("NA", dim(rgSet)[2])
    if(is.null(rgSet$Sex))
        rgSet$Sex<-rep("NA", dim(rgSet)[2])
    commoncolumn<-intersect(names(colData(rgSet)), 
                            names(colData(referenceRGset)))
    colData(referenceRGset)[commoncolumn] <- mapply(FUN = as,
                                        colData(referenceRGset)[commoncolumn],
                                        vapply(colData(rgSet)[commoncolumn],
                                                class, FUN.VALUE=character(1)),
                                                    SIMPLIFY = FALSE)
    colData(referenceRGset)<-colData(referenceRGset)[commoncolumn]
    colData(rgSet)<-colData(rgSet)[commoncolumn]
    referencePd <- colData(referenceRGset)
    combinedRGset <- combineArrays(rgSet, referenceRGset, 
                                    outType = referencePlatform)
    colData(combinedRGset) <- newpd
    colnames(combinedRGset) <- newpd$sampleNames
    rm(referenceRGset)
    if (verbose) 
        message(strwrap("[estimateCellCounts2] Processing user and reference 
                        data together.\n", width = 80, prefix = " ",
                        initial = ""))
    if (compositeCellType == "CordBlood") {
        if (!is(combinedRGset, "RGChannelSet"))
            combinedRGset@preprocessMethod["rg.norm"]<-
                "Raw (no normalization or bg correction)"
        combinedMset <- processMethod(combinedRGset, verbose = subverbose)
        rm(combinedRGset)
        gc()
        compTable <- get(paste0(referencePkg, ".compTable"))
        combinedMset <- combinedMset[which(rownames(combinedMset) %in% 
                                                rownames(compTable)), ]
    } else {
        if (!is(combinedRGset, "RGChannelSet"))
            combinedRGset@preprocessMethod["rg.norm"]<-
                "Raw (no normalization or bg correction)"
        combinedMset <- processMethod(combinedRGset)
        rm(combinedRGset)
        gc()
    }
    referenceMset <- combinedMset[, combinedMset$studyIndex == "reference"]
    colData(referenceMset) <- as(referencePd, "DataFrame")
    mSet <- combinedMset[, combinedMset$studyIndex == "user"]
    colData(mSet) <- as(colData(rgSet), "DataFrame")
    rm(combinedMset)
    if (probeSelect != "IDOL") {
        if (verbose) 
            message(strwrap("[estimateCellCounts2] Picking probes for 
                            composition estimation.\n", width = 80, 
                            prefix = " ", initial = ""))
        compData <- pickCompProbes(referenceMset, cellTypes = cellTypes, 
                                    compositeCellType = compositeCellType, 
                                    probeSelect = probeSelect)
        coefs <- compData$coefEsts
        if (verbose) 
            message("[estimateCellCounts2] Estimating composition.\n")
        counts <- projectCellType(getBeta(mSet)[rownames(coefs), ], coefs,
                                lessThanOne = lessThanOne)
        rownames(counts) <- colnames(rgSet)
        if (meanPlot) {
            smeans <- compData$sampleMeans
            smeans <- smeans[order(names(smeans))]
            sampleMeans <- c(colMeans(minfi::getBeta(mSet)[rownames(coefs), ]), 
                            smeans)
            sampleColors <- c(rep(1, ncol(mSet)), 1 + 
                                as.numeric(factor(names(smeans))))
            plot(sampleMeans, pch = 21, bg = sampleColors)
            legend("bottomleft", c("blood", levels(factor(names(smeans)))), 
                    col = seq_len(7), pch = 15)
        }
        if (returnAll) {
            list(counts = counts, compTable = compData$compTable, 
                normalizedData = mSet)
        } else {
            list(counts = counts)  
        }
    } else {
        if (verbose) 
            message(strwrap("[estimateCellCounts2] Using IDOL L-DMR probes for 
                            composition estimation.\n", width = 80, 
                            prefix = " ", initial = ""))
        p <- getBeta(referenceMset)
        pd <- as.data.frame(colData(referenceMset))
        rm(referenceMset)
        if (!is.null(cellTypes)) {
            if (!all(cellTypes %in% pd$CellType)) 
                stop(strwrap("elements of argument 'cellTypes' is not part of 
                            'referenceMset$CellType'", width = 80, 
                            prefix = " ", initial = ""))
            keep <- which(pd$CellType %in% cellTypes)
            pd <- pd[keep, ]
            p <- p[, keep]
        }
        pd$CellType <- factor(pd$CellType, levels = cellTypes)
        ffComp <- rowFtests(p, pd$CellType)
        tIndexes <- split(seq(along=pd$CellType), pd$CellType)
        prof <- vapply(tIndexes, function(i) rowMeans(p[,i]),
                        FUN.VALUE=numeric(dim(p)[1]))
        r <- rowRanges(p)
        compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 2]))
        names(compTable)[1] <- "Fstat"
        names(compTable)[c(-2, -1, 0) + ncol(compTable)] <- c("low", 
                                                            "high", "range")
        tstatList <- lapply(tIndexes, function(i) {
            x <- rep(0, ncol(p))
            x[i] <- 1
            return(rowttests(p, factor(x)))
        })
        trainingProbes <- IDOLOptimizedCpGs
        trainingProbes<-trainingProbes[trainingProbes%in%rownames(p)]
        p <- p[trainingProbes, ]
        pMeans <- colMeans(p)
        names(pMeans) <- pd$CellType
        form <- as.formula(sprintf("y ~ %s - 1", paste(levels(pd$CellType), 
                                                        collapse = "+")))
        phenoDF <- as.data.frame(model.matrix(~pd$CellType - 1))
        colnames(phenoDF) <- sub("^pd\\$CellType", "", colnames(phenoDF))
        if (ncol(phenoDF) == 2) {
            X <- as.matrix(phenoDF)
            coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
            coefs <- coefEsts
        } else {
            tmp <- validationCellType(Y = p, pheno = phenoDF, modelFix = form)
            coefEsts <- tmp$coefEsts
            coefs <- coefEsts
        }
        compData<-list(coefEsts = coefEsts, compTable = compTable,
                        sampleMeans = pMeans)
        if (verbose) 
            message("[estimateCellCounts2] Estimating composition.\n")
        counts <- projectCellType(getBeta(mSet)[rownames(coefs), ], coefs,
                                lessThanOne = lessThanOne)
        rownames(counts) <- colnames(rgSet)
        if (meanPlot) {
            smeans <- compData$sampleMeans
            smeans <- smeans[order(names(smeans))]
            sampleMeans <- c(colMeans(getBeta(mSet)[rownames(coefs), ]), smeans)
            sampleColors <- c(rep(1, ncol(mSet)), 1 + 
                                as.numeric(factor(names(smeans))))
            plot(sampleMeans, pch = 21, bg = sampleColors)
            legend("bottomleft", c("blood", levels(factor(names(smeans)))), 
                    col = seq_len(7), pch = 15)
        }
        if (returnAll) {
            list(counts = counts, compTable = compTable, normalizedData = mSet)
        } else {
            list(counts = counts)  
        }
    }
}

#These are minfi internal functions here they are called to keep the function 
#above running for the alternative selection ("auto" and "both")
pickCompProbes <- function(mSet, cellTypes = NULL, numProbes = 50, 
                            compositeCellType = compositeCellType, 
                            probeSelect = probeSelect) {
    p <- getBeta(mSet)
    pd <- as.data.frame(colData(mSet))
    if(!is.null(cellTypes)) {
        if(!all(cellTypes %in% pd$CellType))
            stop(strwrap("elements of argument 'cellTypes' are not part of 
                        'mSet$CellType'", width = 80, prefix = " ", 
                        initial = ""))
        keep <- which(pd$CellType %in% cellTypes)
        pd <- pd[keep,]
        p <- p[,keep]
    }
    ## make cell type a factor 
    pd$CellType <- factor(pd$CellType, levels = cellTypes)
    ffComp <- rowFtests(p, pd$CellType)
    tIndexes <- split(seq(along=pd$CellType), pd$CellType)
    prof <- vapply(tIndexes, function(i) rowMeans(p[,i]), 
                    FUN.VALUE=numeric(dim(p)[1]))
    r <- rowRanges(p)
    compTable <- cbind(ffComp, prof, r, abs(r[,1] - r[,2]))
    names(compTable)[1] <- "Fstat"
    names(compTable)[c(-2,-1,0) + ncol(compTable)] <- c("low", "high", "range") 
    tstatList <- lapply(tIndexes, function(i) {
        x <- rep(0,ncol(p))
        x[i] <- 1
        return(rowttests(p, factor(x)))
    })
    if (probeSelect == "any"){
        probeList <- lapply(tstatList, function(x) {
            y <- x[x[,"p.value"] < 1e-8,]
            yAny <- y[order(abs(y[,"dm"]), decreasing=TRUE),]      
            c(rownames(yAny)[seq_len(numProbes*2)])
        })
    } else {
        probeList <- lapply(tstatList, function(x) {
            y <- x[x[,"p.value"] < 1e-8,]
            yUp <- y[order(y[,"dm"], decreasing=TRUE),]
            yDown <- y[order(y[,"dm"], decreasing=FALSE),]
            c(rownames(yUp)[seq_len(numProbes)], 
                rownames(yDown)[seq_len(numProbes)])
        })
    }
    trainingProbes <- unique(unlist(probeList))
    p <- p[trainingProbes,]
    
    pMeans <- colMeans(p)
    names(pMeans) <- pd$CellType
    form <- as.formula(sprintf("y ~ %s - 1", paste(levels(pd$CellType), 
                                                    collapse="+")))
    phenoDF <- as.data.frame(model.matrix(~pd$CellType-1))
    colnames(phenoDF) <- sub("^pd\\$CellType", "", colnames(phenoDF))
    if(ncol(phenoDF) == 2) { # two group solution
        X <- as.matrix(phenoDF)
        coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
    } else { # > 2 group solution
        tmp <- validationCellType(Y = p, pheno = phenoDF, modelFix = form)
        coefEsts <- tmp$coefEsts
    }
    
    out <- list(coefEsts = coefEsts, compTable = compTable,
                sampleMeans = pMeans)
    return(out)
}
#This is the modified version of the QP/CP from EA Houseman used by minfi
#By default assumes non-negative and a less than one constrain to keep unity
###############################################################################
# FUNCTION:  projectCellType
#    This function predicts the underlying cellular composition of heterogeneous
#    tissue types (i.e., WB) using the constrained projection procedure 
#    described by Houseman et al., (2012) and modified in minfi.  
#
# REQUIRES: quadprog
#
# ARGUMENTS:
#      
# Y:    A J x N matrix of methylation beta-values collected from 
#       mixed/ heterogeneous biospecimen (i.e., WB).  Target set.
#
# coefWBC:    A J x K projection matrix;, i.e., within-cell type mean  
#             methylation matrix across J DMLs and K many cell types
#
# contrastWBC: Contrast for cell composition predictions.  The user  
#               needn't modify this 
#
# nonnegative: Should cell predictions be nonnegative.  Defaults to TRUE
#
# lessThanOne: Should the predictions be constrained to exactly one, 
#                in minfi default is FALSE
#
# RETURNS:   A N x K matrix of cell proportion estimates across the K cell 
#            types for each of the N subjects contained in the Target Set.
#    
###############################################################################
projectCellType <- function(Y, coefCellType, contrastCellType=NULL, 
                    nonnegative=TRUE, lessThanOne=lessThanOne){ 
    if(is.null(contrastCellType))
        Xmat <- coefCellType
    else
        Xmat <- tcrossprod(coefCellType, contrastCellType) 
    nCol <- dim(Xmat)[2]
    if(nCol == 2) {
        Dmat <- crossprod(Xmat)
        mixCoef <- t(apply(Y, 2, function(x) {solve(Dmat, crossprod(Xmat, x))}))
        colnames(mixCoef) <- colnames(Xmat)
        return(mixCoef)
    } else {
        nSubj <- dim(Y)[2]
        mixCoef <- matrix(0, nSubj, nCol)
        rownames(mixCoef) <- colnames(Y)
        colnames(mixCoef) <- colnames(Xmat)
        if(nonnegative){
            if(lessThanOne) {
                Amat <- cbind(rep(-1, nCol), diag(nCol))
                b0vec <- c(-1, rep(0, nCol))
            } else {
                Amat <- diag(nCol)
                b0vec <- rep(0, nCol)
            }
            for(i in seq_len(nSubj)) {
                obs <- which(!is.na(Y[,i])) 
                Dmat <- crossprod(Xmat[obs,])
                mixCoef[i,] <- solve.QP(Dmat, crossprod(Xmat[obs,], Y[obs,i]), 
                                        Amat, b0vec)$sol
            }
        } else {
            for(i in seq_len(nSubj)) {
                obs <- which(!is.na(Y[,i])) 
                Dmat <- crossprod(Xmat[obs,])
                mixCoef[i,] <- solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
            }
        }
        return(mixCoef)
    }
}
validationCellType <- function(Y, pheno, modelFix, modelBatch=NULL,
                                L.forFstat = NULL, verbose = FALSE){
    N <- dim(pheno)[1]
    pheno$y <- rep(0, N)
    xTest <- model.matrix(modelFix, pheno)
    sizeModel <- dim(xTest)[2]
    M <- dim(Y)[1]
    if(is.null(L.forFstat)) {
        L.forFstat <- diag(sizeModel)[-1,] # All non-intercept coefficients
        colnames(L.forFstat) <- colnames(xTest) 
        rownames(L.forFstat) <- colnames(xTest)[-1] 
    }
    # Initialize various containers
    sigmaResid <- sigmaIcept <- nObserved <- nClusters <- Fstat <- rep(NA, M)
    coefEsts <- matrix(NA, M, sizeModel)
    coefVcovs <- list()
    if(verbose)
        cat("[validationCellType] ")
    for(j in seq_len(M)) { # For each CpG
        ## Remove missing methylation values
        ii <- !is.na(Y[j,])
        nObserved[j] <- sum(ii)
        pheno$y <- Y[j,]

        if(j%%round(M/10)==0 && verbose)
            cat(".") # Report progress
        
        try({ # Try to fit a mixed model to adjust for plate
            if(!is.null(modelBatch)) {
                fit <- try(lme(modelFix, random=modelBatch, data=pheno[ii,]))
                OLS <- inherits(fit,"try-error") 
                # If LME can't be fit, just use OLS
            } else
                OLS <- TRUE
            
            if(OLS) {
                fit <- lm(modelFix, data=pheno[ii,])
                fitCoef <- fit$coef
                sigmaResid[j] <- summary(fit)$sigma
                sigmaIcept[j] <- 0
                nClusters[j] <- 0
            } else { 
                fitCoef <- fit$coef$fixed
                sigmaResid[j] <- fit$sigma
                sigmaIcept[j] <- sqrt(getVarCov(fit)[1])
                nClusters[j] <- length(fit$coef$random[[1]])
            }
            coefEsts[j,] <- fitCoef
            coefVcovs[[j]] <- vcov(fit)
            
            useCoef <- L.forFstat %*% fitCoef
            useV <- L.forFstat %*% coefVcovs[[j]] %*% t(L.forFstat)
            Fstat[j] <- (t(useCoef) %*% solve(useV, useCoef))/sizeModel
        })
    }
    if(verbose)
        cat(" done\n")
    ## Name the rows so that they can be easily matched to the target data set
    rownames(coefEsts) <- rownames(Y)
    colnames(coefEsts) <- names(fitCoef)
    degFree <- nObserved - nClusters - sizeModel + 1
    
    ## Get P values corresponding to F statistics
    Pval <- 1-pf(Fstat, sizeModel, degFree)
    
    out <- list(coefEsts=coefEsts, coefVcovs=coefVcovs, modelFix=modelFix, 
                modelBatch=modelBatch,
                sigmaIcept=sigmaIcept, sigmaResid=sigmaResid, 
                L.forFstat=L.forFstat, Pval=Pval,
                orderFstat=order(-Fstat), Fstat=Fstat, nClusters=nClusters, 
                nObserved=nObserved,
                degFree=degFree)
    
    out
}