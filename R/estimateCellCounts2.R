#' estimateCellCounts2 function allows the use of customized reference datasets 
#' and IDOL probes L-DMR lists
#' @import minfi
#' @examples
#' # Step 1: Load the library(FlowSorted.Blood.EPIC)
#' data(FlowSorted.Blood.EPIC)
#' # Step 2 separate the reference from the testing dataset if you want to run 
#' # examples for estimations for this function example
#' cell.Mix <- which(FlowSorted.Blood.EPIC$CellType == "MIX")
#' RGsetTargets <- FlowSorted.Blood.EPIC[, cell.Mix]
#' sampleNames(RGsetTargets) <- paste(RGsetTargets$CellType,
#'                             seq(along = cell.Mix), sep = "_")
#' # Step 3: use your favorite package for deconvolution.
#' # You can also read in the IDOL optimized DMR library based on the EPIC 
#' # array.  Thisobject is nothing more than a vector of length 450 consisting 
#' # of the names of the IDOL optimized CpGs.  These CpGs are used as the 
#' # backbone for deconvolution and were selected because their methylation 
#' # signature differs across the six normal leukocyte subtypes.
#' data(IDOLOptimizedCpGs)
#' # Step 4: Deconvolute a target data set consisting of EPIC DNA methlation 
#' # data profiled in blood, using your prefered method estimateCellCounts 
#' # (minfi), or similar.
#' # We recommend using Noob processMethod = "preprocessNoob" in minfi
#' # for the target and reference datasets Cell types included are 
#' # cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")
#' # We also provide an IDOL optimized list of CpGs (IDOLOptimizedCpGs)
#' # that can be used for optimized cell estimations
#' # We offer an adaptation of the popular estimateCellCounts in minfi to allow 
#' # the inclusion of customized reference arrays to select the custom IDOL use 
#' # probeSelect = "IDOL"
#' countsEPIC<-estimateCellCounts2(RGsetTargets, compositeCellType = "Blood", 
#'                                 processMethod = "preprocessNoob",
#'                                 probeSelect = "IDOL", 
#'                                 cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
#'                                 "Mono", "Neu"), 
#'                                 referencePlatform = 
#'                                 "IlluminaHumanMethylationEPIC",
#'                                 referenceset=FlowSorted.Blood.EPIC, 
#'                                 IDOLOptimizedCpGs =IDOLOptimizedCpGs, 
#'                                 returnAll = FALSE)
#' 
#' @export
estimateCellCounts2 <- function(rgSet, compositeCellType = "Blood", 
                                processMethod = c("auto", "preprocessNoob"), 
                                probeSelect = c("auto", "any", "IDOL"), 
                                cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                              "Mono", "Neu"), 
                                referencePlatform = 
                                    c("IlluminaHumanMethylation450k", 
                                      "IlluminaHumanMethylationEPIC", 
                                      "IlluminaHumanMethylation27k"), 
                                referenceset = NULL, IDOLOptimizedCpGs = NULL, 
                                returnAll = FALSE, meanPlot = FALSE, 
                                verbose = TRUE, 
    ...) {
    minfi:::.isRGOrStop(rgSet)
    rgSet <- as(rgSet, "RGChannelSet")
    referencePlatform <- match.arg(referencePlatform)
    rgPlatform <- sub("IlluminaHumanMethylation", "", 
            annotation(rgSet)[which(names(annotation(rgSet)) == "array")])
    platform <- sub("IlluminaHumanMethylation", "", referencePlatform)
    if ((compositeCellType == "CordBlood") && (!"nRBC" %in% cellTypes)) 
message("[estimateCellCounts] Consider including 'nRBC' in argument 'cellTypes' for cord blood estimation.\n")
    #if ((compositeCellType == "Blood") && (referencePlatform == "IlluminaHumanMethylationEPIC") &&
    #    ("Gran" %in% cellTypes)) 
    #    message("[estimateCellCounts] Replace 'Gran' for 'Neu' in argument 'cellTypes' for EPIC blood estimation.\n")
    referencePkg <- sprintf("FlowSorted.%s.%s", compositeCellType, platform)
    subverbose <- max(as.integer(verbose) - 1L, 0L)
    if (!require(referencePkg, character.only = TRUE)) {
        referenceRGset <- referenceset
        if (!is(referenceset, "RGChannelSet")) 
            stop(sprintf("Could not find reference data package for compositeCellType '%s' and referencePlatform '%s' (inferred package name is '%s')", compositeCellType, platform, referencePkg))
    } else {
        data(list = referencePkg)
        referenceRGset <- get(referencePkg)
    }
    if (rgPlatform != platform) {
        rgSet <- convertArray(rgSet, outType = referencePlatform, verbose = TRUE)
    }
    if (!"CellType" %in% names(colData(referenceRGset))) 
        stop(sprintf("the reference sorted dataset (in this case '%s') needs to have a phenoData column called 'CellType'"), names(referencePkg))
    if (sum(colnames(rgSet) %in% colnames(referenceRGset)) > 0) 
        stop("the sample/column names in the user set must not be in the reference data ")
    if (!all(cellTypes %in% referenceRGset$CellType)) 
        stop(sprintf("all elements of argument 'cellTypes' needs to be part of the reference phenoData columns 'CellType' (containg the following elements: '%s')", 
                     paste(unique(referenceRGset$cellType), collapse = "', '")))
    if (length(unique(cellTypes)) < 2) 
        stop("At least 2 cell types must be provided.")
    if ((processMethod == "auto") && (compositeCellType %in% c("Blood", "DLPFC"))) 
        processMethod <- "preprocessQuantile"
    if ((processMethod == "auto") && (!compositeCellType %in% c("Blood", "DLPFC"))) 
        processMethod <- "preprocessNoob"
    processMethod <- get(processMethod)
    if ((probeSelect == "auto") && (compositeCellType == "CordBlood")) {
        probeSelect <- "any"
    }
    if ((probeSelect == "auto") && (compositeCellType != "CordBlood")) {
        probeSelect <- "both"
    }
    if (verbose) 
        message("[estimateCellCounts] Combining user data with reference (flow sorted) data.\n")
    newpd <- DataFrame(sampleNames = c(colnames(rgSet), colnames(referenceRGset)), studyIndex = rep(c("user", "reference"), times = c(ncol(rgSet), ncol(referenceRGset))), stringsAsFactors = FALSE)
    referencePd <- colData(referenceRGset)
    combinedRGset <- combineArrays(rgSet, referenceRGset, outType = referencePlatform)
    colData(combinedRGset) <- newpd
    colnames(combinedRGset) <- newpd$sampleNames
    rm(referenceRGset)
    if (verbose) 
        message("[estimateCellCounts] Processing user and reference data together.\n")
    if (compositeCellType == "CordBlood") {
        combinedMset <- processMethod(combinedRGset, verbose = subverbose)
        compTable <- get(paste0(referencePkg, ".compTable"))
        combinedMset <- combinedMset[which(rownames(combinedMset) %in% rownames(compTable)), ]
    } else {
        combinedMset <- processMethod(combinedRGset)
    }
    rm(combinedRGset)
    referenceMset <- combinedMset[, combinedMset$studyIndex == "reference"]
    colData(referenceMset) <- as(referencePd, "DataFrame")
    mSet <- combinedMset[, combinedMset$studyIndex == "user"]
    colData(mSet) <- as(colData(rgSet), "DataFrame")
    rm(combinedMset)
    if (verbose) 
        message("[estimateCellCounts] Picking probes for composition estimation.\n")
    if (probeSelect != "IDOL") {
        compData <- minfi:::pickCompProbes(referenceMset, cellTypes = cellTypes, compositeCellType = compositeCellType, probeSelect = probeSelect)
        coefs <- compData$coefEsts
        if (verbose) 
            message("[estimateCellCounts] Estimating composition.\n")
        counts <- minfi:::projectCellType(minfi::getBeta(mSet)[rownames(coefs), ], coefs)
        rownames(counts) <- colnames(rgSet)
        if (meanPlot) {
            smeans <- compData$sampleMeans
            smeans <- smeans[order(names(smeans))]
            sampleMeans <- c(colMeans(minfi::getBeta(mSet)[rownames(coefs), ]), smeans)
            sampleColors <- c(rep(1, ncol(mSet)), 1 + as.numeric(factor(names(smeans))))
            plot(sampleMeans, pch = 21, bg = sampleColors)
            legend("bottomleft", c("blood", levels(factor(names(smeans)))), col = 1:7, pch = 15)
        }
        if (returnAll) {
            list(counts = counts, compTable = compData$compTable, normalizedData = mSet)
        } else {
            list(counts = counts)  #,
            # normalizedData = mSet)
        }
    } else {
        referenceMset <- referenceset
        p <- minfi::getBeta(referenceMset)
        pd <- as.data.frame(colData(referenceMset))
        if (!is.null(cellTypes)) {
            if (!all(cellTypes %in% pd$CellType)) 
                stop("elements of argument 'cellTypes' is not part of 'referenceMset$CellType'")
            keep <- which(pd$CellType %in% cellTypes)
            pd <- pd[keep, ]
            p <- p[, keep]
        }
        pd$CellType <- factor(pd$CellType, levels = cellTypes)
        trainingProbes <- IDOLOptimizedCpGs
        p <- p[trainingProbes, ]
        pMeans <- colMeans(p)
        names(pMeans) <- pd$CellType
        form <- as.formula(sprintf("y ~ %s - 1", paste(levels(pd$CellType), collapse = "+")))
        phenoDF <- as.data.frame(model.matrix(~pd$CellType - 1))
        colnames(phenoDF) <- sub("^pd\\$CellType", "", colnames(phenoDF))
        if (ncol(phenoDF) == 2) {
            X <- as.matrix(phenoDF)
            coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
            
            coefs <- coefEsts
        } else {
            tmp <- minfi:::validationCellType(Y = p, pheno = phenoDF, modelFix = form)
            coefEsts <- tmp$coefEsts
            coefs <- coefEsts
        }
        # return(coefs)
        rm(referenceMset)
        if (verbose) 
            message("[estimateCellCounts] Estimating composition.\n")
        counts <- minfi:::projectCellType(getBeta(mSet)[rownames(coefs), ], coefs)
        rownames(counts) <- colnames(rgSet)
        if (meanPlot) {
            smeans <- compData$sampleMeans
            smeans <- smeans[order(names(smeans))]
            sampleMeans <- c(colMeans(getBeta(mSet)[rownames(coefs), ]), smeans)
            sampleColors <- c(rep(1, ncol(mSet)), 1 + as.numeric(factor(names(smeans))))
            plot(sampleMeans, pch = 21, bg = sampleColors)
            legend("bottomleft", c("blood", levels(factor(names(smeans)))), col = 1:7, pch = 15)
        }
        if (returnAll) {
            list(counts = counts, compTable = compData$compTable, normalizedData = mSet)
        } else {
            list(counts = counts)  #,, compTable = coefs
            # normalizedData = mSet)
        }
    }
}
