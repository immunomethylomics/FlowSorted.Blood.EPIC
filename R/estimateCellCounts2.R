#' estimateCellCounts2 function allows the use of customized reference datasets 
#' and IDOL probes L-DMR lists
#' @import minfi
#' @importFrom  graphics legend
#' @importFrom  graphics plot
#' @importFrom  stats as.formula
#' @importFrom  stats model.matrix
#' @importFrom  stats  lm
#' @importFrom  stats  pf
#' @importFrom  stats  vcov
#' @importFrom  utils data
#' @importFrom genefilter rowFtests
#' @importFrom genefilter rowttests
#' @importFrom matrixStats rowRanges
#' @importFrom quadprog solve.QP
#' @importFrom nlme lme
#' @importFrom nlme getVarCov
#' @importClassesFrom S4Vectors DataFrame
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
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
#' # array.  This object is nothing more than a vector of length 450 consisting 
#' # of the names of the IDOL optimized CpGs.  These CpGs are used as the 
#' # backbone for deconvolution and were selected because their methylation 
#' # signature differs across the six normal leukocyte subtypes.
#' data(IDOLOptimizedCpGs)
#' # Step 4: Deconvolute a target data set consisting of EPIC DNA methylation 
#' # data profiled in blood, using your prefered method estimateCellCounts 
#' # (minfi), or similar.
#' # We recommend using Noob processMethod = "preprocessNoob" in minfi
#' # for the target and reference datasets Cell types included are 
#' # cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")
#' # We also provide an IDOL optimized list of CpGs (IDOLOptimizedCpGs)
#' # that can be used for optimized cell estimations
#' # We offer an adaptation of the popular estimateCellCounts in minfi to allow 
#' # the inclusion of customized reference arrays. 
#'  countsEPIC<-estimateCellCounts2(RGsetTargets, compositeCellType = "Blood", 
#'                                 processMethod = "preprocessNoob",
#'                                 probeSelect = "IDOL", 
#'                                 cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
#'                                 "Mono", "Neu"), 
#'                                 referencePlatform = 
#'                                 "IlluminaHumanMethylationEPIC",
#'                                 IDOLOptimizedCpGs =IDOLOptimizedCpGs, 
#'                                 returnAll = FALSE)
#' @param 
#' rgSet           The input RGChannelSet for the procedure.
#' @param
#' compositeCellType   Which composite cell type is being deconvoluted. 
#'                      Should be one of "Blood", "CordBlood", or "DLPFC".
#'                      See details.
#' @param
#' processMethod How should the user and reference data be processed together? 
#'                Default input, "preprocessNoob" in minfi, you can use "auto" 
#'                for preprocessQuantile for Blood and DLPFC in 450K datasets 
#'                and preprocessNoob otherwise, according to existing 
#'                literature. Set it to any preprocessing function as a 
#'                character if you want to override it, like "preprocessFunnorm"
#' @param 
#' probeSelect	 How should probes be selected to distinguish cell types? 
#'                Options include "IDOL", for using a customized set of probes 
#'                obtained from IDOL optimization, "both", which selects an 
#'                equal number (50) of probes (with F-stat p-value < 1E-8) 
#'                with the greatest magnitude of effect from the hyper- and 
#'                hypo-methylated sides, and "any", which selects the 100 probes 
#'                (with F-stat p-value < 1E-8) with the greatest magnitude of 
#'                difference regardless of direction of effect. This according 
#'                to minfi algorithm. Default input "auto" in minfi will use 
#'                "any" for cord blood and "both" otherwise. Please see minfi 
#'                estimateCellCounts references for more details.
#' @param 
#' cellTypes     A vector of length K that contains the cell type names.  
#'                Default: c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu").
#'                Please notice that this platform use Neutrophils instead 
#'                of Granulocytes. See details.
#' @param                 
#' referencePlatform The platform for the reference dataset; if the input rgSet 
#'                    belongs to another platform, it will be converted using 
#'                    minfi function convertArray.
#' @param
#' referenceset A custom reference rgset if it is not a package installed
#' @param
#' IDOLOptimizedCpGs a vector of probe names for cell deconvolution
#' @param 
#' returnAll	Should the composition table and the normalized user supplied 
#'              data be return? Default is False.
#' @param              
#' verbose Should the function be verbose?
#' @param
#' meanPlot	Whether to plots the average DNA methylation across the cell-type 
#'          discriminating probes within the mixed and sorted samples.
#' @param 
#' ...  Other arguments for preprocessquantile
#'@return
#' This function will return a list containing matrix ofcell count (counts), if
#' returnAll=FALSE, or a list containing the counts, mean methylation per 
#' cellType, and the normalized betas (if returnAll is set to TRUE). These 
#' objects are important if you decide to use a different deconvolution 
#' algorithm such as CIBERSORT or RPM.
#' 
#' @export
estimateCellCounts2 <- function(rgSet, compositeCellType = "Blood", 
                                processMethod = "preprocessNoob", 
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
    .isRGOrStop(rgSet)
    rgSet <- as(rgSet, "RGChannelSet")
    referencePlatform <- match.arg(referencePlatform)
    rgPlatform <- sub("IlluminaHumanMethylation", "", 
            annotation(rgSet)[which(names(annotation(rgSet)) == "array")])
    platform <- sub("IlluminaHumanMethylation", "", referencePlatform)
    if ((compositeCellType == "CordBlood") && (!"nRBC" %in% cellTypes)) 
message("[estimateCellCounts] Consider including 'nRBC' in argument 'cellTypes' for cord blood estimation.\n")
    if ((compositeCellType == "Blood") && (referencePlatform == "IlluminaHumanMethylationEPIC") && ("Gran" %in% cellTypes)) 
     message("[estimateCellCounts] Replace 'Gran' for 'Neu' in argument 'cellTypes' for EPIC blood estimation.\n")
    referencePkg <- sprintf("FlowSorted.%s.%s", compositeCellType, platform)
    subverbose <- max(as.integer(verbose) - 1L, 0L)
    if (!is.null(referenceset)){
        referenceRGset<- get(referenceset)
        .isRGOrStop(referenceRGset)} else{ 
    if (!require(referencePkg, character.only = TRUE)) 
            stop(sprintf("Could not find reference data package for compositeCellType '%s' and referencePlatform '%s' (inferred package name is '%s')", compositeCellType, platform, referencePkg))
        #data(list = referencePkg)
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
        compData <- pickCompProbes(referenceMset, cellTypes = cellTypes, compositeCellType = compositeCellType, probeSelect = probeSelect)
        coefs <- compData$coefEsts
        if (verbose) 
            message("[estimateCellCounts] Estimating composition.\n")
        counts <- projectCellType(getBeta(mSet)[rownames(coefs), ], coefs)
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
            list(counts = counts)  
        }
    } else {
        p <- getBeta(referenceMset)
        pd <- as.data.frame(colData(referenceMset))
        if (!is.null(cellTypes)) {
            if (!all(cellTypes %in% pd$CellType)) 
                stop("elements of argument 'cellTypes' is not part of 'referenceMset$CellType'")
            keep <- which(pd$CellType %in% cellTypes)
            pd <- pd[keep, ]
            p <- p[, keep]
        }
        pd$CellType <- factor(pd$CellType, levels = cellTypes)
        ffComp <- rowFtests(p, pd$CellType)
        prof <- sapply(splitit(pd$CellType), function(i) rowMeans(p[, 
                                                                    i]))
        r <- rowRanges(p)
        compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 2]))
        names(compTable)[1] <- "Fstat"
        names(compTable)[c(-2, -1, 0) + ncol(compTable)] <- c("low", 
                                                              "high", "range")
        tIndexes <- splitit(pd$CellType)
        tstatList <- lapply(tIndexes, function(i) {
            x <- rep(0, ncol(p))
            x[i] <- 1
            return(rowttests(p, factor(x)))
        })
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
            tmp <- validationCellType(Y = p, pheno = phenoDF, modelFix = form)
            coefEsts <- tmp$coefEsts
            coefs <- coefEsts
        }
        rm(referenceMset)
        if (verbose) 
            message("[estimateCellCounts] Estimating composition.\n")
        counts <- projectCellType(getBeta(mSet)[rownames(coefs), ], coefs)
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
            list(counts = counts, compTable = compTable, normalizedData = mSet)
        } else {
            list(counts = counts)  
        }
    }
}

#These are a minfi internal functions here they are called to keep the function above running
pickCompProbes <- function(mSet, cellTypes = NULL, numProbes = 50, compositeCellType = compositeCellType, probeSelect = probeSelect) {
    splitit <- function(x) {
        split(seq(along=x), x)
    }
    
    p <- getBeta(mSet)
    pd <- as.data.frame(colData(mSet))
    if(!is.null(cellTypes)) {
        if(!all(cellTypes %in% pd$CellType))
            stop("elements of argument 'cellTypes' is not part of 'mSet$CellType'")
        keep <- which(pd$CellType %in% cellTypes)
        pd <- pd[keep,]
        p <- p[,keep]
    }
    ## make cell type a factor 
    pd$CellType <- factor(pd$CellType, levels = cellTypes)
    ffComp <- rowFtests(p, pd$CellType)
    prof <- sapply(splitit(pd$CellType), function(i) rowMeans(p[,i]))
    r <- matrixStats::rowRanges(p)
    compTable <- cbind(ffComp, prof, r, abs(r[,1] - r[,2]))
    names(compTable)[1] <- "Fstat"
    names(compTable)[c(-2,-1,0) + ncol(compTable)] <- c("low", "high", "range") 
    tIndexes <- splitit(pd$CellType)
    tstatList <- lapply(tIndexes, function(i) {
        x <- rep(0,ncol(p))
        x[i] <- 1
        return(rowttests(p, factor(x)))
    })
    
    if (probeSelect == "any"){
        probeList <- lapply(tstatList, function(x) {
            y <- x[x[,"p.value"] < 1e-8,]
            yAny <- y[order(abs(y[,"dm"]), decreasing=TRUE),]      
            c(rownames(yAny)[1:(numProbes*2)])
        })
    } else {
        probeList <- lapply(tstatList, function(x) {
            y <- x[x[,"p.value"] < 1e-8,]
            yUp <- y[order(y[,"dm"], decreasing=TRUE),]
            yDown <- y[order(y[,"dm"], decreasing=FALSE),]
            c(rownames(yUp)[1:numProbes], rownames(yDown)[1:numProbes])
        })
    }
    
    trainingProbes <- unique(unlist(probeList))
    p <- p[trainingProbes,]
    
    pMeans <- colMeans(p)
    names(pMeans) <- pd$CellType
    
    form <- as.formula(sprintf("y ~ %s - 1", paste(levels(pd$CellType), collapse="+")))
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

projectCellType <- function(Y, coefCellType, contrastCellType=NULL, nonnegative=TRUE, lessThanOne=FALSE){ 
    if(is.null(contrastCellType))
        Xmat <- coefCellType
    else
        Xmat <- tcrossprod(coefCellType, contrastCellType) 
    
    nCol <- dim(Xmat)[2]
    if(nCol == 2) {
        Dmat <- crossprod(Xmat)
        mixCoef <- t(apply(Y, 2, function(x) { solve(Dmat, crossprod(Xmat, x)) }))
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
            for(i in 1:nSubj) {
                obs <- which(!is.na(Y[,i])) 
                Dmat <- crossprod(Xmat[obs,])
                mixCoef[i,] <- solve.QP(Dmat, crossprod(Xmat[obs,], Y[obs,i]), Amat, b0vec)$sol
            }
        } else {
            for(i in 1:nSubj) {
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
    
    ## Initialize various containers
    sigmaResid <- sigmaIcept <- nObserved <- nClusters <- Fstat <- rep(NA, M)
    coefEsts <- matrix(NA, M, sizeModel)
    coefVcovs <- list()
    
    if(verbose)
        cat("[validationCellType] ")
    for(j in 1:M) { # For each CpG
        ## Remove missing methylation values
        ii <- !is.na(Y[j,])
        nObserved[j] <- sum(ii)
        pheno$y <- Y[j,]
        
        if(j%%round(M/10)==0 && verbose)
            cat(".") # Report progress
        
        try({ # Try to fit a mixed model to adjust for plate
            if(!is.null(modelBatch)) {
                fit <- try(lme(modelFix, random=modelBatch, data=pheno[ii,]))
                OLS <- inherits(fit,"try-error") # If LME can't be fit, just use OLS
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
    
    out <- list(coefEsts=coefEsts, coefVcovs=coefVcovs, modelFix=modelFix, modelBatch=modelBatch,
                sigmaIcept=sigmaIcept, sigmaResid=sigmaResid, L.forFstat=L.forFstat, Pval=Pval,
                orderFstat=order(-Fstat), Fstat=Fstat, nClusters=nClusters, nObserved=nObserved,
                degFree=degFree)
    
    out
}

.isRGOrStop<-function (object) {
    if (!is(object, "RGChannelSet")) 
        stop(sprintf("object is of class '%s', but needs to be of class 'RGChannelSet' or 'RGChannelSetExtended'", 
                     class(object)))
}

splitit <- function(x) {
    split(seq(along=x), x)
}