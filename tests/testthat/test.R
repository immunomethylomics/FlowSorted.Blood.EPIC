
test_that("errors if bad parameters", {
    library(FlowSorted.Blood.EPIC)
    library(ExperimentHub)
    hub <- ExperimentHub()
    query(hub, "FlowSorted.Blood.EPIC")
    FlowSorted.Blood.EPIC1 <- hub[["EH1136"]]
    RGsetTargets <- NULL
    RGsetTargets <- FlowSorted.Blood.EPIC1[, FlowSorted.Blood.EPIC1$CellType == "MIX"]
    RGsetTargets <- RGsetTargets[, 1:3]
    sampleNames(RGsetTargets) <- paste(RGsetTargets$CellType,
        seq_len(dim(RGsetTargets)[2]),
        sep = "_"
    )
    expect_error(expect_message(estimateCellCounts2(RGsetTargets,
        compositeCellType = "Blood",
        processMethod = "preprocessNoob",
        probeSelect = "IDOL",
        cellTypes = c(
            "CD8T", "CD4T", "NK", "Bcell",
            "Mono", "Gran"
        ),
        referencePlatform =
            "IlluminaHumanMethylationEPIC",
        referenceset = NULL,
        IDOLOptimizedCpGs = IDOLOptimizedCpGs,
        returnAll = FALSE
    )))
    expect_warning(expect_error(expect_message(estimateCellCounts2(RGsetTargets,
        compositeCellType = "CordBlood",
        processMethod = "preprocessNoob",
        probeSelect = "IDOL",
        cellTypes = c(
            "CD8T", "CD4T", "NK", "Bcell",
            "Mono", "Gran"
        ),
        referencePlatform =
            "IlluminaHumanMethylationEPIC",
        IDOLOptimizedCpGs = IDOLOptimizedCpGs,
        returnAll = FALSE
    ))))
})

test_that("errors if bad parameters", {
    library(FlowSorted.Blood.EPIC)
    library(ExperimentHub)
    hub <- ExperimentHub()
    query(hub, "FlowSorted.Blood.EPIC")
    FlowSorted.Blood.EPIC1 <- hub[["EH1136"]]
    RGsetTargets <- NULL
    RGsetTargets <- FlowSorted.Blood.EPIC1[, FlowSorted.Blood.EPIC1$CellType == "MIX"]
    RGsetTargets <- RGsetTargets[, 1:3]
    sampleNames(RGsetTargets) <- paste(RGsetTargets$CellType,
        seq_len(dim(RGsetTargets)[2]),
        sep = "_"
    )
    RGsetTargets2 <- preprocessRaw(RGsetTargets)
    expect_error(estimateCellCounts2(RGsetTargets2,
        compositeCellType = "Blood",
        processMethod = "preprocessNoob",
        probeSelect = "IDOL",
        cellTypes = c(
            "CD8T", "CD4T", "NK", "Bcell",
            "Mono", "Neu"
        ),
        referencePlatform =
            "IlluminaHumanMethylationEPIC",
        referenceset = NULL,
        IDOLOptimizedCpGs = IDOLOptimizedCpGs,
        returnAll = FALSE
    ), "object is of class 'MethylSet', but needs to be of class 'RGChannelSet' or 'RGChannelSetExtended' to use other methods different to 'preprocessQuantile'")
})
