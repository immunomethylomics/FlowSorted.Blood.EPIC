cell.Mix <- which(FlowSorted.Blood.EPIC$CellType == "MIX")
RGsetTargets <- FlowSorted.Blood.EPIC[, cell.Mix]
sampleNames(RGsetTargets) <- paste(RGsetTargets$CellType,
                                   seq(along = cell.Mix), sep = "_")

countsEPIC0<-estimateCellCounts2(RGsetTargets, compositeCellType = "Blood", 
                                 processMethod = "preprocessNoob",
                                 probeSelect = "IDOL", 
                                 cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                               "Mono", "Neu"), 
                                 referencePlatform = 
                                     "IlluminaHumanMethylationEPIC",
                                 IDOLOptimizedCpGs =IDOLOptimizedCpGs, 
                                 returnAll = FALSE)

expect_is(countsEPIC0, "list")

countsEPIC<-estimateCellCounts2(RGsetTargets, compositeCellType = "Blood", 
                                processMethod = "preprocessNoob",
                                probeSelect = "IDOL", 
                                cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                              "Mono", "Neu"), 
                                referencePlatform = 
                                    "IlluminaHumanMethylationEPIC",
                                IDOLOptimizedCpGs =IDOLOptimizedCpGs, 
                                returnAll = TRUE)
expect_s4_class(countsEPIC$normalizedData, "MethylSet")
expect_is(countsEPIC$counts, "matrix")
expect_is(countsEPIC$compTable, "data.frame")

expect_error(estimateCellCounts2(RGsetTargets, compositeCellType = "Blood", 
                              processMethod = "preprocessNoob",
                              probeSelect = "IDOL", 
                              cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                            "Mono", "Gran"), 
                              referencePlatform = 
                                  "IlluminaHumanMethylationEPIC",
                              IDOLOptimizedCpGs =IDOLOptimizedCpGs, 
                              returnAll = FALSE))
expect_warning(expect_error(estimateCellCounts2(RGsetTargets, compositeCellType = "CordBlood", 
                                   processMethod = "preprocessNoob",
                                   probeSelect = "IDOL", 
                                   cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                                 "Mono", "Gran"), 
                                   referencePlatform = 
                                       "IlluminaHumanMethylationEPIC",
                                   IDOLOptimizedCpGs =IDOLOptimizedCpGs, 
                                   returnAll = FALSE)))
RGsetTargets2<-preprocessRaw(RGsetTargets)
expect_error(estimateCellCounts2(RGsetTargets2, compositeCellType = "Blood", 
                              processMethod = "preprocessNoob",
                              probeSelect = "IDOL", 
                              cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                            "Mono", "Neu"), 
                              referencePlatform = 
                                  "IlluminaHumanMethylationEPIC",
                              IDOLOptimizedCpGs =IDOLOptimizedCpGs, 
                              returnAll = FALSE), "object is of class 'MethylSet', but needs to be of class 'RGChannelSet' or 'RGChannelSetExtended' to use other methods different to 'preprocessQuantile'")
countsEPIC2<-estimateCellCounts2(RGsetTargets2, compositeCellType = "Blood", 
                    processMethod = "preprocessQuantile",
                    probeSelect = "IDOL", 
                    cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                  "Mono", "Neu"), 
                    referencePlatform = 
                        "IlluminaHumanMethylationEPIC",
                    IDOLOptimizedCpGs =IDOLOptimizedCpGs, 
                    returnAll = FALSE)
expect_is(countsEPIC2, "list")
