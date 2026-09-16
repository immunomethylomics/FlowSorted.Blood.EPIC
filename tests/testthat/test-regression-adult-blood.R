test_that("adult EPIC IDOL estimates remain numerically stable", {
  reference_epic <- libraryDataGet("FlowSorted.Blood.EPIC")
  
  target_mixtures <- reference_epic[
    ,
    reference_epic$CellType == "MIX"
  ]
  
  target_mixtures <- target_mixtures[, seq_len(3)]
  
  minfi::sampleNames(target_mixtures) <- paste0(
    "MIX_",
    seq_len(ncol(target_mixtures))
  )
  
  observed <- estimateCellCounts2(
    target_mixtures,
    compositeCellType = "Blood",
    processMethod = "preprocessNoob",
    probeSelect = "IDOL",
    cellTypes = c(
      "CD8T",
      "CD4T",
      "NK",
      "Bcell",
      "Mono",
      "Neu"
    ),
    referencePlatform = "IlluminaHumanMethylationEPIC",
    returnAll = FALSE,
    verbose = FALSE
  )
  
  expected <- structure(
    c(
      0.1915, 0.0464, 0.0672,
      0.0704, 0.1757, 0.1010,
      0.1517, 0.0181, 0.0047,
      0.1906, 0.0444, 0.0225,
      0.1907, 0.0584, 0.1093,
      0.2114, 0.6699, 0.7034
    ),
    dim = c(3L, 6L),
    dimnames = list(
      paste0("MIX_", seq_len(3)),
      c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")
    )
  )
  
  expect_equal(
    observed$prop,
    expected,
    tolerance = 1e-4
  )
  
  expect_equal(
    rowSums(observed$prop),
    c(MIX_1 = 1.0063, MIX_2 = 1.0129, MIX_3 = 1.0081),
    tolerance = 1e-4
  )
})