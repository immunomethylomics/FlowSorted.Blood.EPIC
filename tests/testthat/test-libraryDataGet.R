test_that("libraryDataGet validates title", {
  expect_error(
    libraryDataGet(character()),
    "single, non-empty character string",
    fixed = TRUE
  )
  
  expect_error(
    libraryDataGet(c("first", "second")),
    "single, non-empty character string",
    fixed = TRUE
  )
  
  expect_error(
    libraryDataGet(NA_character_),
    "single, non-empty character string",
    fixed = TRUE
  )
  
  expect_error(
    libraryDataGet(""),
    "single, non-empty character string",
    fixed = TRUE
  )
})

test_that("libraryDataGet returns the adult blood reference", {
  observed <- libraryDataGet("FlowSorted.Blood.EPIC")
  
  expect_s4_class(observed, "RGChannelSet")
  expect_identical(ncol(observed), 49L)
  expect_true("CellType" %in% names(SummarizedExperiment::colData(observed)))
})