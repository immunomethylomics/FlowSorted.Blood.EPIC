test_that("reference loader validates package names", {
  expect_error(
    .loadReferenceObject(NULL),
    "single, non-empty character string",
    fixed = TRUE
  )

  expect_error(
    .loadReferenceObject(""),
    "single, non-empty character string",
    fixed = TRUE
  )
})

test_that("reference loader reports unavailable packages", {
  expect_error(
    .loadReferenceObject(
      "FlowSorted.Package.That.Does.Not.Exist"
    ),
    "Could not find reference data package",
    fixed = TRUE
  )
})

test_that("reference loader retrieves exported reference objects", {
  skip_if_not_installed("FlowSorted.Blood.450k")

  reference <- .loadReferenceObject(
    "FlowSorted.Blood.450k"
  )

  expect_s4_class(reference, "RGChannelSet")
})

test_that("reference loader retrieves traditional data objects", {
  skip_if_not_installed(
    "FlowSorted.CordBloodNorway.450k"
  )

  reference <- .loadReferenceObject(
    "FlowSorted.CordBloodNorway.450k"
  )

  expect_s4_class(reference, "RGChannelSet")
})
