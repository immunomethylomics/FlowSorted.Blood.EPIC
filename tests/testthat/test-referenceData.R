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

test_that("reference resolver handles custom reference objects", {
  skip_if_not_installed("FlowSorted.Blood.450k")

  reference <- .loadReferenceObject(
    "FlowSorted.Blood.450k"
  )
  expect_identical(
    .resolveReferenceSet(reference, environment()),
    reference
  )
  named_reference <- reference
  expect_identical(
    .resolveReferenceSet("named_reference", environment()),
    named_reference
  )

  expect_error(
    .resolveReferenceSet("missing_reference", environment()),
    "was not found",
    fixed = TRUE
  )

  expect_error(
    .resolveReferenceSet(1, environment()),
    "must be an RGChannelSet",
    fixed = TRUE
  )
})
