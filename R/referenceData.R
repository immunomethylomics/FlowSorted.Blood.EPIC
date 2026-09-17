.loadReferenceObject <- function(packageName) {
  if (
    !is.character(packageName) ||
    length(packageName) != 1L ||
    is.na(packageName) ||
    !nzchar(packageName)
  ) {
    stop(
      "`packageName` must be a single, non-empty character string.",
      call. = FALSE
    )
  }

  if (!requireNamespace(packageName, quietly = TRUE)) {
    message <- sprintf(
      paste0(
        "Could not find reference data package '%s'. ",
        "Install the package before using this reference."
      ),
      packageName
    )

    if (identical(
      packageName,
      "FlowSorted.BloodExtended.EPIC"
    )) {
      message <- paste(
        message,
        paste0(
          "Contact Technology.Transfer@dartmouth.edu ",
          "for access."
        )
      )
    }

    stop(message, call. = FALSE)
  }

  experimentHubReferences <- c(
    "FlowSorted.Blood.EPIC",
    "FlowSorted.CordBloodCombined.450k"
  )

  if (packageName %in% experimentHubReferences) {
    return(libraryDataGet(packageName))
  }

  if (packageName %in% getNamespaceExports(packageName)) {
    return(getExportedValue(packageName, packageName))
  }

  dataEnvironment <- new.env(parent = emptyenv())

  suppressWarnings(
    utils::data(
      list = packageName,
      package = packageName,
      envir = dataEnvironment
    )
  )

  if (exists(
    packageName,
    envir = dataEnvironment,
    inherits = FALSE
  )) {
    return(
      get(
        packageName,
        envir = dataEnvironment,
        inherits = FALSE
      )
    )
  }

  stop(
    sprintf(
      paste0(
        "Reference object '%s' was not found ",
        "in package '%s'."
      ),
      packageName,
      packageName
    ),
    call. = FALSE
  )
}

.resolveReferenceSet <- function(referenceset, envir) {
    if (
        is.character(referenceset) &&
        length(referenceset) == 1L &&
        !is.na(referenceset) &&
        nzchar(referenceset)
    ) {
        if (!exists(
            referenceset,
            envir = envir,
            inherits = TRUE
        )) {
            stop(
                sprintf(
                    "Custom reference object '%s' was not found.",
                    referenceset
                ),
                call. = FALSE
            )
        }

        referenceset <- get(
            referenceset,
            envir = envir,
            inherits = TRUE
        )
    }

    if (
        !is(referenceset, "RGChannelSet") &&
        !is(referenceset, "RGChannelSetExtended")
    ) {
        stop(
            paste0(
                "`referenceset` must be an RGChannelSet, ",
                "an RGChannelSetExtended, or a single character ",
                "string naming one of these objects."
            ),
            call. = FALSE
        )
    }

    referenceset
}
