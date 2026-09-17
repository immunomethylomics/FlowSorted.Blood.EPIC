test_that("validationCellType fits ordinary least squares rows", {
    pheno <- data.frame(
        group = factor(rep(c("A", "B", "C"), each = 3))
    )
    modelFix <- y ~ group - 1
    Y <- rbind(
        c(1.0, 1.2, 0.8, 2.0, 2.1, 1.9, 3.0, 3.2, 2.8),
        c(1.5, 1.7, NA, 2.5, 2.6, 2.4, 3.5, 3.7, 3.3),
        c(0.5, 0.7, 0.6, 1.5, 1.4, 1.6, 2.5, 2.7, 2.6)
    )

    observed <- validationCellType(
        Y,
        pheno,
        modelFix,
        verbose = FALSE
    )
    fits <- lapply(seq_len(nrow(Y)), function(ii) {
        data <- pheno
        data$y <- Y[ii, ]
        stats::lm(modelFix, data = data[!is.na(data$y), ])
    })

    expect_equal(
        observed$coefEsts,
        do.call(rbind, lapply(fits, coef))
    )
    expect_equal(
        observed$sigmaResid,
        vapply(fits, function(fit) summary(fit)$sigma, numeric(1))
    )
    expect_true(all(observed$sigmaIcept == 0))
    expect_true(all(observed$nClusters == 0))
    expect_null(observed$modelBatch)
    expect_equal(
        observed$nObserved,
        rowSums(!is.na(Y))
    )
    expect_length(observed$coefVcovs, nrow(Y))
    expect_length(observed$Fstat, nrow(Y))
    expect_length(observed$Pval, nrow(Y))
    expect_length(observed$degFree, nrow(Y))
})
