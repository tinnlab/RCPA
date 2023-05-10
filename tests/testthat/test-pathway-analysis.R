library(testthat)
library(hgu133plus2.db)
library(AnnotationDbi)
library(SummarizedExperiment)
library(limma)
library(RCPA)

getTestDataPA <- function() {
    # generate a random gene expression matrix
    set.seed(123)
    exprs <- round(matrix(2^abs(rnorm(10000, sd = 4)), nrow = 1000, ncol = 10))
    rownames(exprs) <- sample(keys(hgu133plus2.db, keytype = "PROBEID"), nrow(exprs), replace = FALSE)
    colnames(exprs) <- paste0("sample", 1:10)

    controlSamples <- paste0("sample", 1:5)
    conditionSamples <- paste0("sample", 6:10)

    colData <- data.frame(
        row.names = colnames(exprs),
        group = factor(c(rep("control", length(controlSamples)), rep("condition", length(conditionSamples)))),
        pair = factor(c(seq_along(controlSamples), seq_along(conditionSamples)))
    )

    summarizedExperiment <- SummarizedExperiment(
        assays = list(counts = exprs),
        colData = colData
    )

    # control vs condition
    design <- model.matrix(~0 + group, data = colData)
    contrast <- makeContrasts("groupcondition-groupcontrol", levels = design)


    annotation <- RCPA:::.getIDMappingAnnotation("GPL570")
    DERes <- runDEAnalysis(summarizedExperiment, method = "DESeq2", design, contrast, annotation = annotation)
    DERes
}

test_that('SPIA ', {
    skip_on_cran()
    skip_if_offline()
    DERes <- getTestDataPA()
    spiaNetwork <- getSPIAKEGGNetwork("hsa", FALSE)
    spiaRes <- RCPA:::.runSPIA(summarizedExperiment = DERes, network = spiaNetwork[["network"]], pThreshold = 0.05, all = NULL)
    expect_true(all(c("ID", "p.value", "score", "normalizedScore") %in% colnames(spiaRes)))
    expect_true(all(spiaRes$p.value <= 1))
    expect_true(all(spiaRes$p.value >= 0))
})

test_that('CePaORA ', {
    skip_on_cran()
    skip_if_offline()
    DERes <- getTestDataPA()
    cepaNetwork <- getCePaPathwayCatalogue("hsa", FALSE)
    cepaOraRes <- RCPA:::.runCePaORA(DERes, cepaNetwork[["network"]], bk = NULL, iter = 1000, pThreshold = 0.05)
    expect_true(all(c("ID", "p.value", "score", "normalizedScore") %in% colnames(cepaOraRes)))
    expect_true(all(cepaOraRes$p.value <= 1))
    expect_true(all(cepaOraRes$p.value >= 0))
})

test_that('CePaGSA ', {
    skip_on_cran()
    skip_if_offline()
    DERes <- getTestDataPA()
    cepaNetwork <- getCePaPathwayCatalogue("hsa", FALSE)
    cepaGsaRes <- RCPA:::.runCePaGSA(DERes, cepaNetwork[["network"]], nlevel = "tvalue_abs", plevel = "mean", iter = 1000)
    expect_true(all(c("ID", "p.value", "score", "normalizedScore") %in% colnames(cepaGsaRes)))
    expect_true(all(cepaGsaRes$p.value <= 1))
    expect_true(all(cepaGsaRes$p.value >= 0))
})

test_that('Pathway Enrichment Analysis using SPIA', {
    skip_on_cran()
    skip_if_offline()
    DERes <- getTestDataPA()
    spiaNetwork <- getSPIAKEGGNetwork("hsa", FALSE)
    result <- RCPA::runPathwayAnalysis(DERes, spiaNetwork, method = "spia")
    expect_true(all(c("ID", "p.value", "score", "normalizedScore") %in% colnames(result)))
    expect_true(all(result$p.value <= 1))
    expect_true(all(result$p.value >= 0))
})

test_that('Pathway Enrichment Analysis using CePaORA ', {
    skip_on_cran()
    skip_if_offline()
    DERes <- getTestDataPA()
    cepaNetwork <- getCePaPathwayCatalogue("hsa", FALSE)
    result <- RCPA::runPathwayAnalysis(DERes, cepaNetwork, method = "cepaORA")
    expect_true(all(c("ID", "p.value", "score", "normalizedScore") %in% colnames(result)))
    expect_true(all(result$p.value <= 1))
    expect_true(all(result$p.value >= 0))
})

test_that('Pathway Enrichment Analysis using CePaGSA ', {
    skip_on_cran()
    skip_if_offline()
    DERes <- getTestDataPA()
    cepaNetwork <- getCePaPathwayCatalogue("hsa", FALSE)
    result <- RCPA::runPathwayAnalysis(DERes, cepaNetwork, method = "cepaGSA")
    expect_true(all(c("ID", "p.value", "score", "normalizedScore") %in% colnames(result)))
    expect_true(all(result$p.value <= 1))
    expect_true(all(result$p.value >= 0))
})