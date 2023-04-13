library(testthat)
library(SummarizedExperiment)

devtools::load_all()

destDir <- file.path(getwd(), ".tmp")

test_that('Download GEO GSE20153  GPL570', {
    geoDat <- .downloadGEOObject("GSE20153", "GPL570", destDir)
    expect_true("ExpressionSet" %in% class(geoDat))

    assay <- exprs(geoDat)
    colDat <- pData(geoDat)
    expect_true(dim(assay)[1] == 54675)
    expect_true(dim(colDat)[1] == 16)
    expect_true(dim(assay)[2] == dim(colDat)[1])
})

test_that('Download GEO GSE20153 with random GPL', {
    expect_error(.downloadGEOObject("GSE20153", "GPL12345", destDir))
})

test_that('Download GEO with random GEO', {
    expect_error(.downloadGEOObject("GSE123456789", "GPL570", destDir))
})

test_that('Download GEO GSE20153 with non existed destination', {
    destDir <- file.path(getwd(), ".tmp", "non_existed")
    expect_error(.downloadGEOObject("GSE20153", "GPL570", destDir), "The destination directory does not exist.")
})


test_that("Download and Normalize affymetrix GEO dataset", {
    destDir <- file.path(getwd(), ".tmp")
    summarizedExperimentObject <- downloadGEO("GSE20153", "GPL570", "affymetrix", destDir)
    assay <- SummarizedExperiment::assay(summarizedExperimentObject)
    colDat <- SummarizedExperiment::colData(summarizedExperimentObject)
    expect_true(dim(assay)[1] == 54675)
    expect_true(dim(colDat)[1] == 16)
    expect_true(dim(assay)[2] == dim(colDat)[1])
})


test_that("Download and Normalize agilent GEO dataset", {
    summarizedExperimentObject <- downloadGEO("GSE22491", "GPL6480", "agilent", getwd())
    assay <- SummarizedExperiment::assay(summarizedExperimentObject)
    colDat <- SummarizedExperiment::colData(summarizedExperimentObject)
    expect_true(dim(assay)[1] == 45015)
    expect_true(dim(colDat)[1] == 18)
    expect_true(dim(assay)[2] == dim(colDat)[1])
})

test_that("Download and Normalize RNASeq GEO dataset", {
    summarizedExperimentObject <- downloadGEO("GSE165082", "GPL11154", "RNASeq", getwd())
    assay <- SummarizedExperiment::assay(summarizedExperimentObject)
    colDat <- SummarizedExperiment::colData(summarizedExperimentObject)
    expect_true(dim(assay)[1] == 63677)
    expect_true(dim(colDat)[1] == 26)
    expect_true(dim(assay)[2] == dim(colDat)[1])
})