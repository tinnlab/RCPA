library(testthat)
library(SummarizedExperiment)
library(GEOquery)
library(RCPA)

destDir <- file.path(tempdir(), "dataFolder")
if(!dir.exists(destDir)) dir.create(destDir)

test_that('Download GEO GSE20153  GPL570', {
    skip_if_offline()
    geoDat <- RCPA:::.downloadGEOObject("GSE20153", "GPL570", destDir)
    expect_true("ExpressionSet" %in% class(geoDat))

    assay <- exprs(geoDat)
    colDat <- pData(geoDat)
    expect_true(dim(assay)[1] == 54675)
    expect_true(dim(colDat)[1] == 16)
    expect_true(dim(assay)[2] == dim(colDat)[1])
})

test_that('Download GEO GSE20153 with random GPL', {
    skip_if_offline()
    expect_error(RCPA:::.downloadGEOObject("GSE20153", "GPL12345", destDir))
})

test_that('Download GEO with random GEO', {
    skip_if_offline()
    expect_error(RCPA:::.downloadGEOObject("GSE123456789", "GPL570", destDir))
})

test_that('Download GEO GSE20153 with non existed destination', {
    skip_if_offline()
    destDir <- file.path(tempdir(), "GSE20153", "non_existed")
    expect_error(RCPA:::.downloadGEOObject("GSE20153", "GPL570", destDir), "The destination directory does not exist.")
})

test_that('Download samples of GEO dataset GSE20153 ', {
    skip_if_offline()
    result <- RCPA:::.downloadSamples(c("GSM505297", "GSM505298", "GSM505299", "GSM505300", "GSM505301"), "affymetrix", destDir)
    expect_true(result == TRUE)
})

test_that('Download not valid sample IDs', {
    skip_if_offline()
    expect_error(RCPA:::.downloadSamples(c("GSM505297123", "GSM505298111"), "affymetrix", destDir), "Check the specified samples IDs to be valid. No file is found.")
})

test_that('Test .downloadSamples with invalid protocol', {
    skip_if_offline()
    expect_error(RCPA:::.downloadSamples(c("GSM505297", "GSM505298"), "affymrix", destDir), "The specified protocol is not valid.")
})

test_that('Process affymetrix GEO dataset', {
    skip_if_offline()
    geoDat <- RCPA:::.downloadGEOObject("GSE20153", "GPL570", destDir)
    samples <- c("GSM505297", "GSM505298")
    geoMetadata <- geoDat %>% pData()
    geoMetadata <- geoMetadata[geoMetadata$geo_accession %in% samples,]
    RCPA:::.downloadSamples(samples, "affymetrix", destDir)
    processedDat <- RCPA:::.processAffymetrix(geoMetadata, samples, destDir)
    expect_true("SummarizedExperiment" %in% class(processedDat))
})

test_that('Process affymetrix GEO dataset with not matched IDs in metadata', {
    skip_if_offline()
    geoDat <- RCPA:::.downloadGEOObject("GSE20153", "GPL570", destDir)
    samples <- c("GSM505297111", "GSM505298111")
    geoMetadata <- geoDat %>% pData()
    expect_error(RCPA:::.processAffymetrix(geoMetadata, samples, destDir), "The input metadata and sampleIDs do not match. Make sure the sample IDs match with geo_accession IDs from dataset.")
})

test_that('Process agilent GEO dataset', {
    skip_if_offline()
    geoDat <- RCPA:::.downloadGEOObject("GSE22491", "GPL6480", destDir)
    samples <- c("GSM558679", "GSM558680")
    geoMetadata <- geoDat %>% pData()
    geoMetadata <- geoMetadata[geoMetadata$geo_accession %in% samples,]
    RCPA:::.downloadSamples(samples, "agilent", destDir)
    processedDat <- RCPA:::.processAgilent(geoMetadata, samples, destDir)
    expect_true("SummarizedExperiment" %in% class(processedDat))
})

test_that('Process agilent GEO dataset with not matched IDs in metadata', {
    skip_if_offline()
    geoDat <- RCPA:::.downloadGEOObject("GSE22491", "GPL6480", destDir)
    samples <- c("GSM505297111", "GSM505298111")
    geoMetadata <- geoDat %>% pData()
    expect_error(RCPA:::.processAgilent(geoMetadata, samples, destDir), "The input metadata and sampleIDs do not match. Make sure the sample IDs match with geo_accession IDs from dataset.")
})

test_that("Download and Normalize affymetrix GEO dataset", {
    skip_on_cran()
    skip_if_offline()
    summarizedExperimentObject <- downloadGEO("GSE20153", "GPL570", "affymetrix", destDir)
    assay <- SummarizedExperiment::assay(summarizedExperimentObject)
    colDat <- SummarizedExperiment::colData(summarizedExperimentObject)
    expect_true(dim(assay)[1] == 54675)
    expect_true(dim(colDat)[1] == 16)
    expect_true(dim(assay)[2] == dim(colDat)[1])
})


test_that("Download and Normalize agilent GEO dataset", {
    skip_on_cran()
    skip_if_offline()
    summarizedExperimentObject <- downloadGEO("GSE22491", "GPL6480", "agilent", destDir)
    assay <- SummarizedExperiment::assay(summarizedExperimentObject)
    colDat <- SummarizedExperiment::colData(summarizedExperimentObject)
    expect_true(dim(assay)[1] == 45015)
    expect_true(dim(colDat)[1] == 18)
    expect_true(dim(assay)[2] == dim(colDat)[1])
})

test_that("Download and Normalize agilent GEO dataset with invalid protocol", {
    expect_error(downloadGEO("GSE22491", "GPL6480", "aglnt", destDir))
})