library(testthat)
library(SummarizedExperiment)

test_that("Download and Normalize affymetrix GEO dataset", {
  summarizedExperimentObject <- downloadGEO("GSE20153", "GPL570", "affymetrix", getwd())
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