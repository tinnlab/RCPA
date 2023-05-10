library(testthat)
# library(ggplot2)
library(hgu133plus2.db)
library(AnnotationDbi)
library(RCPA)

getTestVolcanoData <- function(){
  set.seed(1)
  exprs <- round(matrix(abs(rnorm(2000 * 10, sd = 4)), nrow = 2000, ncol = 10))
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

  DEResult <- runDEAnalysis(summarizedExperiment, method = "limma", design, contrast, annotation = annotation) %>% rowData()

  DEResult
}

test_that("Plot volcano plot default", {
  skip_if_offline()
  DEResult <- getTestVolcanoData()
  pl <- plotVolcanoDE(DEResult)
  expect_true(is.ggplot(pl))
  expect_equal(pl$labels$y, "-log10 pFDR")
  expect_equal(pl$labels$x, "log2 fold change")
})

test_that("Plot volcano plot with p-value", {
  skip_if_offline()
  DEResult <- getTestVolcanoData()
  pl <- plotVolcanoDE(DEResult, useFDR = FALSE)
  expect_true(is.ggplot(pl))
  expect_equal(pl$labels$y, "-log10 p-value")
})

test_that("Plot volcano plot with pFDR", {
  skip_if_offline()
  DEResult <- getTestVolcanoData()
  pl <- plotVolcanoDE(DEResult, useFDR = TRUE)
  expect_true(is.ggplot(pl))
  expect_equal(pl$labels$y, "-log10 pFDR")
})