library(testthat)
library(ggplot2)
library(hgu133plus2.db)
library(AnnotationDbi)

# devtools::load_all()

set.seed(1)
exprs <- round(matrix(abs(rnorm(20000 * 10, sd = 4)), nrow = 20000, ncol = 10))
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

annotation <- .getIDMappingAnnotation("GPL570")

DEResult <- runDEAnalysis(summarizedExperiment, method = "limma", design, contrast, annotation = annotation) %>% rowData()


test_that("Plot volcano plot default", {
  pl <- plotVolcanoDE(DEResult)
  expect_true(is.ggplot(pl))
  expect_equal(pl$labels$y, "-log10 pFDR")
  expect_equal(pl$labels$x, "logFC")
})

test_that("Plot volcano plot with statistic", {
  pl <- plotVolcanoDE(DEResult, xAxis = "statistic")
  expect_true(is.ggplot(pl))
  expect_equal(pl$labels$x, "statistic")
})

test_that("Plot volcano plot with p-value", {
  pl <- plotVolcanoDE(DEResult, yAxis = "-log10(p.value)")
  expect_true(is.ggplot(pl))
  expect_equal(pl$labels$y, "-log10 p-value")
})

test_that("Plot volcano plot with pFDR", {
  pl <- plotVolcanoDE(DEResult, yAxis = "-log10(pFDR)")
  expect_true(is.ggplot(pl))
  expect_equal(pl$labels$y, "-log10 pFDR")
})