library(testthat)
library(ggplot2)
library(hgu133plus2.db)
library(AnnotationDbi)

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

test_that("Plot MA default", {
    pl <- RCPA::plotMA(DEResult)
    expect_true(is.ggplot(pl))
    expect_equal(pl$labels$y, "Log2 fold change")
    expect_equal(pl$labels$x, "Average expression")
})

test_that("Plot MA with non FDR", {
    pl <- RCPA::plotMA(DEResult, useFDR = F)
    expect_true(is.ggplot(pl))
    expect_equal(pl$labels$y, "Log2 fold change")
    expect_equal(pl$labels$x, "Average expression")
})

test_that("Plot MA with FDR", {
    pl <- RCPA::plotMA(DEResult, useFDR = T)
    expect_true(is.ggplot(pl))
    expect_equal(pl$labels$y, "Log2 fold change")
    expect_equal(pl$labels$x, "Average expression")
})

test_that("Plot MA with labels", {
    labels <- c("label1", "label2", "label3", "label4", "label5", "label6", "label7", "label8", "label9", "label10")
    names(labels) <- sample(DEResult$ID, length(labels), replace = F)
    pl <- RCPA::plotMA(DEResult, labels = labels)
    expect_true(is.ggplot(pl))
    expect_equal(pl$labels$y, "Log2 fold change")
    expect_equal(pl$labels$x, "Average expression")
})