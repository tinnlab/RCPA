library(testthat)
library(hgu133plus2.db)
library(AnnotationDbi)
library(SummarizedExperiment)
library(limma)
library(RCPA)

# generate a random gene expression matrix
set.seed(123)
exprs <- round(matrix(2^abs(rnorm(1000, sd = 4)), nrow = 100, ncol = 10))
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

getTestAnnotation <- function(){
    annotation <- RCPA:::.getIDMappingAnnotation("GPL570")
    annotation
}

test_that("DE analysis with ID mapping", {
    skip_if_offline()
    annotation <- getTestAnnotation()
    DERes <- runDEAnalysis(summarizedExperiment, method = "DESeq2", design, contrast, annotation = annotation)

    expect_true(all(c("PROBEID", "ID", "p.value", "pFDR", "statistic", "logFC", "logFCSE", "sampleSize") %in% colnames(rowData(DERes))))
    expect_true(all(
        c("DEAnalysis.method", "DEAnalysis.design", "DEAnalysis.contrast", "DEAnalysis.mapping") %in% names(metadata(DERes))
    ))
    expect_true(all(metadata(DERes)$DEAnalysis.method == "DESeq2"))
    expect_true(all(metadata(DERes)$DEAnalysis.design == design))
    expect_true(all(metadata(DERes)$DEAnalysis.contrast == contrast))
})

test_that("DE analysis without ID mapping and without platform", {
    expect_error(runDEAnalysis(summarizedExperiment, method = "limma", design, contrast, annotation = NULL))
})

test_that("DE analysis with platform and without ID mapping", {
    metadata(summarizedExperiment)$platform <- "GPL570"
    DERes <- runDEAnalysis(summarizedExperiment, method = "limma", design, contrast, annotation = NULL)
    expect_true(all(c("PROBEID", "ID", "p.value", "pFDR", "statistic", "logFC", "logFCSE", "sampleSize") %in% colnames(rowData(DERes))))
})


