library(testthat)
library(hgu133plus2.db)
library(AnnotationDbi)
library(SummarizedExperiment)
library(limma)
library(ggplot2)
library(RCPA)

getTestKEGGMapData <- function(){
    DEResults <- lapply(1:3, function(seed) {
        set.seed(seed)
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

        runDEAnalysis(summarizedExperiment, method = "limma", design, contrast, annotation = annotation) %>% rowData()
    })

    DEResults
}

test_that("Plot KEGG map of path:hsa05200", {
    skip_if_offline()
    DEResults <- getTestKEGGMapData()
    res <- plotKEGGMap(DEResults, "path:hsa05200", statistic = "logFC", useFDR = TRUE, pThreshold = 0.05)
    expect_true(all(names(res) %in% c("plot", "width", "height")))
    expect_true(is.ggplot(res$plot))
    expect_true(res$width > 0)
    expect_true(res$height > 0)
})

test_that("Plot KEGG map of path:hsa05200 using p-value", {
    skip_if_offline()
    DEResults <- getTestKEGGMapData()
    res <- plotKEGGMap(DEResults, "path:hsa05200", statistic = "p.value", useFDR = FALSE, pThreshold = 0.05)
    expect_true(all(names(res) %in% c("plot", "width", "height")))
    expect_true(is.ggplot(res$plot))
    expect_true(res$width > 0)
    expect_true(res$height > 0)
})

test_that("Plot KEGG map of random pathway", {
    skip_if_offline()
    DEResults <- getTestKEGGMapData()
    expect_error(plotKEGGMap(DEResults, "path:hsa05201", statistic = "logFC", useFDR = TRUE, pThreshold = 0.05))
})