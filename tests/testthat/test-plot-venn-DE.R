library(testthat)
library(hgu133plus2.db)
library(AnnotationDbi)
library(SummarizedExperiment)
library(limma)
library(ggplot2)
library(RCPA)

getTestVennData <- function(){
    # generate a random gene expression matrix
    DEResults <- lapply(1:3, function(seed) {
        set.seed(seed)
        exprs <- round(matrix(abs(rnorm(20000 * 10, sd = 4)), nrow = 20000, ncol = 10))
        rownames(exprs) <- sample(keys(hgu133plus2.db, keytype = "PROBEID"), nrow(exprs), replace = FALSE)
        colnames(exprs) <- paste0("sample", 1:10)

        controlSamples <- paste0("sample", 1:5)
        conditionSamples <- paste0("sample", 6:10)

        exprs[, conditionSamples] <- exprs[, conditionSamples] + 2*sample(c(1,-1), nrow(exprs), replace = TRUE)

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

test_that("Plot venn diagram default", {
    skip_if_offline()
    DEResults <- getTestVennData()
    pl <- plotVennDE(DEResults)
    expect_true(is.ggplot(pl))

    plb <- ggplot_build(pl)
    expect_true(length(plb$data) == 4)
})

test_that("Plot venn diagram with logFC", {
    skip_if_offline()
    DEResults <- getTestVennData()
    pl <- plotVennDE(DEResults, stat = "logFC", statThreshold = 10)
    expect_true(is.ggplot(pl))

    plb <- ggplot_build(pl)
    expect_true(length(plb$data) == 4)
})

test_that("Plot venn diagram with p-value", {
    skip_if_offline()
    DEResults <- getTestVennData()
    pl <- plotVennDE(DEResults, useFDR = F)
    expect_true(is.ggplot(pl))

    plb <- ggplot_build(pl)
    expect_true(length(plb$data) == 4)
})

test_that("Plot venn diagram with FDR", {
    skip_if_offline()
    DEResults <- getTestVennData()
    pl <- plotVennDE(DEResults, useFDR = T)
    expect_true(is.ggplot(pl))

    plb <- ggplot_build(pl)
    expect_true(length(plb$data) == 4)
})

test_that("Plot venn diagram with random stat", {
    skip_if_offline()
    DEResults <- getTestVennData()
    expect_error(plotVennDE(DEResults, stat = "logFC2"))
})