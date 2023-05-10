library(testthat)
library(hgu133plus2.db)
library(AnnotationDbi)
library(SummarizedExperiment)
library(limma)
library(RCPA)

getTestGSData <- function(){
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

    # # control vs condition
    design <- model.matrix(~0 + group, data = colData)
    contrast <- makeContrasts("groupcondition-groupcontrol", levels = design)

    # two class paired
    # design <- model.matrix(~0 + group + pair, data = colData)
    # contrast <- makeContrasts("groupcondition-groupcontrol", levels = design.paired)
    annotation <- RCPA:::.getIDMappingAnnotation("GPL570")
    DERes <- runDEAnalysis(summarizedExperiment, method = "DESeq2", design, contrast, annotation = annotation)
    DERes
}

getTestGenesetsGA <- function(DERes){
    gs <- lapply(1:100, function(x) {
        sample(rownames(DERes), runif(1, 100, 500))
    })
    names(gs) <- paste0("geneset", 1:100)

    gs_fullNames <- paste0("path:", 1:100)
    names(gs_fullNames) <- names(gs)

    genesets <- list(
        database = "TEST",
        genesets = gs,
        names = gs_fullNames
    )
    genesets
}

test_that('ORA', {
    skip_if_offline()
    DERes <- getTestGSData()
    genesets <- getTestGenesetsGA(DERes)
    oraRes <- RCPA:::.runORA(summarizedExperiment = DERes, genesets = genesets[["genesets"]], pThreshold = 0.05)
    expect_true(all(c("ID", "p.value", "score", "normalizedScore") %in% colnames(oraRes)))
    expect_true(all(oraRes$p.value <= 1))
    expect_true(all(oraRes$p.value >= 0))
})

test_that('fgsea', {
    skip_if_offline()
    DERes <- getTestGSData()
    genesets <- getTestGenesetsGA(DERes)
    fgseaRes <- RCPA:::.runFgsea(summarizedExperiment = DERes, genesets = genesets[["genesets"]], eps = 1e-50, scoreType = "std")
    expect_true(all(c("ID", "p.value", "score", "normalizedScore") %in% colnames(fgseaRes)))
    expect_true(all(fgseaRes$p.value <= 1))
    expect_true(all(fgseaRes$p.value >= 0))
})

test_that('GSA', {
    skip_if_offline()
    DERes <- getTestGSData()
    genesets <- getTestGenesetsGA(DERes)
    gsaRes <- RCPA:::.runGSA(DERes, genesets[["genesets"]], method = "maxmean")
    expect_true(all(c("ID", "p.value", "score", "normalizedScore") %in% colnames(gsaRes)))
    expect_true(all(gsaRes$p.value <= 1))
    expect_true(all(gsaRes$p.value >= 0))
})

test_that('KS ', {
    skip_if_offline()
    DERes <- getTestGSData()
    genesets <- getTestGenesetsGA(DERes)
    ksRes <- RCPA:::.runKsWilcox(DERes, genesets[["genesets"]], sTest = "ks")
    expect_true(all(c("ID", "p.value", "score", "normalizedScore") %in% colnames(ksRes)))
    expect_true(all(ksRes$p.value <= 1))
    expect_true(all(ksRes$p.value >= 0))
})

test_that('Wilcox ', {
    skip_if_offline()
    DERes <- getTestGSData()
    genesets <- getTestGenesetsGA(DERes)
    wilcoxRes <- RCPA:::.runKsWilcox(DERes, genesets[["genesets"]], sTest = "wilcox")
    expect_true(all(c("ID", "p.value", "score", "normalizedScore") %in% colnames(wilcoxRes)))
    expect_true(all(wilcoxRes$p.value <= 1))
    expect_true(all(wilcoxRes$p.value >= 0))
})

test_that('GeneSet Enrichment Analysis with wilcox ', {
    skip_if_offline()
    DERes <- getTestGSData()
    genesets <- getTestGenesetsGA(DERes)
    result <- runGeneSetAnalysis(DERes, genesets, method = "wilcox")
    wilcoxRes <- RCPA:::.runKsWilcox(DERes, genesets[["genesets"]], sTest = "wilcox")
    expect_true(all(c("ID", "p.value", "score", "normalizedScore", "sampleSize") %in% colnames(result)))
    expect_true(all(result$p.value == wilcoxRes$p.value))
    expect_true(all(result$p.value <= 1))
    expect_true(all(result$p.value >= 0))
})

test_that('GeneSet Enrichment Analysis with fgsea ', {
    skip_if_offline()
    DERes <- getTestGSData()
    genesets <- getTestGenesetsGA(DERes)
    result <- runGeneSetAnalysis(DERes, genesets, method = "fgsea", FgseaArgs = list(eps = 1e-50, scoreType = "std", nPermSimple = 1000))
    expect_true(all(c("ID", "p.value", "score", "normalizedScore", "sampleSize") %in% colnames(result)))
    expect_true(all(result$p.value <= 1))
    expect_true(all(result$p.value >= 0))
})

test_that('GeneSet Enrichment Analysis with fgsea with wrong arguments ', {
    skip_if_offline()
    DERes <- getTestGSData()
    genesets <- getTestGenesetsGA(DERes)
    expect_error(runGeneSetAnalysis(DERes, genesets, method = "fgsea", FgseaArgs = list(scoretype = "std", nPermSimple = 1000)))
    expect_error(runGeneSetAnalysis(DERes, genesets = NULL, method = "fgsea", FgseaArgs = list(scoreType = "std", nPermSimple = 1000)))
})