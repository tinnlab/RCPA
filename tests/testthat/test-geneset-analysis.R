library(testthat)
library(hgu133plus2.db)
library(AnnotationDbi)
library(SummarizedExperiment)
library(limma)

devtools::load_all()

# generate a random gene expression matrix
set.seed(123)
exprs <- round(matrix(2^abs(rnorm(100000, sd = 4)), nrow = 10000, ncol = 10))
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

annotation <- .getIDMappingAnnotation("GPL570")
DERes <- runDEAnalysis(summarizedExperiment, method = "DESeq2", design, contrast, annotation = annotation)

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

test_that('ORA', {
    oraRes <- .runORA(summarizedExperiment = DERes, genesets = genesets[["genesets"]], pThreshold = 0.05)
    expect_true(all(c("ID", "p.value", "score", "normalizedScore") %in% colnames(oraRes)))
    expect_true(all(oraRes$p.value <= 1))
    expect_true(all(oraRes$p.value >= 0))
})

test_that('fgsea', {
    fgseaRes <- .runFgsea(summarizedExperiment = DERes, genesets = genesets[["genesets"]], eps = 1e-50, scoreType = "std")
    expect_true(all(c("ID", "p.value", "score", "normalizedScore") %in% colnames(fgseaRes)))
    expect_true(all(fgseaRes$p.value <= 1))
    expect_true(all(fgseaRes$p.value >= 0))
})

test_that('GSA', {
    gsaRes <- .runGSA(DERes, genesets[["genesets"]], method = "maxmean")
    expect_true(all(c("ID", "p.value", "score", "normalizedScore") %in% colnames(gsaRes)))
    expect_true(all(gsaRes$p.value <= 1))
    expect_true(all(gsaRes$p.value >= 0))
})

test_that('KS ', {
    ksRes <- .runKsWilcox(DERes, genesets[["genesets"]], sTest = "ks")
    expect_true(all(c("ID", "p.value", "score", "normalizedScore") %in% colnames(ksRes)))
    expect_true(all(ksRes$p.value <= 1))
    expect_true(all(ksRes$p.value >= 0))
})

test_that('Wilcox ', {
    wilcoxRes <- .runKsWilcox(DERes, genesets[["genesets"]], sTest = "wilcox")
    expect_true(all(c("ID", "p.value", "score", "normalizedScore") %in% colnames(wilcoxRes)))
    expect_true(all(wilcoxRes$p.value <= 1))
    expect_true(all(wilcoxRes$p.value >= 0))
})

test_that('GeneSet Enrichment Analysis with wilcox ', {
    result <- runGeneSetAnalysis(DERes, genesets, method = "wilcox")
    wilcoxRes <- .runKsWilcox(DERes, genesets[["genesets"]], sTest = "wilcox")
    expect_true(all(c("ID", "p.value", "score", "normalizedScore", "sampleSize") %in% colnames(result)))
    expect_true(all(result$p.value == wilcoxRes$p.value))
    expect_true(all(result$p.value <= 1))
    expect_true(all(result$p.value >= 0))
})

test_that('GeneSet Enrichment Analysis with fgsea ', {
    result <- runGeneSetAnalysis(DERes, genesets, method = "fgsea", FgseaArgs = list(eps = 1e-50, scoreType = "std", nPermSimple = 1000))
    expect_true(all(c("ID", "p.value", "score", "normalizedScore", "sampleSize") %in% colnames(result)))
    expect_true(all(result$p.value <= 1))
    expect_true(all(result$p.value >= 0))
})

test_that('GeneSet Enrichment Analysis with fgsea with wrong arguments ', {
    expect_error(runGeneSetAnalysis(DERes, genesets, method = "fgsea", FgseaArgs = list(scoretype = "std", nPermSimple = 1000)))
    expect_error(runGeneSetAnalysis(DERes, genesets = NULL, method = "fgsea", FgseaArgs = list(scoreType = "std", nPermSimple = 1000)))
})