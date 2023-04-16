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

# control vs condition
design <- model.matrix(~0 + group, data = colData)
contrast <- makeContrasts("groupcondition-groupcontrol", levels = design)

annotation <- .getIDMappingAnnotation("GPL570")
DERes <- runDEAnalysis(summarizedExperiment, method = "DESeq2", design, contrast, annotation = annotation)

genesets <- lapply(1:100, function(x) {
    sample(rownames(DERes), runif(1, 100, 500))
})
names(genesets) <- paste0("geneset", 1:100)

test_that('ORA', {
    oraRes <- .runORA(DERes, genesets, pThreshold = 0.05)
    expect_true(all(c("pathway", "p.value", "ES", "NES") %in% colnames(oraRes)))
    expect_true(all(oraRes$p.value <= 1))
    expect_true(all(oraRes$p.value >= 0))
})

test_that('fgsea', {
    fgseaRes <- .runFgsea(DERes, genesets, nperm = 1000)
    expect_true(all(c("pathway", "p.value", "ES", "NES") %in% colnames(fgseaRes)))
    expect_true(all(fgseaRes$p.value <= 1))
    expect_true(all(fgseaRes$p.value >= 0))
})

test_that('GSA', {
    gsaRes <- .runGSA(DERes, genesets, nperm = 1000)
    expect_true(all(c("pathway", "p.value", "ES", "NES") %in% colnames(gsaRes)))
    expect_true(all(gsaRes$p.value <= 1))
    expect_true(all(gsaRes$p.value >= 0))
})