library(testthat)
library(limma)
library(RCPA)

# generate a random gene expression matrix
set.seed(123)
exprs <- round(matrix(2^abs(rnorm(1000, sd = 4)), nrow = 100, ncol = 10))
rownames(exprs) <- paste0("gene", 1:100)
colnames(exprs) <- paste0("sample", 1:10)

controlSamples <- paste0("sample", 1:5)
conditionSamples <- paste0("sample", 6:10)

colData <- data.frame(
    row.names = colnames(exprs),
    group = factor(c(rep("control", length(controlSamples)), rep("condition", length(conditionSamples)))),
    pair = factor(c(seq_along(controlSamples), seq_along(conditionSamples)))
)

# control vs condition
design <- model.matrix(~0 + group, data = colData)
contrast <- makeContrasts("groupcondition-groupcontrol", levels = design)

design.paired <- model.matrix(~0 + group + pair, data = colData)
contrast.paired <- makeContrasts("groupcondition-groupcontrol", levels = design.paired)

# condition vs control
design.rev <- model.matrix(~0 + group, data = colData)
contrast.rev <- makeContrasts("groupcontrol-groupcondition", levels = design.rev)

design.paired.rev <- model.matrix(~0 + group + pair, data = colData)
contrast.paired.rev <- makeContrasts("groupcontrol-groupcondition", levels = design.paired.rev)

# Limma

test_that("Limma unpaired", {
    limmaRes <- RCPA:::.runLimma(exprs, design, contrast)

    expect_true(all(c("ID", "p.value", "statistic", "logFC", "logFCSE") %in% colnames(limmaRes)))
    expect_true(all(limmaRes$p.value <= 1))
    expect_true(all(limmaRes$p.value >= 0))
    expect_true(all(limmaRes$ID %in% rownames(exprs)))
})

test_that("Limma paired", {
    limmaRes <- RCPA:::.runLimma(exprs, design.paired, contrast.paired)

    expect_true(all(c("ID", "p.value", "statistic", "logFC", "avgExpr", "logFCSE") %in% colnames(limmaRes)))
    expect_true(all(limmaRes$p.value <= 1))
    expect_true(all(limmaRes$p.value >= 0))
    expect_true(all(limmaRes$ID %in% rownames(exprs)))
})

test_that("Limma unpaired reverse", {
    limmaRes <- RCPA:::.runLimma(exprs, design, contrast)
    limmaRes.rev <- RCPA:::.runLimma(exprs, design.rev, contrast.rev)

    expect_true(all(limmaRes$p.value == limmaRes.rev$p.value))
    expect_true(all(limmaRes$statistic == -limmaRes.rev$statistic))
})

test_that("Limma paired reverse", {
    limmaRes <- RCPA:::.runLimma(exprs, design.paired, contrast.paired)
    limmaRes.rev <- RCPA:::.runLimma(exprs, design.paired.rev, contrast.paired.rev)

    expect_true(all(limmaRes$p.value == limmaRes.rev$p.value))
    expect_true(all(limmaRes$statistic == -limmaRes.rev$statistic))
})

test_that("Limma unpaired vs paired", {
    limmaRes <- RCPA:::.runLimma(exprs, design, contrast)
    limmaRes.paired <- RCPA:::.runLimma(exprs, design.paired, contrast.paired)

    expect_true(!all(limmaRes$p.value == limmaRes.paired$p.value))
    expect_true(!all(limmaRes$statistic == limmaRes.paired$statistic))
})

# DESeq2

test_that("DESeq2 unpaired", {
    deseq2Res <- RCPA:::.runDESeq2(exprs, design, contrast)

    expect_true(all(c("ID", "p.value", "logFC", "logFCSE") %in% colnames(deseq2Res)))
    expect_true(all(deseq2Res$p.value <= 1))
    expect_true(all(deseq2Res$p.value >= 0))
    expect_true(all(deseq2Res$ID %in% rownames(exprs)))
})

test_that("DESeq2 paired", {
    deseq2Res <- RCPA:::.runDESeq2(exprs, design.paired, contrast.paired)

    expect_true(all(c("ID", "p.value", "logFC", "logFCSE") %in% colnames(deseq2Res)))
    expect_true(all(deseq2Res$p.value <= 1))
    expect_true(all(deseq2Res$p.value >= 0))
    expect_true(all(deseq2Res$ID %in% rownames(exprs)))
})

test_that("DESeq2 unpaired reverse", {
    deseq2Res <- RCPA:::.runDESeq2(exprs, design, contrast)
    deseq2Res.rev <- RCPA:::.runDESeq2(exprs, design.rev, contrast.rev)

    expect_true(all(deseq2Res$p.value == deseq2Res.rev$p.value))
    expect_true(all(deseq2Res$statistic == -deseq2Res.rev$statistic))
})

test_that("DESeq2 paired reverse", {
    deseq2Res <- RCPA:::.runDESeq2(exprs, design.paired, contrast.paired)
    deseq2Res.rev <- RCPA:::.runDESeq2(exprs, design.paired.rev, contrast.paired.rev)

    expect_true(all(deseq2Res$p.value == deseq2Res.rev$p.value))
    expect_true(all(deseq2Res$statistic == -deseq2Res.rev$statistic))
})

test_that("DESeq2 unpaired vs paired", {
    deseq2Res <- RCPA:::.runDESeq2(exprs, design, contrast)
    deseq2Res.paired <- RCPA:::.runDESeq2(exprs, design.paired, contrast.paired)

    expect_true(!all(deseq2Res$p.value == deseq2Res.paired$p.value))
    expect_true(!all(deseq2Res$statistic == deseq2Res.paired$statistic))
})

# edgeR

test_that("edgeR unpaired", {
    edgeRRes <- RCPA:::.runEdgeR(exprs, design, contrast)

    expect_true(all(c("ID", "p.value", "logFC", "logFC", "logFCSE") %in% colnames(edgeRRes)))
    expect_true(all(edgeRRes$p.value <= 1))
    expect_true(all(edgeRRes$p.value >= 0))
    expect_true(all(edgeRRes$ID %in% rownames(exprs)))
})

test_that("edgeR paired", {
    edgeRRes <- RCPA:::.runEdgeR(exprs, design.paired, contrast.paired)

    expect_true(all(c("ID", "p.value", "logFC", "logFC", "logFCSE") %in% colnames(edgeRRes)))
    expect_true(all(edgeRRes$p.value <= 1))
    expect_true(all(edgeRRes$p.value >= 0))
    expect_true(all(edgeRRes$ID %in% rownames(exprs)))
})

test_that("edgeR unpaired reverse", {
    edgeRRes <- RCPA:::.runEdgeR(exprs, design, contrast)
    edgeRRes.rev <- RCPA:::.runEdgeR(exprs, design.rev, contrast.rev)

    expect_true(all(edgeRRes$p.value == edgeRRes.rev$p.value))
    expect_true(all(edgeRRes$statistic == -edgeRRes.rev$statistic))
})

test_that("edgeR paired reverse", {
    edgeRRes <- RCPA:::.runEdgeR(exprs, design.paired, contrast.paired)
    edgeRRes.rev <- RCPA:::.runEdgeR(exprs, design.paired.rev, contrast.paired.rev)

    expect_true(all(edgeRRes$p.value == edgeRRes.rev$p.value))
    expect_true(all(edgeRRes$statistic == -edgeRRes.rev$statistic))
})

test_that("edgeR unpaired vs paired", {
    edgeRRes <- RCPA:::.runEdgeR(exprs, design, contrast)
    edgeRRes.paired <- RCPA:::.runEdgeR(exprs, design.paired, contrast.paired)

    expect_true(!all(edgeRRes$p.value == edgeRRes.paired$p.value))
    expect_true(!all(edgeRRes$statistic == edgeRRes.paired$statistic))
})



