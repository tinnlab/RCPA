library(testthat)
library(hgu133plus2.db)
library(AnnotationDbi)
library(RCPA)

# generate a random gene expression matrix
set.seed(123)
exprs <- round(matrix(2^abs(rnorm(1000, sd = 4)), nrow = 100, ncol = 10))
rownames(exprs) <- sample(keys(hgu133plus2.db, keytype = "PROBEID"), nrow(exprs), replace = FALSE)

# simulate DE results

DEResults <- data.frame(
    ID = rownames(exprs),
    p.value = runif(nrow(exprs), 0, 1),
    logFC = rnorm(nrow(exprs)),
    dispersion = runif(nrow(exprs), 0, 1),
    statistic = rnorm(nrow(exprs))
)

platform <- "GPL570"
annotationDB <- hgu133plus2.db

test_that("Get annotation from GPL570", {
    skip_if_offline()
    mapping <- RCPA:::.getIDMappingAnnotation("GPL570")

    expect_true(all(colnames(mapping) == c("FROM", "TO")))
    expect_true(all(mapping$FROM %in% keys(hgu133plus2.db, keytype = "PROBEID")))
    expect_true(all(mapping$TO %in% keys(hgu133plus2.db, keytype = "ENTREZID")))
})

test_that("Get annotation from NULL", {
    skip_if_offline()
    expect_error(RCPA:::.getIDMappingAnnotation(NULL))
})

test_that("Get annotation from non-existing platform", {
    skip_if_offline()
    expect_error(RCPA:::.getIDMappingAnnotation("GPL123"))
})

test_that("Map expresison from GPL570 to ENTREZ", {
    skip_if_offline()
    annotation <- RCPA:::.getIDMappingAnnotation("GPL570")
    mappedExprs <- .mapIDs(exprs, annotation, DEResults)

    expect_true(all(colnames(mappedExprs$exprs) == colnames(exprs)))
    expect_true(all(rownames(mappedExprs$exprs) %in% keys(hgu133plus2.db, keytype = "ENTREZID")))
    expect_true(all(mappedExprs$DEResults$ID %in% keys(hgu133plus2.db, keytype = "ENTREZID")))
    expect_true(nrow(mappedExprs$DEResults) <= nrow(DEResults))
    expect_true(all(rownames(mappedExprs$exprs) == rownames(mappedExprs$mapping$TO)))
    expect_true(all(rownames(mappedExprs$exprs) == rownames(mappedExprs$DEResults)))
})

test_that("Map expresison from GPL570 to ENTREZ with NULL annotation", {
    mappedExprs <- RCPA:::.mapIDs(exprs, NULL, DEResults)

    expect_true(all(colnames(mappedExprs$exprs) == colnames(exprs)))
    expect_true(all(rownames(mappedExprs$exprs) %in% keys(hgu133plus2.db, keytype = "PROBEID")))
    expect_true(all(mappedExprs$DEResults$ID %in% keys(hgu133plus2.db, keytype = "PROBEID")))
    expect_true(nrow(mappedExprs$DEResults) == nrow(DEResults))
})

test_that("Map expresison from GPL570 to ENTREZ with NULL annotation and NULL DE results", {
    mappedExprs <- RCPA:::.mapIDs(exprs, NULL, NULL)

    expect_true(all(colnames(mappedExprs$exprs) == colnames(exprs)))
    expect_true(all(rownames(mappedExprs$exprs) %in% keys(hgu133plus2.db, keytype = "PROBEID")))
    expect_true(is.null(mappedExprs$DEResults))
})

test_that("Map expresison from GPL570 to ENTREZ with NULL DE results", {
    skip_if_offline()
    annotation <- RCPA:::.getIDMappingAnnotation("GPL570")
    expect_error(RCPA:::.mapIDs(exprs, annotation, NULL))
})

test_that("Map expresison from GPL570 to ENTREZ with NULL expression", {
    annotation <- RCPA:::.getIDMappingAnnotation("GPL570")
    expect_error(RCPA:::.mapIDs(NULL, annotation, DEResults))
})

