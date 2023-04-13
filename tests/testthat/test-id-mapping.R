library(testthat)

# generate a random gene expression matrix
set.seed(123)
exprs <- round(matrix(2^abs(rnorm(1000, sd = 4)), nrow = 100, ncol = 10))
rownames(exprs) <- sample(AnnotationDbi::keys(hgu133plus2.db::hgu133plus2.db, keytype = "PROBEID"), nrow(exprs), replace = FALSE)

# simulate DE results

DE <- data.frame(
    PROBEID = rownames(exprs),
    p.value = runif(nrow(exprs), 0, 1),
    statistic = rnorm(nrow(exprs))
)

platform <- "GPL570"
annotationDB <- hgu133plus2.db::hgu133plus2.db
annotationDF