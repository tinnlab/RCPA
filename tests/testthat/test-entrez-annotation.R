library(testthat)

test_that("Get Entrez Annotation for a single ID", {
    entrezIds <- "127"
    entrezAnnotation <- getEntrezAnnotation(entrezIds)
    expect_equal(entrezAnnotation$ID, entrezIds)
    expect_equal(entrezAnnotation$Symbol, "ADH4")
    expect_equal(entrezAnnotation$Description, "alcohol dehydrogenase 4 (class II), pi polypeptide")
    expect_equal(entrezAnnotation$OtherDesignations, "all-trans-retinol dehydrogenase [NAD(+)] ADH4|alcohol dehydrogenase class II pi chain|aldehyde reductase|epididymis secretory protein Li 4")
    expect_equal(entrezAnnotation$OtherAliases, "ADH-2, HEL-S-4")
    expect_equal(entrezAnnotation$Chromosome, "4")
})

test_that("Entrez Annotation for multiple IDs", {
    entrezIds <- c("127", "128")
    entrezAnnotation <- getEntrezAnnotation(entrezIds)
    expect_equal(entrezAnnotation$ID, entrezIds)
    expect_equal(entrezAnnotation$Symbol, c("ADH4", "ADH5"))
    expect_equal(entrezAnnotation$Description, c("alcohol dehydrogenase 4 (class II), pi polypeptide", "alcohol dehydrogenase 5 (class III), chi polypeptide"))
    expect_equal(entrezAnnotation$Chromosome, c("4", "4"))
})

test_that("Entrez Annotation for multiple IDs with duplicates", {
    entrezIds <- c("127", "128", "127")
    entrezAnnotation <- getEntrezAnnotation(entrezIds)
    expect_equal(entrezAnnotation$ID, c("127", "128"))
    expect_equal(entrezAnnotation$Symbol, c("ADH4", "ADH5"))
    expect_equal(entrezAnnotation$Description, c("alcohol dehydrogenase 4 (class II), pi polypeptide", "alcohol dehydrogenase 5 (class III), chi polypeptide"))
    expect_equal(entrezAnnotation$Chromosome, c("4", "4"))
})

test_that("Entrez Annotation for multiple IDs with duplicates and NA", {
    entrezIds <- c("127", "128", "127", NA)
    entrezAnnotation <- getEntrezAnnotation(entrezIds)
    expect_equal(entrezAnnotation$ID, c("127", "128"))
    expect_equal(entrezAnnotation$Symbol, c("ADH4", "ADH5"))
    expect_equal(entrezAnnotation$Description, c("alcohol dehydrogenase 4 (class II), pi polypeptide", "alcohol dehydrogenase 5 (class III), chi polypeptide"))
    expect_equal(entrezAnnotation$Chromosome, c("4", "4"))
})

test_that("Entrez Annotation for multiple IDs with duplicates and NA and empty string", {
    entrezIds <- c("127", "128", "127", NA, "")
    entrezAnnotation <- getEntrezAnnotation(entrezIds)
    expect_equal(entrezAnnotation$ID, c("127", "128"))
    expect_equal(entrezAnnotation$Symbol, c("ADH4", "ADH5"))
    expect_equal(entrezAnnotation$Description, c("alcohol dehydrogenase 4 (class II), pi polypeptide", "alcohol dehydrogenase 5 (class III), chi polypeptide"))
    expect_equal(entrezAnnotation$Chromosome, c("4", "4"))
})

test_that("Entrez Annotation for multiple IDs with duplicates and NA and empty string and invalid IDs", {
    entrezIds <- c("127", "128", "127", NA, "", "invalid")
    entrezAnnotation <- getEntrezAnnotation(entrezIds)
    expect_equal(entrezAnnotation$ID, c("127", "128"))
    expect_equal(entrezAnnotation$Symbol, c("ADH4", "ADH5"))
    expect_equal(entrezAnnotation$Description, c("alcohol dehydrogenase 4 (class II), pi polypeptide", "alcohol dehydrogenase 5 (class III), chi polypeptide"))
    expect_equal(entrezAnnotation$Chromosome, c("4", "4"))
})

