library(testthat)
library(RCPA)

test_that("Get KEGG names for hsa", {
    skip_if_offline()
    keggNames <- RCPA:::.getKEGGPathwayNames("hsa")
    expect_true(length(keggNames) > 300)
    expect_true("path:hsa00010" %in% names(keggNames))
})

test_that("Get KEGG names for random org", {
    skip_if_offline()
    expect_error(RCPA:::.getKEGGPathwayNames("abcd"))
})

test_that("Get KEGG gene sets for hsa", {
    skip_if_offline()
    keggGeneSets <- RCPA:::.getKEGGGeneSets("hsa")

    expect_true(all(names(keggGeneSets) %in% c("database", "genesets", "names")))
    expect_true(keggGeneSets$database == "KEGG")
    expect_true(length(keggGeneSets$genesets) > 300)
    expect_true("path:hsa00010" %in% names(keggGeneSets$genesets))
    expect_true(length(keggGeneSets$genesets[["path:hsa00010"]]) > 10)
    expect_true(sum(is.na(as.numeric(unique(unlist(keggGeneSets$genesets))))) == 0)
})

test_that("Get KEGG gene sets for random org", {
    skip_if_offline()
    expect_warning(expect_error(RCPA:::.getKEGGGeneSets("abcd")))
})

test_that("Get GO names for namespace biological_process", {
    skip_if_offline()
    goNames <- RCPA:::.getGOTermNames("biological_process")
    expect_true(length(goNames) > 1000)
    expect_true("GO:0000003" %in% names(goNames))
})

test_that("Get GO names for namespace molecular_function", {
    skip_if_offline()
    goNames <- RCPA:::.getGOTermNames("molecular_function")
    expect_true(length(goNames) > 1000)
    expect_true("GO:0000016" %in% names(goNames))
})

test_that("Get GO names for namespace cellular_component", {
    skip_if_offline()
    goNames <- RCPA:::.getGOTermNames("cellular_component")
    expect_true(length(goNames) > 1000)
    expect_true("GO:0000015" %in% names(goNames))
})

test_that("Get GO names for random namespace", {
    skip_if_offline()
    expect_error(RCPA:::.getGOTermNames("abcd"))
})

test_that("Get GO terms for human, namespace biological_process", {
    skip_if_offline()
    goTerms <- RCPA:::.getGOTerms(9606, "biological_process")

    expect_true(all(names(goTerms) %in% c("database", "genesets", "names")))
    expect_true(goTerms$database == "GO")
    expect_true(length(goTerms$genesets) > 1000)
    expect_true("GO:0000723" %in% names(goTerms$genesets))
    expect_true(length(goTerms$genesets[["GO:0000723"]]) > 10)
    expect_true(sum(is.na(as.numeric(unique(unlist(goTerms$genesets))))) == 0)
})

test_that("Get GO terms for human, namespace molecular_function", {
    skip_if_offline()
    goTerms <- RCPA:::.getGOTerms(9606, "molecular_function")

    expect_true(all(names(goTerms) %in% c("database", "genesets", "names")))
    expect_true(goTerms$database == "GO")
    expect_true(length(goTerms$genesets) > 1000)
    expect_true("GO:0003674" %in% names(goTerms$genesets))
    expect_true(length(goTerms$genesets[["GO:0003674"]]) > 10)
    expect_true(sum(is.na(as.numeric(unique(unlist(goTerms$genesets))))) == 0)
})

test_that("Get GO terms for human, namespace cellular_component", {
    skip_if_offline()
    goTerms <- RCPA:::.getGOTerms(9606, "cellular_component")

    expect_true(all(names(goTerms) %in% c("database", "genesets", "names")))
    expect_true(goTerms$database == "GO")
    expect_true(length(goTerms$genesets) > 1000)
    expect_true("GO:0005575" %in% names(goTerms$genesets))
    expect_true(length(goTerms$genesets[["GO:0005575"]]) > 10)
    expect_true(sum(is.na(as.numeric(unique(unlist(goTerms$genesets))))) == 0)
})

test_that("Get GO terms for random species, namespace biological_process", {
    skip_if_offline()
    expect_error(RCPA:::.getGOTerms(1234567, "biological_process"))
})

test_that("Get GO terms for human, random namespace", {
    skip_if_offline()
    expect_error(RCPA:::.getGOTerms(9606, "abcd"))
})

test_that("Exported function: KEGG and hsa", {
    skip_if_offline()
    gensets <- getGeneSets("KEGG", "hsa", minSize = 10, maxSize = 1000)
    expect_true(all(names(gensets) %in% c("database", "genesets", "names")))
    expect_true(gensets$database == "KEGG")
    expect_true(all(sapply(gensets$genesets, length) >= 10))
    expect_true(all(sapply(gensets$genesets, length) <= 1000))
})

test_that("Exported function: KEGG and NULL org", {
    expect_error(getGeneSets("KEGG", NULL))
})

test_that("Exported function: KEGG and random org", {
    skip_if_offline()
    expect_warning(expect_error(getGeneSets("KEGG", "abcd")))
})

test_that("Exported function: GO and human, namespace biological_process", {
    skip_if_offline()
    gensets <- getGeneSets("GO", taxid = 9606, namespace = "biological_process", minSize = 10, maxSize = 1000)
    expect_true(all(names(gensets) %in% c("database", "genesets", "names")))
    expect_true(gensets$database == "GO")
    expect_true(all(sapply(gensets$genesets, length) >= 10))
    expect_true(all(sapply(gensets$genesets, length) <= 1000))
})

test_that("Exported function: GO and NULL taxId, namespace molecular_function", {
    expect_error(getGeneSets("GO", taxid = NULL, namespace = "molecular_function"))
})

test_that("Exported function: GO and human, namespace random", {
    expect_error(getGeneSets("GO", taxid = 9606, namespace = "abcd"))
})