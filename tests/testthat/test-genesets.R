library(testthat)

test_that("Get KEGG names", {
    keggNames <- .getKEGGPathwayNames("hsa")
    expect_true(length(keggNames) > 300)
    expect_true("hsa00010" %in% names(keggNames))
})

