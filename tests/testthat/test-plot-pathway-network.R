library(testthat)
library(RCyjs)
library(RCPA)

styleFile <- system.file(package="RCPA", "extdata", "pieStyle.js")

gs <- getGeneSets("KEGG")
genesets <- gs$genesets

results <- lapply(1:3, function(i) {
    set.seed(i)

    data.frame(
        ID = names(genesets),
        name = gs$names,
        p.value = runif(length(genesets))/10,
        pFDR = runif(length(genesets))/10,
        ES = runif(length(genesets)),
        NES = runif(length(genesets)),
        nDE = sample(1:1000, length(genesets))
    )
})

names(results) <- c("ORA", "FGSEA", "GSA")

IDs <- sample(names(gs$genesets), 10)

genesets <- gs$genesets[IDs]
pThreshold <- 0.05
useFDR <- FALSE
edgeThreshold <- 0.1

test_that('plot pathway network with default params', {
    rcy <- plotPathwayNetwork(results, genesets, pThreshold = pThreshold, useFDR = useFDR)
    expect_true("RCyjs" %in% class(rcy))
    expect_true(rcy@port > 10000)
    expect_true(rcy@port < 20000)
    expect_true(all(rcy@graph@nodes %in% names(genesets)))
    expect_true(all(names(genesets) %in% rcy@graph@nodes))
})

test_that('plot pathway network without FDR', {
    rcy <- plotPathwayNetwork(results, genesets, pThreshold = pThreshold, useFDR = FALSE)
    expect_true("RCyjs" %in% class(rcy))
    expect_true(rcy@port > 10000)
    expect_true(rcy@port < 20000)
    expect_true(all(rcy@graph@nodes %in% names(genesets)))
    expect_true(all(names(genesets) %in% rcy@graph@nodes))
})

test_that('plot pathway network with edgeThreshold is 0', {
    rcy <- plotPathwayNetwork(results, genesets, pThreshold = pThreshold, useFDR = TRUE, edgeThreshold = 0)
    expect_true("RCyjs" %in% class(rcy))
    expect_true(rcy@port > 10000)
    expect_true(rcy@port < 20000)
    expect_true(all(rcy@graph@nodes %in% names(genesets)))
    expect_true(all(names(genesets) %in% rcy@graph@nodes))
})

test_that('plot pathway network with custom labels', {
    rcy <- plotPathwayNetwork(results, genesets, pThreshold = pThreshold, labels = gs$names[IDs])
    expect_true("RCyjs" %in% class(rcy))
    expect_true(rcy@port > 10000)
    expect_true(rcy@port < 20000)
    expect_true(all(rcy@graph@nodes %in% names(genesets)))
    expect_true(all(names(genesets) %in% rcy@graph@nodes))
})
