library(testthat)
library(RCyjs)
library(RCPA)
library(RCPA)

styleFile <- system.file(package="RCPA", "extdata", "pieStyle.js")

pThreshold <- 0.05
useFDR <- FALSE
edgeThreshold <- 0.1

getTestPathwayNetworkGeneset <- function (){
    gs <- getGeneSets("KEGG")
    gs
}

getTestPathwayNetworkData <- function(gs){
    results <- lapply(1:3, function(i) {
        set.seed(i)

        data.frame(
            ID = names(gs$genesets),
            name = gs$names,
            p.value = runif(length(gs$genesets))/10,
            pFDR = runif(length(gs$genesets))/10,
            score = runif(length(gs$genesets), -2, 2)
        )
    })

    names(results) <- c("ORA", "FGSEA", "GSA")
    
    results
}

test_that('plot pathway network with default params', {
    skip_on_cran()
    skip_if_offline()
    gs <- getTestPathwayNetworkGeneset()
    results <- getTestPathwayNetworkData(gs)
    IDs <- sample(names(gs$genesets), 10)
    genesets <- gs$genesets[IDs]
    rcy <- plotPathwayNetwork(results, genesets, pThreshold = pThreshold, useFDR = useFDR)
    expect_true("RCyjs" %in% class(rcy))
    expect_true(rcy@port > 10000)
    expect_true(rcy@port < 20000)
    expect_true(all(rcy@graph@nodes %in% names(genesets)))
    expect_true(all(names(genesets) %in% rcy@graph@nodes))
})

test_that('plot pathway network without FDR', {
    skip_on_cran()
    skip_if_offline()
    gs <- getTestPathwayNetworkGeneset()
    results <- getTestPathwayNetworkData(gs)
    IDs <- sample(names(gs$genesets), 10)
    genesets <- gs$genesets[IDs]
    rcy <- plotPathwayNetwork(results, genesets, pThreshold = pThreshold, useFDR = FALSE)
    expect_true("RCyjs" %in% class(rcy))
    expect_true(rcy@port > 10000)
    expect_true(rcy@port < 20000)
    expect_true(all(rcy@graph@nodes %in% names(genesets)))
    expect_true(all(names(genesets) %in% rcy@graph@nodes))
})

test_that('plot pathway network with edgeThreshold is 0', {
    skip_on_cran()
    skip_if_offline()
    gs <- getTestPathwayNetworkGeneset()
    results <- getTestPathwayNetworkData(gs)
    IDs <- sample(names(gs$genesets), 10)
    genesets <- gs$genesets[IDs]
    rcy <- plotPathwayNetwork(results, genesets, pThreshold = pThreshold, useFDR = TRUE, edgeThreshold = 0)
    expect_true("RCyjs" %in% class(rcy))
    expect_true(rcy@port > 10000)
    expect_true(rcy@port < 20000)
    expect_true(all(rcy@graph@nodes %in% names(genesets)))
    expect_true(all(names(genesets) %in% rcy@graph@nodes))
})

test_that('plot pathway network with custom labels', {
    skip_on_cran()
    skip_if_offline()
    gs <- getTestPathwayNetworkGeneset()
    results <- getTestPathwayNetworkData(gs)
    IDs <- sample(names(gs$genesets), 10)
    genesets <- gs$genesets[IDs]
    rcy <- plotPathwayNetwork(results, genesets, pThreshold = pThreshold, labels = gs$names[IDs])
    expect_true("RCyjs" %in% class(rcy))
    expect_true(rcy@port > 10000)
    expect_true(rcy@port < 20000)
    expect_true(all(rcy@graph@nodes %in% names(genesets)))
    expect_true(all(names(genesets) %in% rcy@graph@nodes))
})

test_that('plot pathway network discrete mode', {
    skip_on_cran()
    skip_if_offline()
    gs <- getTestPathwayNetworkGeneset()
    results <- getTestPathwayNetworkData(gs)
    IDs <- sample(names(gs$genesets), 10)
    genesets <- gs$genesets[IDs]
    rcy <- plotPathwayNetwork(results, genesets, mode = "discrete")
    expect_true("RCyjs" %in% class(rcy))
    expect_true(rcy@port > 10000)
    expect_true(rcy@port < 20000)
    expect_true(all(rcy@graph@nodes %in% names(genesets)))
    expect_true(all(names(genesets) %in% rcy@graph@nodes))
})