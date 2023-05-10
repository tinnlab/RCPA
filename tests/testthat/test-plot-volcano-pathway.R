library(testthat)
library(ggplot2)
library(RCPA)

results <- data.frame(
    ID = paste0("geneset", 1:100),
    name = paste0("Pathway ", 1:100),
    description = paste0("Description ", 1:100),
    p.value = runif(100) / 10,
    pFDR = runif(100) / 5,
    pathwaySize = runif(100, 100, 500),
    nDE = runif(100, 10, 100),
    score = runif(100, -2, 2),
    normalizedScore = runif(100)
)

test_that("Plot volcano plot default", {
    pl <- plotVolcanoPathway(results)
    expect_true(is.ggplot(pl))
    expect_equal(pl$labels$y, "-log10 pFDR")
    expect_equal(pl$labels$x, "Normalized score")
})

test_that("Plot volcano plot with score", {
    pl <- plotVolcanoPathway(results, xAxis = "score")
    expect_true(is.ggplot(pl))
    expect_equal(pl$labels$x, "Score")
})

test_that("Plot volcano plot with normalizedScore", {
    pl <- plotVolcanoPathway(results, xAxis = "normalizedScore")
    expect_true(is.ggplot(pl))
    expect_equal(pl$labels$x, "Normalized score")
})

test_that("Plot volcano plot with p-value", {
    pl <- plotVolcanoPathway(results, yAxis = "-log10(p.value)")
    expect_true(is.ggplot(pl))
    expect_equal(pl$labels$y, "-log10 p-value")
})

test_that("Plot volcano plot with pFDR", {
    pl <- plotVolcanoPathway(results, yAxis = "-log10(pFDR)")
    expect_true(is.ggplot(pl))
    expect_equal(pl$labels$y, "-log10 pFDR")
})