library(testthat)
library(ggplot2)
library(RCPA)

results <- lapply(1:3, function(i) {
    set.seed(i)

    data.frame(
        ID = paste0("geneset", 1:100),
        name = paste0("Pathway ", 1:100),
        description = paste0("Description ", 1:100),
        p.value = runif(100) / 10,
        pFDR = runif(100) / 5,
        size = runif(100, 100, 500),
        nDE = runif(100, 10, 100),
        score = runif(100, -2, 2),
        normalizedScore = runif(100, -2, 2)
    )
})

# devtools::document()
# devtools::load_all()

test_that("Plot bar plot default", {
    pl <- plotBarChart(results)
    expect_true(is.ggplot(pl))
    expect_equal(length(pl$labels), 9)
})

test_that("Plot bar plot with non FDR", {
    pl <- plotBarChart(results, useFDR = FALSE)
    expect_true(is.ggplot(pl))
    expect_equal(length(pl$labels), 9)
})

test_that("Plot bar plot with FDR", {
    pl <- plotBarChart(results, useFDR = TRUE)
    expect_true(is.ggplot(pl))
    expect_equal(length(pl$labels), 9)
})

test_that("Plot bar plot with p-value", {
    pl <- plotBarChart(results, by = "p.value")
    expect_true(is.ggplot(pl))
    expect_equal(pl$labels$y, "-log10 p.value")
})

test_that("Plot bar plot by p-value without FDR", {
    pl <- plotBarChart(results, by = "p.value", useFDR = FALSE)
    expect_true(is.ggplot(pl))
    expect_equal(pl$labels$y, "-log10 p.value")
})

test_that("Plot bar plot with pFDR", {
    pl <- plotBarChart(results, by = "pFDR", useFDR = TRUE)
    expect_true(is.ggplot(pl))
    expect_equal(pl$labels$y, "-log10 pFDR")
})

test_that("Plot bar plot by pFDR with useFDR false", {
    pl <- plotBarChart(results, by = "pFDR", useFDR = FALSE)
    expect_true(is.ggplot(pl))
    expect_equal(pl$labels$y, "-log10 pFDR")
})

test_that("Plot bar plot by size", {
    expect_error(plotBarChart(results, by = "size"))
})

test_that("Plot bar plot by nDE", {
    expect_error(plotBarChart(results, by = "nDE"))
})

test_that("Plot bar plot by score", {
    pl <- plotBarChart(results, by = "score")
    expect_true(is.ggplot(pl))
    expect_equal(pl$labels$y, "Score")
})

test_that("Plot bar plot by normalizedScore", {
    pl <- plotBarChart(results, by = "normalizedScore")
    expect_true(is.ggplot(pl))
    expect_equal(pl$labels$y, "Normalized score")
})

