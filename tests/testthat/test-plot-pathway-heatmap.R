library(testthat)
library(ggplot2)
library(RCPA)

set.seed(123)

result1 <- data.frame(
  ID = paste0("geneset", 1:10),
  name = paste0("Pathway ", 1:10),
  description = paste0("Description ", 1:10),
  p.value = runif(10) / 10,
  pFDR = runif(10) / 5,
  pathwaySize = runif(10, 100, 500),
  score = runif(10, -2, 2),
  normalizedScore = runif(10)
)
rownames(result1) <- result1$ID

result2 <- data.frame(
  ID = paste0("geneset", 1:10),
  name = paste0("Pathway ", 1:10),
  description = paste0("Description ", 1:10),
  p.value = runif(10) / 10,
  pFDR = runif(10) / 5,
  pathwaySize = runif(10, 100, 500),
  score = runif(10, -2, 2),
  normalizedScore = runif(10)
)
rownames(result2) <- result2$ID

resultsLst <- list(
  "study1" = result1,
  "study2" = result2
)

test_that("Plot pathway heatmap plot default", {
  pl <- plotPathwayHeatmap(resultsLst, yAxis = "ID")
  expect_true(is.ggplot(pl))
  expect_equal(pl$labels$size, "Normalized score")
  expect_equal(pl$labels$fill, "Direction")
  expect_equal(pl$labels$fill_new, "p.value.scaled")
})

test_that("Plot pathway heatmap plot with 'name' as yAxis", {
  pl <- plotPathwayHeatmap(resultsLst, yAxis = "name")
  expect_true(is.ggplot(pl))
  expect_equal(pl$labels$size, "Normalized score")
  expect_equal(pl$labels$fill, "Direction")
  expect_equal(pl$labels$fill_new, "p.value.scaled")
})



