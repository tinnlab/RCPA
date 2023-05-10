library(testthat)
library(ggplot2)
library(RCPA)

set.seed(123)

result1 <- data.frame(
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
rownames(result1) <- result1$ID

result2 <- data.frame(
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
rownames(result2) <- result2$ID

resultsLst <- list(
  "study1" = result1,
  "study2" = result2
)

test_that("Plot forest plot default", {
  pl <- plotForest(resultsLst, yAxis = "ID")
  expect_true(is.ggplot(pl))
  expect_equal(length(pl$labels), 8)
})

test_that("Plot forest plot name as yAxis", {
  pl <- plotForest(resultsLst, yAxis = "name")
  expect_true(is.ggplot(pl))
  expect_equal(length(pl$labels), 8)
})

