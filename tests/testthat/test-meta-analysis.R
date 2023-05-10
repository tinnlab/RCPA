library(testthat)
library(dplyr)
library(RCPA)

#generate three simulated dataframes of pathway analysis results
genes <- paste0("gene", seq(1:500))
genesets <- lapply(1:10, function(x) {
  sample(genes, runif(1, 50, 300))
})
names(genesets) <- paste0("geneset", 1:10)
genesets.size <- genesets %>% lapply(function (gs) length(gs)) %>% unlist() %>% as.vector()

p.values <- runif(30)
scores <- runif(30, -2, 2)
normalizedScores <- runif(30, -2, 2)
nDEs <- ceiling(runif(30, 0, 100))
sample.sizes <- c(rep(10, 10), rep(12, 10), rep(50, 10))

DF1 <- data.frame(
  ID = names(genesets),
  p.value = p.values[1:10],
  score = scores[1:10],
  normalizedScore = normalizedScores[1:10],
  sampleSize = sample.sizes[1:10],
  name = names(genesets),
  pFDR = p.adjust(p.values[1:10], method = "fdr"),
  pathwaySize = genesets.size,
  stringsAsFactors = FALSE
)

DF2 <- data.frame(
  ID = names(genesets),
  p.value = p.values[11:20],
  score = scores[11:20],
  normalizedScore = normalizedScores[11:20],
  sampleSize = sample.sizes[11:20],
  name = names(genesets),
  pFDR = p.adjust(p.values[11:20], method = "fdr"),
  pathwaySize = genesets.size,
  stringsAsFactors = FALSE
)

DF3 <- data.frame(
  ID = names(genesets),
  p.value = p.values[21:30],
  score = scores[21:30],
  normalizedScore = normalizedScores[21:30],
  sampleSize = sample.sizes[21:30],
  name = names(genesets),
  pFDR = p.adjust(p.values[21:30], method = "fdr"),
  pathwaySize = genesets.size,
  stringsAsFactors = FALSE
)

allDFs <- list(DF1, DF2, DF3)
allData <- allDFs %>% bind_rows()

test_that('Combine Pvalues with Fisher ', {
  metaRes <- RCPA:::.runFisher(allData$p.value)
  expect_true(metaRes <= 1)
  expect_true(metaRes >= 0)
})

test_that('Combine Pvalues with Stouffer ', {
  metaRes <- RCPA:::.runStouffer(allData$p.value)
  expect_true(metaRes <= 1)
  expect_true(metaRes >= 0)
})

test_that('Combine Pvalues with addCLT ', {
  metaRes <- RCPA:::.runAddCLT(allData$p.value)
  expect_true(metaRes <= 1)
  expect_true(metaRes >= 0)
})

test_that('Combine Pvalues with geoMean ', {
  metaRes <- RCPA:::.runGeoMean(allData$p.value)
  expect_true(metaRes <= 1)
  expect_true(metaRes >= 0)
})

test_that('Combine Pathway Analysis Results using REML ', {
  metaRes <- combineEnrichmentAnalysisResults(allDFs, method = "REML")
  selected_pval <- metaRes$p.value[metaRes$ID == "geneset1"]
  DFs_pvals <- c(DF1$p.value[DF1$ID == "geneset1"], DF2$p.value[DF2$ID == "geneset1"], DF3$p.value[DF3$ID == "geneset1"])
  expect_true(all(c("ID", "p.value", "score", "pFDR", "normalizedScore") %in% colnames(metaRes)))
  expect_true(! selected_pval %in% DFs_pvals)
  expect_true(all(metaRes$p.value <= 1))
  expect_true(all(metaRes$p.value >= 0))
  expect_true(all(metaRes$score >= -2))
  expect_true(all(metaRes$score <= 2))
})

test_that('Combine Pathway Analysis Results using Stouffer ', {
  metaRes <- combineEnrichmentAnalysisResults(allDFs, method = "stouffer")
  selected_pval <- metaRes$p.value[metaRes$ID == "geneset1"]
  DFs_pvals <- c(DF1$p.value[DF1$ID == "geneset1"], DF2$p.value[DF2$ID == "geneset1"], DF3$p.value[DF3$ID == "geneset1"])
  expect_true(all(c("ID", "p.value", "score", "pFDR", "normalizedScore") %in% colnames(metaRes)))
  expect_true(! selected_pval %in% DFs_pvals)
  expect_true(all(metaRes$p.value <= 1))
  expect_true(all(metaRes$p.value >= 0))
  expect_true(all(metaRes$score >= -2))
  expect_true(all(metaRes$score <= 2))
})

test_that('Combine Pathway Analysis Results check errors ', {
  singleList <- list(DF1)
  notCompatibleCols <- list(DF1, DF2[,1:3])
  DF4 <- NULL
  nullList <- list(DF1, DF2, DF4, DF3)
  expect_error(combineEnrichmentAnalysisResults(singleList, method = "REML"), "Meta analysis is valid for two or more studies.")
  expect_error(combineEnrichmentAnalysisResults(notCompatibleCols, method = "REML"), "All the dataframes in the input list must have ID, p.value, normalizedScore and sampleSize columns.")
  expect_error(combineEnrichmentAnalysisResults(nullList, method = "REML"), "There is null object in the input list.")
})



