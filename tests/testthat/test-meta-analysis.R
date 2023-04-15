library(testthat)
library(dplyr)

devtools::load_all()
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
  pathway = names(genesets),
  pathway.name = names(genesets),
  pathway.size = genesets.size,
  p.value = p.values[1:10],
  score = scores[1:10],
  normalizedScore = normalizedScores[1:10],
  nDE = nDEs[1:10],
  sample.size = sample.sizes[1:10],
  stringsAsFactors = FALSE
)

DF2 <- data.frame(
  pathway = names(genesets),
  pathway.name = names(genesets),
  pathway.size = genesets.size,
  p.value = p.values[11:20],
  score = scores[11:20],
  normalizedScore = normalizedScores[11:20],
  nDE = nDEs[11:20],
  sample.size = sample.sizes[11:20],
  stringsAsFactors = FALSE
)

DF3 <- data.frame(
  pathway = names(genesets),
  pathway.name = names(genesets),
  pathway.size = genesets.size,
  p.value = p.values[21:30],
  score = scores[21:30],
  normalizedScore = normalizedScores[21:30],
  nDE = nDEs[21:30],
  sample.size = sample.sizes[21:30],
  stringsAsFactors = FALSE
)

allDFs <- list(DF1, DF2, DF3)
allData <- allDFs %>% bind_rows()

test_that('Combine Pvalues with Fisher ', {
  metaRes <- .combinePvalues(allData, method = "fisher")
  selected_pval <- metaRes$p.value[metaRes$pathway == "geneset1"]
  DFs_pvals <- c(DF1$p.value[DF1$pathway == "geneset1"], DF2$p.value[DF2$pathway == "geneset1"], DF3$p.value[DF3$pathway == "geneset1"])
  expect_true(all(c("pathway", "p.value", "score", "score.sd", "count") %in% colnames(metaRes)))
  expect_true(! selected_pval %in% c())
  expect_true(all(metaRes$p.value <= 1))
  expect_true(all(metaRes$p.value >= 0))
})

test_that('Combine Pvalues with Stouffer ', {
  metaRes <- .combinePvalues(allData, method = "stouffer")
  selected_pval <- metaRes$p.value[metaRes$pathway == "geneset1"]
  DFs_pvals <- c(DF1$p.value[DF1$pathway == "geneset1"], DF2$p.value[DF2$pathway == "geneset1"], DF3$p.value[DF3$pathway == "geneset1"])
  expect_true(all(c("pathway", "p.value", "score", "score.sd", "count") %in% colnames(metaRes)))
  expect_true(! selected_pval %in% DFs_pvals)
  expect_true(all(metaRes$p.value <= 1))
  expect_true(all(metaRes$p.value >= 0))
})

test_that('Combine Pvalues with addCLT ', {
  metaRes <- .combinePvalues(allData, method = "addCLT")
  selected_pval <- metaRes$p.value[metaRes$pathway == "geneset1"]
  DFs_pvals <- c(DF1$p.value[DF1$pathway == "geneset1"], DF2$p.value[DF2$pathway == "geneset1"], DF3$p.value[DF3$pathway == "geneset1"])
  expect_true(all(c("pathway", "p.value", "score", "score.sd", "count") %in% colnames(metaRes)))
  expect_true(! selected_pval %in% DFs_pvals)
  expect_true(all(metaRes$p.value <= 1))
  expect_true(all(metaRes$p.value >= 0))
})

test_that('Combine Pvalues with geoMean ', {
  metaRes <- .combinePvalues(allData, method = "geoMean")
  selected_pval <- metaRes$p.value[metaRes$pathway == "geneset1"]
  DFs_pvals <- c(DF1$p.value[DF1$pathway == "geneset1"], DF2$p.value[DF2$pathway == "geneset1"], DF3$p.value[DF3$pathway == "geneset1"])
  expect_true(all(c("pathway", "p.value", "score", "score.sd", "count") %in% colnames(metaRes)))
  expect_true(! selected_pval %in% DFs_pvals)
  expect_true(all(metaRes$p.value <= 1))
  expect_true(all(metaRes$p.value >= 0))
})

test_that('Combine Pvalues with minP ', {
  metaRes <- .combinePvalues(allData, method = "minP")
  selected_pval <- metaRes$p.value[metaRes$pathway == "geneset1"]
  DFs_pvals <- c(DF1$p.value[DF1$pathway == "geneset1"], DF2$p.value[DF2$pathway == "geneset1"], DF3$p.value[DF3$pathway == "geneset1"])
  expect_true(all(c("pathway", "p.value", "score", "score.sd", "count") %in% colnames(metaRes)))
  expect_true(selected_pval == min(DFs_pvals))
  expect_true(all(metaRes$p.value <= 1))
  expect_true(all(metaRes$p.value >= 0))
})

test_that('Combine Enrichment Scores ', {
  metaRes <- .combineEnrichmentScores(allData)
  selected_pval <- metaRes$p.value[metaRes$pathway == "geneset1"]
  DFs_pvals <- c(DF1$p.value[DF1$pathway == "geneset1"], DF2$p.value[DF2$pathway == "geneset1"], DF3$p.value[DF3$pathway == "geneset1"])
  expect_true(all(c("pathway", "p.value", "score", "score.sd", "count") %in% colnames(metaRes)))
  expect_true(! selected_pval %in% DFs_pvals)
  expect_true(all(metaRes$p.value <= 1))
  expect_true(all(metaRes$p.value >= 0))
  expect_true(all(metaRes$score >= -2))
  expect_true(all(metaRes$score <= 2))
})

test_that('Combine Pathway Analysis Results ', {
  metaRes <- combinePathwayAnalysisResults(allDFs, method = "ES")
  selected_pval <- metaRes$p.value[metaRes$pathway == "geneset1"]
  DFs_pvals <- c(DF1$p.value[DF1$pathway == "geneset1"], DF2$p.value[DF2$pathway == "geneset1"], DF3$p.value[DF3$pathway == "geneset1"])
  expect_true(all(c("pathway", "p.value", "score", "score.sd", "count") %in% colnames(metaRes)))
  expect_true(! selected_pval %in% DFs_pvals)
  expect_true(all(metaRes$p.value <= 1))
  expect_true(all(metaRes$p.value >= 0))
  expect_true(all(metaRes$score >= -2))
  expect_true(all(metaRes$score <= 2))
})

test_that('Combine Pathway Analysis Results ', {
  metaRes <- combinePathwayAnalysisResults(allDFs, method = "ES")
  selected_pval <- metaRes$p.value[metaRes$pathway == "geneset1"]
  DFs_pvals <- c(DF1$p.value[DF1$pathway == "geneset1"], DF2$p.value[DF2$pathway == "geneset1"], DF3$p.value[DF3$pathway == "geneset1"])
  expect_true(all(c("pathway", "p.value", "score", "score.sd", "count") %in% colnames(metaRes)))
  expect_true(! selected_pval %in% DFs_pvals)
  expect_true(all(metaRes$p.value <= 1))
  expect_true(all(metaRes$p.value >= 0))
  expect_true(all(metaRes$score >= -2))
  expect_true(all(metaRes$score <= 2))
})

test_that('Combine Pathway Analysis Results ', {
  metaRes <- combinePathwayAnalysisResults(allDFs, method = "stouffer")
  selected_pval <- metaRes$p.value[metaRes$pathway == "geneset1"]
  DFs_pvals <- c(DF1$p.value[DF1$pathway == "geneset1"], DF2$p.value[DF2$pathway == "geneset1"], DF3$p.value[DF3$pathway == "geneset1"])
  expect_true(all(c("pathway", "p.value", "score", "score.sd", "count") %in% colnames(metaRes)))
  expect_true(! selected_pval %in% DFs_pvals)
  expect_true(all(metaRes$p.value <= 1))
  expect_true(all(metaRes$p.value >= 0))
  expect_true(all(metaRes$score >= -2))
  expect_true(all(metaRes$score <= 2))
})

test_that('Combine Pathway Analysis Results ', {
  singleList <- list(DF1)
  notCompatibleDim <- list(DF1, DF2[1:5,])
  DF2.tmp <- DF2
  rownames(DF2.tmp)[1:10] <- LETTERS[1:10]
  notCompatibleRownames <- list(DF1, DF2.tmp)
  DF4 <- NULL
  nullList <- list(DF1, DF2, DF4, DF3)
  expect_error(combinePathwayAnalysisResults(singleList, method = "ES"), "Meta analysis is valid for two or more studies.")
  expect_error(combinePathwayAnalysisResults(notCompatibleDim, method = "ES"), "All the dataframes in the input list must have the same dimension.")
  expect_error(combinePathwayAnalysisResults(notCompatibleRownames, method = "ES"), "All the dataframes in the input list must have the same row names.")
  expect_error(combinePathwayAnalysisResults(nullList, method = "ES"), "There is null object in the input list.")
})



