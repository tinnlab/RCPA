library(testthat)
library(RCPA)

set.seed(123)

result1 <- data.frame(
    ID = paste0("geneset", 1:100),
    name = paste0("Pathway ", 1:100),
    description = paste0("Description ", 1:100),
    p.value = runif(100) / 10,
    pFDR = runif(100) / 5,
    pathwaySize = runif(100, 100, 500),
    score = runif(100, -2, 2),
    normalizedScore = runif(100)
)

result2 <- data.frame(
    ID = paste0("geneset", 1:100),
    name = paste0("Pathway ", 1:100),
    description = paste0("Description ", 1:100),
    p.value = runif(100) / 10,
    pFDR = runif(100) / 5,
    pathwaySize = runif(100, 100, 500),
    score = runif(100, -2, 2),
    normalizedScore = runif(100)
)

result3 <- data.frame(
    ID = paste0("geneset", 1:100),
    name = paste0("Pathway ", 1:100),
    description = paste0("Description ", 1:100),
    p.value = runif(100) / 10,
    pFDR = runif(100) / 5,
    pathwaySize = runif(100, 100, 500),
    score = runif(100, -2, 2),
    normalizedScore = runif(100)
)

resultsLst <- list(
    "study1" = result1,
    "study2" = result2,
    "study3" = result3
)

test_that(".runWeightedMean ", {
    result <- .runWeightedMean(resultsLst, c(1, 1, 1), useFDR = TRUE)
    expect_true(any(result$p.value >=0) & any(result$p.value <= 1))
    expect_true(nrow(result) == 100)
    expect_true(all(c("ID", "p.value", "name", "pathwaySize") %in% colnames(result)))
    expect_error(.runWeightedMean(resultsLst, c(1, 2)))
})

test_that(".runRankPathways ", {
    result <- .runRankPathways(resultsLst, rankParam = "pFDR")
    expect_true(length(result) == 3)
})

test_that("runConsensusAnalysis with weighted.mean as method ", {
    result <- runConsensusAnalysis(resultsLst, method = "weightedZMean", weightsList = c(2,1,1))
    expect_true(nrow(result) == 100 & ncol(result) == 5)
    expect_true(all(result$p.value >= 0 & result$p.value <= 1))
    expect_error(runConsensusAnalysis(resultsLst, method = "weightedZMean", weightsList = c(2,1)))

    resultsLst2 <- list(
        "study1" = result1[1:10,],
        "study2" = result2[1:20,],
        "study3" = result3
    )
    resultt <- runConsensusAnalysis(resultsLst2, method = "weightedZMean", weightsList = c(2,1,1))
    expect_true(nrow(resultt) == 10 & ncol(resultt) == 5)

    resultsLst3 <- list(
        "study1" = result1[1:10,],
        "study2" = result2[11:20,],
        "study3" = result3[21:30,]
    )
    expect_error(runConsensusAnalysis(resultsLst3, method = "weightedZMean", weightsList = c(2,1,1)), "There is no common pathways among input data!")
})

test_that("runConsensusAnalysis with RRA as method ", {
    result <- runConsensusAnalysis(resultsLst, method = "RRA", rank.by = "pFDR")
    expect_true(nrow(result) == 100 & ncol(result) == 5)
    expect_true(all(result$p.value >= 0 & result$p.value <= 1))

    space1 <- resultsLst[[1]]$ID %>% as.vector()
    space2 <- resultsLst[[2]]$ID %>% as.vector()
    space3 <- resultsLst[[3]]$ID %>% as.vector()
    resultt <- runConsensusAnalysis(resultsLst, method = "RRA", rank.by = "pFDR", backgroundSpace = list(space1, space2, space3))
    expect_true(nrow(resultt) == 100 & ncol(resultt) == 5)
    expect_true(all(resultt$p.value >= 0 & resultt$p.value <= 1))
    expect_error(runConsensusAnalysis(resultsLst, method = "RRA", rank.by = "pFDR", backgroundSpace = list(space2, space3)))
})
