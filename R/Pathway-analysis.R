#' @title Pathway Enrichment Analysis using ORA
#' @description This function performs pathway analysis based on ORA (Over Representation Analysis).
#' This function is used internally by runPathwayAnalysis.
#' @param summarizedExprimentObj The generated SummarizedExpriment object from DE analysis result.
#' @param genesets The genesets definition, ex. KEGG genesets.
#' @return A dataframe of pathway analysis results
#' @details This function is used internally by runPathwayAnalysis.
#' @importFrom SummarizedExperiment SummarizedExperiment rowData assay
#' @importFrom dplyr %>% filter
#' @importFrom tidyr drop_na
.runORA <- function(summarizedExprimentObj, gs, nperm){
    genesets.names <- gs$names
    genesets <- gs$genesets
    assay <- assay(summarizedExprimentObj)
    DE.Analysis.res <- rowData(summarizedExprimentObj)
    DE.genes <- DE.Analysis.res %>% filter(p.value <= 0.05) %>% rownames(.)
    background.genes <- DE.Analysis.res %>% rownames(.)

    res <- lapply(1:length(genesets), function(i){
      gs <- genesets[[i]]
      wBallDraw <- intersect(gs, DE.genes) %>% length() - 1
      if (wBallDraw < 0) return(
                                data.frame(
                                          pathway = names(genesets)[i],
                                          p.value = 1,
                                          nDE = 0,
                                          stringsAsFactors = FALSE
                                          )
                                )
      wBall <- length(DE.genes)
      bBall <- length(background.genes) - length(DEGenes)
      ballDraw <- length(intersect(gs, background.genes))
      p.value <- 1 - phyper(wBallDraw, wBall, bBall, ballDraw)
      data.frame(
        pathway = names(genesets)[i],
        p.value = p.value,
        nDE = wBallDraw + 1,
        stringsAsFactors = FALSE
      )
    }) %>% do.call(what = rbind)

    finalRes <- data.frame(pathway = names(genesets) %>% unique(), stringsAsFactors = F)
    rownames(finalRes) <- finalRes$pathway
    finalRes$p.value <- 1
    finalRes$nDE <- 0

    finalRes[res$pathway %>% as.character(), 'p.value'] <- res$p.value
    finalRes[res$pathway %>% as.character(), 'nDE'] <- res$nDE
    finalRes$pFDR <- fdr
    colnames(res) <- c("p.value","pathway")

    res <- res %>% drop_na()
    res
}

#' @title Pathway Enrichment Analysis using fgsea
#' @description This function performs pathway analysis using fgsea (fast geneset enrichment analysis).
#' This function is used internally by runPathwayAnalysis.
#' @param summarizedExprimentObj The generated SummarizedExpriment object from DE analysis result.
#' @param genesets The genesets definition, ex. KEGG genesets.
#' @param nperm The number of permutations to run pathway analysis.
#' @return A dataframe of pathway analysis results
#' @details This function is used internally by runPathwayAnalysis.
#' @importFrom SummarizedExperiment SummarizedExperiment rowData
#' @importFrom dplyr %>% select
#' @importFrom fgsea fgsea
#' @importFrom tidyr drop_na
.runFgsea <- function(summarizedExprimentObj, gs, nperm){
    genesets.names <- gs$names
    genesets <- gs$genesets
    DE.Analysis.res <- rowData(summarizedExprimentObj)
    statistic <- DE.Analysis.res %>% select(statistic) %>% unlist() %>% as.vector()
    names(statistic) <- rownames(DE.Analysis.res)
    res <- fgsea(pathways = genesets,
                 stats = statistic, nperm = nperm)

    res <- res$pvals[,c('pathway', 'pval')] %>% drop_na()
    colnames(res) <- c("pathway", "p.value")

    return(res)
}

#' @title Pathway Enrichment Analysis using GSA
#' @description This function performs patwhay analysis using GSA method (GeneSet Analysis).
#' This function is used internally by runPathwayAnalysis.
#' @param summarizedExprimentObj The generated SummarizedExpriment object from DE analysis result.
#' @param genesets The genesets definition, ex. KEGG genesets.
#' @param nperm The number of permutations to run pathway analysis.
#' @return A dataframe of pathway analysis results
#' @details This function is used internally by runPathwayAnalysis.
#' @importFrom SummarizedExperiment SummarizedExperiment rowData
#' @importFrom dplyr %>% select
#' @importFrom GSA GSA
.runGSA <- function(summarizedExprimentObj, gs, nperm){ #Fix this
    genesets.names <- gs$names
    genesets <- gs$genesets
    res <- GSA(x = exprs, y = (group == "d") + 1, nperms=nperm, genesets=genesets, resp.type = "Two class unpaired",
               genenames = rownames(exprs), random.seed = 1,
               method = "maxmean")

    p.values <- cbind(res$pvalues.lo, res$pvalues.hi) %>% apply(1, min)
    p.values <- (p.values * 2) %>% data.frame(pathway = names(genesets))
    colnames(p.values) <- c('p.value', 'pathway')

    return(drop_na(p.values))
}

#' @title Pathway Enrichment Analysis
#' @description This function performs patwhay analysis using GSA method (GeneSet Analysis).
#' This function is used internally by runPathwayAnalysis.
#' @param summarizedExprimentObj The generated SummarizedExpriment object from DE analysis result.
#' @param genesets The genesets definition, ex. KEGG genesets from getGeneSets function.
#' @param method The pathway analsyis method, including ORA, fgsea, and GSA.
#' @param nperm The number of permutations to run pathway analysis.
#' @return A dataframe of pathway analysis result
#' @details This function is used internally by runPathwayAnalysis.
#' @importFrom SummarizedExperiment SummarizedExperiment rowData
#' @importFrom dplyr %>% select
runPathwayAnalysis <- function(summarizedExprimentObj, genesets, method = "ORA", nperm = 1000){
    analysisFunc <- switch(method,
                           ORA = .runORA,
                           fgsea = .runFgsea,
                           GSA = .runGSA
                          )

    analysisResult <- analysisFunc(summarizedExprimentObj, genesets$genesets, nperm)
    return(analysisResult)
}