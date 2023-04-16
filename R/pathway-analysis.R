#' @title Pathway Enrichment Analysis using ORA
#' @description This function performs pathway analysis based on ORA (Over Representation Analysis).
#' This function is used internally by runPathwayAnalysis.
#' @param summarizedExperiment The generated SummarizedExpriment object from runDEAnalysis.
#' @param genesets The genesets definition, ex. KEGG genesets.
#' @return A dataframe of pathway analysis results
#' @details This function is used internally by runPathwayAnalysis.
#' @importFrom SummarizedExperiment SummarizedExperiment rowData assay
#' @importFrom dplyr %>% filter
#' @importFrom tidyr drop_na
.runORA <- function(summarizedExperiment, genesets, pThreshold = 0.05) {
  DE.Analysis.res <- rowData(summarizedExperiment) %>% as.data.frame()
  DE.genes <- DE.Analysis.res %>%
    filter(.$p.value <= pThreshold) %>%
    rownames()
  background.genes <- DE.Analysis.res %>% rownames(.)
  
  GSOverlap <- sapply(genesets, function(gs) length(intersect(gs, background.genes)))
  DEOverlap <- sapply(genesets, function(gs) length(intersect(gs, DE.genes)))
  NoneDEInBackground <- length(background.genes) - length(DE.genes)
  Expected <- GSOverlap * length(DE.genes) / length(background.genes)
  
  pvals <- 1 - phyper(DEOverlap - 1, length(DE.genes), NoneDEInBackground, GSOverlap)
  ES <- DEOverlap / Expected
  
  data.frame(
    pathway = names(genesets),
    p.value = pvals,
    ES = ES,
    NES = ES
  )
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
.runFgsea <- function(summarizedExprimentObj, gs, nperm) {
  genesets <- gs$genesets
  DE.Analysis.res <- rowData(summarizedExprimentObj)
  statistic <- DE.Analysis.res %>%
    select(statistic) %>%
    unlist() %>%
    as.vector()
  names(statistic) <- rownames(DE.Analysis.res)
  res <- fgsea(pathways = genesets,
               stats = statistic, nperm = nperm)
  
  res <- res$pvals[, c('pathway', 'pval')] %>% drop_na()
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
.runGSA <- function(summarizedExprimentObj, gs, nperm) { #Fix this
  genesets.names <- gs$names
  genesets <- gs$genesets
  res <- GSA(x = exprs, y = (group == "d") + 1, nperms = nperm, genesets = genesets, resp.type = "Two class unpaired",
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
runPathwayAnalysis <- function(summarizedExprimentObj, genesets, method = "ORA", nperm = 1000) {
  analysisFunc <- switch(method,
                         ORA = .runORA,
                         fgsea = .runFgsea,
                         GSA = .runGSA
  )
  
  analysisResult <- analysisFunc(summarizedExprimentObj, genesets$genesets, nperm)
  return(analysisResult)
}