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

#' @title Pathway Enrichment Analysisbiological_process
#' @description This function performs patwhay analysis using GSA method (GeneSet Analysis).
#' This function is used internally by runPathwayAnalysis.
#' @param summarizedExprimentObj The generated SummarizedExpriment object from DE analysis result.
#' @param genesets The genesets definition, ex. KEGG genesets from getGeneSets function.
#' @param method The pathway analsyis method, including ORA, fgsea, and GSA.
#' @param nperm The number of permutations to run pathway analysis.
#' @return A dataframe of pathway analysis result
#' @examples
#' \dontrun{
#' # Loading libraries
#' library(hgu133plus2.db)
#' library(AnnotationDbi)
#' library(SummarizedExperiment)
#' library(limma)
#' library(RCPA)
#' # Set seed for reproducibilty
#' set.seed(123)
#' exprs <- round(matrix(2^abs(rnorm(100000, sd = 4)), nrow = 10000, ncol = 10))
#' # Assign gene names
#' rownames(exprs) <- sample(keys(hgu133plus2.db, keytype = "PROBEID"), nrow(exprs), replace = FALSE)
#' # Assign sample names
#' colnames(exprs) <- paste0("sample", 1:10)
#' # Generate control and condition samples
#' controlSamples <- paste0("sample", 1:5)
#' conditionSamples <- paste0("sample", 6:10)
#' # Get colData
#' colData <- data.frame(
#'     row.names = colnames(exprs),
#'     group = factor(c(rep("control", length(controlSamples)), rep("condition", length(conditionSamples)))),
#'     pair = factor(c(seq_along(controlSamples), seq_along(conditionSamples)))
#' )
#' 
#' # Construct summarizedExperiment object
#' summarizedExperiment <- SummarizedExperiment(
#'     assays = list(counts = exprs),
#'     colData = colData
#' )
#' # Construct design and contrast tables
#' # control vs condition
#' # design <- model.matrix(~0 + group, data = colData)
#' # contrast <- makeContrasts("groupcondition-groupcontrol", levels = design)
#' 
#' # Perform DE analysis
#' DERes <- runDEAnalysis(summarizedExperiment, method = "limma", design, contrast, annotation = 'GPL570')
#' 
#' # Simulate genesets
#' gs <- lapply(1:100, function(x) {
#' sample(rownames(DERes), runif(1, 100, 500))
#' })
#' # Set names to genesets
#' names(gs) <- paste0("geneset", 1:100)
#' # Set descriptions to genesets
#' gs_fullNames <- paste0("path:", 1:100)
#' names(gs_fullNames) <- names(gs)
#' # Store genesets information in a list
#' genesets <- list(
#'     database = "TEST",
#'     genesets = gs,
#'     names = gs_fullNames
#' )
#' result <- runGeneSetEnrichmentAnalysis(DERes, genesets, method = "ora")
#' }
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