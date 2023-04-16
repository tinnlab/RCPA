#' @title Pathway Enrichment Analysis using ORA
#' @description This function performs pathway analysis based on ORA (Over Representation Analysis).
#' This function is used internally by runPathwayAnalysis.
#' @param DESummarizedExperiment The generated SummarizedExpriment object from runDEAnalysis.
#' @param genesets The genesets definition, ex. KEGG genesets.
#' @return A dataframe of pathway analysis results
#' @details This function is used internally by runPathwayAnalysis.
#' @importFrom SummarizedExperiment SummarizedExperiment rowData assay
#' @importFrom dplyr %>% filter
#' @importFrom tidyr drop_na
.runORA <- function(DESummarizedExperiment, genesets, pThreshold = 0.05) {
    DE.Analysis.res <- rowData(DESummarizedExperiment) %>% as.data.frame()
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
        score = ES,
        normalizedScore = ES,
        stringsAsFactors = FALSE
    )
}

#' @title Pathway Enrichment Analysis using fgsea
#' @description This function performs pathway analysis using fgsea (fast geneset enrichment analysis).
#' This function is used internally by runPathwayAnalysis.
#' @param DESummarizedExperiment The generated SummarizedExpriment object from DE analysis result.
#' @param genesets The genesets definition, ex. KEGG genesets.
#' @param nperm The number of permutations to run pathway analysis.
#' @return A dataframe of pathway analysis results
#' @details This function is used internally by runPathwayAnalysis.
#' @importFrom SummarizedExperiment SummarizedExperiment rowData
#' @importFrom dplyr %>% select
#' @importFrom fgsea fgsea
#' @importFrom tidyr drop_na
.runFgsea <- function(DESummarizedExperiment, genesets, nperm = 1000) {
    DE.Analysis.res <- rowData(DESummarizedExperiment) %>% as.data.frame()
    statistic <- DE.Analysis.res %>%
        .$statistic %>%
        unlist() %>%
        as.vector()
    names(statistic) <- rownames(DE.Analysis.res)

    res <- fgsea(pathways = genesets,
                 stats = statistic, nperm = nperm)

    res <- res %>% drop_na()

    data.frame(
      pathway = res$pathway,
      p.value = res$pval,
      ES = res$ES,
      NES = res$NES,
      stringsAsFactors = FALSE
    )
}

#' @title Pathway Enrichment Analysis using GSA
#' @description This function performs patwhay analysis using GSA method (GeneSet Analysis).
#' This function is used internally by runPathwayAnalysis.
#' @param DESummarizedExperiment The generated SummarizedExpriment object from DE analysis result.
#' @param genesets The genesets definition, ex. KEGG genesets.
#' @param nperm The number of permutations to run pathway analysis.
#' @return A dataframe of pathway analysis results
#' @details This function is used internally by runPathwayAnalysis.
#' @importFrom SummarizedExperiment SummarizedExperiment rowData
#' @importFrom dplyr %>% select
#' @importFrom GSA GSA
.runGSA <- function(assay, group_data, genesets, gsa.args = list()) {

    if(is.null(group_data)){
        stop("The group data in colData cannot be empty.")
    }

    DE_metadata <- group_data$DEAnaysis
    design_mtx <- DE_metadata$design
    if(is.null(design_mtx) | dim(design_mtx)[1] == 0 | dim(design_mtx)[2] == 0){
        stop("The design matrix cannot be empty.")
    }

    contrast_mtx <- DE_metadata$contrast
    if(is.null(contrast_mtx) | dim(contrast_mtx)[1] == 0 | dim(contrast_mtx)[2] == 0){
        stop("The contrast matrix cannot be empty.")
    }

    pair_data <- .extractPairInfo(design_mtx, contrast_mtx)
    if(is.na())

    gsa_res <- GSA(x = assay, y = (group == "d") + 1, nperms = gsa.args$nperms, genesets = genesets, resp.type = gsa.args$resp.type,
               genenames = rownames(assay), random.seed = gsa.args$random.seed, knn.neighbors = gsa.args$knn.neighbors,
               method = gsa.args$method, s0 = gsa.args$s0, s0.perc = gsa.args$s0.perc, minsize = gsa.args$minsize, maxsize = gsa.args$maxsize, restand = gsa.args$restand,
               restand.basis = gsa.args$restand.basis, xl.mode = gsa.args$xl.mode, xl.time = gsa.args$xl.time, xl.prevfit = gsa.args$xl.prevfit)

    p.values <- cbind(gsa_res$pvalues.lo, gsa_res$pvalues.hi) %>% apply(1, min)
    p.values <- (p.values * 2) %>% data.frame(pathway = names(genesets))
}

#' @title Pathway Enrichment Analysis using GSA
#' @description This function performs patwhay analysis using GSA method (GeneSet Analysis).
#' This function is used internally by runPathwayAnalysis.
#' @param DESummarizedExperiment The generated SummarizedExpriment object from DE analysis result.
#' @param genesets The genesets definition, ex. KEGG genesets.
#' @param nperm The number of permutations to run pathway analysis.
#' @return A dataframe of pathway analysis results
#' @details This function is used internally by runPathwayAnalysis.
#' @importFrom SummarizedExperiment SummarizedExperiment rowData
#' @importFrom dplyr %>% select
#' @importFrom GSA GSA
.runKSWilcox <- function(DE_data, genesets, method){
    ranks <- DE_data$statistic
    names(ranks) <- rownames(DE_data)

    test <- if (method == "ks") ks.test else wilcox.test

    background.genes <- names(DE_genes)

    PA_res <- lapply(genesets, function(gs){

        DEhit <- ranks[background.genes[background.genes %in% gs]]
        DEmiss <- ranks[background.genes[!background.genes %in% gs]]

        if (length(DEhit) == 0 | length(DEmiss) == 0) return(NA)

        test(DEhit, DEmiss)$p.value
    }) %>% unlist() %>% data.frame(pathway=names(genesets))

    DEhit <- sapply(genesets, function(gs) ranks[background.genes[background.genes %in% gs]])
    DEmiss <- sapply(genesets, function(gs) ranks[background.genes[!background.genes %in% gs]])

    PA_res <- lapply(1:length(DEhit), function (i){
        cur.DEhit <- DEhit[[i]]
        cur.DEmiss <- DEmiss[[i]]

        if(length(cur.DEhit) == 0 | length(cur.DEmiss) == 0) return(1)
        test(cur.DEhit, cur.DEmiss)$p.value
    }) %>% unlist() %>% data.frame(pathway=names(genesets))

    PA_res <- PA_res %>% drop_na()

    data.frame(
      pathway = PA_res$pathway,
      p.value = PA_res$p.value,
      score = rep(0, nrow(PA_res)),
      normalizedScore = rep(0, nrow(PA_res)),
      stringsAsFactors = FALSE
    )
}

#' @title Pathway Enrichment Analysis using GSA
#' @description This function performs patwhay analysis using GSA method (GeneSet Analysis).
#' This function is used internally by runPathwayAnalysis.
#' @param DESummarizedExperiment The generated SummarizedExpriment object from DE analysis result.
#' @param genesets The genesets definition, ex. KEGG genesets.
#' @param nperm The number of permutations to run pathway analysis.
#' @return A dataframe of pathway analysis results
#' @details This function is used internally by runPathwayAnalysis.
#' @importFrom SummarizedExperiment SummarizedExperiment rowData
#' @importFrom dplyr %>% select
#' @importFrom GSA GSA
.runSPIA <- function(){
    .SPIAMod()
}

#' @title Pathway Enrichment Analysis using GSA
#' @description This function performs patwhay analysis using GSA method (GeneSet Analysis).
#' This function is used internally by runPathwayAnalysis.
#' @param DESummarizedExperiment The generated SummarizedExpriment object from DE analysis result.
#' @param genesets The genesets definition, ex. KEGG genesets.
#' @param nperm The number of permutations to run pathway analysis.
#' @return A dataframe of pathway analysis results
#' @details This function is used internally by runPathwayAnalysis.
#' @importFrom SummarizedExperiment SummarizedExperiment rowData
#' @importFrom dplyr %>% select
#' @importFrom GSA GSA
.runCePaORA <- function(){

}

#' @title Pathway Enrichment Analysis using GSA
#' @description This function performs patwhay analysis using GSA method (GeneSet Analysis).
#' This function is used internally by runPathwayAnalysis.
#' @param DESummarizedExperiment The generated SummarizedExpriment object from DE analysis result.
#' @param genesets The genesets definition, ex. KEGG genesets.
#' @param nperm The number of permutations to run pathway analysis.
#' @return A dataframe of pathway analysis results
#' @details This function is used internally by runPathwayAnalysis.
#' @importFrom SummarizedExperiment SummarizedExperiment rowData
#' @importFrom dplyr %>% select
#' @importFrom GSA GSA
.runCePaGSA <- function(){

}

#' @title Pathway Enrichment Analysis
#' @description This function performs patwhay analysis using GSA method (GeneSet Analysis).
#' This function is used internally by runPathwayAnalysis.
#' @param DESummarizedExperiment The generated SummarizedExpriment object from DE analysis result.
#' @param genesets The genesets definition, ex. KEGG genesets from getGeneSets function.
#' @param method The pathway analsyis method, including ORA, fgsea, and GSA.
#' @param nperm The number of permutations to run pathway analysis.
#' @return A dataframe of pathway analysis result
#' @details This function is used internally by runPathwayAnalysis.
#' @importFrom SummarizedExperiment SummarizedExperiment rowData
#' @importFrom dplyr %>% select
runGeneSetEnrichmentAnalysis <- function(summarizedExperiment, genesets, method = c("ORA", "fgsea", "GSA", "KS", "wilcox"),
                             ora.args = list(pThreshold = 0.05),
                             fgsea.args = list(sampleSize = 101, minSize = 1, maxSize = Inf, eps = 1e-50, scoreType = "std",
                                nproc = 0, gseaParam = 1, BPPARAM = NULL, nPermSimple = 1000, absEps = NULL),
                             gsa.args = list(method = "maxmean", resp.type = "Two class unpaired", random.seed = NULL, knn.neighbors = 10,
                                s0 = NULL, s0.perc = NULL, minsize = 15, maxsize = 500, restand = TRUE, restand.basis = "catalog",
                                nperms = 200, xl.mode = "regular", xl.time = NULL,xl.prevfit = NULL)
){
    method <- match.arg(method)

    assay <- assay(summarizedExperiment)
    if(is.null(assay) | dim(assay)[1] == 0 | dim(assay)[2] == 0){
        stop("No expression data is in input data.")
    }

    DE_data <- rowData(summarizedExperiment)
    if(is.null(DE_data) | dim(DE_data)[1] == 0 | dim(DE_data)[2] == 0){
        stop("No differential analysis data is in input data.")
    }

    group_data <- colData(summarizedExperiment)

    result <- switch(method,
           ora = .runORA(DE_data, genesets, ora.args),
           fgsea = .runFgsea(DE_data, genesets, fgsea.args),
           gsa = .runGSA(assay, group_data, genesets, gsa.args),
           ks = .runKSWilcox(DE_data, genesets, method = "ks"),
           wilcox = .runKSWilcox(DE_data, genesets, method = "wilcox")
    )

    if(is.null(result)){
        stop("There is an error in geneset analysis procedure.")
    }

    return(result)
}

runPathwayAnalysis <- function(summarizedExperiment, network, method = c("SPIA", "cepaORA", "cepaGSA"),
                               spia.args = list(de=NULL, all=NULL, organism="hsa", data.dir=NULL, pathids=NULL,
                                 nB=2000, verbose=TRUE, beta=NULL, combine="fisher"),
                              cepaORA.args = list(dif = NULL, bk = NULL, mat = NULL, label = NULL, pc, cen = default.centralities,
                                cen.name = sapply(cen, function(x) ifelse(mode(x) == "name", deparse(x), x)),
                                nlevel = "tvalue_abs", plevel = "mean", iter = 1000),
                              cepaGSA.args = list(
                                ...
                              )
){
    method <- match.arg(method)

    assay <- assay(summarizedExperiment)
    if(is.null(assay) | dim(assay)[1] == 0 | dim(assay)[2] == 0){
        stop("No expression data is in input data.")
    }

    if(is.null(network)){
        stop("The network needs to be prepared before running pathway enrichment analysis.")
    }

    pc <- makeCePaPathwayCat("hsa")

    result <- switch(method,
                     SPIA = .runSPIA(),
                     cepaORA = .runCePaORA(),
                     cepaGSA = .runCePaGSA()
                     )

    if(is.null(result)){
        stop("There is an error in pathway analysis procedure.")
    }

    return(result)
}