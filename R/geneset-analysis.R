#' @title Geneset Enrichment Analysis using ORA
#' @description This function performs geneset analysis based on ORA (Over Representation Analysis).
#' This function is used internally by runGeneSetEnrichmentAnalysis.
#' @param summarizedExperiment The generated SummarizedExpriment object from DE analysis result.
#' @param genesets The genesets definition, ex. KEGG genesets from getGeneSets function.
#' @param pThreshold The p.value cutoff threashold.
#' @return A dataframe of geneset analysis results
#' @details This function is used internally by runGeneSetEnrichmentAnalysis.
#' @importFrom SummarizedExperiment SummarizedExperiment rowData assay
#' @importFrom dplyr %>% filter
#' @importFrom tidyr drop_na
#' @importFrom stats phyper
.runORA <- function(summarizedExperiment, genesets, pThreshold) {

    DE_data <- rowData(summarizedExperiment)
    if (is.null(DE_data) |
        dim(DE_data)[1] == 0 |
        dim(DE_data)[2] == 0) {
        stop("No differential analysis data is in input data.")
    }

    DE.genes <- DE_data[DE_data$p.value <= pThreshold,] %>%
        rownames()

    background.genes <- DE_data %>% rownames(.)

    GSOverlap <- sapply(genesets, function(gs) length(intersect(gs, background.genes)))
    DEOverlap <- sapply(genesets, function(gs) length(intersect(gs, DE.genes)))
    NoneDEInBackground <- length(background.genes) - length(DE.genes)
    Expected <- GSOverlap * length(DE.genes) / length(background.genes)

    pvals <- 1 - phyper(DEOverlap - 1, length(DE.genes), NoneDEInBackground, GSOverlap)
    ES <- DEOverlap / Expected

   data.frame(
        ID = names(genesets),
        p.value = pvals,
        score = ES,
        normalizedScore = ES,
        stringsAsFactors = FALSE
    )
}

#' @title Geneset Enrichment Analysis using fgsea
#' @description This function performs geneset analysis using fgsea (fast geneset enrichment analysis).
#' This function is used internally by runGeneSetEnrichmentAnalysis.
#' @param summarizedExperiment The generated SummarizedExpriment object from DE analysis result.
#' @param genesets The genesets definition, ex. KEGG genesets from getGeneSets function.
#' @param ... A list of other passed arguments to fgsea. See fgsea function.
#' @return A dataframe of geneset analysis results
#' @details This function is used internally by runGeneSetEnrichmentAnalysis.
#' @importFrom SummarizedExperiment SummarizedExperiment rowData assay
#' @importFrom dplyr %>%
#' @importFrom fgsea fgsea
#' @importFrom tidyr drop_na
.runFgsea <- function(summarizedExperiment, genesets,...) {

    DE_data <- rowData(summarizedExperiment)
    if (is.null(DE_data) |
        dim(DE_data)[1] == 0 |
        dim(DE_data)[2] == 0) {
        stop("No differential analysis data is in input data.")
    }

    statistic <- DE_data$statistic %>%
        unlist() %>%
        as.vector()
    names(statistic) <- rownames(DE_data)

    res <- fgsea(pathways = genesets, stats = statistic, ...)

    res <- res %>% drop_na()

    data.frame(
        ID = res$pathway,
        p.value = res$pval,
        score = res$ES,
        normalizedScore = res$NES,
        stringsAsFactors = FALSE
    )
}

#' @title Geneset Enrichment Analysis using GSA
#' @description This function performs geneset analysis using GSA method (GeneSet Analysis).
#' This function is used internally by runGeneSetEnrichmentAnalysis.
#' @param summarizedExperiment The generated SummarizedExpriment object from DE analysis result.
#' @param genesets The genesets definition, ex. KEGG genesets from getGeneSets function.
#' @param ... A list of other passed arguments to GSA. See GSA function.
#' @return A dataframe of geneset analysis results
#' @details This function is used internally by runGeneSetEnrichmentAnalysis.
#' @importFrom dplyr %>%
#' @importFrom SummarizedExperiment SummarizedExperiment rowData assay metadata
#' @importFrom GSA GSA
.runGSA <- function(summarizedExperiment, genesets, ...) {

    assay <- assay(summarizedExperiment)
    if (is.null(assay) |
        dim(assay)[1] == 0 |
        dim(assay)[2] == 0) {
        stop("No expression data is in input data.")
    }

    group_data <- metadata(summarizedExperiment)
    if (is.null(group_data)) {
        stop("The group data in colData cannot be empty.")
    }

    design_mtx <- group_data$DEAnalysis.design
    if (is.null(design_mtx) |
        dim(design_mtx)[1] == 0 |
        dim(design_mtx)[2] == 0) {
        stop("The design matrix cannot be empty.")
    }

    contrast_mtx <- group_data$DEAnalysis.contrast
    if (is.null(contrast_mtx) |
        dim(contrast_mtx)[1] == 0 |
        dim(contrast_mtx)[2] == 0) {
        stop("The contrast matrix cannot be empty.")
    }

    pair_data <- .extractPairInfo(design_mtx, contrast_mtx)
    group <- pair_data$group
    pairs <- pair_data$pair

    num_classes <- length(unique(group))
    resp_type = NULL

    if(num_classes == 1){
        stop("The group classes in design matrix must be at least two.")
    }else if(num_classes == 2 & !is.null(pairs)){
        resp_type <- "Two class paired"
        pairs_count <- unique(pairs)
        new_group <- NULL
        for(p in pairs_count){
            samples <- pairs[pairs == p] %>% names()
            if(group[[samples[1]]] == 1){
                tmp <- c(p, -p)
            }else{
                tmp <- c(-p, p)
            }
            names(tmp) <- samples
            new_group <- c(new_group, tmp)
        }
        group <- new_group
    }else if(num_classes == 2){
        resp_type <- "Two class unpaired"
    }else{
        resp_type <- "Multiclass"
    }

    gsa_res <- GSA(x = assay, y = group, genesets = genesets, resp.type = resp_type, genenames = rownames(assay), ...)

    pvalues <- cbind(gsa_res$pvalues.lo, gsa_res$pvalues.hi) %>% apply(1, min)
    res <- (pvalues * 2) %>% data.frame(ID = names(genesets), p.value = ., stringsAsFactors = FALSE)

    data.frame(
        ID = res$ID,
        p.value = res$p.value,
        score = gsa_res$GSA.scores,
        normalizedScore = gsa_res$GSA.scores,
        stringsAsFactors = FALSE
    )
}

#' @title Geneset Enrichment Analysis using KS or Wilcox
#' @description This function performs geneset analysis using KS or Wilcox test.
#' This function is used internally by runGeneSetEnrichmentAnalysis.
#' @param summarizedExperiment The generated SummarizedExpriment object from DE analysis result.
#' @param genesets The genesets definition, ex. KEGG genesets from getGeneSets function.
#' @param sTest The selected test to perform geneset analysis, either 'ks' or 'wilcox'.
#' @return A dataframe of geneset analysis results
#' @details This function is used internally by runGeneSetEnrichmentAnalysis.
#' @importFrom dplyr %>%
#' @importFrom SummarizedExperiment SummarizedExperiment rowData assay
#' @importFrom tidyr drop_na
#' @importFrom stats ks.test wilcox.test
.runKsWilcox <- function(summarizedExperiment, genesets, sTest) {

    DE_data <- rowData(summarizedExperiment)
    if (is.null(DE_data) |
        dim(DE_data)[1] == 0 |
        dim(DE_data)[2] == 0) {
        stop("No differential analysis data is in input data.")
    }

    ranks <- DE_data$statistic
    names(ranks) <- rownames(DE_data)

    test <- if (sTest == "ks") ks.test else wilcox.test

    background.genes <- rownames(DE_data)

    DEhit <- sapply(genesets, function(gs) ranks[background.genes[background.genes %in% gs]])
    DEmiss <- sapply(genesets, function(gs) ranks[background.genes[!background.genes %in% gs]])

    PA_res <- lapply(1:length(DEhit), function(i) {
        cur.DEhit <- DEhit[[i]]
        cur.DEmiss <- DEmiss[[i]]

        if (length(cur.DEhit) == 0 | length(cur.DEmiss) == 0) return(1)
        test(cur.DEhit, cur.DEmiss)$p.value
    }) %>%
        unlist() %>%
        data.frame(ID = names(genesets), p.value = ., stringsAsFactors = FALSE)

    PA_res <- PA_res %>% drop_na() %>% `rownames<-`(.$ID)

    oraRes <- .runORA(summarizedExperiment, genesets, pThreshold = 0.05)
    oraRes <- oraRes[oraRes$ID %in% PA_res$ID,] %>% `rownames<-`(.$ID)
    oraRes <- oraRes[PA_res$ID,]

    data.frame(
        ID = PA_res$ID,
        p.value = PA_res$p.value,
        score = oraRes$score,
        normalizedScore = oraRes$normalizedScore,
        stringsAsFactors = FALSE
    )
}

#' @title Geneset Enrichment Analysis
#' @description This function performs geneset analysis using either ORA, fgsea, GSA, ks, or wilcox approaches.
#' @param summarizedExperiment The generated SummarizedExpriment object from DE analysis result.
#' @param genesets The genesets definition, ex. KEGG genesets from getGeneSets function.
#' @param method The geneset analsyis method, including ORA, fgsea, GSA, ks, and wilcox.
#' @param ora.args A list of other passed arguments to ORA. pThreshold is used as p.value cutoff to pick DE genes.
#' @param fgsea.args A list of other passed arguments to fgsea. See fgsea function.
#' @param gsa.args A list of other passed arguments to GSA. See GSA function.
#' @return A dataframe of geneset analysis result
#' @importFrom SummarizedExperiment SummarizedExperiment rowData assay metadata
#' @importFrom dplyr %>%
#' @export
runGeneSetEnrichmentAnalysis <- function(summarizedExperiment, genesets, method = c("ora", "fgsea", "gsa", "ks", "wilcox"),
                                         ora.args = list(pThreshold = 0.05),
                                         fgsea.args = list(sampleSize = 101, minSize = 1, maxSize = Inf, eps = 1e-50, scoreType = "std",
                                                           nproc = 0, gseaParam = 1, BPPARAM = NULL, nPermSimple = 1000, absEps = NULL),
                                         gsa.args = list(method = "maxmean", random.seed = NULL, knn.neighbors = 10,
                                                         s0 = NULL, s0.perc = NULL, minsize = 15, maxsize = 500, restand = TRUE, restand.basis = "catalog",
                                                         nperms = 200, xl.mode = "regular", xl.time = NULL, xl.prevfit = NULL)
) {
    method <- match.arg(method)

    if(length(genesets) < 3 | !any(names(genesets) %in% c("database", "genesets", "names"))){
        stop("The genesets should be an object similar to returend object from getGeneSets().")
    }

    if(length(genesets[["genesets"]]) == 0){
        stop("There is no genesets.")
    }

    if(length(genesets[["genesets"]]) != length(genesets[["names"]])){
        stop("The length of genesets and their names do not match.")
    }

    gs.args <- NULL

    if (method == "fgsea") {

        fgsea.args.default <- list(sampleSize = 101, minSize = 1, maxSize = Inf, eps = 1e-50, scoreType = "std",
                                   nproc = 0, gseaParam = 1, BPPARAM = NULL, nPermSimple = 1000, absEps = NULL)

        if (any(!names(fgsea.args) %in% names(fgsea.args.default))) {
            stop("The names of arguments should be matched with fgsea definition.")
        }

        tmp <- fgsea.args
        fgsea.args <- fgsea.args.default
        fgsea.args[names(tmp)] <- fgsea.args.default[names(tmp)]

        if (!fgsea.args$scoreType %in% c("std", "pos", "neg")) {
            stop("The scoreType value is not valid.")
        }

        if (fgsea.args$minSize < 0 | fgsea.args$maxSize < 0) {
            stop("The minSize/maxSize cannot be negative.")
        }

        if (fgsea.args$nPermSimple <= 0) {
            stop("The number of permutations cannot be zero/negative.")
        }

        if (fgsea.args$nproc != 0 & is.null(fgsea.args$BPPARAM)) {
            stop("Set BPPARAM to use nproc workers.")
        }

        gs.args <- fgsea.args
    }

    if (method == "gsa") {
        gsa.args.default <- list(method = "maxmean", random.seed = NULL, knn.neighbors = 10,
                                 s0 = NULL, s0.perc = NULL, minsize = 15, maxsize = 500, restand = TRUE, restand.basis = "catalog",
                                 nperms = 200, xl.mode = "regular", xl.time = NULL, xl.prevfit = NULL)

        if (any(!names(gsa.args) %in% names(gsa.args.default))) {
            stop("The names of arguments should be matched with GSA definition.")
        }

        tmp <- gsa.args
        gsa.args <- gsa.args.default
        gsa.args[names(tmp)] <- gsa.args.default[names(tmp)]

        if (!gsa.args$method %in% c("maxmean", "mean", "absmean")) {
            stop("The method value in GSA is not valid.")
        }

        if (!gsa.args$restand.basis %in% c("catalog", "data")) {
            stop("The restand.basis value is not valid.")
        }

        if (!gsa.args$xl.mode %in% c("regular", "firsttime", "next20", "lasttime")) {
            stop("The xl.mode value is not valid.")
        }

        if (gsa.args$minsize < 0 | gsa.args$maxsize < 0) {
            stop("The minsize/maxsize cannot be negative.")
        }

        if (gsa.args$knn.neighbors < 0) {
            stop("The knn.neighbors cannot be negative.")
        }

        if (gsa.args$nperms <= 0) {
            stop("The number of permutations cannot be zero/negative.")
        }

        gs.args <- gsa.args
    }

    if (method == "ora") {
        ora.args.default <- list(pThreshold = 0.05)

        if (any(!names(ora.args) %in% names(ora.args.default))) {
            stop("The names of arguments should be matched with ORA definition.")
        }

        tmp <- ora.args
        ora.args <- ora.args.default
        ora.args[names(tmp)] <- ora.args.default[names(tmp)]

        if (ora.args$pThreshold < 0 | ora.args$pThreshold > 1) {
            stop("The pThreshold must be between zero and one.")
        }

        gs.args <- ora.args
    }

    if (method %in% c("ks", "wilcox")) {
        gs.args$sTest = method
    }

    gs.args$summarizedExperiment <- summarizedExperiment
    gs.args$genesets <- genesets[["genesets"]]

    methodFnc <- switch(method,
                        ora = .runORA,
                        fgsea = .runFgsea,
                        gsa = .runGSA,
                        ks = .runKsWilcox,
                        wilcox = .runKsWilcox)

    result <- do.call(methodFnc, gs.args)

    if (is.null(result)) {
        stop("There is an error in geneset analysis procedure.")
    }

    result$sample.size = ncol(assay(summarizedExperiment))
    genesets_names <- genesets[["names"]]
    result$name = genesets_names[as.character(result$ID)]

    return(result)
}