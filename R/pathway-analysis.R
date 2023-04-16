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
        ES = ES,
        NES = ES,
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

.est.s0 <- function(tt, sd, s0.perc = seq(0, 1, by = .05)) {
    ## estimate s0 (exchangeability) factor for denominator.
    ## returns the actual estimate s0 (not a percentile)

    br = unique(quantile(sd, seq(0, 1, len = 101)))
    nbr = length(br)

    a <- cut(sd, br, labels = F)
    a[is.na(a)] <- 1
    cv.sd <- rep(0, length(s0.perc))

    for (j in 1:length(s0.perc)) {
        w <- quantile(sd, s0.perc[j])
        w[j == 1] <- 0
        tt2 <- tt * sd / (sd + w)
        tt2[tt2 == Inf] = NA
        sds <- rep(0, nbr - 1)

        for (i in 1:(nbr - 1)) {
            sds[i] <- mad(tt2[a == i], na.rm = TRUE)
        }

        cv.sd[j] <- sqrt(var(sds)) / mean(sds)
    }

    o = (1:length(s0.perc))[cv.sd == min(cv.sd)]

    s0.hat = quantile(sd[sd != 0], s0.perc[o])

    return(list(s0.perc = s0.perc, cv.sd = cv.sd, s0.hat = s0.hat))
}

.GSA.func = function(tt,
                     gs.mat,
                     method = c("maxmean", "mean", "absmean"),
                     resp.type = c("Quantitative", "Two class unpaired", "Multiclass", "Two class paired")) {

    resp.type <- match.arg(resp.type)
    method <- match.arg(method)

    gs.ind <- (1:dim(gs.mat)[1])

    ngenes <- gs.mat %>% apply(1, function (cur.row) length(unique(cur.row)))

    mean.all = mean(tt)
    sd.all = sqrt(var(tt))
    mean.pos = mean(tt * (tt > 0))
    sd.pos = sqrt(var(tt * (tt > 0)))
    mean.neg = -mean(tt * (tt < 0))
    sd.neg = sqrt(var(tt * (tt < 0)))

    tt2 = c(tt, 0)
    ttt = matrix(tt2[gs.mat], nrow = nrow(gs.mat))

    if (method == "maxmean") {
        s2 = abs(ttt)
        rpos = rowSums((ttt + s2) / 2) / ngenes
        rneg = rowSums((s2 - ttt) / 2) / ngenes

        rpos[is.na(rpos)] = 0
        rneg[is.na(rneg)] = 0

        rpos = (rpos - mean.pos) / (sd.pos)

        if (resp.type != "Multiclass") {
            rneg = (rneg - mean.neg) / (sd.neg)
        }
        rr <- pmax(rpos, rneg)
        rr[rneg > rpos] = -1 * rr[rneg > rpos]
    }

    if (method == "mean") {
        rr = rowSums(ttt) / ngenes
        rr <- (rr - mean.all) / (sd.all / sqrt(ngenes))
    }

    if (method == "absmean") {
        rr <- rowSums(abs(ttt)) / ngenes
        rr <- (rr - mean.all) / (sd.all / sqrt(ngenes))
    }

    rrr = rep(NA, length(genesets))
    rrr[gs.ind] = rr
    rrr
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
.runGSA <- function(DESummarizedExperiment, genesets, nperm = 1000,
                    method = c("maxmean", "mean", "absmean"),
                    resp.type = c("Quantitative", "Two class unpaired", "Multiclass", "Two class paired"), maxsize = 500) {
    assay <- DESummarizedExperiment %>% assay()

    k = length(genesets)
    gs.ind = (1:k)

    ngenes = rep(NA, length(gs.ind))

    gs.mat = matrix(nrow(assay) + 1, nrow = length(gs.ind), ncol = maxsize)
    for (i in gs.ind) {
        gene.set = match(genesets[[i]], rownames(assay))
        gene.set = gene.set[!is.na(gene.set)]
        if (length(gene.set) > 0) {
            gs.mat[i, 1:length(gene.set)] = gene.set
            ngenes[i] = length(gene.set)
        }
    }
    gs.mat = gs.mat[, 1:max(ngenes, na.rm = TRUE)]

    tt <- rowData(DESummarizedExperiment)$statistic
    sd <- rowData(DESummarizedExperiment)$dispersion
    s0 <- .est.s0(tt, sd, s0.perc = seq(0, 1, by = .05))

    observerd.tt <- tt * sd / (sd + s0$s0.hat)
    observerd.score <- .GSA.func(observerd.tt, gs.mat, method = "maxmean", resp.type = "Two class unpaired")
    names(observerd.score) <- names(genesets)

    DEMethod <- metadata(DESummarizedExperiment)$DEAnalysis.method
    DEFunc <- switch(
      DEMethod,
      "DESeq2" = .runDESeq2,
      "edgeR" = .runEdgeR,
      "limma" = .runLimma
    )

    design <- metadata(DESummarizedExperiment)$DEAnalysis.design
    contrast <- metadata(DESummarizedExperiment)$DEAnalysis.contrast
    toPermuteGroup <- names(which(contrast[, 1] !=0))
    observerd.score.normalized = -qnorm(1 - pt(observerd.score, df = nrow(design) - 2))

    permRes <- lapply(1:1000, function(i) {
        set.seed(i)
        message(i)
        designPerm <- design[sample(rownames(design)),] %>% `rownames<-`(rownames(design))
        sePerm <- DEFunc(assay, designPerm, contrast)
        perm.tt <- sePerm$statistic
        perm.sd <- sePerm$dispersion

        perm.tt.adj <- perm.tt * perm.sd / (perm.sd + s0$s0.hat)
        perm.score <- .GSA.func(perm.tt.adj, gs.mat, method = "maxmean")
        perm.score.normalized = -qnorm(1 - pt(perm.score, df = nrow(design) - 2))
        data.frame(
          pathway = names(genesets),
          perm.score = perm.score,
          perm.norm.score = perm.score.normalized
        )
    }) %>% do.call(what = rbind)

    pvalues.Df <- permRes %>% group_by(pathway) %>% group_split() %>% lapply(function(data){
        null_scores <- data$perm.score
        original.score <- observerd.score[names(observerd.score) == data$pathway[1]]
        pval.hi <- length(null_scores[null_scores > original.score])/nperm
        pval.lo <- length(null_scores[null_scores < original.score])/nperm
        #pval <- (length(null_scores[null_scores >= original.score]) + 1) / (nperm + 1)
        data.frame(
          pathway = data$pathway[1],
          p.value.hi = pval.hi,
          p.value.lo = pval.lo,
          stringsAsFactors = FALSE
        )
    }) %>% do.call(what = rbind)

    pvalues.Df$p.value <- cbind(pvalues.Df$p.value.lo, pvalues.Df$p.value.hi) %>% apply(1, min)
    pvalues.Df$p.value <- pvalues.Df$p.value *2

    exprs <- assay
    group <- c(rep("c",5), rep("d",5))
    names(group) <- paste0("sample", seq(1:10))
    res <- GSA::GSA(x = exprs, y = (group == "d") + 1, nperms = 1000, genesets = genesets, resp.type = "Two class unpaired",
               genenames = rownames(exprs), random.seed = 1,
               method = "maxmean")

    gsa.res <- cbind(res$pvalues.lo, res$pvalues.hi) %>% apply(1, min)
    gsa.res <- (gsa.res * 2) %>% data.frame(pathway = names(genesets))
    colnames(gsa.res) <- c('p.value', 'pathway')
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
runPathwayAnalysis <- function(DESummarizedExperiment, genesets, method = "ORA", nperm = 1000) {
    analysisFunc <- switch(method,
                           ORA = .runORA,
                           fgsea = .runFgsea,
                           GSA = .runGSA
    )

    analysisResult <- analysisFunc(summarizedExprimentObj, genesets$genesets, nperm)
    return(analysisResult)
}