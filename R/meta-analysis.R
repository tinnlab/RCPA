#' @title Combine PValues
#' @description This function combines P-Values based on the selected method.
#' This function is used internally by combinePathwayAnalysisResults.
#' @param inputData A dataframe of pathway analysis results from multiple studies.
#' @param method The selected method to combine PValues, including Fisher, Stouffer, addCLT, GeoMean, and minP.
#' @return A dataframe of meta analysis results.
#' @details This function is used internally by combinePathwayAnalysisResults.
#' @importFrom dplyr %>%
.combinePvalues <- function(inputData, method) {

    combineFunc <- switch(method,
                          fisher = .runFisher,
                          stouffer = .runStouffer,
                          minP = min,
                          addCLT = .runAddCLT,
                          geoMean = .runGeoMean
    )

    metaPvalRes <- inputData %>%
        group_by(ID) %>%
        group_split() %>%
        lapply(function(data) {
            pValues <- data$p.value
            meta.Pval <- combineFunc(pValues)
            data.frame(
                ID = data$ID[1],
                p.value = meta.Pval,
                stringsAsFactors = FALSE
            )
        }) %>%
        do.call(what = rbind)

    metaScoreRes <- .combineEnrichmentScores(inputData)

    metaRes <- metaPvalRes
    metaRes$score <- metaRes$count <- 0
    metaRes$score <- metaScoreRes$score[match(metaScoreRes$ID, metaRes$ID)]
    #metaRes$score.sd <- metaScoreRes$score.sd[match(metaScoreRes$ID, metaRes$ID)]
    metaRes$count <- metaScoreRes$count[match(metaScoreRes$ID, metaRes$ID)]

    metaRes <- metaRes[, c("ID", "p.value", "score", "count")]

    return(metaRes)
}

#' @title Combine P-Valuse using Fisher
#' @description This function combines P-Values based on Fisher method.
#' This function is used internally by .combinePvalues.
#' @param pvals The vector of P-Values to be combined.
#' @return A combined P-Value
#' @details This function is used internally by .combinePvalues.
#' @importFrom stats qnorm pchisq
.runFisher <- function(pvals) {
    pvals[pvals == 0] <- .Machine$double.eps
    p.value <- pchisq(-2 * sum(log(pvals)), df = 2 * length(pvals), lower.tail = FALSE)
    return(p.value)
}

#' @title Combine P-Valuse using Stouffer
#' @description This function combines P-Values based on Stouffer method.
#' This function is used internally by .combinePvalues.
#' @param pvals The vector of P-Values to be combined.
#' @return A combined P-Value
#' @details This function is used internally by .combinePvalues.
#' @importFrom stats pnorm qnorm
.runStouffer <- function(pvals) {
    pvals[pvals == 0] <- .Machine$double.eps
    p.value <- pnorm(sum(qnorm(pvals)) / sqrt(length(pvals)))
    return(p.value)
}

#' @title Combine P-Valuse using addCLT
#' @description This function combines P-Values based on addCLT method.
#' This function is used internally by .combinePvalues.
#' @param pvals The vector of P-Values to be combined.
#' @return A combined P-Value
#' @details This function is used internally by .combinePvalues.
#' @importFrom stats pnorm
.runAddCLT <- function(pvals) {
    pvals <- pvals[!is.na(pvals)]
    pvals[pvals == 0] <- .Machine$double.eps
    n <- length(pvals)
    p.value <- 1
    if (n <= 20) {
        x <- sum(pvals)
        p.value <- 1 / factorial(n) * sum(sapply(0:floor(x), function(k) (-1)^k * choose(n, k) * (x - k)^(n)))
    }else {
        p.value <- pnorm(sum(pvals), n / 2, sqrt(n / 12), lower.tail = TRUE)
    }
    return(p.value)
}

#' @title Combine P-Valuse using Geometric Mean
#' @description This function combines P-Values by computing their geometric mean.
#' This function is used internally by .combinePvalues.
#' @param pvals The vector of P-Values to be combined.
#' @return A combined P-Value
#' @details This function is used internally by .combinePvalues.
.runGeoMean <- function(pvals) {
    pvals[pvals == 0] <- .Machine$double.eps
    p.value <- exp(mean(log(pvals)))
    return(p.value)
}

#' @title Combine Normalized Enrichment Scores
#' @description This function combines normalized enrichment scores based on REML method.
#' This function is used internally by combinePathwayAnalysisResults and .combinePvalues.
#' @param inputData A dataframe of pathway analysis results from multiple studies.
#' @return A dataframe of meta analysis results.
#' @details This function is used internally by combinePathwayAnalysisResults and .combinePvalues.
#' @importFrom meta metagen
#' @importFrom dplyr %>%
.combineEnrichmentScores <- function(inputData) {

    metaRes <- inputData %>%
        group_by(ID) %>%
        group_split() %>%
        lapply(function(data) {
            data$normalizedScore.sd <- abs((data$normalizedScore - ifelse(data$normalizedScore > 0, 1, -1)) / qnorm(data$p.value))
            meta.res <- metagen(data = data,
                                studlab = ID,
                                TE = normalizedScore,
                                seTE = normalizedScore.sd,
                                sm = "SMD",
                                n.e = sampleSize,
                                method.tau = "REML",
                                hakn = TRUE)

            normalizedScore.combined <- meta.res$TE.fixed
            normalizedScore.combined.sd <- meta.res$seTE.fixed

            pval <- pnorm((ifelse(normalizedScore.combined > 0, 1, -1) - normalizedScore.combined) / normalizedScore.combined.sd)
            if (normalizedScore.combined < 0) pval <- 1 - pval

            data.frame(
                ID = data$ID[1],
                p.value = pval,
                score = normalizedScore.combined,
                #score.sd = normalizedScore.combined.sd,
                count = nrow(data),
                stringsAsFactors = F
            )
        }) %>%
        do.call(what = rbind) %>%
        as.data.frame()

    metaRes <- metaRes[, c("ID", "p.value", "score", "count")]

    return(metaRes)
}

#' @title Perform Meta Analysis
#' @description This function performs meta analysis on multiple pathway analysis results.
#' @param PAResults A list of data frames obtained from runPathwayAnalysis.
#' @param method A method used to combine pathway analysis results
#' @return A dataframe of meta analysis results including combined normalized score and combined p-value for each pathway.
#' @examples
#' \dontrun{
#' #' # Loading libraries
#' library(SummarizedExperiment)
#' library(limma)
#' library(RCPA)
#' data("data")
#' # Get affymetrix dataset
#' affyDataset <- data$affyDataset
#' # Create the analysis design
#' affyDesign <- model.matrix(~0 + condition + region, data = colData(affyDataset))
#' affyContrast <- limma::makeContrasts("conditionalzheimer-conditionnormal", levels=affyDesign)
#' # Perform DE analysis affymetrix dataset
#' affyDEExperiment <- RCPA::runDEAnalysis(affyDataset, 
#'                                         method = "limma", 
#'                                         design = affyDesign, 
#'                                         contrast = affyContrast, 
#'                                         annotation = "GPL570")
#' 
#' agilDataset <- data$agilDataset
#' # Create the analysis design
#' agilDesign <- model.matrix(~0 + condition, data = colData(agilDataset))
#' agilContrast <- limma::makeContrasts(conditionalzheimer-conditionnormal, levels=agilDesign)
#' # Perform genID mapping
#' GPL4133Anno <- GEOquery::dataTable(GEOquery::getGEO("GPL4133"))@table
#' GPL4133GeneMapping <- data.frame(FROM = GPL4133Anno$SPOT_ID, 
#'                                  TO = as.character(GPL4133Anno$GENE), 
#'                                  stringsAsFactors = F)
#' GPL4133GeneMapping <- GPL4133GeneMapping[!is.na(GPL4133GeneMapping$TO), ]
#' # Perform DE analysis for agilent dataset
#' agilDEExperiment <- RCPA::runDEAnalysis(agilDataset, 
#'                                         method = "limma", 
#'                                         design = agilDesign, 
#'                                         contrast = agilContrast, 
#'                                         annotation = GPL4133GeneMapping)
#' # Get RNA-Seq dataset
#' RNASeqDataset <- data$RNASeqDataset
#' # perform geneID mapping
#' if (!require("org.Hs.eg.db", quietly = TRUE)) {
#'   BiocManager::install("org.Hs.eg.db")
#' }
#' library(org.Hs.eg.db)
#' ENSEMBLMapping <- AnnotationDbi::select(org.Hs.eg.db, 
#'                                         keys = rownames(RNASeqDataset), 
#'                                         columns = c("SYMBOL", "ENTREZID"), 
#'                                         keytype = "SYMBOL")
#' colnames(ENSEMBLMapping) <- c("FROM", "TO")
#' # Create the analysis design
#' RNASeqDesign <- model.matrix(~0 + condition, data = colData(RNASeqDataset))
#' RNASeqContrast <- limma::makeContrasts(conditionalzheimer-conditionnormal, 
#'                                        levels=RNASeqDesign)
#' # Perform DE analysis for RNA-Seq dataset
#' RNASeqDEExperiment <- RCPA::runDEAnalysis(RNASeqDataset, 
#'                                           method = "DESeq2", 
#'                                           design = RNASeqDesign, 
#'                                           contrast = RNASeqContrast, 
#'                                           annotation = ENSEMBLMapping)
#' # Get genesets                                           
#' genesets <- getGeneSets("KEGG", "hsa", minSize = 10, maxSize = 1000)
#' 
#' # Get enrichment analysis for three datasets
#' DF1 <- runGeneSetEnrichmentAnalysis(affyDEExperiment, genesets, method = "fgsea")
#' rownames(DF1) <- DF1$ID
#' DF2 <- runGeneSetEnrichmentAnalysis(agilDEExperiment, genesets, method = "fgsea")
#' rownames(DF2) <- DF2$ID
#' DF3 <- runGeneSetEnrichmentAnalysis(RNASeqDEExperiment, genesets, method = "fgsea")
#' rownames(DF3) <- DF3$ID
#' # Get common pathways
#' common_rownames <- intersect(intersect(rownames(DF1), rownames(DF2)), rownames(DF3))
#' 
#' DF1 <- DF1[common_rownames,]
#' DF2 <- DF2[common_rownames,]
#' DF3 <- DF3[common_rownames,]
#' # Merge into a list
#' allDFs <- list(DF1, DF2, DF3)
#' allData <- allDFs %>% bind_rows()
#' 
#' # Combine result using Enrichment Score
#' metaRes <- combinePathwayAnalysisResults(allDFs, method = "score")
#' }
#' @details This function performs mata analysis on multiple pathway analysis results.
#' @importFrom dplyr %>% bind_rows
#' @export
combinePathwayAnalysisResults <- function(PAResults, method = c("stouffer", "fisher", "addCLT", "geoMean", "minP", "REML")) {

    method <- match.arg(method)

    if (length(PAResults) == 1) {
        stop("Meta analysis is valid for two or more studies.")
    }

    for (PAResult in PAResults) {
        if (is.null(PAResult)) {
            stop("There is null object in the input list.")
        }

        if (!all(c("ID", "p.value", "normalizedScore", "sampleSize") %in% colnames(PAResult))) {
            stop("All the dataframes in the input list must have ID, p.value, normalizedScore and sampleSize columns.")
        }
    }

    if (method != "REML"){
        combinePFunc <- switch(method,
                           fisher = .runFisher,
                           stouffer = .runStouffer,
                           minP = min,
                           addCLT = .runAddCLT,
                           geoMean = .runGeoMean)

        pvalRes <- PAResults %>%
            lapply(function(df) {
                df[, c("ID", "p.value", "normalizedScore")] %>% as.data.frame()
            }) %>%
            bind_rows() %>%
            drop_na() %>%
            mutate(
                left.p = ifelse(.$normalizedScore < 0, .$p.value, 1 - .$p.value),
                right.p = ifelse(.$normalizedScore > 0, .$p.value, 1 - .$p.value)
            ) %>%
            group_by(ID) %>%
            summarise(
                left.p = combinePFunc(.data$left.p),
                right.p = combinePFunc(.data$right.p),
                n = length(.data$ID)
            ) %>%
            filter(n == length(PAResults))
    }

    metagenRes <- PAResults %>%
        lapply(function(df) {
            df[, c("ID", "normalizedScore", "p.value", "sampleSize")] %>% as.data.frame()
        }) %>%
        bind_rows() %>%
        group_by(ID) %>%
        group_split() %>%
        lapply(function(dat) {
            if (nrow(dat) < length(PAResults)) {
                return(NULL)
            }

            res <- try({
                metagen(data = dat,
                        studlab = ID,
                        TE = normalizedScore,
                        # seTE = logFCSE,
                        pval = p.value,
                        sm = "SMD",
                        method.tau = "REML",
                        hakn = TRUE,
                        n.e = sampleSize
                ) }, silent = TRUE)
            if (inherits(res, "try-error")) {
                res <- NULL
            }

            return(res)
        }) %>%
        do.call(what = rbind)


    metaResult <- metagenRes[, c("studlab", "TE.fixed", "seTE.fixed", "pval.fixed")] %>%
        as.data.frame() %>%
        mutate(
            ID = .data$studlab %>% sapply(`[`, 1),
            p.value = unlist(.data$pval.fixed),
            score = unlist(.data$TE.fixed),
            normalizedScore = .data$score
        ) %>%
        mutate(
            pFDR = p.adjust(.data$p.value, method = "fdr")
        )

    if (method != "REML"){
        metaResult <- metaResult %>%
            select("ID", "score", "normalizedScore") %>%
            inner_join(pvalRes, by = "ID") %>%
            mutate(
                p.value = ifelse(.$normalizedScore < 0, .$left.p, .$right.p),
                pFDR = p.adjust(.data$p.value, method = "fdr")
            )
    }

    metaResult$name <- PAResults[[1]][match(metaAnalysisResult$ID, PAResults[[1]]$ID), "name"]
    metaResult$pathwaySize <- PAResults[[1]][match(metaAnalysisResult$ID, PAResults[[1]]$ID), "pathwaySize"]
    metaResult[, c("ID", "name", "p.value", "pFDR", "score", "normalizedScore", "pathwaySize")]
}


#' @title Combine DE analysis results
#' @description This function performs mata analysis on multiple DE analysis results.
#' @param DEResults A list of dataframes containing DE analysis results.
#' Each dataframe must have ID, p.value, logFC and logFCSE columns.
#' @param method The method to combine p-values. It can be one of "fisher", "stouffer", "geoMean", "addCLT", "minP".
#' @return A dataframe containing combined DE analysis results.
#' The dataframe has ID, p.value, pDFR, logFC, and logFCSE columns.
#' @importFrom dplyr %>% bind_rows group_by summarize mutate
#' @importFrom meta metagen
#' @importFrom stats p.adjust
#' @importFrom tidyr drop_na
#' @export
combineDEAnalysisResults <- function(DEResults, method = c("stouffer", "fisher", "addCLT", "geoMean", "minP", "REML")) {

    method <- match.arg(method)

    if (length(DEResults) == 1) {
        stop("Meta analysis is valid for two or more studies.")
    }

    for (DEResult in DEResults) {
        if (is.null(DEResult)) {
            stop("There is null object in the input list.")
        }

        if (!all(c("ID", "p.value", "logFC", "logFCSE", "sampleSize") %in% colnames(DEResult))) {
            stop("All the dataframes in the input list must have p.value, logFC, logFCSE, and sampleSize columns.")
        }
    }

    if (method != "REML"){
        combinePFunc <- switch(method,
                           meanZ = .runMeanZ,
                           fisher = .runFisher,
                           stouffer = .runStouffer,
                           minP = min,
                           addCLT = .runAddCLT,
                           geoMean = .runGeoMean
    )

    pvalRes <- DEResults %>%
        lapply(function(df) {
            df[, c("ID", "p.value", "logFC")] %>% as.data.frame()
        }) %>%
        bind_rows() %>%
        drop_na() %>%
        mutate(
            left.p = ifelse(.$logFC < 0, .$p.value, 1 - .$p.value),
            right.p = ifelse(.$logFC > 0, .$p.value, 1 - .$p.value)
        ) %>%
        group_by(ID) %>%
        summarise(
            left.p = combinePFunc(.data$left.p),
            right.p = combinePFunc(.data$right.p),
            n = length(.data$ID)
        ) %>%
        filter(n == length(DEResults))
    }

    metagenRes <- DEResults %>%
        lapply(function(df) {
            df[, c("ID", "logFC", "logFCSE", "p.value", "sampleSize")] %>% as.data.frame()
        }) %>%
        bind_rows() %>%
        group_by(ID) %>%
        group_split() %>%
        lapply(function(dat) {
            if (nrow(dat) < length(DEResults)) {
                return(NULL)
            }

            res <- try({
                metagen(data = dat,
                        studlab = ID,
                        TE = logFC,
                        # seTE = logFCSE,
                        pval = p.value,
                        sm = "SMD",
                        method.tau = "REML",
                        hakn = TRUE,
                        n.e = sampleSize
                ) }, silent = TRUE)
            if (inherits(res, "try-error")) {
                res <- NULL
            }

            return(res)
        }) %>%
        do.call(what = rbind)


    metaResult <- metagenRes[, c("studlab", "TE.fixed", "seTE.fixed", "pval.fixed")] %>%
        as.data.frame() %>%
        mutate(
            ID = .data$studlab %>% sapply(`[`, 1),
            p.value = unlist(.data$pval.fixed),
            logFC = unlist(.data$TE.fixed),
            logFCSE = unlist(.data$seTE.fixed)
        ) %>%
        mutate(
            pFDR = p.adjust(.data$p.value, method = "fdr")
        ) %>%
        select("ID", "p.value", "pFDR", "logFC", "logFCSE") %>%
        drop_na() %>%
        as.data.frame()

    if (method != "REML"){
        metaResult <- metaResult %>%
            select("ID", "logFC", "logFCSE") %>%
            inner_join(pvalRes, by = "ID") %>%
            mutate(
                p.value = ifelse(.$logFC < 0, .$left.p, .$right.p),
                pFDR = p.adjust(.data$p.value, method = "fdr")
            ) %>%
            select("ID", "p.value", "pFDR", "logFC", "logFCSE")
    }

    return(metaResult)
}