#' Differential expression analysis
#' @description Functions to perform differential expression analysis
#' These functions are used internally by runDEAnalysis.
#' @param exprs Normalized expression matrix. Rows are genes and columns are samples.
#' @param design A design model output by model.matrix.
#' @param contrast A contrast matrix output by limma::makeContrasts.
#' @return A data frame with DE analysis results.
#' Must contain the following columns: PROBEID, p.value, logFC, statistic, dispersion.
#' @importFrom limma lmFit contrasts.fit eBayes topTable makeContrasts
#' @importFrom stats model.matrix
#' @importFrom dplyr %>% mutate
#' @name runDEInternal
.runLimma <- function(exprs, design, contrast) {
    DERes <- exprs %>%
        lmFit(design) %>%
        contrasts.fit(contrast) %>%
        eBayes()

    resTable <- DERes %>%
        topTable(coef = 1, number = nrow(exprs)) %>%
        mutate(PROBEID = rownames(.),
               p.value = .$P.Value,
               statistic = .$t
        )

    resTable$dispersion <- (DERes$stdev.unscaled * sqrt(DERes$s2.post))[rownames(resTable),1]
    resTable[rownames(exprs), c("PROBEID", "p.value", "statistic", "logFC", "dispersion")]
}

#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom dplyr %>% mutate
#' @rdname runDEInternal
.runDESeq2 <- function(exprs, design, contrast) {

    DERes <- DESeqDataSetFromMatrix(
        countData = exprs,
        colData = data.frame(row.names = colnames(exprs), dummy.var = rep(1, ncol(exprs))),
        design = design
    ) %>%
        DESeq() %>%
        results(contrast)

    resTable <- DERes %>%
        as.data.frame() %>%
        mutate(PROBEID = rownames(.), statistic = .$stat, p.value = .$pvalue, logFC = .$log2FoldChange)

    resTable$p.value[is.na(resTable$p.value)] <- 1
    resTable$dispersion <- resTable$lfcSE
    resTable[rownames(exprs), c("PROBEID", "p.value", "statistic", "logFC", "dispersion")]
}

#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmQLFit glmQLFTest topTags
#' @importFrom dplyr %>% mutate
#' @rdname runDEInternal
.runEdgeR <- function(exprs, design, contrast) {

    DERes <- DGEList(counts = exprs) %>%
        calcNormFactors() %>%
        estimateDisp(design = design) %>%
        glmQLFit(design = design, contrast = contrast) %>%
        glmQLFTest(contrast = contrast)

    resTable <- DERes$table %>%
        mutate(PROBEID = rownames(.), p.value = .$PValue, statistic = .$logFC, logFC = .$logFC)

    resTable$p.value[is.na(resTable$p.value)] <- 1
    resTable$dispersion <- DERes$dispersion
    resTable[rownames(exprs), c("PROBEID", "p.value", "statistic", "logFC", "dispersion")]
}

#' @title Differential expression analysis
#' @description This function performs differential expression analysis using either limma, DESeq2 or edgeR.
#' @param summarizedExperiment SummarizedExperiment object
#' @param controlSamples Vector of control samples. Each sample is a column in the SummarizedExperiment object.
#' @param conditionSamples Vector of condition samples. Each sample is a column in the SummarizedExperiment object.
#' @param method Method to use for differential expression analysis. Can be "limma", "DESeq2" or "edgeR".
#' @param annotation A data frame, ChipDb, or OrgDb object with mapping between probe IDs and gene symbols.
#' If not provided, the function will try to get the mapping from the platform annotation in the SummarizedExperiment object.
#' If the annotation is not available, the function will return the probe IDs.
#' Regardless of the type of annotation, it must contains two columns: PROBEID and SYMBOL.
#' @return A SummarizedExperiment object with DE analysis results appended to the rowData slot with the following columns:
#' \itemize{
#' \item{logFC}{log2 fold change}
#' \item{p.value}{p-value}
#' \item{pFDR}{p-value adjusted for multiple testing using Benjamini-Hochberg method}
#' \item{control.mean}{mean expression in control samples}
#' \item{condition.mean}{mean expression in condition samples}
#' }
#' The assay slot will contain only samples from the control and condition vectors,
#' and the rownames will be the gene symbols.
#' Other slots will be the same as in the input SummarizedExperiment object.
#' @examples
#'
#' @importFrom SummarizedExperiment SummarizedExperiment rowData assay
#' @importFrom S4Vectors SimpleList
#' @importFrom dplyr %>%
#' @export
runDEAnalysis <- function(summarizedExperiment, method = c("limma", "DESeq2", "edgeR"), design, contrast, annotation = NULL) {
    method <- match.arg(method)

    if (is.null(design)) {
        stop("Design matrix must be provided")
    }

    if (is.null(contrast)) {
        stop("Contrast matrix must be provided")
    }

    # get ID mapping annotation
    if (is.null(annotation)) {
        annotation <- .getIDMappingAnnotation(platform = summarizedExperiment$platform)
    }

    if ("data.frame" %in% class(annotation) & !all(c("PROBEID", "SYMBOL") %in% colnames(annotation))) {
        stop("Annotation data frame must have columns PROBEID and SYMBOL")
    } else if (any(class(annotation) %in% c("ChipDb", "OrgDb"))) {
        annotation <- AnnotationDbi::select(annotation, keys = rownames(summarizedExperiment), columns = c("PROBEID", "SYMBOL"), keytype = "PROBEID")
    } else {
        stop("Annotation must be a data frame, ChipDb or OrgDb object")
    }

    # get expression matrix
    exprs <- assay(summarizedExperiment)

    DEFunc <- switch(method,
                     limma = .runLimma,
                     DESeq2 = .runDESeq2,
                     edgeR = .runEdgeR
    )

    # run DE analysis
    DEResult <- DEFunc(exprs, design, contrast)

    # map probe IDs to gene symbols
    mappedResults <- .mapProbeIDsToGeneSymbols(exprs, annotation, DEResult)

    # create a new SummarizedExperiment object
    newSummarizedExperiment <- SummarizedExperiment(
        assays = SimpleList(counts = mappedResults$exprs),
        rowData = cbind(
            rowData(summarizedExperiment),
            mappedResults$DEResult[, c("logFC", "p.value", "statistic")]
        ),
        colData = colData(summarizedExperiment)
    )

    newSummarizedExperiment
}


