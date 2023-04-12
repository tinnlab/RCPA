#' Differential expression analysis
#' @description Functions to perform differential expression analysis
#' These functions are used internally by runDEAnalysis.
#' @param exprs Normalized expression matrix. Rows are genes and columns are samples.
#' @param design A design model output by model.matrix.
#' @param contrast A contrast matrix output by limma::makeContrasts.
#' @return A data frame with DE analysis results.
#' Must contain the following columns: PROBEID, p.value, logFC, statistic.
#' @importFrom limma lmFit contrasts.fit eBayes topTable makeContrasts
#' @importFrom stats model.matrix
#' @importFrom dplyr %>% mutate rename
#' @name runDEInternal
.runLimma <- function(exprs, design, contrast) {
    exprs %>%
        lmFit(design) %>%
        contrasts.fit(contrast) %>%
        eBayes() %>%
        topTable(coef = 1, number = nrow(exprs)) %>%
        mutate(PROBEID = rownames(.), p.value = .$P.Value, statistic = .$t)
}

#' @rdname runDEInternal
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
.runDESeq2 <- function(exprs, design, contrast) {

    DERes <- DESeqDataSetFromMatrix(
        countData = exprs,
        colData = data.frame(row.names = colnames(exprs), dummy.var = rep(1, ncol(exprs))),
        design = design
    ) %>%
        DESeq() %>%
        results(contrast) %>%
        as.data.frame() %>%
        mutate(PROBEID = rownames(.), statistic = .$stat, p.value = .$pvalue, logFC = .$log2FoldChange)

    DERes$p.value[is.na(DERes$p.value)] <- 1
    DERes
}

#' @rdname runDEInternal
#' @importFrom edgeR DGEList calcNormFactors estimateDisp exactTest topTags
.runEdgeR <- function(exprs, design, contrast) {

    DERes <- DGEList(counts = exprs) %>%
        calcNormFactors() %>%
        estimateDisp(design = design) %>%
        glmQLFit(design = design, contrast = contrast) %>%
        glmQLFTest(contrast = contrast) %>%
        topTags(n = nrow(exprs)) %>%
        `$`("table") %>%
        mutate(PROBEID = rownames(.), p.value = .$PValue, statistic = .$logFC, logFC = .$logFC)

    DERes$p.value[is.na(DERes$p.value)] <- 1
    DERes
}

#' @title Get ID mapping annotation from GEO platform
#' @description This function gets ID mapping annotation from GEO platform.
#' This function is used internally by runDEAnalysis.
#' @param GEO platform ID. E.g., GPL570
#' @param dbObj A database object from AnnotationDbi package.
#' @return A data frame with ID mapping annotation. The first column is the probe ID and the second column is the gene symbol.
.getIDMappingAnnotation <- function(platform) {

    annotations <- list(
        GPL96 = "hgu133a.db",
        GPL97 = "hgu133b.db",
        GPL570 = "hgu133plus2.db",
        GPL6244 = "hugene10sttranscriptcluster.db",
        GPL4866 = "hgu133plus2.db",
        GPL16311 = "hgu133plus2.db",
        GPL571 = "hgu133a2.db",
        GPL10739 = "hugene10sttranscriptcluster.db",
        GPL201 = "hgfocus.db",
        GPL8300 = "hgu95av2.db",
        GPL80 = "hu6800.db",
        GPL1261 = "mouse4302.db",
        GPL6246 = "mogene10sttranscriptcluster.db",
        GPL81 = "mgu74a.db",
        GPL339 = "moe430a.db",
        GPL10740 = "mogene10sttranscriptcluster.db",
        GPL16570 = "mogene20sttranscriptcluster.db",
        GPL11044 = "mouse4302.db",
        GPL6193 = "moex10sttranscriptcluster.db",
        GPL13667 = "hgu219.db",
        GPL23126 = "clariomdhumantranscriptcluster.db",
        GPL6102 = "illuminaHumanv2.db",
        GPL6947 = "illuminaHumanv3.db"
    )

    if (grep("GPL", platform)) {
        if (!is.null(annotations[[annotation]])) {
            anno <- annotations[[annotation]]
            .requirePackage(anno)
            annotation <- get(anno)
        } else {
            stop(paste0("Platform ", platform, " is not supported. Please pass an AnnotationDbi object instead"))
        }
    }

    annotation
}

#' @title Map probe IDs to gene symbols
#' @description This function maps probe IDs to gene symbols.
#' This function is used internally by runDEAnalysis.
#' @param exprs Expression matrix. Rows are genes and columns are samples.
#' @param annotation Annotation data frame with mapping between probe IDs and gene symbols. The first column is the probe ID and the second column is the gene symbol.
#' @param DEResults DE analysis results data frame. Must have a column named ID and p.value.
#' @param mappingMethod Method to use for mapping probe IDs to gene symbols. Can be "mean", "max", "median", "mostSignificant".
#' @return A list with two elements:
#' \itemize{
#'  \item exprs: Expression matrix with probe IDs mapped to gene symbols. Rows are genes and columns are samples.
#'  \item DEResults: DE analysis results data frame with probe IDs mapped to gene symbols.
#' }
#' @importFrom dplyr %>% arrange select group_by mutate first rename left_join summarize_all
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr drop_na
.mapProbeIDsToGeneSymbols <- function(exprs, annotation, DEResults) {

    if (is.null(annotation)) {
        return(list(exprs = exprs, DEResults = DEResults))
    }

    top <- DEResults %>% arrange(.data$p.value)

    mapping <- top %>%
        select(PROBEID) %>%
        left_join(annotation, by = "PROBEID") %>%
        group_by(PROBEID) %>%
        mutate(SYMBOL.FIRST = first(SYMBOL)) %>%
        select(PROBEID, SYMBOL.FIRST) %>%
        unique() %>%
        group_by(SYMBOL.FIRST) %>%
        mutate(PROBEID.FIRST = first(PROBEID)) %>%
        select(PROBEID.FIRST, SYMBOL.FIRST) %>%
        unique() %>%
        rename(PROBEID = PROBEID.FIRST, SYMBOL = SYMBOL.FIRST)

    mappedExprs <- exprs[mapping$PROBEID,] %>%
        `rownames<-`(mapping$SYMBOL)

    mappedDEResults <- DEResults %>%
        left_join(mapping, by = "PROBEID") %>%
        drop_na()

    list(
        exprs = mappedExprs,
        DEResults = mappedDEResults
    )
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


