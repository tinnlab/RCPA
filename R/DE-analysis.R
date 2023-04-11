#' @title Differential expression analysis using limma
#' @description This function performs differential expression analysis using limma.
#' This function is used internally by runDEAnalysis.
#' @param exprs Normalized expression matrix. Rows are genes and columns are samples.
#' @param controlSamples Vector of control samples.
#' @param conditionSamples Vector of condition samples.
#' @return A data frame with DE analysis results
#' @importFrom limma lmFit contrasts.fit eBayes topTable makeContrasts
#' @importFrom stats model.matrix
#' @importFrom dplyr %>% mutate rename
.runLimma <- function(exprs, controlSamples, conditionSamples, isPaired) {
    group <- c(rep("c", length(controlSamples)), rep("d", length(conditionSamples)))
    block <- c(rep(1, length(controlSamples)), rep(2, length(conditionSamples)))

    if (isPaired) {
        block <- factor(block)
        design <- model.matrix(~0 + group + block)
        colnames(design) <- substr(colnames(design), 2, 100)
    }  else {
        design <- model.matrix(~0 + group)
        colnames(design) <- levels(group)
    }

    exprs[, c(controlSamples, conditionSamples)] %>%
        lmFit(design) %>%
        contrasts.fit(makeContrasts(contrasts = "d-c", levels = design)) %>%
        eBayes() %>%
        topTable(coef = 1, number = nrow(exprs)) %>%
        mutate(PROBEID = rownames(.)) %>%
        rename(p.value = P.Value, statistic = t)
}

#' @title Differential expression analysis using t-test
#' @description This function performs differential expression analysis using t-test.
#' This function is used internally by runDEAnalysis.
#' @param exprs Normalized expression matrix. Rows are genes and columns are samples.
#' @param controlSamples Vector of control samples.
#' @param conditionSamples Vector of condition samples.
#' @return A data frame with DE analysis results
#' @importFrom matrixTests row_t_equalvar
.runTTest <- function(exprs, controlSamples, conditionSamples, isPaired) {

    X <- exprs[, controlSamples]
    Y <- exprs[, conditionSamples]

    DEfnc <- if (isPaired) matrixTests::row_t_paired else matrixTests::row_t_equalvar
    DEfnc(X, Y, alternative = "two.sided") %>%
        mutate(PROBEID = rownames(.)) %>%
        rename(p.value = pvalue)
}

#' @title Differential expression analysis using DESeq2
#' @description This function performs differential expression analysis using limma.
#' This function is used internally by runDEAnalysis.
#' @param exprs Raw counts matrix. Rows are genes and columns are samples.
#' @param controlSamples Vector of control samples.
#' @param conditionSamples Vector of condition samples.
#' @return A data frame with DE analysis results
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom dplyr %>% mutate rename
.runDESeq2 <- function(exprs, controlSamples, conditionSamples) {

    group <- c(rep("c", length(controlSamples)), rep("d", length(conditionSamples)))
    exprs <- exprs[, c(controlSamples, conditionSamples)]

    colData <- data.frame(
        row.names = colnames(exprs),
        condition = group
    )

    DESeqDataSetFromMatrix(
        countData = exprs,
        colData = colData,
        design = ~condition
    ) %>%
        DESeq() %>%
        results(contrast = c("condition", "d", "c")) %>%
        as.data.frame() %>%
        mutate(PROBEID = rownames(.)) %>%
        rename(p.value = pvalue)
}

#' @title Differential expression analysis using DESeq2
#' @description This function performs differential expression analysis using limma.
#' This function is used internally by runDEAnalysis.
#' @param exprs Raw counts matrix. Rows are genes and columns are samples.
#' @param controlSamples Vector of control samples.
#' @param conditionSamples Vector of condition samples.
#' @return A data frame with DE analysis results
#' @importFrom edgeR DGEList calcNormFactors estimateDisp exactTest topTags
#' @importFrom dplyr %>% mutate rename
.runEdgeR <- function(exprs, controlSamples, conditionSamples) {

    group <- c(rep("c", length(controlSamples)), rep("d", length(conditionSamples)))
    exprs <- exprs[, c(controlSamples, conditionSamples)]

    DGEList(counts = exprs, group = factor(group)) %>%
        calcNormFactors() %>%
        estimateDisp(design = model.matrix(~group), robust = TRUE) %>%
        exactTest() %>%
        topTags(n = nrow(exprs)) %>%
        `$`("table") %>%
        mutate(PROBEID = rownames(.)) %>%
        rename(p.value = pvalue)
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

    top <- DEResults %>% arrange(p.value)

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
runDEAnalysis <- function(summarizedExperiment, controlSamples, conditionSamples, isPaired = FALSE, method = c("limma", "DESeq2", "edgeR"), annotation = NULL) {
    method <- match.arg(method)
    mappingMethod <- match.arg(mappingMethod)

    # check if all samples are in the SummarizedExperiment object
    if (any(!controlSamples %in% colnames(summarizedExperiment))) {
        stop("Some control samples are not in the SummarizedExperiment object")
    }

    if (any(!conditionSamples %in% colnames(summarizedExperiment))) {
        stop("Some condition samples are not in the SummarizedExperiment object")
    }

    # check if there are any samples in both control and condition vectors
    if (any(intersect(controlSamples, conditionSamples))) {
        stop("Some samples are in both control and condition vectors")
    }

    if (isPaired) {
        if (length(controlSamples) != length(conditionSamples)) {
            stop("Number of control and condition samples must be the same for paired analysis")
        }
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

    # filter out samples that are not in the control or condition vectors

    summarizedExperiment <- summarizedExperiment[, intersect(c(controlSamples, conditionSamples), colnames(summarizedExperiment))]

    # get expression matrix
    exprs <- assay(summarizedExperiment)

    DEFunc <- switch(method,
                     limma = .runLimma,
                     DESeq2 = .runDESeq2,
                     edgeR = .runEdgeR
    )

    # run DE analysis
    DEResult <- DEFunc(exprs, controlSamples, conditionSamples, isPaired)

    # map probe IDs to gene symbols
    mappedResults <- .mapProbeIDsToGeneSymbols(exprs, annotation, DEResult)

    # create a new SummarizedExperiment object
    newSummarizedExperiment <- SummarizedExperiment(
        assays = SimpleList(counts = mappedResults$exprs),
        rowData = cbind(rowData(summarizedExperiment), mappedResults$DEResult),
        colData = colData(summarizedExperiment)
    )

    newSummarizedExperiment
}


