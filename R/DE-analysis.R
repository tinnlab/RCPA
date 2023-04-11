#' @title Differential expression analysis using limma
#' @description This function performs differential expression analysis using limma.
#' This function is used internally by runDEAnalysis.
#' @param exprs Normalized expression matrix. Rows are genes and columns are samples.
#' @param controlSamples Vector of control samples.
#' @param conditionSamples Vector of condition samples.
#' @return A data frame with DE analysis results
#' @details This function is used internally by runDEAnalysis.
.runLimma <- function(exprs, controlSamples, conditionSamples) {

}

#' @title Differential expression analysis using DESeq2
#' @description This function performs differential expression analysis using limma.
#' This function is used internally by runDEAnalysis.
#' @param exprs Raw counts matrix. Rows are genes and columns are samples.
#' @param controlSamples Vector of control samples.
#' @param conditionSamples Vector of condition samples.
#' @return A data frame with DE analysis results
#' @details This function is used internally by runDEAnalysis.
.runDESeq2 <- function(exprs, controlSamples, conditionSamples) {


}

#' @title Differential expression analysis using DESeq2
#' @description This function performs differential expression analysis using limma.
#' This function is used internally by runDEAnalysis.
#' @param exprs Raw counts matrix. Rows are genes and columns are samples.
#' @param controlSamples Vector of control samples.
#' @param conditionSamples Vector of condition samples.
#' @return A data frame with DE analysis results
#' @details This function is used internally by runDEAnalysis.
.runEdgeR <- function(exprs, controlSamples, conditionSamples) {


}

#' @title Get ID mapping annotation from GEO platform
#' @description This function gets ID mapping annotation from GEO platform.
#' This function is used internally by runDEAnalysis.
#' @param GEO platform ID. E.g., GPL570
#' @return A data frame with ID mapping annotation. The first column is the probe ID and the second column is the gene symbol.
.getIDMappingAnnotation <- function(platform) {

}

#' @title Map probe IDs to gene symbols
#' @description This function maps probe IDs to gene symbols.
#' This function is used internally by runDEAnalysis.
#' @param exprs Expression matrix. Rows are genes and columns are samples.
#' @param annotation Annotation data frame with mapping between probe IDs and gene symbols. The first column is the probe ID and the second column is the gene symbol.
#' @param DEResults DE analysis results data frame. Must have a column named ID and p.value.
#' @return A list with two elements:
#' \itemize{
#'  \item exprs: Expression matrix with probe IDs mapped to gene symbols. Rows are genes and columns are samples.
#'  \item DEResults: DE analysis results data frame with probe IDs mapped to gene symbols.
#' }
.mapProbeIDsToGeneSymbols <- function(exprs, annotation, DEResults) {

}

#' @title Differential expression analysis
#' @description This function performs differential expression analysis using either limma, DESeq2 or edgeR.
#' @param summarizedExperiment SummarizedExperiment object
#' @param controlSamples Vector of control samples. Each sample is a column in the SummarizedExperiment object.
#' @param conditionSamples Vector of condition samples. Each sample is a column in the SummarizedExperiment object.
#' @param method Method to use for differential expression analysis. Can be "limma", "DESeq2" or "edgeR".
#' @param IDMapping A data frame with mapping between probe IDs and gene symbols.
#' If not provided, the function will try to get the mapping from the platform annotation in the SummarizedExperiment object.
#' If the mapping is not available, the function will return the probe IDs.
#' The data frame should have two columns: PROBEID and SYMBOL.
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
#' @export
runDEAnalysis <- function(summarizedExperiment, controlSamples, conditionSamples, method = c("limma", "DESeq2", "edgeR"), annotation = NULL) {
    method <- match.arg(method)

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

    # get ID mapping annotation
    if (is.null(annotation)) {
        annotation <- .getIDMappingAnnotation(summarizedExperiment$platform)
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
    DEResult <- DEFunc(exprs, controlSamples, conditionSamples)

    # map probe IDs to gene symbols
    mappedResults <- .mapProbeIDsToGeneSymbols(exprs, annotation, DEResult)

    # create a new SummarizedExperiment object
    newSummarizedExperiment <- SummarizedExperiment(
        assays = SimpleList(counts = mappedResults$exprs),
        rowData = rowData(summarizedExperiment),
        colData = colData(summarizedExperiment)
    )

    # add DE results to the rowData slot
    rowData(newSummarizedExperiment) <- cbind(rowData(newSummarizedExperiment), mappedResults$DEResult)
}


