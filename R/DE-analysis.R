#' @title Differential expression analysis
#' @description Functions to perform differential expression analysis
#' These functions are used internally by runDEAnalysis.
#' @param exprs Normalized expression matrix. Rows are genes and columns are samples.
#' @param design A design model output by model.matrix.
#' @param contrast A contrast matrix output by limma::makeContrasts.
#' @return A data frame with DE analysis results.
#' Must contain the following columns: ID, p.value, logFC, statistic, dispersion.
#' @importFrom limma lmFit contrasts.fit eBayes topTable makeContrasts
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmQLFit glmQLFTest topTags
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
        mutate(ID = rownames(.),
               p.value = .$P.Value,
               statistic = .$t
        )

    resTable$dispersion <- (DERes$stdev.unscaled * sqrt(DERes$s2.post))[rownames(resTable), 1]
    resTable[rownames(exprs), c("ID", "p.value", "statistic", "logFC", "dispersion")]
}

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
        mutate(ID = rownames(.), statistic = .$stat, p.value = .$pvalue, logFC = .$log2FoldChange)

    resTable$p.value[is.na(resTable$p.value)] <- 1
    resTable$dispersion <- resTable$lfcSE
    resTable[rownames(exprs), c("ID", "p.value", "statistic", "logFC", "dispersion")]
}

#' @rdname runDEInternal
.runEdgeR <- function(exprs, design, contrast) {

    DERes <- DGEList(counts = exprs) %>%
        calcNormFactors() %>%
        estimateDisp(design = design) %>%
        glmQLFit(design = design, contrast = contrast) %>%
        glmQLFTest(contrast = contrast)

    resTable <- DERes$table %>%
        mutate(ID = rownames(.), p.value = .$PValue, statistic = .$logFC, logFC = .$logFC)

    resTable$p.value[is.na(resTable$p.value)] <- 1
    resTable$dispersion <- DERes$dispersion
    resTable[rownames(exprs), c("ID", "p.value", "statistic", "logFC", "dispersion")]
}

#' @title Differential expression analysis
#' @description This function performs differential expression analysis using either limma, DESeq2 or edgeR.
#' @param summarizedExperiment SummarizedExperiment object
#' @param method Method to use for differential expression analysis. Can be "limma", "DESeq2" or "edgeR".
#' @param design A design model output by model.matrix.
#' @param contrast A contrast matrix. See limma::makeContrasts.
#' @param annotation A data frame mapping between probe IDs and entrez gene IDs.
#' If not provided, the function will try to get the mapping from the platform annotation in the SummarizedExperiment object.
#' If the annotation is not available, the function will return the probe IDs.
#' Regardless of the type of annotation, it must contains two columns: FROM and TO,
#' where FROM is the probe ID and TO is the entrez gene ID.
#' @return A SummarizedExperiment object with DE analysis results appended to the rowData slot with the following columns:
#' \itemize{
#' \item{ID}{gene ID. If annotation is provided, this will be the entrez gene ID. Otherwise, it will be the probe ID.}
#' \item{logFC}{log2 fold change}
#' \item{p.value}{p-value from the DE analysis using the specified method}
#' \item{pFDR}{p-value adjusted for multiple testing using Benjamini-Hochberg method}
#' \item{statistic}{statistic from the DE analysis using the specified method.
#' \item{dispersion}{dispersion from the DE analysis using the specified method}
#' For limma, this is the t-statistic.
#' For DESeq2, this is the Wald statistic.
#' For edgeR, this is the log fold change.}
#' \item{dispersion}{dispersion from the DE analysis using the specified method}
#' }
#' The assay slot will contain the input expression/count matrix,
#' and the rownames will be mapped to the gene IDs if annotation is found in the input SummarizedExperiment object
#' or in the annotation parameter.
#' Other slots will be the same as in the input SummarizedExperiment object.
#' @examples
#' #TODO add examples
#' @importFrom SummarizedExperiment SummarizedExperiment rowData assay colData
#' @importFrom S4Vectors SimpleList
#' @importFrom dplyr %>%
#' @seealso \code{\link{limma::makeContrasts}}
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
        if (is.null(metadata(summarizedExperiment)$platform)) {
            stop("Platform annotation is not found in the meta data of the SummarizedExperiment object. Please provide an annotation data frame instead.")
        }
        annotation <- .getIDMappingAnnotation(platform = metadata(summarizedExperiment)$platform)
    }

    if ("character" %in% class(annotation)) {
        annotation <- .getIDMappingAnnotation(platform = annotation)
    }

    if ("data.frame" %in% class(annotation)) {
        if (!all(c("FROM", "TO") %in% colnames(annotation))) {
            stop("Annotation data frame must have columns FROM and TO")
        }
    } else {
        stop("Annotation must be a data frame or a string")
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
    mappedResults <- .mapIDs(exprs, annotation, DEResult)
    mappedResults$DEResult$pFDR <- p.adjust(mappedResults$DEResult$p.value, method = "BH")

    # create a new SummarizedExperiment object
    newSummarizedExperiment <- SummarizedExperiment(
        assays = SimpleList(counts = mappedResults$exprs),
        rowData = data.frame(
            cbind(
                rowData(summarizedExperiment)[mappedResults$mapping$ID,],
                PROBEID = mappedResults$mapping$ID,
                mappedResults$DEResult
            ),
            row.names = mappedResults$mapping$TO
        ),
        colData = colData(summarizedExperiment),
        metadata = c(
            metadata(summarizedExperiment),
            DEAnalysis = list(
                method = method,
                design = design,
                contrast = contrast,
                mapping = mappedResults$mapping
            )
        )
    )

    newSummarizedExperiment
}


