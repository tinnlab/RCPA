#' @title Differential expression analysis
#' @description Functions to perform differential expression analysis
#' These functions are used internally by runDEAnalysis.
#' @param exprs Normalized expression matrix. Rows are genes and columns are samples.
#' @param design A design model output by model.matrix.
#' @param contrast A contrast matrix output by makeContrasts from limma package.
#' @return A data frame with DE analysis results.
#' Must contain the following columns: ID, p.value, logFC, statistic, avgExpr
#' @examples
#' \dontrun{
#' #Loading necessary libraries
#' library(AnnotationDbi)
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
#' rowData(affyDEExperiment)
#' 
#' # Get Agilent dataset
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
#' rowData(agilDEExperiment)
#' 
#' # Get RNA-Seq dataset
#' RNASeqDataset <- data$RNASeqDataset
#' # Oerform geneID mapping
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
#' rowData(RNASeqDEExperiment)
#' }
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
               statistic = .$t,
               avgExpr = .$AveExpr
        )

    # resTable$dispersion <- (DERes$stdev.unscaled * sqrt(DERes$s2.post))[rownames(resTable), 1]
    resTable[rownames(exprs), c("ID", "p.value", "statistic", "logFC", "avgExpr")]
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
        mutate(ID = rownames(.), statistic = .$stat, p.value = .$pvalue, logFC = .$log2FoldChange, avgExpr = log2(.$baseMean + 1))

    resTable$p.value[is.na(resTable$p.value)] <- 1
    # resTable$dispersion <- resTable$lfcSE
    resTable[rownames(exprs), c("ID", "p.value", "statistic", "logFC", "avgExpr")]
}

#' @rdname runDEInternal
.runEdgeR <- function(exprs, design, contrast) {

    DERes <- DGEList(counts = exprs) %>%
        calcNormFactors() %>%
        estimateDisp(design = design) %>%
        glmQLFit(design = design, contrast = contrast) %>%
        glmQLFTest(contrast = contrast)

    resTable <- DERes$table %>%
        mutate(ID = rownames(.), p.value = .$PValue, statistic = .$logFC, logFC = .$logFC, avgExpr = .$logCPM)

    resTable$p.value[is.na(resTable$p.value)] <- 1
    # resTable$dispersion <- DERes$dispersion
    resTable[rownames(exprs), c("ID", "p.value", "statistic", "logFC", "avgExpr")]
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
#' For limma, this is the t-statistic.
#' For DESeq2, this is the Wald statistic.
#' For edgeR, this is the log fold change.}
#' \item{avgExpr}{
#' For limma, it is the average expression.
#' For DESeq2, it is the log base mean.
#' For edgeR, it is the log CPM.
#' }
#' }
#' The assay slot will contain the input expression/count matrix,
#' and the rownames will be mapped to the gene IDs if annotation is found in the input SummarizedExperiment object
#' or in the annotation parameter.
#' Other slots will be the same as in the input SummarizedExperiment object.
#' @examples
#' \dontrun{
#' #Loading necessary libraries
#' library(hgu133plus2.db)
#' library(AnnotationDbi)
#' library(SummarizedExperiment)
#' library(limma)
#' library(RCPA)
#' # generate a random gene expression matrix
#' set.seed(123)
#' exprs <- round(matrix(2^abs(rnorm(1000, sd = 4)), nrow = 100, ncol = 10))
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
#'     group = factor(c(rep("control", 
#'     length(controlSamples)), 
#'     rep("condition", 
#'     length(conditionSamples)))),
#'     pair = factor(c(seq_along(controlSamples), seq_along(conditionSamples)))
#' )
#' # Construct summarizedExperiment object
#' summarizedExperiment <- SummarizedExperiment(
#'     assays = list(counts = exprs),
#'     colData = colData
#' )
# # Construct design and contrast tables
#' # control vs condition
#' design <- model.matrix(~0 + group, data = colData)
#' contrast <- makeContrasts("groupcondition-groupcontrol", levels = design)
#' # Perform DE analysis
#' DERes <- runDEAnalysis(summarizedExperiment, 
#' method = "DESeq2", 
#' design, contrast, 
#' annotation = "GPL570")
#' 
#' # View DE analysis data frame
#' rowData(DERes)
#' }
#' @importFrom SummarizedExperiment SummarizedExperiment rowData assay colData
#' @importFrom S4Vectors SimpleList metadata
#' @importFrom dplyr %>%
#' @importFrom stats p.adjust
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

    if (!is.null(annotation)){
        exprs <- exprs[intersect(unique(annotation$FROM), rownames(exprs)),]
    }

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


