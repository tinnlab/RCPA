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
