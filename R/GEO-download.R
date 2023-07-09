#' @title Download GEO object
#' @description This function downloads the GEO object using getGEO function.
#' This function is used internally by downloadGEO.
#' @param GEOID The GEO dataset under query, ex. GSE14762.
#' @param platform The platform of selected GEO dataset, ex. GPL4866.
#' @param destDir The user path to save the downloaded files.
#' @return A GEO object with class ExpressionSet.
#' @details This function is used internally by downloadGEO.
#' @importFrom Biobase annotation
#' @importFrom GEOquery getGEO
#' @noRd
.downloadGEOObject <- function(GEOID, platform, destDir) {

    if (!dir.exists(destDir)) {
        stop("The destination directory does not exist.")
    }

    gsets <- getGEO(GEOID, GSEMatrix = TRUE, getGPL = TRUE, destdir = destDir)
    platforms <- sapply(gsets, annotation)

    if (!platform %in% platforms) {
        stop("The platform is not available in the GEO dataset.")
    }

    gsets[[which(platforms == platform)]]
}

#' @title Download samples from GEO object
#' @description This function downloads the corresponding samples from queried GEO object. This function is only used if the protocol is affymetrix or agilent.
#' This function is used internally by downloadGEO.
#' @param sampleIDs A vector of GEO accession sample IDs to be downloaded.
#' @param protocol The protocol of the selected GEO dataset, including affymetrix and agilent.
#' @param destDir The user path to save the downloaded files.
#' @return TRUE
#' @details This function is used internally by downloadGEO.
#' @importFrom GEOquery getGEOSuppFiles
#' @importFrom dplyr %>%
#' @importFrom utils URLdecode
#' @importFrom httr HEAD
#' @noRd
.downloadSamples <- function(sampleIDs, protocol, destDir) {

    if (!dir.exists(destDir)) {
        stop("The destination directory does not exist.")
    }

    if(!(protocol %in% c("affymetrix", "agilent"))) {
        stop("The specified protocol is not valid.")
    }

    for (id in sampleIDs) {

        if(protocol == "affymetrix"){

            if (file.exists(file.path(destDir, paste0(id, ".CEL.gz")))) next()

            filesURL <- getGEOSuppFiles(id, baseDir = destDir, makeDirectory = FALSE, fetch_files = FALSE)

            actualInfo <- lapply(1:nrow(filesURL), function(curRow){
                actualFileLength <- httr::HEAD(filesURL[curRow, 2])$headers$`content-length`
                data.frame(
                    fname = filesURL[curRow, 1],
                    fsize = actualFileLength,
                    stringsAsFactors = FALSE
                )
            }) %>% do.call(what = rbind)

            downloadedInfo <- getGEOSuppFiles(id, baseDir = destDir, makeDirectory = FALSE)
            downloadedFiles <- downloadedInfo %>% rownames()

            isDeletedFlag <- FALSE
            for(i in 1:2){
                 lapply(nrow(downloadedInfo), function(curRow){
                    if(downloadedInfo[curRow, "size"] != actualInfo[curRow, "fsize"]){
                        file.remove(downloadedFiles[curRow])
                        isDeletedFlag <- TRUE
                    }
                })

                if(i == 2 & isDeletedFlag == TRUE){
                    stop(paste0("There is an error in downloading sample ", id))
                }else if(isDeletedFlag){
                    downloadedInfo <- getGEOSuppFiles(id, baseDir = destDir, makeDirectory = FALSE)
                    isDeletedFlag <- FALSE
                }else break
            }

            #downloadedFiles <- getGEOSuppFiles(id, baseDir = destDir, makeDirectory = FALSE) %>% rownames()

            if(is.null(downloadedFiles)){
                stop("Check the specified samples IDs to be valid. No file is found.")
            }

            downloadedFiles <- sapply(downloadedFiles, function(fileName){
                if (!file.exists(fileName)) {
                    URLdecode(fileName)
                } else {
                    fileName
                }
            }) %>% as.vector()

            downloadedFiles[grep(".cel.gz", downloadedFiles, ignore.case = TRUE)] %>% file.rename(paste0(destDir, "/", id, ".CEL.gz"))

        }else{

            if (file.exists(file.path(destDir, paste0(id, ".TXT.gz")))) next()

            downloadedFiles <- getGEOSuppFiles(id, baseDir = destDir, makeDirectory = FALSE) %>% rownames()

            if(is.null(downloadedFiles)){
                stop("Check the specified samples IDs to be valid. No file is found.")
            }

            downloadedFiles <- sapply(downloadedFiles, function(fileName){
                if (!file.exists(fileName)) {
                    URLdecode(fileName)
                } else {
                    fileName
                }
            }) %>% as.vector()

            downloadedFiles[grep(".txt.gz", downloadedFiles, ignore.case = TRUE)] %>% file.rename(paste0(destDir, "/", id, ".TXT.gz"))
        }
    }

    return(TRUE)
}

#' @title Process and normalize affymetrix-based dataset
#' @description This function normalize expression data of the downloaded GEO object.
#' The normalization process includes background correction using RMA and data normalization using quantile method.
#' This function is only used if the protocol is affymetrix or agilent.
#' This function is used internally by downloadGEO.
#' @param metadata The metadat of downloaded GEO object.
#' @param sampleIDs A vector of downloaded samples IDs from the queried GEO object.
#' @param destDir The user path to save the downloaded files.
#' @return A SummarizedExperiment object with appended expression values as assay and metadata as colData.
#' @details This function is used internally by downloadGEO.
#' @importFrom SummarizedExperiment SummarizedExperiment colData assay
#' @importFrom dplyr %>%
#' @importFrom Biobase exprs
#' @noRd
.processAffymetrix <- function(metadata, sampleIDs, destDir) {

    if (!dir.exists(destDir)) {
        stop("The destination directory does not exist.")
    }

    metadata <- metadata[metadata$geo_accession %in% sampleIDs,]

    if(nrow(metadata) < length(sampleIDs)){
        stop("The input metadata and sampleIDs do not match. Make sure the sample IDs match with geo_accession IDs from dataset.")
    }

    if (!.requirePackage("oligo")){
        return(NULL)
    }

    #Normalize expression data based on RMA method
    # expression <- try({
    #     ReadAffy(verbose = TRUE, celfile.path = destDir, sampleNames = sampleIDs, filenames = paste0(sampleIDs, '.CEL.gz')) %>%
    #         threestep(background.method = "RMA.2", normalize.method = "quantile", summary.method = "median.polish") %>%
    #         exprs() %>%
    #         as.data.frame()
    # })

    # if ("try-error" %in% class(expression)) {
    expression <- oligo::read.celfiles(file.path(destDir, paste0(sampleIDs, '.CEL.gz'))) %>%
        oligo::rma(normalize = TRUE) %>%
        Biobase::exprs() %>%
        as.data.frame()
    if (sum(is.na(expression)) > 0) stop("There is NA in expression data.")
    colnames(expression) <- sampleIDs
    # }

    if(dim(expression)[1] == 0 | dim(expression)[2] == 0){
        stop("The expression matrix is empty.")
    }

    summarizedExperiment <- SummarizedExperiment(
        assays = expression %>% log2() %>% as.matrix(),
        colData = data.frame(metadata)
    )

    return(summarizedExperiment)
}

#' @title Process and normalize agilent-based dataset
#' @description This function normalize expression data of the downloaded GEO object.
#' In the normalization process the limma package functions are utilized, including background correction using normexp and data normalization using loess and quantile methods.
#' This function is only used if the protocol is agilent.
#' This function is used internally by downloadGEO.
#' @param metadata The metadat of downloaded GEO object.
#' @param sampleIDs A vector of downloaded samples IDs from the queried GEO object.
#' @param destDir The user path to save the downloaded files.
#' @param greenOnly Logical, for use with source, should the green (Cy3) channel only be read, or are both red and green required.
#' @return A SummarizedExperiment object with appended expression values as assay and metadata as colData.
#' @details This function is used internally by downloadGEO.
#' @importFrom SummarizedExperiment SummarizedExperiment colData assay
#' @importFrom limma read.maimages backgroundCorrect normalizeWithinArrays normalizeBetweenArrays
#' @importFrom dplyr %>%
#' @noRd
.processAgilent <- function(metadata, sampleIDs, destDir, greenOnly) {

    if (!dir.exists(destDir)) {
        stop("The destination directory does not exist.")
    }

    metadata <- metadata[metadata$geo_accession %in% sampleIDs,]

    if(nrow(metadata) < length(sampleIDs)){
        stop("The input metadata and sampleIDs do not match. Make sure the sample IDs match with geo_accession IDs from dataset.")
    }

    raw.data <- read.maimages(file.path(destDir, paste0(sampleIDs, ".TXT.gz")),
                              source = "agilent",
                              green.only = greenOnly,
                              names = sampleIDs
    )

    #Correct expression for background using the normexp method
    background_corrected_data <- backgroundCorrect(raw.data, method = "normexp")

    if(inherits(background_corrected_data, "RGList")){
        # Normalize background-corrected data using the loess method for two color array
        norm1.data <- normalizeWithinArrays(background_corrected_data, method = "loess")

        # Normalize background-corrected data using the quantile method
        norm2.data <- normalizeBetweenArrays(norm1.data, method = "quantile")

        expression <- norm2.data$A
        rownames(expression) <- norm2.data$genes[, "ProbeName"]
    }else{
        # Normalize background-corrected data using the quantile method
        norm.data <- normalizeBetweenArrays(background_corrected_data, method = "quantile")

        expression <- norm.data$E
        rownames(expression) <- norm.data$genes[, "ProbeName"]
    }

    # if (!is.null(background_corrected_data$G) & !is.null(background_corrected_data$R)) {
    #     # Normalize background-corrected data using the loess method for two color array
    #     norm1.data <- normalizeWithinArrays(background_corrected_data, method = "loess")
    #
    #     # Normalize background-corrected data using the quantile method
    #     norm2.data <- normalizeBetweenArrays(norm1.data, method = "quantile")
    #
    #     expression <- norm2.data$G
    #     rownames(expression) <- norm2.data$genes[, "ProbeName"]
    # } else {
    #     # Normalize background-corrected data using the quantile method
    #     norm.data <- normalizeBetweenArrays(background_corrected_data, method = "quantile")
    #
    #     expression <- norm.data$E
    #     rownames(expression) <- norm.data$genes[, "ProbeName"]
    # }

    if(dim(expression)[1] == 0 | dim(expression)[2] == 0){
        stop("The expression matrix is empty.")
    }

    summarizedExperiment <- SummarizedExperiment(
        assays = expression %>% as.matrix(),
        colData = data.frame(metadata)
    )

    return(summarizedExperiment)
}

#' @title Download GEO data
#' @description This function download and process data from GEO for microarray and RNASeq data.
#' @param GEOID The ID of the GEO dataset.
#' @param platform The platform of selected GEO dataset.
#' @param destDir A path to save downloaded data.
#' @param protocol The protocol of the selected GEO dataset, including affymetrix and agilent.
#' @param greenOnly Logical, for use with source, should the green (Cy3) channel only be read, or are both red and green required.
#' @return A SummarizedExperiment object including the processed data.
#' @examples
#' \donttest{
#' library(RCPA)
#' # Affymetrix
#' downloadPath <- file.path(tempdir(), "GSE230534")
#' if(!dir.exists(downloadPath)) dir.create(downloadPath)
#' affyDataset <- RCPA::downloadGEO(GEOID = "GSE230534", platform ="GPL23126",
#'                                  protocol ="affymetrix", destDir = downloadPath)
#' }
#' @importFrom SummarizedExperiment SummarizedExperiment colData assay
#' @importFrom dplyr %>%
#' @importFrom Biobase pData
#' @export
downloadGEO <- function(GEOID, platform, protocol = c("affymetrix", "agilent"), destDir, greenOnly = TRUE) {
    protocol <- match.arg(protocol)
    protocol <- protocol %>% tolower()

    #Download data with the specified platform from GEO
    GEOObject <- .downloadGEOObject(GEOID, platform, destDir)

    #Extract metadata from GEOObject
    GEOObject.metadata <- GEOObject %>% pData()

    #Object to store processed and normalized data
    summarizedExperimentObject <- NULL

    #Download samples from the current GEO object
    sampleIDs <- GEOObject.metadata["geo_accession"] %>% apply(MARGIN = 1, FUN = function(x) x)
    downloadRes <- .downloadSamples(sampleIDs, protocol, destDir)

    if(downloadRes != TRUE){
        stop("There is an error in downloading samples.")
    }

    if(!is.logical(greenOnly)){
        stop("The greenOnly parameter can only be TRUE or FALSE! Please refer to document.")
    }

    if (protocol == "affymetrix"){
        #Normalize expression data for affymetrix using RMA method
        summarizedExperimentObject <- .processAffymetrix(GEOObject.metadata, sampleIDs, destDir)
    }else{
        #Normalize expression data for affymetrix using limma normexp, loess, and quantile methods
        summarizedExperimentObject <- .processAgilent(GEOObject.metadata, sampleIDs, destDir, greenOnly)
    }

    if(is.null(summarizedExperimentObject)){
        if (pkgEnv$isMissingDependency){
            return(NULL)
        }
        stop("There is an error in processing and normalizing data.")
    }

    return(summarizedExperimentObject)
}