#' @title Download GEO object
#' @description This function downloads the GEO object using getGEO function.
#' This function is used internally by downloadGEO.
#' @param GEOID The GEO dataset under query, ex. GSE14762.
#' @param platform The platform of selected GEO dataset, ex. GPL4866.
#' @param destDir The user path to save the downloaded files.
#' @return A GEO object
#' @details This function is used internally by downloadGEO.
#' @importFrom Biobase annotation
#' @importFrom GEOquery getGEO
.downloadGEOObject <- function(GEOID, platform, destDir){
  data <- NULL
  #Download GEO dataset
  gsets <- getGEO(GEOID, GSEMatrix = T, getGPL = T, destdir = destDir)
  #Filter GEO objects to keep object comaptible with the selected platform
  for (gset in gsets){
    if (annotation(gset) == platform)
      data <- gset
  }
  return(data)
}

#' @title Download samples from GEO object
#' @description This function downloads the corresponding samples from queried GEO object. This function is only used if the protocol is affymetrix or agilent.
#' This function is used internally by downloadGEO.
#' @param GEOObject The downloaded GEO object.
#' @param destDir The user path to save the downloaded files.
#' @return A vector of samples IDs of the queried GEO object.
#' @details This function is used internally by downloadGEO.
#' @importFrom GEOquery getGEOSuppFiles
#' @importFrom Biobase pData
#' @importFrom dplyr %>%
.downloadSamples <- function(GEOObject, destDir){
  samples =  pData(GEOObject)["geo_accession"] %>% apply(MARGIN = 1, FUN = function(x) x)
  for (sample in samples){
    if (file.exists(paste0(destDir, "/", sample, ".CEL.gz"))) next()
    else {
      downloadedFiles = getGEOSuppFiles(sample, baseDir = destDir, makeDirectory = F) %>% rownames()
      if(protocol == "affymetrix")
        downloadedFiles[grep(".cel.gz", downloadedFiles, ignore.case = T)] %>% file.rename(paste0(destDir, "/", sample, ".CEL.gz"))
      else
        downloadedFiles[grep(".txt.gz", downloadedFiles, ignore.case = T)] %>% file.rename(paste0(destDir, "/", sample, ".TXT.gz"))
    }
  }
  return(samples)
}

#' @title Process and normalize affymetrix-based dataset
#' @description This function normalize expression data of the downloaded GEO object.
#' The normalization process includes background correction using RMA and data normalization using quantile method.
#' This function is only used if the protocol is affymetrix or agilent.
#' This function is used internally by downloadGEO.
#' @param GEOObject The downloaded GEO object.
#' @param samples The downloaded samples from the queried GEO object.
#' @param destDir The user path to save the downloaded files.
#' @return A SummarizedExperiment object with appended expression values as assay and metadata as colData.
#' @details This function is used internally by downloadGEO.
#' @importFrom affy ReadAffy
#' @importFrom affyPLM threestep
#' @importFrom oligo read.celfiles rma
#' @importFrom Biobase pData
#' @importFrom SummarizedExperiment SummarizedExperiment colData assay
#' @importFrom dplyr %>%
.processAffymetrix <- function(GEOObject, samples, destDir){
  #Normalize expression data based on RMA method
  expression <-  try({
    ReadAffy(verbose = T, celfile.path = destDir, sampleNames = samples, filenames = paste0(samples, '.CEL.gz')) %>%
      threestep(background.method="RMA.2", normalize.method="quantile", summary.method="median.polish") %>%
      exprs() %>%
      as.data.frame()
  })

  if (class(expression) == "try-error"){
    expression <- read.celfiles(paste0(destDir, "/", samples, '.CEL.gz')) %>% rma(normalize=F) %>% exprs() %>% as.data.frame()
    if (sum(is.na(expression)) > 0) stop("NA in expression")
    colnames(expression) <- samples
  }

  metadata <- GEOObject %>% pData()

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
#' @param GEOObject The downloaded GEO object.
#' @param samples The downloaded samples from the queried GEO object.
#' @param destDir The user path to save the downloaded files.
#' @return A SummarizedExperiment object with appended expression values as assay and metadata as colData.
#' @details This function is used internally by downloadGEO.
#' @importFrom SummarizedExperiment SummarizedExperiment colData assay
#' @importFrom Biobase pData
#' @importFrom limma read.maimages backgroundCorrect normalizeWithinArrays normalizeBetweenArrays
#' @importFrom dplyr %>%
.processAgilent <- function(GEOObject, samples, destDir){
  raw.data <- read.maimages(paste0(destDir, "/", samples, ".TXT.gz"),
                            source = "agilent",
                            green.only = TRUE,
                            names = samples
  )
  #Correct expression for background using the normexp method
  background_corrected_data <- backgroundCorrect(raw.data, method = "normexp")

  if(!is.null(background_corrected_data$G) & !is.null(background_corrected_data$R)){
    # Normalize background-corrected data using the loess method for two color array
    norm1.data <- normalizeWithinArrays(background_corrected_data, method = "loess")

    # Normalize background-corrected data using the quantile method
    norm2.data <- normalizeBetweenArrays(norm1.data, method = "quantile")

    expression <- norm2.data$G
    rownames(expression) <- norm2.data$genes[,"ProbeName"]
  } else{
    # Normalize background-corrected data using the quantile method
    norm.data <- normalizeBetweenArrays(background_corrected_data, method = "quantile")

    expression <- norm.data$E
    rownames(expression) <- norm.data$genes[,"ProbeName"]
  }

  metadata <- GEOObject %>% pData()

  summarizedExperiment <- SummarizedExperiment(
    assays = expression %>% as.matrix(),
    colData = data.frame(metadata)
  )
  return(summarizedExperiment)
}

#' @title Download and process RNASeq-based dataset
#' @description This function downlaods the counts matrix from supplementary files of a specific GEO object.
#' This function is only used if the protocol is RNASeq
#' This function is used internally by downloadGEO.
#' @param GEOID The GEO dataset under query, ex. GSE165082.
#' @param GEOObject The downloaded GEO object.
#' @param destDir The user path to save the downloaded files.
#' @return A SummarizedExperiment object with appended counts data as assay and metadata as colData.
#' @details This function is used internally by downloadGEO.
#' @importFrom SummarizedExperiment SummarizedExperiment colData assay
#' @importFrom GEOquery getGEOSuppFiles
#' @importFrom Biobase pData
#' @importFrom dplyr %>% select
#' @importFrom utils read.table
.processRNASeq <- function(GEOID, GEOObject, destDir){
  #Download supplementary files
  supplementaryFiles <- getGEOSuppFiles(GEOID, fetch_files = TRUE, baseDir = destDir) %>% rownames()
  #Get the counts file
  countsFile <- supplementaryFiles[grepl("counts", supplementaryFiles)]
  #Normalize counts data
  countsData <- read.table(countsFile, header = TRUE, sep = "\t", fill = 0)
  expression <- (countsData[,-1] + 1) %>% log2()
  rownames(expression) <- countsData[,1]
  colnames(expression) <- lapply(colnames(countsData[,-1]), function (name) strsplit(name, "X")[[1]][2])

  metadata <- GEOObject %>% pData()
  rownames(metadata) <- GEOObject %>% pData() %>% select(title) %>% unlist() %>% as.vector()
  expression <- expression[, metadata$title]

  summarizedExperiment <- SummarizedExperiment(
    assays = expression %>% as.matrix(),
    colData = data.frame(metadata)
  )
  return(summarizedExperiment)
}

#' @title Download GEO data
#' @description This function download and process data from GEO for microarray and RNASeq data.
#' @param GEOID The ID of the GEO dataset
#' @param platform The platform of selected GEO dataset
#' @param protocol The protocol of the selected GEO dataset, including affymetrix, agilent, and RNASeq.
#' @return A SummarizedExperiment object including the processed data
#' @examples
#' downloadGEO("GSE20153", "GPL570", "affymetrix", getwd())
#' @importFrom SummarizedExperiment SummarizedExperiment colData assay
#' @importFrom dplyr %>%
#' @export
downloadGEO <- function(GEOID, platform, protocol, destDir){
  #Download data with the specified platform from GEO
  GEOObject <- .downloadGEOObject(GEOID, platform, destDir)

  protocol <- protocol %>% tolower()
  summarizedExperimentObject <- NULL
  if(protocol %in% c("affymetrix", "agilent")){
    #Download samples from the current GEO object
    samples <- .downloadSamples(GEOObject, destDir)

    if(protocol == "affymetrix")
      #Normalize expression data for affymetrix using RMA method
      summarizedExperimentObject <- .processAffymetrix(GEOObject, samples, destDir)
    else
      #Normalize expression data for affymetrix using limma normexp, loess, and quantile methods
      summarizedExperimentObject <- .processAgilent(GEOObject, samples, destDir)
  }else{
    #Download and normalize read counts data for RNASeq
    summarizedExperimentObject <- .processRNASeq(GEOID, GEOObject, destDir)
  }

  return(summarizedExperimentObject)
}