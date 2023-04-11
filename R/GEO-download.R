#' @title Download GEO object
#' @description This function downloads the GEO object using getGEO function.
#' This function is used internally by downloadGEO.
#' @param GEOID The GEO dataset under query, ex. GSE14762.
#' @param platform The platform of selected GEO dataset, ex. GPL4866.
#' @param destDir The user path to save the downloaded files.
#' @return A GEO object
#' @details This function is used internally by downloadGEO.
.downloadGEOObject <- function(GEOID, platform, destDir){
  data <- NULL
  destdir <- destDir
  #Download GEO dataset
  gsets <- getGEO(GEOID, GSEMatrix = T, getGPL = T, destdir = getwd())
  #Filter GEO objects to keep object comaptible with the selected platform
  for (gset in gsets){
    if (Biobase::annotation(gset) == platform)
      data <- gset
  }
  return(data)
}

#' @title Extract applied protocol on GEO dataset
#' @description This function extracts the applied protocol to prepare dataset. The protocols include agilent, affymetrix, and RNASeq. 
#' This function is used internally by downloadGEO.
#' @param GEOObject The downloaded GEO object.
#' @return A string containing protocol information. 
#' @details This function is used internally by downloadGEO.
.getProtocol <- function(GEOObject){
  #Select one row of metadata of GEO object
  cur.row <- pData(GEOObject)[1,] %>% as.vector() %>% lapply(tolower)
  
  #Search through all values of selected row to find protocol information
  protocol <- ifelse(any(grepl("agilent", cur.row)), "agilent", ifelse(any(grepl("affymetrix", cur.row)), "affymetrix", "RNASeq"))
  
  return(protocol)
}

#' @title Download samples from GEO object
#' @description This function downloads the corresponding samples from queried GEO object. This function is only used if the protocol is affymetrix or agilent.
#' This function is used internally by downloadGEO.
#' @param GEOObject The downloaded GEO object.
#' @param destDir The user path to save the downloaded files.
#' @return A vector of samples IDs of the queried GEO object.
#' @details This function is used internally by downloadGEO.
.downlaodSamples <- function(GEOObject, destDir){
  samples =  pData(GEOObject)["geo_accession"] %>% apply(MARGIN = 1, FUN = function(x) x)
  destdir <- destDir
  for (sample in samples){
    if (file.exists(paste0(destDir, "/", sample, ".CEL.gz"))) next()
    else {
      downloadedFiles = getGEOSuppFiles(sample, baseDir = destDir, makeDirectory = F) %>% rownames()
      downloadedFiles[grep(".cel.gz", downloadedFiles, ignore.case = T)] %>% file.rename(paste0(destDir, "/", sample, ".CEL.gz"))
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
.processAffymetrix <- function(GEOObject, samples, destDir){
  #Normalize expression data based on RMA method
  expression <-  try({
    ReadAffy(verbose = T, celfile.path = destDir, sampleNames = samples, filenames = paste0(samples, '.CEL.gz')) %>%
      threestep(background.method="RMA.2", normalize.method="quantile", summary.method="median.polish") %>%
      exprs() %>%
      as.data.frame()
  })
  
  if (class(expression) == "try-error"){
    expression <- oligo::read.celfiles(paste0(destDir, "/", samples, '.CEL.gz')) %>% oligo::rma(normalize=F) %>% exprs() %>% as.data.frame()
    if (sum(is.na(expression)) > 0) stop("NA in expression")
    colnames(expression) <- samples
  }
  
  metadata <- GEOObject %>% pData()
  
  summarizedExperiment <- SummarizedExperiment::SummarizedExperiment(
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
.processAgilent <- function(samples, destDir){
  raw.data <- limma::read.maimages(paste0(destDir, "/", samples), 
                                   source = "agilent",
                                   green.only = TRUE, 
                                   names = samples
                                   )
  #Correct expression for background using the normexp method
  background_corrected_data <- limma::backgroundCorrect(raw.data, method = "normexp")
  background_corrected_data$R <- log2(background_corrected_data$R)
  background_corrected_data$G <- log2(background_corrected_data$G)
  
  # Normalize background-corrected data using the loess method
  expression <- limma::normalizeWithinArrays(background_corrected_data, method = "loess")
  # Normalize background-corrected data using the quantile method
  expression <- limma::normalizeBetweenArrays(background_corrected_data, method = "quantile")
  
  metadata <- GEOObject %>% pData()
  
  summarizedExperiment <- SummarizedExperiment::SummarizedExperiment(
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
#' @param destDir The user path to save the downloaded files.
#' @return A SummarizedExperiment object with appended counts data as assay and metadata as colData.
#' @details This function is used internally by downloadGEO.
.processRNASeq <- function(GEOID, destDir){
  #Download supplementary files
  supplementary.files <- getGEOSuppFiles(GEOID, fetch_files = TRUE, baseDir = destDir) %>% rownames()
  #Get the counts file
  counts.file <- supplementary.files[grepl("counts", supplementary.files)]
  #Normalize counts data
  expression <- read.table(counts.file, header = TRUE, sep = "\t", fill = NA) %>% log2()
  
  metadata <- .downloadGEOObject(GEOID, platform, destDir) %>% pData()
  
  summarizedExperiment <- SummarizedExperiment::SummarizedExperiment(
    assays = expression %>% as.matrix(),
    colData = data.frame(metadata)
  )
  return(summarizedExperiment)
}

#' @title Download GEO data
#' @description This function download and process data from GEO for microarray and RNASeq data.
#' @param GEOID The ID of the GEO dataset
#' @param platform The platform of selected GEO dataset
#' @return A SummarizedExperiment object including the processed data
#' @examples
#'
#' @importFrom SummarizedExperiment SummarizedExperiment rowData assay
#' @importFrom GEOQuery getGEO
#' @importFrom S4Vectors SimpleList
#' @export
downloadGEO <- function(GEOID, platform, destDir){
  #Download data with the specified platform from GEO  
  GEOObject <- .downloadGEOObject(GEOID, platform, destDir)
  
  #Extract protocol information from GEO object
  protocol <- .getProtocol(GEOObject)
   
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
    summarizedExperimentObject <- .processRNASeq(GEOID, destDir)
  }
    
  return(summarizedExperimentObject)
}