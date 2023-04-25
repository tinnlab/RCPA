#' @title Combine PValues
#' @description This function combines P-Values based on the selected method.
#' This function is used internally by combinePathwayAnalysisResults.
#' @param inputData A dataframe of pathway analysis results from multiple studies.
#' @param method The selected method to combine PValues, including Fisher, Stouffer, addCLT, GeoMean, and minP.
#' @return A dataframe of meta analysis results.
#' @details This function is used internally by combinePathwayAnalysisResults.
#' @importFrom dplyr %>%
.combinePvalues <- function(inputData, method){

    combineFunc <- switch(method,
           fisher = .runFisher,
           stouffer = .runStouffer,
           minP = min,
           addCLT = .runAddCLT,
           geoMean = .runGeoMean
          )

    metaPvalRes <- inputData %>% group_by(ID) %>% group_split() %>% lapply(function (data){
                pValues <- data$p.value
                meta.Pval <- combineFunc(pValues)
                data.frame(
                  ID = data$ID[1],
                  p.value = meta.Pval,
                  stringsAsFactors = FALSE
                )
              }) %>% do.call(what = rbind)

    metaScoreRes <- .combineEnrichmentScores(inputData)

    metaRes <- metaPvalRes
    metaRes$score <- metaRes$score.sd <- metaRes$count <- 0
    metaRes$score <- metaScoreRes$score[match(metaScoreRes$ID, metaRes$ID)]
    metaRes$score.sd <- metaScoreRes$score.sd[match(metaScoreRes$ID, metaRes$ID)]
    metaRes$count <- metaScoreRes$count[match(metaScoreRes$ID, metaRes$ID)]

    metaRes <- metaRes[, c("ID", "p.value", "score", "score.sd", "count")]

    return(metaRes)
}

#' @title Combine P-Valuse using Fisher
#' @description This function combines P-Values based on Fisher method.
#' This function is used internally by .combinePvalues.
#' @param pvals The vector of P-Values to be combined.
#' @return A combined P-Value
#' @details This function is used internally by .combinePvalues.
#' @importFrom stats qnorm pchisq
.runFisher <- function(pvals){
    p.value <- pchisq(-2 * sum(log(pvals)), df=2*length(pvals), lower.tail=FALSE)
    return(p.value)
}

#' @title Combine P-Valuse using Stouffer
#' @description This function combines P-Values based on Stouffer method.
#' This function is used internally by .combinePvalues.
#' @param pvals The vector of P-Values to be combined.
#' @return A combined P-Value
#' @details This function is used internally by .combinePvalues.
#' @importFrom stats pnorm qnorm
.runStouffer <- function(pvals){
    p.value <- pnorm(sum(qnorm(pvals)) / sqrt(length(pvals)))
    return(p.value)
}

#' @title Combine P-Valuse using addCLT
#' @description This function combines P-Values based on addCLT method.
#' This function is used internally by .combinePvalues.
#' @param pvals The vector of P-Values to be combined.
#' @return A combined P-Value
#' @details This function is used internally by .combinePvalues.
#' @importFrom stats pnorm
.runAddCLT <- function(pvals){
    n <- length(pvals)
    p.value <- 1
    if(n <= 20){
        x <- sum(pvals)
        p.value <- 1/factorial(n) * sum(sapply(0:floor(x), function(k) (-1)^k * choose(n,k) * (x-k)^(n)))
    }else{
        p.value <- pnorm(sum(pvals),n/2,sqrt(n/12),lower.tail=TRUE)
    }
    return(p.value)
}

#' @title Combine P-Valuse using Geometric Mean
#' @description This function combines P-Values by computing their geometric mean.
#' This function is used internally by .combinePvalues.
#' @param pvals The vector of P-Values to be combined.
#' @return A combined P-Value
#' @details This function is used internally by .combinePvalues.
.runGeoMean <- function(pvals){
    p.value <- exp(mean(log(pvals)))
    return(p.value)
}

#' @title Combine Normalized Enrichment Scores
#' @description This function combines normalized enrichment scores based on REML method.
#' This function is used internally by combinePathwayAnalysisResults and .combinePvalues.
#' @param inputData A dataframe of pathway analysis results from multiple studies.
#' @return A dataframe of meta analysis results.
#' @details This function is used internally by combinePathwayAnalysisResults and .combinePvalues.
#' @importFrom meta metagen
#' @importFrom dplyr %>%
.combineEnrichmentScores <- function(inputData){

    metaRes <- inputData %>% group_by(ID) %>% group_split() %>% lapply(function(data){
        data$normalizedScore.sd <- abs((data$normalizedScore - ifelse(data$normalizedScore > 0, 1, -1))/qnorm(data$p.value))
        meta.res <- metagen(data = data,
                            studlab = ID,
                            TE = normalizedScore ,
                            seTE = normalizedScore.sd,
                            sm = "SMD",
                            n.e = sample.size,
                            method.tau = "REML",
                            hakn = TRUE)

        normalizedScore.combined <- meta.res$TE.fixed
        normalizedScore.combined.sd <- meta.res$seTE.fixed

        pval <- pnorm((ifelse(normalizedScore.combined > 0, 1, -1) - normalizedScore.combined)/normalizedScore.combined.sd)
        if (normalizedScore.combined < 0) pval <- 1 - pval

        data.frame(
          ID = data$ID[1],
          p.value = pval,
          score = normalizedScore.combined,
          score.sd = normalizedScore.combined.sd,
          count = nrow(data),
          stringsAsFactors = F
        )
      }) %>% do.call(what = rbind) %>% as.data.frame()

    metaRes <- metaRes[, c("ID", "p.value", "score", "score.sd", "count")]

    return(metaRes)
}

#' @title Perform Meta Analysis
#' @description This function performs mata analysis on multiple pathway analysis results.
#' @param DFsList A list of dataframes obtained from runPathwayAnalysis.
#' @param method A method used to combine pathway analysis results
#' @return A dataframe of meta analysis results including combined normalized score and combined p-value for each pathway.
#' @details This function performs mata analysis on multiple pathway analysis results.
#' @importFrom dplyr %>% bind_rows
#' @export
combinePathwayAnalysisResults <- function(DFsList, method = c("fisher", "stouffer", "min", "score", "geoMean", "addCLT", "minP")){

    method <- match.arg(method)

    if(length(DFsList) == 1){
      stop("Meta analysis is valid for two or more studies.")
    }

    if(any(sapply(DFsList, is.null))) {
      stop("There is null object in the input list.")
    }

    dim.lst <- DFsList %>% lapply(function(df) dim(df))
    if(length(unique(dim.lst)) != 1){
      stop("All the dataframes in the input list must have the same dimension.")
    }

    colnames.lst <- DFsList %>% lapply(function(df) colnames(df))
    if(length(unique(colnames.lst)) != 1){
      stop("All the dataframes in the input list must have the same column names.")
    }

    rownames.lst <- DFsList %>% lapply(function(df) rownames(df))
    if(length(unique(rownames.lst)) != 1){
      stop("All the dataframes in the input list must have the same row names.")
    }

    allData <- bind_rows(DFsList)

    metaResult <- NULL
    if(method == "score")
      metaResult <- .combineEnrichmentScores(allData)
    else metaResult <- .combinePvalues(allData)

    if(is.null(metaResult)){
        stop("There is an error in meta analysis.")
    }

    metaResult$name <- allData$name[match(metaResult$ID, allData$ID)]
    metaResult$pathway.size <- allData$pathway.size[match(metaResult$ID, allData$ID)]

    return(metaResult)
}