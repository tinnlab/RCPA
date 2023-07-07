#' @title Plot a pathway network
#' @description This function plots a pathway network.
#' This function needs an interactive environment with browser view support to work.
#' @param PAResults A named list of data frame of Pathway analysis results.
#' The columns of each data frame should be at least ID, name, p.value and pFDR.
#' An optional column "color" can be used to specify the color of the nodes.
#' If the column "color" is not specified, the color of the nodes will be determined by the mode and the statistic.
#' @param genesets A named list of character vectors of gene sets.
#' @param statistic A character value of the statistic to use for the network.
#' The statistic should be one of the columns of the data frame in the results.
#' @param mode A character value of the mode to use to color the nodes.
#' The mode should be one of "discrete", "continous".
#' If the mode is "discrete", the color of the nodes are determined by whether the p-value is significant or not.
#' If the mode is "continous", the color of the nodes are proportional to the statistic.
#' @param labels A named character vector of labels for the results.
#' @param pThreshold A numeric value of p-value threshold.
#' @param useFDR A logical value indicating whether to use FDR or not.
#' This parameter is independent of the pThreshold.
#' @param edgeThreshold A numeric from 0 to 1 indicating the threshold to draw edges.
#' edgeThreshold of 0.1 means that edges are drawn if the number of genes in common is greater than 1\% of the smaller gene set.
#' @param statLimit A numeric value of the maximum absolute value of the statistic.
#' If statistic is p.value or pFDR, this parameter is the limit of -log10(p-value).
#' @param discreteColors A character vector of colors to use for the discrete mode.
#' The length of the vector must be the same as the number of results.
#' @param continuousScaleFunc A function that takes a numeric value from -1 to 1 and returns a color.
#' @param NAColor A character value of the color to use for NA values.
#' @param borderColor A character value of the color to use for the border of the nodes.
#' @param nodeSizeFnc A function that takes a character value of the ID of the node and returns a numeric value of the size of the node.
#' @param borderWidthFnc A function that takes a character value of the ID of the node and returns a numeric value of the width of the border of the node.
#' @param edgeWidthFnc A function that takes a character value of the ID of the from node and a character value of the ID of the to node and returns a numeric value of the width of the edge.
#' @param styleFile A character value of the path to the style file.
#' If NULL, the default style file will be used, which is located at system.file(package="RCPA", "extdata", "pieStyle.js")
#' @return A RCyjs object.
#' @details
#' The function will plot a pathway network using the results of pathway analysis.
#' The nodes of the network are the pathways and the edges are the pathways that have at least a certain number of genes in common defined by the edgeThreshold.
#' The size of the nodes are proportional to the number of genes in the pathway.
#' The color of the nodes are proportional to the statistic used if the mode is "continuous".
#' If the mode is "discrete", the color of the nodes are determined by whether the p-value is significant or not.
#' The width of the edges are proportional to the number of genes in common.
#' @export
#' @examples
#' \donttest{
#' # Loading libraries
#' #if (!require("RCyjs", quietly = TRUE)) {
#'  #BiocManager::install("RCyjs")
#' #}
#' #library(RCyjs)
#' library(RCPA)
#' affyFgseaResult <- loadData("affyFgseaResult")
#' agilFgseaResult <- loadData("agilFgseaResult")
#' RNASeqFgseaResult <- loadData("RNASeqFgseaResult")
#' metaPAResult <- loadData("metaPAResult")
#' genesets <- loadData("genesets")
#'
#' PAResults <- list(
#'     "Affymetrix - GSE5281" = affyFgseaResult,
#'     "Agilent - GSE61196" = agilFgseaResult,
#'     "RNASeq - GSE153873" = RNASeqFgseaResult,
#'     "Meta-analysis" = metaPAResult
#' )
#'
#' genesetsToPlot <- metaPAResult$ID[order(metaPAResult$pFDR)][1:30]
#'if (require("RCyjs", quietly = TRUE)){
#' pltObj <- RCPA::plotPathwayNetwork(
#'     PAResults,
#'     genesets = genesets$genesets[genesetsToPlot],
#'     labels = genesets$names[genesetsToPlot],
#'     edgeThreshold = 0.75,
#'     mode = "continuous",
#'     statistic = "normalizedScore"
#' )
#' }
#'
#' # pltObj$plot() This function requires a browser to work
#'
#' }
#'
#' @importFrom graph graphNEL addEdge `nodeDataDefaults<-` `nodeData<-`
#' @importFrom grDevices colorRampPalette
plotPathwayNetwork <- function(PAResults, genesets,
                               statistic = "pFDR",
                               mode = c("continuous", "discrete"),
                               labels = NULL,
                               pThreshold = 0.05,
                               useFDR = TRUE,
                               edgeThreshold = 0.5,
                               statLimit = 4,
                               discreteColors = NULL,
                               continuousScaleFunc = NULL,
                               NAColor = "#dddddd",
                               borderColor = "#333333",
                               nodeSizeFnc = function(id) length(genesets[[id]])^.75,
                               borderWidthFnc = function(id) 1,
                               edgeWidthFnc = function(from, to) 1,
                               styleFile = system.file(package = "RCPA", "extdata", "pieStyle.js")) {

    mode <- match.arg(mode)

    if (!.requirePackage("BrowserViz")){
        return(NULL)
    }
    if (!.requirePackage("RCyjs")){
        return(NULL)
    }

    cyjsQueryFnc <- function(queryString)
    {
        ampersand.loc <- as.integer(regexpr("&", queryString, fixed = TRUE))
        if (ampersand.loc > 0) {
            queryString <- substring(queryString, 1, ampersand.loc - 1);
        }
        questionMark.loc <- as.integer(regexpr("?", queryString, fixed = TRUE));

        if (questionMark.loc == 1)
            queryString <- substring(queryString, 2, nchar(queryString))

        filename <- queryString

        if (!file.exists(filename)) {
            return(list(contentType = "text/plain", body = sprintf("file not found: %s", filename)))
        }
        text <- paste(scan(filename, what = character(0), sep = "\n", quiet = TRUE), collapse = "\n")
        return(list(contentType = "text/plain", body = text));
    }

    if (is.null(styleFile)) {
        styleFile <- system.file(package = "RCPA", "extdata", "pieStyle.js")
    }

    for (res in PAResults) {
        if (!all(c("ID", "name", "p.value", "pFDR") %in% colnames(res))) {
            stop("The columns of the data frame in the results should be at least ID, name, p.value, and pFDR")
        }

        if (!statistic %in% colnames(res)) {
            stop("The statistic should be one of the columns of the data frame in the results.")
        }
    }

    if (is.null(discreteColors)) {
        discretePieColors <- c(
            "#316b9d",
            "#fce397",
            "#99cc83",
            "#f77a65",
            "#a6a1d0",
            "#fea9c4",
            "#74e7bc",
            "#febb73",
            "#1db4db",
            "#ffc5a6",
            "#b6c9fa",
            "#ee5437"
        )
    } else if (mode == "discrete") {
        if (length(discreteColors) != length(PAResults)) {
            stop("The length of the discreteColors should be the same as the number of results.")
        }
    }

    if (is.null(continuousScaleFunc)) {

        continuousScaleFunc <- function(x) {
            x <- ifelse(x >= 1, 0.999999, x)
            x <- ifelse(x <= -1, -0.999999, x)

            colorRampPalette(c("blue", "white", "red"))(1000)[ceiling((x + 1) / 2 * 1000)]
        }

    }

    pathwayInfo <- data.frame(
        ID = names(genesets),
        label = if (!is.null(labels)) labels else names(genesets),
        size = sapply(names(genesets), nodeSizeFnc),
        borderWidth = sapply(names(genesets), borderWidthFnc),
        borderColor = borderColor,
        nResult = length(PAResults),
        legendTitle = (
            if (statistic == "p.value") "-log10 p-value"
            else if (statistic == "pFDR") "-log10 FDR"
            else statistic
        ),
        statLimitUpper = statLimit,
        statLimitLower = (
            if (statistic == "p.value") 0
            else if (statistic == "pFDR") 0
            else -statLimit
        )
    )

    for (i in seq_along(PAResults)) {

        pathwayInfo[[paste0("stat", i)]] <- NA
        pathwayInfo[[paste0("isSig", i)]] <- FALSE
        # pathwayInfo[[paste0("nDE", i)]] <- NA

        statValues <- PAResults[[i]][[statistic]]
        if (statistic == "p.value" || statistic == "pFDR") {
            statValues <- -log10(statValues)
            statValues[statValues > statLimit] <- statLimit
            statValues <- statValues / statLimit
        } else {
            statValues[statValues > statLimit] <- statLimit
            statValues[statValues < -statLimit] <- -statLimit
            statValues <- statValues / statLimit
        }

        idx <- match(pathwayInfo$ID, PAResults[[i]]$ID)

        pathwayInfo[[paste0("stat", i)]] <- statValues[idx]
        # pathwayInfo[[paste0("nDE", i)]] <- results[[i]]$nDE[idx]

        if (useFDR) {
            pathwayInfo[[paste0("isSig", i)]] <- PAResults[[i]]$pFDR[idx] < pThreshold
        } else {
            pathwayInfo[[paste0("isSig", i)]] <- PAResults[[i]]$p.value[idx] < pThreshold
        }

        if (mode == "discrete") {
            pathwayInfo[[paste0("color", i)]] <- ifelse(pathwayInfo[[paste0("isSig", i)]], discretePieColors[i], NAColor)
        } else {
            pathwayInfo[[paste0("color", i)]] <- continuousScaleFunc(pathwayInfo[[paste0("stat", i)]])
            pathwayInfo[[paste0("color", i)]][is.na(pathwayInfo[[paste0("stat", i)]])] <- NAColor
        }
    }

    gsIDs <- pathwayInfo$ID
    graphEdges <- lapply(1:(length(gsIDs) - 1), function(i) {
        lapply((i + 1):length(gsIDs), function(j) {
            data.frame(
                from = gsIDs[i],
                to = gsIDs[j],
                from.size = length(genesets[[i]]),
                to.size = length(genesets[[j]]),
                common = length(intersect(genesets[[i]], genesets[[j]])),
                stringsAsFactors = FALSE
            )
        }) %>% do.call(what = rbind)
    }) %>%
        do.call(what = rbind) %>%
        filter(.$common / min(.$from.size, .$to.size) > edgeThreshold)

    graphEdges$weight <- apply(graphEdges[, c("from", "to")], 1, function(x) {
        edgeWidthFnc(x[1], x[2])
    })

    graphObj <- graphNEL(pathwayInfo$ID, edgemode = "undirected")
    graphObj <- graph::addEdge(graphEdges$from, graphEdges$to, graphObj, graphEdges$weight)

    if (is.null(names(PAResults))) {
        names(PAResults) <- paste0("Result ", seq_along(PAResults))
    }

    for (i in seq_along(PAResults)) {
        pathwayInfo[[paste0("resultName", i)]] <- names(PAResults)[i]
    }

    for (attr in colnames(pathwayInfo)) {
        nodeDataDefaults(graphObj, attr = attr) <- NA
        nodeData(graphObj, pathwayInfo$ID, attr) <- pathwayInfo[[attr]]
    }

    rCy <- BrowserViz::BrowserViz(portRange = 10000:20000,
                      title = "Pathway network",
                      quiet = TRUE,
                      browserFile = system.file(package = "RCyjs", "browserCode", "dist", "rcyjs.html"),
                      httpQueryProcessingFunction = cyjsQueryFnc
    )

    class(rCy) <- "RCyjs"

    list(
        rCy = rCy,
        graphObj = graphObj,
        pathwayInfo = pathwayInfo,
        plot = function(layout = "cose") {
            RCyjs::setGraph(rCy, graph = graphObj)
            RCyjs::setBackgroundColor(rCy, "#ffffff")
            try({RCyjs::loadStyleFile(rCy, styleFile)})
            RCyjs::redraw(rCy)
            Sys.sleep(2)
            RCyjs::layout(rCy, layout)
        }
    )
}

