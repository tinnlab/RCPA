.RCyjs <- setClass("RCyjs",
                   representation = representation(graph = "graph"),
                   contains = "BrowserViz"
)

#' Plot a pathway network
#' @description This function plots a pathway network.
#' @param results A named list of data frame of Pathway analysis results.
#' The columns of each data frame should be at least ID, name, p.value and pFDR, ES, NES, and nDE.
#' @param genesets A named list of character vectors of gene sets.
#' @param labels A named character vector of labels for the results.
#' @param pThreshold A numeric value of p-value threshold.
#' @param useFDR A logical value indicating whether to use FDR or not.
#' @param edgeThreshold A numeric from 0 to 1 indicating the threshold to draw edges.
#' edgeThreshold of 0.1 means that edges are drawn if the number of genes in common is greater than 1% of the smaller gene set.
#' @param styleFile A character value of the path to the style file.
#' If NULL, the default style file will be used, which is located at system.file(package="RCPA", "extdata", "pieStyle.js")
#' @return A RCyjs object.
#' @export
#' @examples
#' #TODO add examples
#' @importFrom RCyjs setGraph redraw layout loadStyleFile
#' @importFrom BrowserViz BrowserViz
#' @importFrom graph graphNEL addEdge nodeDataDefaults nodeData
#' @importFrom dplyr filter mutate %>%
plotPathwayNetwork <- function(results, genesets, labels = NULL, pThreshold = 0.05, useFDR = TRUE, edgeThreshold = 0.1, styleFile = NULL) {

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

    pathwayInfo <- data.frame(
        ID = names(genesets),
        label = if (!is.null(labels)) labels else names(genesets),
        size = sapply(names(genesets), function(x) length(genesets[[x]])),
        nResult = length(results)
    )

    for (i in seq_along(results)) {
        pathwayInfo[[paste0("pvalue", i)]] <- NA
        pathwayInfo[[paste0("isSig", i)]] <- NA
        pathwayInfo[[paste0("nDE", i)]] <- NA

        idx <- match(pathwayInfo$ID, results[[i]]$ID)

        pathwayInfo[[paste0("nDE", i)]] <- results[[i]]$nDE[idx]

        if (useFDR) {
            pathwayInfo[[paste0("pvalue", i)]] <- results[[i]]$pFDR[idx]
            pathwayInfo[[paste0("isSig", i)]] <- results[[i]]$pFDR[idx] < pThreshold
        } else {
            pathwayInfo[[paste0("pvalue", i)]] <- results[[i]]$p.value[idx]
            pathwayInfo[[paste0("isSig", i)]] <- results[[i]]$p.value[idx] < pThreshold
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
                nCommon = length(intersect(genesets[[i]], genesets[[j]])),
                stringsAsFactors = FALSE
            )
        }) %>% do.call(what = rbind)
    }) %>%
        do.call(what = rbind) %>%
        filter(.$nCommon / min(.$from.size, .$to.size) > edgeThreshold) %>%
        mutate(
            weight = nCommon
        )

    graphObj <- graphNEL(pathwayInfo$ID, edgemode = "undirected")
    graphObj <- graph::addEdge(graphEdges$from, graphEdges$to, graphObj, graphEdges$weight)

    for (attr in colnames(pathwayInfo)) {
        nodeDataDefaults(graphObj, attr = attr) <- NA
        nodeData(graphObj, pathwayInfo$ID, attr) <- pathwayInfo[[attr]]
    }

    rCy <- BrowserViz(portRange = 10000:20000,
                      title = "Pathway network",
                      quiet = T,
                      browserFile = system.file(package="RCyjs", "browserCode", "dist", "rcyjs.html"),
                      httpQueryProcessingFunction = cyjsQueryFnc,
    ) %>% .RCyjs(graph = graphObj)

    RCyjs::setGraph(rCy, graph = graphObj)
    RCyjs::loadStyleFile(rCy, styleFile)
    RCyjs::redraw(rCy)
    Sys.sleep(2)
    RCyjs::layout(rCy, "cose")

    rCy
}

