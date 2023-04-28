#' @title Plot forest plot from pathway/geneset/meta analysis results
#' @description Plot forest plot from pathway/geneset/meta analysis results.
#' @param resultsList A named list of dataframes from pathway analysis, geneset analysis, and/or meta analysis results.
#' The columns are ID, name, description, p.value, pFDR, size, nDE, score and normalizedScore.
#' @param yAxis The column to use for the y-axis.
#' @param statLims A vector of length 2 specifying the limits for score to use in the x axis.
#' @param useFDR A boolean parameter to specify if adjusted p-value should be considered.
#' @return A list of ggplot objects.
#' @examples
#' \dontrun{
#' # Loading the library
#' library(RCPA)
#' library(ggplot2)
#' # Set seed for repro
#' set.seed(123)
#' # Simulate the enrichment analysis for the study 1
#' result1 <- data.frame(
#'   ID = paste0("geneset", 1:100),
#'   name = paste0("Pathway ", 1:100),
#'   description = paste0("Description ", 1:100),
#'   p.value = runif(100) / 10,
#'   pFDR = runif(100) / 5,
#'   size = runif(100, 100, 500),
#'   nDE = runif(100, 10, 100),
#'   score = runif(100, -2, 2),
#'   normalizedScore = runif(100)
#' )
#' rownames(result1) <- result1$ID
#' 
#' # Simulate the enrichment analysis for the study 2
#' result2 <- data.frame(
#'   ID = paste0("geneset", 1:100),
#'   name = paste0("Pathway ", 1:100),
#'   description = paste0("Description ", 1:100),
#'   p.value = runif(100) / 10,
#'   pFDR = runif(100) / 5,
#'   size = runif(100, 100, 500),
#'   nDE = runif(100, 10, 100),
#'   score = runif(100, -2, 2),
#'   normalizedScore = runif(100)
#' )
#' rownames(result2) <- result2$ID
#' 
#' # Store results from two studies in a list
#' resultsLst <- list(
#'   "study1" = result1,
#'   "study2" = result2
#' )
#' 
#' # Forest plot
#' plotForest(resultsLst, yAxis = "ID")
#' }
#' @importFrom ggplot2 ggplot aes geom_point geom_hline theme_minimal theme theme_bw geom_vline scale_color_gradient scale_size_continuous labs geom_rect coord_cartesian geom_errorbarh ggtitle unit
#' @importFrom gridExtra grid.arrange
#' @importFrom dplyr %>%
#' @export
plotForest <- function(resultsList, yAxis = c("ID", "name"), statLims = c(-2.5, 2.5), useFDR = TRUE) {
    yAxis <- match.arg(yAxis)

    for (result in resultsList) {
        if (is.null(result$ID) |
            is.null(result$name) |
            is.null(result$normalizedScore)) {
            stop("All dataframes in the input list must have 'ID', 'name', 'normalizedScore' columns.")
        }

        if (useFDR & is.null(result$pFDR)) {
            stop("All dataframes in the input list must have 'pFDR' column.")
        }

        if (!useFDR & is.null(result$p.value)) {
            stop("All dataframes in the input list must have 'p.value' column.")
        }
    }

    commonIds <- Reduce(intersect, lapply(resultsList, function(x) x$ID))
    if (length(commonIds) == 0) {
        stop("All dataframes in the input list must have at least one common ID.")
    }

    if (!is.numeric(statLims[1]) | !is.numeric(statLims[2])) {
        stop("statLims should be a numeric vector of length two.")
    }

    if (statLims[1] >= statLims[2]) {
        stop("statLims parameter must be a vector of two values with minimum as the first value and maximum as the second value. ")
    }

    if (is.null(names(resultsList))) {
        names(resultsList) <- paste0("Dataset", 1:length(resultsList))
    }

    plotData <- lapply(names(resultsList), function(n) {
        data <- resultsList[[n]]
        rownames(data) <- data$ID
        data[commonIds, unique(c("ID", yAxis, "normalizedScore", "p.value", "pFDR"))] %>% mutate(dataset = n)
    }) %>% do.call(what = rbind)

    plotData$dataset <- factor(plotData$dataset, levels = names(resultsList))

    pvalues <- if (useFDR) plotData$pFDR else plotData$p.value

    plotData$sd <- abs((plotData$normalizedScore - ifelse(plotData$normalizedScore > 0, 1, -1)) / qnorm(pvalues))
    plotData$sd[plotData$sd > 0.5] <- 0.5
    plotData$min <- plotData$normalizedScore - plotData$sd * 2
    plotData$max <- plotData$normalizedScore + plotData$sd * 2
    plotData$min[plotData$min < statLims[1]] <- statLims[1]
    plotData$max[plotData$max > statLims[2]] <- statLims[2]
    plotData$label <- factor(plotData[[yAxis]], levels = unique(plotData[[yAxis]]))

    statRange <- statLims[2] - statLims[1]
    gap <- 0.5

    plotData$normalizedScoreShifted <- plotData$normalizedScore +
        (as.numeric(plotData$dataset) - 1) * statRange +
        (as.numeric(plotData$dataset)) * gap
    plotData$minShifted <- plotData$min +
        (as.numeric(plotData$dataset) - 1) * statRange +
        (as.numeric(plotData$dataset)) * gap
    plotData$maxShifted <- plotData$max +
        (as.numeric(plotData$dataset) - 1) * statRange +
        (as.numeric(plotData$dataset)) * gap

    x2_breaks <- (seq_along(resultsList) - 1) * statRange - statRange / 2 + statLims[2]
    x2_labels <- names(resultsList)

    min_values <- (seq_along(resultsList) - 1) * statRange - statRange / 2 + gap / 2
    max_values <- (seq_along(resultsList) - 1) * statRange + statRange / 2 - gap / 2
    x1_breaks <- NULL

    for (i in 1:length(min_values)) {
        x1_breaks <- c(x1_breaks, ceiling(min_values[i]):floor(max_values[i]))
    }

    x1_labels <- rep(c(ceiling(min_values[1]):floor(max_values[1])), length(min_values))

    ggplot() +
        geom_rect(
            plotData,
            mapping = aes(
                xmin = -Inf,
                xmax = Inf,
                ymin = as.numeric(.data$label) - 0.5,
                ymax = as.numeric(.data$label) + 0.5,
                fill = as.character(as.numeric(.data$label) %% 2)
            )
        ) +
        scale_fill_manual(
            values = c("white", "#eeeeee"),
            guide = FALSE
        ) +
        geom_rect(
            mapping = aes(
                ymin = -Inf,
                ymax = Inf,
                xmin = (seq_along(resultsList) - 1) * statRange - statRange / 2 + gap / 2,
                xmax = (seq_along(resultsList) - 1) * statRange + statRange / 2 - gap / 2
            ),
            fill = "transparent",
            color = "black"
        ) +
        geom_vline(
            mapping = aes(
                xintercept = x1_breaks
            ),
            color = "#cccccc",
            size = 0.5,
            linetype = "dashed"
        ) +
        geom_vline(
            mapping = aes(
                xintercept = x1_breaks[x1_labels == 0]
            ),
            color = "#D00000",
            size = 0.5,
            linetype = "dashed"
        ) +
        geom_errorbarh(data = plotData,
                       aes(
                           y = as.numeric(.data$label),
                           xmin = .data$minShifted,
                           xmax = .data$maxShifted),
                       height = 0.2) +
        geom_point(data = plotData,
                   aes(
                       x = .data$normalizedScoreShifted,
                       y = as.numeric(.data$label)),
                   color = "red") +
        scale_y_discrete(limits = unique(plotData$label), name = "") +
        scale_x_continuous(breaks = x1_breaks,
                           labels = x1_labels,
                           name = "Normalized Score",
                           sec.axis = sec_axis(~. / 1, breaks = x2_breaks, labels = x2_labels)) +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        )
}