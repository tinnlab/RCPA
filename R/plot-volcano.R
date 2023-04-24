#' @title Plot volcano plot from Pathway analysis results
#' @description Plot volcano plot from Pathway analysis results
#' @param results A data frame with Pathway analysis results.
#' The columns are ID, name, description, p.value, pFDR, size, nDE, score and normalizedScore.
#' @param xAxis The column to use for the x-axis.
#' @param yAxis The column to use for the y-axis.
#' @param pThreshold The p-value threshold to use for the horizontal line.
#' @param label The column to use for the labels. Default is "name".
#' @param IDsToLabel A vector of IDs to label.
#' When NULL, the top pathways are labeled. Default is NULL.
#' @param topToLabel The number of top pathways to label when IDsToLabels is NULL.
#' @importFrom ggplot2 ggplot aes geom_point geom_hline theme_minimal theme theme_bw geom_vline scale_color_gradient scale_size_continuous labs
#' @importFrom ggrepel geom_label_repel
#' @importFrom dplyr %>% pull
#' @importFrom utils head
#' @export
plotVolcanoPathway <- function(results, xAxis = c("normalizedScore", "score"), yAxis = c("-log10(pFDR)", "-log10(p.value)"), pThreshold = 0.05, label = "name", IDsToLabel = NULL, topToLabel = 20) {

    xAxis <- match.arg(xAxis)
    yAxis <- match.arg(yAxis)

    if (!label %in% colnames(results)) {
        stop("The label column is not in the results data frame.")
    }

    if (!xAxis %in% colnames(results)) {
        stop("The x-axis column is not in the results data frame.")
    }

    if (yAxis == "-log10(pFDR)" && !("pFDR" %in% colnames(results))) {
        stop("The y-axis column is not in the results data frame.")
    }

    if (yAxis == "-log10(p.value)" && !("p.value" %in% colnames(results))) {
        stop("The y-axis column is not in the results data frame.")
    }

    if (is.null(IDsToLabel)) {
        IDsToLabel <- results %>%
            arrange(.data$pFDR) %>%
            head(topToLabel) %>%
            pull(.data$ID)
    }

    plotDat <- data.frame(
        x = results[[xAxis]],
        y = if (yAxis == "-log10(pFDR)") {
            -log10(results$pFDR)
        } else {
            -log10(results$p.value)
        },
        size = results$size,
        label = ifelse(results$ID %in% IDsToLabel, results[[label]], "")
    )

    ggplot(plotDat, aes(x = .data$x, y = .data$y, color = .data$x)) +
        geom_hline(yintercept = -log10(pThreshold), linetype = "dashed", color = "red") +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_point(
            aes(size = .data$size)
        ) +
        geom_label_repel(
            aes(label = .data$label),
            # point.padding = 0.5,
            # segment.alpha = 0.5,
            # force = 0.5,
            # nudge_x = -0.25,
            # nudge_y = 0.01,
            color = "black",

        ) +
        scale_color_gradient(low = "blue", high = "red") +
        scale_size_continuous(range = c(1, 10)) +
        labs(
            x = if (xAxis == "normalizedScore") {
                "Normalized score"
            } else {
                "Score"
            },
            y = if (yAxis == "-log10(pFDR)") {
                "-log10 pFDR"
            } else {
                "-log10 p-value"
            }
        ) +
        theme_bw() +
        theme(
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line.x = element_line(color = "darkgray"),
            axis.line.y = element_line(color = "darkgray"),
            legend.position = "none"
        )
}

#' @title Plot volcano plot from Pathway analysis results
#' @description Plot volcano plot from Pathway analysis results
#' @param results A data frame with Pathway analysis results.
#' The columns are ID, name, description, p.value, pFDR, size, nDE, score and normalizedScore.
#' @param xAxis The column to use for the x-axis.
#' @param yAxis The column to use for the y-axis.
#' @param pThreshold The p-value threshold to use for the horizontal line.
#' @param statThreshold The absolute value of the statistic threshold to use for the vertical line.
#' @importFrom ggplot2 ggplot aes geom_point geom_hline theme_minimal theme theme_bw geom_vline scale_color_gradient scale_size_continuous labs
#' @importFrom ggrepel geom_label_repel
#' @importFrom dplyr %>% pull
#' @importFrom utils head
#' @export
plotVolcanoDE <- function(DEResult, xAxis = "logFC", yAxis = c("-log10(pFDR)", "-log10(p.value)"), pThreshold = 0.05, statThreshold = 2) {

    yAxis <- match.arg(yAxis)

    if (!xAxis %in% colnames(DEResult)) {
        stop("The x-axis column is not in the results data frame.")
    }

    if (yAxis == "-log10(pFDR)" && !("pFDR" %in% colnames(DEResult))) {
        stop("The y-axis column is not in the results data frame.")
    }

    if (yAxis == "-log10(p.value)" && !("p.value" %in% colnames(DEResult))) {
        stop("The y-axis column is not in the results data frame.")
    }

    pvalues <- if (yAxis == "-log10(pFDR)") {
        DEResult$pFDR
    } else {
        DEResult$p.value
    }

    plotDat <- data.frame(
        x = DEResult[[xAxis]],
        y = -log10(pvalues),
        color = ifelse(abs(DEResult[[xAxis]]) > statThreshold & pvalues < pThreshold,  DEResult[[xAxis]], NA)
    )

    if (sum(is.na(plotDat$color)) == nrow(plotDat)) {
        plotDat$color <- plotDat$x
    }

    ggplot(plotDat, aes(x = .data$x, y = .data$y, color = .data$color)) +
        geom_hline(yintercept = -log10(pThreshold), linetype = "dashed", color = "black") +
        geom_vline(xintercept = -statThreshold, linetype = "dashed") +
        geom_vline(xintercept = statThreshold, linetype = "dashed") +
        geom_point() +
        labs(
            x = xAxis,
            y = if (yAxis == "-log10(pFDR)") {
                "-log10 pFDR"
            } else {
                "-log10 p-value"
            }
        ) +
        theme_bw() +
        theme(
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line.x = element_line(color = "darkgray"),
            axis.line.y = element_line(color = "darkgray"),
            legend.position = "none"
        ) +
        scale_color_gradient(low = "blue", high = "red", na.value = "gray")

}