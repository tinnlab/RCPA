#' @title Plot a bar chart of the pathway analysis results
#' @description This function plots a bar chart of the pathway analysis results.
#' @param results A data frame with pathway analysis results.
#' The columns are ID, name, description, p.value, pFDR, size, nDE, score and normalizedScore.
#' @param limit The maximum number of pathways to plot.
#' @param label The column to use for the labels.
#' @param by The column to use for the bar heights.
#' @param maxNegLog10PValue The maximum -log10(p-value) to plot.
#' @param pThreshold The p-value threshold to use for significance.
#' @param useFDR If TRUE, use FDR adjusted p-values. Otherwise, use raw p-values.
#' @importFrom ggplot2 ggplot aes geom_bar theme_minimal theme geom_text coord_flip scale_x_discrete
#' @importFrom ggplot2 scale_y_continuous guide_legend element_blank scale_fill_manual element_line labs geom_hline
#' @importFrom dplyr %>% select mutate arrange
#' @return A ggplot2 object.
#' @examples
#' #TODO add example
#' @export
plotBarChart <- function(results, limit = 50, label = "name", by = c("pFDR", "p.value", "score", "normalizedScore"), maxNegLog10PValue = 5, pThreshold = 0.05, useFDR = TRUE) {

    by <- match.arg(by)

    if (!by %in% colnames(results)) {
        stop(paste0("The column '", by, "' does not exist in the results."))
    }

    if (useFDR && !("pFDR" %in% colnames(results))) {
        stop("The column 'pFDR' does not exist in the results.")
    }

    if (!label %in% colnames(results)) {
        stop(paste0("The column '", label, "' does not exist in the results."))
    }

    plotData <- results %>%
        mutate(
            isSignificant = ifelse(ifelse(rep(useFDR, nrow(.)), .data$pFDR, .data$p.value) <= pThreshold, "yes", "no"),
        ) %>%
        mutate(
            p.value = pmin(-log10(.data$p.value), maxNegLog10PValue),
            pFDR = pmin( -log10(.data$pFDR), maxNegLog10PValue)
        ) %>%
        arrange(desc(abs(.data[[by]]))) %>%
        head(limit) %>%
        arrange(.data[[by]]) %>%
        mutate(
            ID = factor(.data$ID, levels = .data$ID),
        ) %>%
        select(.data$ID, sym(label), sym(by), .data$isSignificant)

    pl <- ggplot(plotData, aes(x = .data$ID, y = .data[[by]], fill= .data$isSignificant)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line.x = element_line(color = "darkgray"),
            axis.line.y = element_line(color = "darkgray")
        ) +
        scale_x_discrete(labels = plotData[[label]], expand = c(-1,-1)) +
        scale_y_continuous(expand = c(0,0)) +
        scale_fill_manual(
            values = c("yes" = "steelblue", "no" = "grey"),
            guide = guide_legend(title = "Significant")
        ) +
        labs(x= element_blank())

    if (by == "p.value" | by == "pFDR"){
        pl <- pl +
            geom_hline(yintercept = -log10(pThreshold), linetype = "dashed", color = "red") +
            labs(y = paste0("-log10 ", by))
    }

    if (by == "score"){
        pl <- pl + labs(y = "Score")
    }

    if (by == "normalizedScore"){
        pl <- pl + labs(y = "Normalized score")
    }

    pl
}