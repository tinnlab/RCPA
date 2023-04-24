#' @title Plot a bar chart of the pathway analysis results
#' @description This function plots a bar chart of the pathway analysis results.
#' @param results A named list of data frame with pathway analysis results.
#' The columns are ID, name, p.value, pFDR, size, nDE, score and normalizedScore.
#' @param limit The maximum number of pathways to plot.
#' The pathway will be sorted by the average absolute value of the -log10(p-value) or -log10(pFDR) depending on the useFDR parameter.
#' @param label The column to use for the labels.
#' @param by The column to use for the bar heights.
#' @param maxNegLog10PValue The maximum -log10(p-value) to plot.
#' @param pThreshold The p-value threshold to use for significance.
#' @param useFDR If TRUE, use FDR adjusted p-values for the significance threshold. Otherwise, use raw p-values.
#' This parameter is used to mark the color of the bars and is independent of the 'by' parameter.
#' @return A ggplot2 object.
#' @examples
#' \dontrun{
#' # Loading libraries
#' library(ggplot2)
#' library(ggpattern)
#' library(RCPA)
#' # Simulating DE analysis result
#' results <- lapply(1:3, function(i) {
#'   set.seed(i)
#'   
#'   data.frame(
#'     ID = paste0("geneset", 1:100),
#'     name = paste0("Pathway ", 1:100),
#'     description = paste0("Description ", 1:100),
#'     p.value = runif(100) / 10,
#'     pFDR = runif(100) / 5,
#'     size = runif(100, 100, 500),
#'     nDE = runif(100, 10, 100),
#'     score = runif(100, -2, 2),
#'     normalizedScore = runif(100)
#'   )
#' })
#' # Plot the barchart
#' plotBarChart(results)
#' 
#' # Plot barchart with p-value
#' plotBarChart(results, by = "p.value")
#' }
#' @export
#' @importFrom ggplot2 ggplot aes geom_bar theme_minimal theme geom_text coord_flip scale_x_discrete scale_colour_discrete scale_fill_discrete scale_color_manual
#' @importFrom ggplot2 scale_y_continuous guide_legend element_blank scale_fill_manual element_line labs geom_hline
#' @importFrom dplyr %>% select mutate arrange desc
#' @importFrom utils head
#' @importFrom rlang sym
#' @importFrom ggpattern geom_bar_pattern scale_pattern_manual
plotBarChart <- function(results, limit = 10, label = "name", by = c("pFDR", "p.value", "score", "normalizedScore"), maxNegLog10PValue = 5, pThreshold = 0.05, useFDR = TRUE) {

    by <- match.arg(by)

    for (result in results) {

        if (!"ID" %in% colnames(result)) {
            stop("The column 'ID' does not exist in the results.")
        }

        if (!by %in% colnames(result)) {
            stop(paste0("The column '", by, "' does not exist in the results."))
        }

        if (useFDR && !("pFDR" %in% colnames(result))) {
            stop("The column 'pFDR' does not exist in the results.")
        }
        if (!label %in% colnames(result)) {
            stop(paste0("The column '", label, "' does not exist in the results."))
        }
    }

    if (is.null(names(results))) {
        names(results) <- paste0("Dataset ", seq_along(results))
    }

    commonColNames <- Reduce(intersect, lapply(results, colnames))

    pathwayOrder <- lapply(results, function(r){
        r[, c("ID", ifelse(useFDR, "pFDR", "p.value"))]
    }) %>%
        do.call(rbind, .) %>%
        mutate(
            logP = if (useFDR) -log10(.data$pFDR) else -log10(.data$p.value)
        ) %>%
        group_by(.data$ID) %>%
        dplyr::summarize(
            avgLogP = mean(.data$logP, na.rm = TRUE)
        ) %>%
        arrange(desc(.data$avgLogP)) %>%
        head(limit) %>%
        pull(.data$ID)

    plotData <- names(results) %>%
        lapply(function(n) {
            results[[n]][, commonColNames] %>%
                mutate(
                    dataset = n
                )
        }) %>%
        do.call(rbind, .) %>%
        filter(.data$ID %in% pathwayOrder) %>%
        mutate(
            isSignificant = ifelse(ifelse(rep(useFDR, nrow(.)), .data$pFDR, .data$p.value) <= pThreshold, "Yes", "No"),
            p.value = pmin(-log10(.data$p.value), maxNegLog10PValue),
            pFDR = pmin(-log10(.data$pFDR), maxNegLog10PValue),
            ID = factor(.data$ID, levels = pathwayOrder),
            dataset = factor(.data$dataset, levels = names(results))
        ) %>%
        select("ID", sym(label), sym(by), "isSignificant", "dataset")

    pl <- ggplot(plotData, aes(x = .data$ID, y = .data[[by]], fill = .data$dataset, pattern = .data$isSignificant)) +
        geom_bar_pattern(
            stat = "identity",
            position = if (length(results) > 1) "dodge" else "fill",
            width = ifelse(length(results) > 1, 0.9, 1),
            pattern_size = 0,
            pattern_alpha = 0.4,
            pattern_fill = "white",
            pattern_spacing = 0.01
        ) +
        coord_flip() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line.x = element_line(color = "darkgray"),
            axis.line.y = element_line(color = "darkgray")
        ) +
        scale_x_discrete(labels = plotData[[label]], expand = expansion(add = c(0.75, 0))) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_fill_discrete(
            guide = guide_legend(title = "Dataset")
        ) +
        scale_pattern_manual(
            values = c("Yes" = "stripe", "No" = "none"),
            guide = guide_legend(title = "Significant")
        ) +
        labs(x = element_blank())

    if (by == "p.value" | by == "pFDR") {
        pl <- pl +
            geom_hline(yintercept = -log10(pThreshold), linetype = "dashed", color = "red") +
            labs(y = paste0("-log10 ", by))
    }

    if (by == "score") {
        pl <- pl + labs(y = "Score")
    }

    if (by == "normalizedScore") {
        pl <- pl + labs(y = "Normalized score")
    }

    pl
}