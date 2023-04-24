#' @title Plot venndiagram from multiple DE Analysis results
#' @description Plot a venndiagram from multiple DE Analysis results.
#' @param DEResults A list of data frames with the results of DE analysis.
#' @param pThreshold The p-value threshold to determine if a gene is differentially expressed.
#' @param useFDR Use the FDR adjusted p-value instead of the raw p-value.
#' @param stat The additional statistis column to use for filtering differentially expressed genes.
#' @param statThreshold The absolute value of the statistic threshold to use for filtering differentially expressed genes.
#' Default is 0, which means no filtering.
#' @importFrom ggvenn ggvenn
#' @importFrom ggplot2 scale_fill_gradient theme
#' @importFrom dplyr %>% filter
#' @importFrom scales trans_new
#' @export
plotVennDE <- function(DEResults, pThreshold = 0.05, useFDR = TRUE, stat = "logFC", statThreshold = 0) {

    if (length(DEResults) < 2) {
        stop("The number of DE results must be at least 2.")
    }

    for (DERes in DEResults) {
        if (useFDR && !("pFDR" %in% colnames(DERes))) {
            stop("The FDR adjusted p-value column is not in the results data frame.")
        } else {
            if (!("p.value" %in% colnames(DERes))) {
                stop("The p.value column is not in the results data frame.")
            }
        }

        if (!stat %in% colnames(DERes)) {
            stop("The statistic column is not in the results data frame.")
        }
    }

    plotDat <- lapply(DEResults, function(DERes) {
        data.frame(DERes) %>%
            filter(
                abs(.data[[stat]]) > statThreshold & (
                    if (useFDR) {
                        .data$pFDR < pThreshold
                    } else {
                        .data$p.value < pThreshold
                    }
                )
            ) %>%
            `[[`("ID")
    })

    if (is.null(names(plotDat))) {
        names(plotDat) <- paste0("Dataset ", seq_along(plotDat))
    }

    ggvenn(plotDat,
           fill_color = c(
               "#316b9d",
               # "#fce397",
               # "#99cc83",
               "#f77a65",
               "#a6a1d0",
               "#fea9c4",
               "#74e7bc",
               "#febb73",
               "#1db4db",
               "#ffc5a6",
               "#b6c9fa",
               "#ee5437"),
           stroke_size = 0.5,
           # set_name_size = 4,
           fill_alpha = 0.75
    )
}

