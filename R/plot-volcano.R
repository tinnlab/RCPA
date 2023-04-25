#' @title Plot volcano plot from Pathway analysis results
#' @description Plot volcano plot from Pathway analysis results
#' @param PAResult A data frame with Pathway analysis results.
#' The columns are ID, name, description, p.value, pFDR, size, nDE, score and normalizedScore.
#' @param xAxis The column to use for the x-axis.
#' @param yAxis The column to use for the y-axis.
#' @param pThreshold The p-value threshold to use for the horizontal line.
#' @param label The column to use for the labels. Default is "name".
#' @param IDsToLabel A vector of IDs to label.
#' When NULL, the top pathways are labeled. Default is NULL.
#' @param topToLabel The number of top pathways to label when IDsToLabels is NULL.
#' @examples
#' \dontrun{
#' library(AnnotationDbi)
#' library(SummarizedExperiment)
#' library(limma)
#' library(RCPA)
#' data("data")
#' # Get affymetrix dataset
#' affyDataset <- data$affyDataset
#' # Create the analysis design
#' affyDesign <- model.matrix(~0 + condition + region, data = colData(affyDataset))
#' affyContrast <- limma::makeContrasts("conditionalzheimer-conditionnormal", levels=affyDesign)
#' # Perform DE analysis affymetrix dataset
#' affyDEExperiment <- RCPA::runDEAnalysis(affyDataset, 
#'                                         method = "limma", 
#'                                         design = affyDesign, 
#'                                         contrast = affyContrast, 
#'                                         annotation = "GPL570")
#' # Get DE data frame
#' DEResult <- rowData(affyDEExperiment)
#' # Plot volcano
#' plotVolcanoDE(DEResult)
#' }
#' 
#' @importFrom ggplot2 ggplot aes geom_point geom_hline theme_minimal theme theme_bw geom_vline scale_color_gradient scale_size_continuous labs
#' @importFrom ggrepel geom_label_repel
#' @importFrom dplyr %>% pull
#' @importFrom utils head
#' @export
plotVolcanoPathway <- function(PAResult, xAxis = c("normalizedScore", "score"), yAxis = c("-log10(pFDR)", "-log10(p.value)"), pThreshold = 0.05, label = "name", IDsToLabel = NULL, topToLabel = 20) {

    xAxis <- match.arg(xAxis)
    yAxis <- match.arg(yAxis)

    if (!label %in% colnames(PAResult)) {
        stop(paste0("The label column '", label, "' is not in the results data frame."))
    }

    if (!xAxis %in% colnames(PAResult)) {
        stop(paste0("The xAxis column '", xAxis, "' is not in the results data frame."))
    }

    if (yAxis == "-log10(pFDR)" && !("pFDR" %in% colnames(PAResult))) {
        stop("The pFDR column is not in the results data frame")
    }

    if (yAxis == "-log10(p.value)" && !("p.value" %in% colnames(PAResult))) {
        stop("The p.value column is not in the results data frame")
    }

    if (!"pathwaySize" %in% colnames(PAResult)) {
        stop("The size column is not in the results data frame")
    }

    if (is.null(IDsToLabel)) {
        IDsToLabel <- PAResult %>%
            arrange(
                if (yAxis == "-log10(pFDR)") {
                    .data$pFDR
                } else {
                    .data$p.value
                }
            ) %>%
            head(topToLabel) %>%
            pull(.data$ID)
    }

    plotDat <- data.frame(
        x = PAResult[[xAxis]],
        y = if (yAxis == "-log10(pFDR)") {
            -log10(PAResult$pFDR)
        } else {
            -log10(PAResult$p.value)
        },
        size = PAResult$pathwaySize,
        label = ifelse(PAResult$ID %in% IDsToLabel, PAResult[[label]], "")
    )

    ggplot(plotDat, aes(x = .data$x, y = .data$y, color = .data$x)) +
        geom_point(
            aes(size = .data$size)
        ) +
        geom_hline(yintercept = -log10(pThreshold), linetype = "dashed", color = "red") +
        geom_vline(xintercept = 0, linetype = "dashed") +
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
#' @param pThreshold The p-value threshold to use for the horizontal line.
#' @param useFDR Whether to use the pFDR column instead of the p.value column.
#' @param logFCThreshold The logFC threshold to use for the vertical line.
#' @importFrom ggplot2 ggplot aes geom_point geom_hline theme_minimal theme theme_bw geom_vline scale_color_gradient scale_size_continuous labs
#' @importFrom ggrepel geom_label_repel
#' @importFrom dplyr %>% pull
#' @importFrom utils head
#' @export
plotVolcanoDE <- function(DEResult, pThreshold = 0.05, useFDR = TRUE, logFCThreshold = 1) {

    if (!"logFC" %in% colnames(DEResult)) {
        stop("The logFC column is not in the results data frame.")
    }

    if (useFDR && !("pFDR" %in% colnames(DEResult))) {
        stop("The pFDR column is not in the results data frame.")
    } else if (!("p.value" %in% colnames(DEResult))) {
        stop("The p.value column is not in the results data frame.")
    }

    pvalues <- if (useFDR) {
        DEResult$pFDR
    } else {
        DEResult$p.value
    }

    plotDat <- data.frame(
        x = DEResult$logFC,
        y = -log10(pvalues),
        color = ifelse(abs(DEResult$logFC) > logFCThreshold & pvalues < pThreshold,  DEResult$logFC, NA)
    )

    isNoSig <- FALSE
    if (sum(is.na(plotDat$color)) == nrow(plotDat)) {
        isNoSig <- TRUE
        plotDat$color <- "gray"
    }

    pl <- ggplot(plotDat, aes(x = .data$x, y = .data$y, color = .data$color)) +
        geom_point() +
        geom_hline(yintercept = -log10(pThreshold), linetype = "dashed", color = "black") +
        geom_vline(xintercept = -logFCThreshold, linetype = "dashed") +
        geom_vline(xintercept = logFCThreshold, linetype = "dashed") +
        labs(
            x = "log2 fold change",
            y = if (useFDR) {
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


    if (!isNoSig) {
        pl <- pl + scale_color_gradient(low = "blue", high = "red", na.value = "gray")
    } else {
        pl <- pl + scale_color_manual(values = "gray")
    }

    pl

}