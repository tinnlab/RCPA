#' @title Plot gene heatmap from a SummarizedExperiment object with DE analysis results
#' @description Plot gene heatmap from a SummarizedExperiment object with DE analysis results.
#' The heatmap contains p-values and log fold changes from the DE analysis.
#' @param summarizedExperiment SummarizedExperiment object
#' @param genes A vector of gene id (e.g. Entrez IDs) to plot.
#' The genes should be in the rownames of the SummarizedExperiment object.
#' @param labels A vector of labels for the genes. If not provided, the gene IDs will be used as labels.
#' @param useFDR If TRUE, use FDR adjusted p-values. Otherwise, use raw p-values.
#' @param logFCLims A vector of length 2 specifying the minimum and maximum log fold change to plot.
#' @param negLog10pValueLims A vector of length 2 specifying the minimum and maximum -log10(p-value) to plot.
#' @return A heatmap of the genes from ggplot2.
#' @examples
#' #TODO add example
#' @importFrom SummarizedExperiment rowData
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient2 theme_minimal theme geom_tile coord_flip scale_x_discrete scale_y_discrete guide_legend element_blank
#' @importFrom dplyr %>% select mutate
#' @importFrom tidyr gather
#' @importFrom ggthemes scale_fill_gradient_tableau scale_fill_gradient2_tableau
#' @importFrom ggnewscale new_scale_fill
#' @export
plotDEGeneHeatmap <- function(summarizedExperiment, genes, useFDR = TRUE, labels = NULL, logFCLims = c(-5, 5), negLog10pValueLims = c(0, 5)){
    commonGenes <- intersect(genes, rownames(summarizedExperiment))

    if (length(commonGenes) == 0) {
        stop("No common genes found between the input genes and the genes in the SummarizedExperiment object")
    }

    DEdf <- as.data.frame(rowData(summarizedExperiment))[commonGenes,]

    if (is.null(labels)) {
        labels <- commonGenes
    } else {
        names(labels) <- genes
        labels <- labels[row.names(DEdf)]
    }

    scaleMinMax <- function(x, minx, maxx) {
        x[x < minx] <- minx
        x[x > maxx] <- maxx
        x
    }

    plotData <- DEdf %>%
        mutate(
            p.value = ifelse(rep(useFDR, nrow(DEdf)), .$pFDR, .$p.value) %>% log10() %>% abs() %>% scaleMinMax(negLog10pValueLims[1], negLog10pValueLims[2]),
            logFC = scaleMinMax(.$logFC, logFCLims[1], logFCLims[2]),
            label = factor(labels, levels = labels)
        ) %>%
        select("label", "logFC", "p.value") %>%
        gather("type", "value", -label)

    ggplot() +
        geom_tile(data = plotData[plotData$type == "logFC", ], aes(x = label, y = type, fill = value)) +
        scale_fill_gradient2_tableau(
            palette =  "Red-Blue Diverging",
            na.value = "white",
            limits = logFCLims,
            guide = guide_legend(title = "log2 FC")
        ) +
        new_scale_fill() +
        geom_tile(data =  plotData[plotData$type == "p.value", ], aes(x = label, y = type, fill = value)) +
        scale_fill_gradient_tableau(
            palette = "Red",
            na.value = "white",
            limits = negLog10pValueLims,
            guide = guide_legend(title = "-log10 p-value")
        ) +
        theme_minimal() +
        coord_flip() +
        theme(
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        scale_x_discrete(
            labels = labels,
            expand = c(0, 0)
        ) +
        scale_y_discrete(
            labels = c(
                "log2 FC",
                paste0("-log10", ifelse(useFDR, " pFDR", " p-value"))
            ),
            expand = c(0, 0)
        )
}