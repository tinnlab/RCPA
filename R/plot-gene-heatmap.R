#' @title Plot gene heatmap from a SummarizedExperiment object with DE analysis results
#' @description Plot gene heatmap from a SummarizedExperiment object with DE analysis results.
#' The heatmap contains p-values and log fold changes from the DE analysis.
#' @param DEResults A named list of data frame of DE analysis results.
#' @param genes A vector of gene id (e.g. Entrez IDs) to plot.
#' The genes must be in the ID column of the data frame in DEResults.
#' @param labels A vector of labels for the genes. If not provided, the gene IDs will be used as labels.
#' @param useFDR If TRUE, use FDR adjusted p-values. Otherwise, use raw p-values.
#' @param logFCLims A vector of length 2 specifying the minimum and maximum log fold change to plot.
#' @param negLog10pValueLims A vector of length 2 specifying the minimum and maximum -log10(p-value) to plot.
#' @return A heatmap of the genes from ggplot2.
#' @examples
#' \dontrun{
#' library(hgu133plus2.db)
#' library(AnnotationDbi)
#' library(SummarizedExperiment)
#' library(limma)
#' library(ggplot2)
#' library(RCPA)
#' 
#' # generate diferential expression analysis using a random gene expression matrix
#' DEResults <- lapply(1:3, function(seed) {
#'   set.seed(seed)
#'   exprs <- round(matrix(abs(rnorm(20000 * 10, sd = 4)), nrow = 20000, ncol = 10))
#'   rownames(exprs) <- sample(keys(hgu133plus2.db, keytype = "PROBEID"), nrow(exprs), replace = FALSE)
#'   colnames(exprs) <- paste0("sample", 1:10)
#'   
#'   controlSamples <- paste0("sample", 1:5)
#'   conditionSamples <- paste0("sample", 6:10)
#'   
#'   exprs[, conditionSamples] <- exprs[, conditionSamples] + 
#'   2*sample(c(1,-1), nrow(exprs), replace = TRUE)
#'   
#'   colData <- data.frame(
#'     row.names = colnames(exprs),
#'     group = factor(c(rep("control", length(controlSamples)), rep("condition", 
#'     length(conditionSamples)))),
#'     pair = factor(c(seq_along(controlSamples), seq_along(conditionSamples)))
#'   )
#'   
#'   summarizedExperiment <- SummarizedExperiment(
#'     assays = list(counts = exprs),
#'    colData = colData
#' )
#'  
#'  # control vs condition
#'  design <- model.matrix(~0 + group, data = colData)
#'  contrast <- makeContrasts("groupcondition-groupcontrol", levels = design)
#'  runDEAnalysis(summarizedExperiment,
#'   method = "limma", 
#'   design, 
#'   contrast, 
#'   annotation = "GPL570") %>% rowData()
#'})
#'
#'# Get common genes from different studies
#'genes <- sample(Reduce(intersect, lapply(DEResults, function(res) res$ID)), 30)
#'# Get gene annotation
#'annotation <- getEntrezAnnotation(genes)
#'# Get gene description
#'labels <- annotation[genes, "Description"]
#'# Plot gene heatmap
#'plotDEGeneHeatmap(DEResults, genes, labels = labels)
#' }
#' @importFrom SummarizedExperiment rowData
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient2 theme_minimal theme geom_tile coord_flip scale_x_discrete scale_y_discrete guide_legend element_blank expansion
#' @importFrom dplyr %>% select mutate
#' @importFrom tidyr gather
#' @importFrom ggthemes scale_fill_gradient_tableau scale_fill_gradient2_tableau
#' @importFrom ggnewscale new_scale_fill
#' @export
plotDEGeneHeatmap <- function(DEResults, genes, useFDR = TRUE, labels = NULL, logFCLims = c(-5, 5), negLog10pValueLims = c(0, 5)) {
    commonGenes <- lapply(DEResults, function(x) intersect(genes, x$ID)) %>%
        Reduce(intersect, .)

    if (any(commonGenes == 0)) {
        stop("No common genes found between the input genes and the genes in the DE results")
    }

    if (length(commonGenes) < length(genes)) {
        warning("Some input genes are not found in the DE results")
    }

    DEdfs <- lapply(DEResults, function(x) x[match(commonGenes, x$ID),])

    if (is.null(labels)) {
        labels <- commonGenes
    } else {
        names(labels) <- genes
        labels <- labels[commonGenes]
    }

    scaleMinMax <- function(x, minx, maxx) {
        x[x < minx] <- minx
        x[x > maxx] <- maxx
        x
    }

    if (is.null(names(DEdfs))) {
        names(DEdfs) <- paste0("Dataset ", seq_along(DEdfs))
    }

    plotData <- names(DEdfs) %>%
        lapply(function(n) {
            DEdf <- as.data.frame(DEdfs[[n]])
            DEdf %>%
                mutate(
                    p.value = ifelse(rep(useFDR, nrow(DEdf)), .$pFDR, .$p.value) %>%
                        log10() %>%
                        abs() %>%
                        scaleMinMax(negLog10pValueLims[1], negLog10pValueLims[2]),
                    logFC = scaleMinMax(.$logFC, logFCLims[1], logFCLims[2]),
                    label = factor(labels, levels = labels),
                    dataset = n
                )
        }) %>%
        do.call(what = rbind) %>%
        select("label", "logFC", "p.value", "dataset") %>%
        gather("type", "value", -"label", -"dataset") %>%
        mutate(
            label = factor(.data$label, levels = labels),
            dataset = factor(.data$dataset, levels = names(DEdfs)),
            type = factor(.data$type, levels = c("p.value", "logFC")),
            colOrder = as.numeric(.data$dataset)*2 + as.numeric(.data$type) + as.numeric(.data$dataset)*0.1 + as.numeric(.data$type)*0.01
        )

    uniqueY <- sort(unique(plotData$colOrder))

    ggplot() +
        geom_tile(data = plotData[plotData$type == "logFC",], aes(x = .data$label, y = .data$colOrder, fill = .data$value, width = 1, height = 1)) +
        scale_fill_gradient2(
            high = "#B80F0A",
            low = "#004F98",
            mid = "white",
            na.value = "white",
            limits = logFCLims,
            guide = guide_colorbar(
                title = "log2 FC"
            )
        ) +
        new_scale_fill() +
        geom_tile(data = plotData[plotData$type == "p.value",], aes(x = .data$label, y = .data$colOrder, fill = .data$value, width = 1, height = 1)) +
        scale_fill_gradient(
            low = "white",
            high = "#B80F0A",
            na.value = "white",
            limits = negLog10pValueLims,
            guide = guide_colorbar(title = "-log10 p-value")
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
            expand = expansion(mult = c(0, 1/length(labels)))
        ) +
        scale_y_continuous(
            breaks = uniqueY,
            labels = rep(c(
                paste0("-log10", ifelse(useFDR, " pFDR", " p-value")),
                "log2 FC"
            ), length(DEdfs)),
            expand = c(0, 0)
        ) +
        annotate(
            "text",
            x = length(labels) + 0.5,
            y = sapply(seq_along(DEdfs), function(i) mean(uniqueY[(i-1)*2 + 1:2])),
            label = names(DEdfs),
            hjust = 0.5,
            vjust = -0.5
        )
}