#' @title Plot pathways heatmap plot from pathway/geneset/meta analysis results
#' @description pathways heatmap forest plot from pathway/geneset/meta analysis results.
#' @param resultsList A named list of dataframes from pathway analysis, geneset analysis, and/or meta analysis results.
#' The columns are ID, name, description, p.value, pFDR, size, nDE, score and normalizedScore.
#' @param yAxis The column to use for the y-axis.
#' @param negLog10pValueLims A vector of length 2 specifying the minimum and maximum -log10(p-value) to plot.
#' @return A ggplot2 object for presenting the heatmap of the pathways.
#' @examples
#' #TODO add example
#' @importFrom ggplot2 ggplot aes geom_point geom_hline theme_minimal theme theme_bw geom_vline scale_color_gradient scale_size_continuous labs scale_fill_continuous scale_size geom_tile
#' @importFrom ggrepel geom_label_repel
#' @importFrom dplyr %>%
#' @export
plotPathwayHeatmap <- function(resultsList, yAxis = c("ID", "name"), negLog10pValueLims = c(0, 5)){

  yAxis = match.arg(yAxis)

  studyIDs <- names(resultsList)

  if(any(sapply(studyIDs, is.null))){
    stop("The names of the input list should be specified.")
  }

  cols_list <- lapply(resultsList, function(data) colnames(data))

  if(!all(sapply(cols_list, function(x) c("ID", "name", "normalizedScore", "p.value") %in% x))){
    stop("All dataframes in the input list must have 'ID', 'name', 'normalizedScore', and 'p.value' columns.")
  }

  rows_list <- lapply(resultsList, function(data) data$ID %>% unlist() %>% as.vector())

  if(!all(lengths(rows_list) == length(rows_list[[1]]))){
    stop("All dataframes in the input list must have the same number of rows.")
  }

  initial_names <- resultsList[[1]] %>% .$ID %>% unlist() %>% as.vector() %>% sort()

  if(!all(sapply(rows_list, function(x) all.equal(sort(x), initial_names)))){
    stop("All dataframes in the input list must have the same set of pathways.")
  }

  plotData <- lapply(1:length(resultsList), function (i){

    data <- resultsList[[i]]
    data$dataset <- studyIDs[i]
    data$Direction <- ifelse(data$normalizedScore <= 0, "Down", "Up")
    data$abs.normalizedScore <- data$normalizedScore %>% abs()
    data %>% data.frame(stringsAsFactors = FALSE)

  }) %>% do.call(what = rbind)

  scaleMinMax <- function(x, minx, maxx) {
    x[x < minx] <- minx
    x[x > maxx] <- maxx
    x
  }

  if(yAxis == "ID"){
    plotData$yLabel <- plotData$ID
  }else{
    plotData$yLabel <- plotData$name
  }

  plotData$p.value.scaled <- plotData$p.value %>% log10() %>% abs() %>% scaleMinMax(negLog10pValueLims[1], negLog10pValueLims[2])

  ggplot(plotData, aes(y = factor(yLabel),  x = factor(dataset))) +
    labs(
      size = "Normalized score",
      fill = "-Log10 P-value",
      x = "",
      y = "") +
    geom_tile(
      aes(fill = p.value.scaled)
    ) +
    scale_fill_continuous(
      low = "#F1F1F1",
      high = "#68BC36",
      limits = c(negLog10pValueLims[1], negLog10pValueLims[2]),
      breaks = c(negLog10pValueLims[1], (negLog10pValueLims[1] + negLog10pValueLims[2])/2, negLog10pValueLims[2])
    ) +
    geom_point(aes(colour = Direction, size = abs.normalizedScore)) +
      scale_color_manual(values = c("Up" = "#f77a65", "Down" = "#316b9d"),
                         guide = guide_legend(override.aes = list(shape = 15, size = 8))) +

      theme_bw()
}