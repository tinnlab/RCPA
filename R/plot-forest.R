#' @title Plot forest plot from pathway/geneset/meta analysis results
#' @description Plot forest plot from pathway/geneset/meta analysis results.
#' @param resultsList A named list of dataframes from pathway analysis, geneset analysis, and/or meta analysis results.
#' The columns are ID, name, description, p.value, pFDR, size, nDE, score and normalizedScore.
#' @param yAxis The column to use for the y-axis.
#' @param axisFontSize The font size of axis labels.
#' @param titleFontSize The font size for plot title.
#' @param pointSize The red point size in the plot.
#' @return A grid of one or more forest plots from grid.arrange.
#' @examples
#' #TODO add example
#' @importFrom ggplot2 ggplot aes geom_point geom_hline theme_minimal theme theme_bw geom_vline scale_color_gradient scale_size_continuous labs geom_rect coord_cartesian geom_errorbarh ggtitle
#' @importFrom gridExtra grid.arrange
#' @importFrom dplyr %>%
#' @export
plotForest <- function(resultsList, yAxis = c("ID", "name"), axisFontSize = 10, titleFontSize = 12, pointSize = 7) {

  yAxis = match.arg(yAxis)

  studyIDs <- names(resultsList)
  if(any(sapply(studyIDs, is.null))){
    stop("The names of the input list should be specified.")
  }

  cols_list <- lapply(resultsList, function(data) colnames(data))

  if(!all(sapply(cols_list, function(x) c("ID", "name", "normalizedScore") %in% x))){
    stop("All dataframes in the input list must have 'ID', 'name', 'normalizedScore' columns.")
  }

  rows_list <- lapply(resultsList, function(data) data$ID %>% unlist() %>% as.vector())

  if(!all(lengths(rows_list) == length(rows_list[[1]]))){
    stop("All dataframes in the input list must have the same number of rows.")
  }

  initial_names <- resultsList[[1]] %>% .$ID %>% unlist() %>% as.vector() %>% sort()

  if(!all(sapply(rows_list, function(x) all.equal(sort(x), initial_names)))){
    stop("All dataframes in the input list must have the same set of pathways.")
  }

  sorted_data <- lapply(resultsList, function (data){
    data <- data[initial_names,] %>% data.frame(stringsAsFactors = FALSE)
  })

  plotData <- lapply(sorted_data, function(data){

    sd <- data$normalizedScore

    sd[sd > 0.5] <- 0.5

    data$min <- data$normalizedScore - sd*2
    data$max <- data$normalizedScore + sd*2

    data$min[data$min < -2.5] <- -2.5
    data$max[data$max > 2.5] <- 2.5

    if(yAxis == "ID")
      data$yLabel <- data$ID
    else
      data$yLabel <- data$name

    data$yLabel <- data$yLabel %>% as.factor()

    data %>% data.frame(stringsAsFactors = FALSE)
  })

  names(plotData) <- studyIDs

  plts <- lapply(1:length(plotData), function(i){
    pltDat <- plotData[[i]]
    ggplot(pltDat, aes(y = yLabel, x = normalizedScore, xmin = min, xmax = max)) +
      geom_point(size = pointSize, color = "red") +
      geom_rect(
        aes(
          xmin = -Inf,
          xmax = Inf,
          ymin = as.numeric(pltDat$yLabel)-0.5,
          ymax = as.numeric(pltDat$yLabel)+0.5
        ),
        fill = ifelse((as.numeric(pltDat$yLabel) %% 2 == 0), "white", "#eeeeee"),
        color = "white"
      ) +
      geom_point(size = pointSize, color = "red") +
      geom_vline(xintercept = c(-1,1), colour="#FA8072", linetype = "longdash") +
      geom_vline(xintercept = c(0), colour="grey", linetype = "solid") +
      coord_cartesian(clip = "off", xlim = c(-2.5,2.5))+
      geom_errorbarh(height=.2, linewidth=1) +
      theme_bw() +
      theme(
        axis.text.y = if (i == 1) {
          element_text(size = axisFontSize, face = "bold")
        } else {
          element_blank()
        },
        axis.title.y = if (i == 1) {
          element_text(size = titleFontSize, face = "bold")
        } else {
          element_blank()
        },
        axis.ticks.y = if (i == 1) {
          NULL
        } else {
          element_blank()
        },
        plot.margin = unit(c(5,5,7,5), "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        text = element_text(size = axisFontSize, face = "bold"),
        plot.title = element_text(size = titleFontSize, face = "bold", hjust = 0.5),
        axis.text = element_text(size = axisFontSize, face = "bold", hjust = 0.5),
        axis.text.x = element_text(vjust = -0.75)
      ) +
      labs(
        x = "Normalized score",
        y = ifelse(yAxis == "ID", "ID", "name")
      ) +
      ggtitle(names(plotData)[i])
  })

  gridExtra::grid.arrange(grobs = plts,
  nrow = 1,
  widths = c(1.3, rep(1, length(plts) - 1))
  )
}