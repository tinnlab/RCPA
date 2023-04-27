#' @title Plot forest plot from pathway/geneset/meta analysis results
#' @description Plot forest plot from pathway/geneset/meta analysis results.
#' @param resultsList A named list of dataframes from pathway analysis, geneset analysis, and/or meta analysis results.
#' The columns are ID, name, description, p.value, pFDR, size, nDE, score and normalizedScore.
#' @param yAxis The column to use for the y-axis.
#' @param statLims A vector of length 2 specifying the limits for score to use in the x axis.
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
plotForest <- function(resultsList, yAxis = c("ID", "name"), statLims = c(-2.5,2.5)) {

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

  if(!is.numeric(statLims[1]) | !is.numeric(statLims[2])){
    stop("statLims should be a numeric vector of length two.")
  }

  if(statLims[1] >= statLims[2]){
    stop("statLims parameter must be a vector of two values with minimum as the first value and maximum as the second value. ")
  }

  sorted_data <- lapply(resultsList, function (data){
    rownames(data) <- data$ID
    data <- data[initial_names,] %>% data.frame(stringsAsFactors = FALSE)
  })

  plotData <- lapply(sorted_data, function(data){

    sd <- abs((data$normalizedScore - ifelse(data$normalizedScore > 0, 1, -1))/qnorm(data$p.value))

    sd[sd > 0.5] <- 0.5

    data$min <- data$normalizedScore - sd*2
    data$max <- data$normalizedScore + sd*2

    data$min[data$min < statLims[1]] <- statLims[1]
    data$max[data$max > statLims[2]] <- statLims[2]

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
    ggplot(pltDat, aes(y = pltDat$yLabel, x = pltDat$normalizedScore, xmin = pltDat$min, xmax = pltDat$max)) +
      geom_point(color = "red") +
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
      geom_vline(xintercept = c(0), colour="grey", linetype = "solid") +
      coord_cartesian(clip = "off", xlim = c(statLims[1],statLims[2]))+
      geom_errorbarh(height=0.2) +
        geom_point(color = "red") +
      theme_bw() +
      theme(
        axis.text.y = if (i == 1) {
          element_text()
        } else {
          element_blank()
        },
        axis.title.y = element_blank(),
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
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(hjust = 0.5),
        axis.text.x = element_text(vjust = -0.75)
      ) +
      labs(
        x = "Normalized score"
      ) +
      ggtitle(names(plotData)[i])
  })

  plts
}