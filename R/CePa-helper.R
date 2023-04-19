#' @title Get KEGG pathway catalouge (network) for CePa.ORA and CePa.GSA methods
#' @description Get KEGG pathway catalouge for CePa.ORA and CePa.GSA methods
#' @param org The organism abbreviation. E.g, hsa, mmu, dme, etc.
#' Visit https://www.genome.jp/kegg/catalog/org_list.html to see the full list of supported organisms.
#' @return A named list of KEGG pathways catalouge for CePa.ORA and CePa.GSA methods.
#' @examples
#' @examples
#' \dontrun{
#' #Loading CPA library
#' library(RCPA)
#' genesets <- getCePaPathwayCatalogue("hsa")
#' }
#' @export
#' @importFrom ROntoTools keggPathwayGraphs
#' @importFrom CePa set.pathway.catalogue
#' @importFrom graph nodes edges
#' @importFrom dplyr %>% filter select group_by group_split
#' @importFrom stringr str_split
getCePaPathwayCatalogue <- function(org = "hsa"){
  keggPathway <- keggPathwayGraphs(org, updateCache = T)

  interactionList <- lapply(keggPathway, function(pathway){
    pathway@edgeData %>% names()
  }) %>% do.call(what=c) %>% unique()

  pathList <- lapply(keggPathway, function(pathway){
    (interactionList %in% (pathway@edgeData %>% names())) %>% which() %>% as.character()
  })

  interactionList <- seq(length(interactionList)) %>% as.character() %>% cbind(
    interactionList %>% strsplit('\\|') %>% do.call(what=rbind)
  ) %>% data.frame(stringsAsFactors = FALSE)

  colnames(interactionList) <- c("interaction.id", "input", "output")
  rownames(interactionList) <- interactionList$interaction.id

  interactionList_modified <- interactionList
  interactionList_modified$input <- interactionList$input %>% strsplit(":") %>% sapply(function(x) x[2]) %>% as.character()
  interactionList_modified$output <- interactionList$output %>% strsplit(":") %>% sapply(function(x) x[2]) %>% as.character()

  mapping <- lapply(keggPathway, function(pathway) pathway@nodes) %>% unlist() %>% unique()
  mapping_modified <- mapping %>% strsplit(":") %>% sapply(function(x) x[2]) %>% as.character()
  mapping_modified <- data.frame(node.id = mapping_modified, symbol = mapping_modified, stringsAsFactors = FALSE)

  cat <- set.pathway.catalogue(pathList, interactionList_modified, mapping_modified, min.node = 2, max.node = 1e+6)

  return(cat)
}