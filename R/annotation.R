#' @title Get Entrez annotation
#' @description This function gets Entrez annotation.
#' @param entrezIds A vector of Entrez IDs.
#' @return A data frame with Entrez annotation. The columns are ID, Symbol, Description, OtherDesignations, OtherAliases and Chromosome.
#' @examples
#' \donttest{
#' library(RCPA)
#' getEntrezAnnotation(c("77267466", "77267467"))
#' }
#' @importFrom httr POST content
#' @importFrom dplyr %>%
#' @export
getEntrezAnnotation <- function(entrezIds) {

    .requirePackage("XML")

    url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    query <- list(
        db = "gene",
        id = as.character(entrezIds) %>%
            unique() %>%
            paste0(collapse = ",")
    )
    res <- POST(url, query = query) %>% content(as = "text")
    xml <- XML::xmlParse(res)

    xpath <- "//DocumentSummary"
    XML::xpathApply(xml, xpath, function(xmlDoc) {
        listDat <- XML::xmlToList(xmlDoc)
        data.frame(
            ID = listDat$.attrs[["uid"]],
            Symbol = ifelse(is.null(listDat$Name), NA, listDat$Name),
            Description = ifelse(is.null(listDat$Description), NA, listDat$Description),
            OtherDesignations = ifelse(is.null(listDat$OtherDesignations), NA, listDat$OtherDesignations),
            OtherAliases = ifelse(is.null(listDat$OtherAliases), NA, listDat$OtherAliases),
            Chromosome = ifelse(is.null(listDat$Chromosome), NA, listDat$Chromosome),
            stringsAsFactors = FALSE
        )
    }) %>%
        do.call(what = rbind) %>%
        `rownames<-`(.$ID)
}
