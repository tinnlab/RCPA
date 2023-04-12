#' @title Get KEGG gene sets names
#' @description This function retrieves KEGG gene sets names for a given organism.
#' @param org The organism abbreviation. E.g, hsa, mmu, dme, etc.
#' Visit https://www.genome.jp/kegg/catalog/org_list.html to see the full list of supported organisms.
#' @return A named vector with KEGG gene sets names.
#' @examples
#' getKEGGPathwayNames("hsa")
#' @importFrom stringr str_split
#' @importFrom dplyr %>%
#' @importFrom stats setNames
#' @importFrom utils read.table
.getKEGGPathwayNames <- function(org = "hsa") {
    gsNames <- read.table(paste0("http://rest.kegg.jp/list/pathway/", org), sep = "\t", header = F, stringsAsFactors = F)
    gsNames[, 2] %>%
        str_split(" - ") %>%
        sapply(function(x) paste0(x[1:(length(x) - 1)], collapse = " - ")) %>%
        setNames(paste0("path:", gsNames[, 1]))
}

#' @title Get KEGG gene sets
#' @description This function retrieves KEGG gene sets for a given organism.
#' @param org The organism abbreviation. E.g, hsa, mmu, dme, etc.
#' Visit https://www.genome.jp/kegg/catalog/org_list.html to see the full list of supported organisms.
#' @return A named list with three elements: database, genesets and names.
#' @examples
#' getKEGGGeneSets("hsa")
#' @importFrom stringr str_split
#' @importFrom dplyr %>% left_join group_by group_split mutate select
#' @importFrom tidyr drop_na
#' @importFrom stats setNames
.getKEGGGeneSets <- function(org = "hsa") {
    geneLink <- read.table(paste0("https://rest.kegg.jp/link/", org, "/pathway"), sep = "\t", header = F, stringsAsFactors = F) %>%
        `colnames<-`(c("geneset", "gene"))

    gensets <- geneLink %>%
        group_by(.$geneset) %>%
        group_split() %>%
        lapply(function(df) {
            list(
                name = df$geneset[1] %>% as.character(),
                gene = df$gene %>% as.character()
            )
        }) %>%
        `names<-`(lapply(., function(gl) gl$name)) %>%
        lapply(function(gl) gl$gene)

    list(
        database = "KEGG",
        genesets = gensets,
        names = .getKEGGPathwayNames(org)[names(gensets)]
    )
}

#' @title Get GO terms names
#' @description This function retrieves GO terms for a given organism.
#' @param namespace The namespace of the GO terms. E.g, biological_process, molecular_function, cellular_component.
#' @return A named vector with GO terms names.
#' @importFrom stringr str_split_1 str_starts str_detect str_match
#' @importFrom dplyr filter
#' @importFrom stats setNames
.getGOTermNames <- function(namespace = c("biological_process", "molecular_function", "cellular_component")) {
    namespace <- match.arg(namespace)

    con <- url("http://purl.obolibrary.org/obo/go.obo")
    rawText <- readLines(con) %>% paste(collapse = "\n")
    close(con)

    terms <- str_split_1(rawText, "\n\n") %>%
        `[`(which(sapply(., function(x) str_starts(x, "\\[Term\\]") && !str_detect(x, "is_obsolete")))) %>%
        lapply(function(x) {
            data.frame(
                id = str_match(x, "\nid: (.*)")[, 2],
                name = str_match(x, "\nname: (.*)")[, 2],
                ns = str_match(x, "\nnamespace: (.*)")[, 2],
                stringsAsFactors = F
            )
        }) %>%
        do.call(what = rbind) %>%
        filter(.$ns == namespace)

    terms$name %>% setNames(terms$id)
}

#' @title Get GO terms
#' @description This function retrieves GO terms for a given organism.
#' @param taxid The NCBI taxonomy ID of the organism.
#' @param namespace The namespace of the GO terms. E.g, biological_process, molecular_function, cellular_component.
#' @return A named list with three elements: database, genesets and names.
#' @importFrom utils read.table
#' @importFrom dplyr %>% filter group_by group_split rename select
#' @importFrom stats setNames
.getGOTerms <- function(taxid = 9606, namespace = c("biological_process", "molecular_function", "cellular_component")) {

    namespace <- match.arg(namespace)

    con <- gzcon(url("https://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz"))
    txt <- readLines(con)
    close(con)

    cat <- switch(namespace,
                  biological_process = "Process",
                  molecular_function = "Function",
                  cellular_component = "Component")

    genesets <- read.table(textConnection(txt), sep = "\t", header = T, stringsAsFactors = F, fill = TRUE, comment.char = "!") %>%
        rename(tax_id = "X.tax_id") %>%
        filter(.$tax_id == taxid & .$Category == cat) %>%
        group_by(.$GO_ID) %>%
        group_split() %>%
        lapply(function(df) {
            list(
                name = df$GO_ID[1] %>% as.character(),
                geneIds = df$GeneID %>% as.character()
            )
        }) %>%
        setNames(lapply(., function(gl) gl$name)) %>%
        lapply(function(gl) gl$geneIds)

    if (length(genesets) == 0) {
        stop("No GO terms found for the given organism and namespace.")
    }

     goTermNames <- .getGOTermNames(namespace)

    list(
        database = "GO",
        genesets = genesets,
        names = goTermNames[names(genesets)]
    )
}

#' @title Get gene sets
#' @description This function retrieves gene sets for a given organism.
#' @param database The database of the gene sets. E.g, KEGG, GO.
#' @param org The organism abbreviation. E.g, hsa, mmu, dme, etc.
#' Visit https://www.genome.jp/kegg/catalog/org_list.html to see the full list of supported organisms.
#' This parameter is only used when database is KEGG.
#' @param taxid The NCBI taxonomy ID of the organism.
#' This parameter is only used when database is GO.
#' @param namespace The namespace of the GO terms. E.g, biological_process, molecular_function, cellular_component.
#' This parameter is only used when database is GO.
#' @return A named list with three elements: database, genesets and names.
#' @examples
#' getGeneSets("KEGG", "abcd")
#' getGeneSets("GO", taxid = 9606, namespace = "biological_process")
#' @export
getGeneSets <- function(database = c("KEGG", "GO"), org = "hsa", taxid = 9606, namespace = c("biological_process", "molecular_function", "cellular_component"), minSize = 10, maxSize = 1000) {
    database <- match.arg(database)

    if (database == "KEGG") {
        if (is.null(org)) {
            stop("Organism must be specified")
        }

        gs <- .getKEGGGeneSets(org)
    } else if (database == "GO") {

        if (is.null(taxid)) {
            stop("Taxonomy ID must be specified")
        }

        namespace <- match.arg(namespace)

        gs <- .getGOTerms(taxid, namespace)
    } else {
        stop("Database not supported")
    }

    gsLength <- sapply(gs$genesets, length)
    keeps <- gsLength >= minSize & gsLength <= maxSize
    gs$genesets <- gs$genesets[keeps]
    gs$names <- gs$names[keeps]

    gs
}
