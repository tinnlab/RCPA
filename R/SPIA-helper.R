#' @title Get KEGG pathway network for SPIA method
#' @description Get KEGG pathway network for SPIA method
#' @param org The organism abbreviation. E.g, hsa, mmu, dme, etc.
#' Visit https://www.genome.jp/kegg/catalog/org_list.html to see the full list of supported organisms.
#' @return A named list of KEGG pathway networks for SPIA method.
#' @examples
#' getSPIAKEGGNetwork("hsa")
#' @export
#' @importFrom ROntoTools keggPathwayGraphs keggPathwayNames
#' @importFrom graph nodes edges
#' @importFrom dplyr %>% filter select group_by group_split
#' @importFrom tidyr spread
#' @importFrom stringr str_split str_length
getSPIAKEGGNetwork <- function(org = "hsa", updateCache = F) {

    keggPathway <- ROntoTools::keggPathwayGraphs(org, updateCache = updateCache)
    pathwayNames <- ROntoTools::keggPathwayNames(org)

    relationships <- c("activation", "compound", "binding/association",
                       "expression", "inhibition", "activation_phosphorylation",
                       "phosphorylation", "inhibition_phosphorylation",
                       "inhibition_dephosphorylation", "dissociation", "dephosphorylation",
                       "activation_dephosphorylation", "state change", "activation_indirect effect",
                       "inhibition_ubiquination", "ubiquination", "expression_indirect effect",
                       "inhibition_indirect effect", "repression", "dissociation_phosphorylation",
                       "indirect effect_phosphorylation", "activation_binding/association",
                       "indirect effect", "activation_compound", "activation_ubiquination")

    replacements <- list(
        c('ubiquitination', 'ubiquination'),
        c(',missing interaction', ''),
        c('missing interaction', ''),
        c('compound,activation', 'activation,compound')
    )

    keggPathwayEntrez <- lapply(keggPathway, function(g) {
        graph::nodes(g) <- graph::nodes(g) %>%
            strsplit(":") %>%
            sapply(function(x) x[2]) %>%
            as.character()
        g
    })

    allRelationships <- lapply(names(keggPathwayEntrez), function(pathway) {
        rels <- keggPathwayEntrez[[pathway]]@edgeData@data %>% sapply(function(s) s$subtype)

        fromTo <- names(rels) %>%
            strsplit("\\|") %>%
            do.call(what = rbind) %>%
            as.data.frame() %>%
            `colnames<-`(c("to", "from"))

        data.frame(
            pathway = pathway,
            relationship = rels,
            from = fromTo$from,
            to = fromTo$to
        )
    }) %>%
        do.call(what = rbind) %>%
        filter(.data$relationship %in% relationships)

    for (r in replacements) {
        allRelationships$relationship <- sub(r[1], r[2], allRelationships$relationship)
    }

    allRelationships$relationship <- sub(pattern = ",", replacement = "_", allRelationships$relationship)

    allRelationships <- allRelationships %>% filter(str_length(.data$relationship) > 0)

    uniqueRelationships <- unique(allRelationships$relationship)

    pathInfo <- allRelationships %>%
        group_by(.data$pathway) %>%
        group_split() %>%
        lapply(function(df) {
            nodes <- graph::nodes(keggPathwayEntrez[[df$pathway[1]]])

            rels <- df %>%
                group_by(.data$relationship) %>%
                group_split() %>%
                lapply(function(dfRel) {
                    reactions <- dfRel %>%
                        select(.data$from, .data$to) %>%
                        tidyr::spread(key = .data$from, value = .data$from)

                    dat <- matrix(NA, nrow = length(nodes), ncol = length(nodes), dimnames = list(nodes, nodes))

                    from <- colnames(reactions)[-1]
                    to <- reactions$to
                    dat[to, from] <- as.numeric(as.matrix(reactions[, -1]))
                    dat[!is.na(dat)] <- 1
                    dat[is.na(dat)] <- 0
                    list(
                        relationship = dfRel$relationship[1],
                        dat = t(dat)
                    )
                }) %>%
                setNames(sapply(., function(x) x$relationship)) %>%
                lapply(function(x) x$dat)

            remainingRels <- setdiff(uniqueRelationships, names(rels))

            for (r in remainingRels) {
                rels[[r]] <- matrix(0, nrow = length(nodes), ncol = length(nodes), dimnames = list(nodes, nodes))
            }

            rels$dissociation_phosphorylation <- matrix(0, nrow = length(nodes), ncol = length(nodes), dimnames = list(nodes, nodes))
            rels <- rels[relationships]
            rels$nodes <- nodes
            rels$NumberOfReactions <- 0

            list(
                pathway = df$pathway[1],
                rels = rels
            )
        }) %>%
        setNames(sapply(., function(x) x$pathway)) %>%
        lapply(function(x) x$rels)

    for (pathwayId in names(pathInfo)) {
        pathInfo[[pathwayId]]$title <- pathwayNames[pathwayId] %>% as.character()
    }

    pathInfo
}

#' @title SPIA combfunc method.
#' @description This function is combfunc from the original SPIA method.
#' @param p1 See SPIA function
#' @param p2 See SPIA function
#' @param combine See SPIA function
#' @return See SPIA function
combfunc <- function (p1 = NULL, p2 = NULL, combine)
{
    tm = na.omit(c(p1, p2))
    if (!all(tm >= 0 & tm <= 1)) {
        stop("values of p1 and p2 have to be >=0 and <=1 or NAs")
    }
    if (combine == "fisher") {
        k = p1 * p2
        comb = k - k * log(k)
        comb[is.na(p1)] <- p2[is.na(p1)]
        comb[is.na(p2)] <- p1[is.na(p2)]
        return(comb)
    }
    if (combine == "norminv") {
        comb = pnorm((qnorm(p1) + qnorm(p2))/sqrt(2))
        comb[is.na(p1)] <- p2[is.na(p1)]
        comb[is.na(p2)] <- p1[is.na(p2)]
        return(comb)
    }
}

#' @title SPIA method modified with inputs are KEGG pathway networks instead of a folder.
#' @description This function is modified from the original SPIA method to accept KEGG pathway networks as inputs.
#' @param de See SPIA function
#' @param all See SPIA function
#' @param path.info pathway information generated by getSPIAKEGGNetwork function
#' @param nB See SPIA function
#' @param verbose See SPIA function
#' @param beta See SPIA function
#' @param combine See SPIA function
#' @return See SPIA function
.SPIAMod <- function(de = NULL, all = NULL, path.info, nB = 2000, verbose = TRUE, beta = NULL, combine = "fisher") {
    if (is.null(de) | is.null(all)) {
        stop("de and all arguments can not be NULL!")
    }

    rel <- c("activation", "compound", "binding/association",
             "expression", "inhibition", "activation_phosphorylation",
             "phosphorylation", "inhibition_phosphorylation", "inhibition_dephosphorylation",
             "dissociation", "dephosphorylation", "activation_dephosphorylation",
             "state change", "activation_indirect effect", "inhibition_ubiquination",
             "ubiquination", "expression_indirect effect", "inhibition_indirect effect",
             "repression", "dissociation_phosphorylation", "indirect effect_phosphorylation",
             "activation_binding/association", "indirect effect",
             "activation_compound", "activation_ubiquination")

    if (is.null(beta)) {
        beta = c(1, 0, 0, 1, -1, 1, 0, -1, -1, 0, 0, 1, 0, 1,
                 -1, 0, 1, -1, -1, 0, 0, 1, 0, 1, 1)
        names(beta) <- rel
    }else {
        if (!all(names(beta) %in% rel) | length(names(beta)) !=
            length(rel)) {
            stop(paste("beta must be a numeric vector of length",
                       length(rel), "with the following names:", "\n",
                       paste(rel, collapse = ",")))
        }
    }

    datpT <- path.info

    datp <- list()
    path.names <- NULL
    hasR <- NULL
    for (jj in 1:length(datpT)) {
        sizem <- dim(datpT[[jj]]$activation)[1]
        s <- 0
        con <- 0
        for (bb in 1:length(rel)) {
            if(is.null(datpT[[jj]][[rel[bb]]])) next()
            con = con + datpT[[jj]][[rel[bb]]] * abs(sign(beta[rel[bb]]))
            s = s + datpT[[jj]][[rel[bb]]] * beta[rel[bb]]
        }
        z = matrix(rep(apply(con, 2, sum), dim(con)[1]), dim(con)[1],
                   dim(con)[1], byrow = TRUE)
        z[z == 0] <- 1
        datp[[jj]] <- s / z
        path.names <- c(path.names, datpT[[jj]]$title)
        hasR <- c(hasR, datpT[[jj]]$NumberOfReactions >= 1)
    }
    names(datp) <- names(datpT)
    names(path.names) <- names(datpT)
    tor <- lapply(datp, function(d) {
        sum(abs(d))
    }) == 0 |
        hasR |
        is.na(path.names)
    datp <- datp[!tor]
    path.names <- path.names[!tor]
    IDsNotP <- names(de)[!names(de) %in% all]
    if (length(IDsNotP) / length(de) > 0.01) {
        stop("More than 1% of your de genes have IDs are not present in the reference array!. Are you sure you use the right reference array?")
    }
    if (!length(IDsNotP) == 0) {
        cat("The following IDs are missing from all vector...:\n")
        cat(paste(IDsNotP, collapse = ","))
        cat("\nThey were added to your universe...")
        all <- c(all, IDsNotP)
    }
    if (length(intersect(names(de), all)) != length(de)) {
        stop("de must be a vector of log2 fold changes. The names of de should be included in the refference array!")
    }
    ph <- pb <- pcomb <- nGP <- pSize <- smPFS <- tA <- tAraw <- KEGGLINK <- NULL
    set.seed(1)

    for (i in 1:length(names(datp))) {
        path <- names(datp)[i]
        M <- datp[[path]]
        diag(M) <- diag(M) - 1
        X <- de[rownames(M)]
        noMy <- sum(!is.na(X))
        nGP[i] <- noMy
        okg <- intersect(rownames(M), all)
        ok <- rownames(M) %in% all
        pSize[i] <- length(okg)
        if ((noMy) > 0 & (abs(det(M)) > 1e-07)) {
            gnns <- paste(names(X)[!is.na(X)], collapse = "+")
            X[is.na(X)] <- 0
            pfs <- solve(M, -X)
            smPFS[i] <- sum(pfs - X)
            tAraw[i] <- smPFS[i]

            ph[i] <- phyper(q = noMy - 1, m = pSize[i], n = length(all) -
                pSize[i], k = length(de), lower.tail = FALSE)
            pfstmp <- NULL
            for (k in 1:nB) {
                x <- rep(0, length(X))
                names(x) <- rownames(M)
                x[ok][sample(1:sum(ok), noMy)] <- as.vector(sample(de,
                                                                   noMy))
                tt <- solve(M, -x)
                pfstmp <- c(pfstmp, sum(tt - x))
            }
            mnn <- median(pfstmp)
            pfstmp <- pfstmp - mnn
            ob <- smPFS[i] - mnn
            tA[i] <- ob
            if (ob > 0) {
                pb[i] <- sum(pfstmp >= ob) / length(pfstmp) *
                    2
                if (pb[i] <= 0) {
                    pb[i] <- 1 / nB / 100
                }
                if (pb[i] > 1) {
                    pb[i] <- 1
                }
            }
            if (ob < 0) {
                pb[i] <- sum(pfstmp <= ob) / length(pfstmp) *
                    2
                if (pb[i] <= 0) {
                    pb[i] <- 1 / nB / 100
                }
                if (pb[i] > 1) {
                    pb[i] <- 1
                }
            }
            if (ob == 0) {
                if (all(pfstmp == 0)) {
                    pb[i] <- NA
                }
                else {
                    pb[i] <- 1
                }
            }

            pcomb[i] <- combfunc(pb[i], ph[i], combine)
        }
        else {
            pb[i] <- ph[i] <- smPFS[i] <- pcomb[i] <- tAraw[i] <- tA[i] <- NA
        }
        if (verbose) {
            cat("\n")
            cat(paste("Done pathway ", i, " : ", substr(path.names[names(datp)[i]],
                                                        1, 30), "..", sep = ""))
        }
    }


    pcombFDR = p.adjust(pcomb, "fdr")
    phFdr = p.adjust(ph, "fdr")
    pcombfwer = p.adjust(pcomb, "bonferroni")
    Name = path.names[names(datp)]
    Status = ifelse(tA > 0, "Activated", "Inhibited")
    res <- data.frame(Name, ID = names(datp), pSize, NDE = nGP,
                      pNDE = ph, tA, pPERT = pb, pG = pcomb, pGFdr = pcombFDR,
                      pGFWER = pcombfwer, Status, stringsAsFactors = FALSE)
    res <- res[!is.na(res$pNDE),]
    res <- res[order(res$pG),]
    rownames(res) <- NULL
    res
}
