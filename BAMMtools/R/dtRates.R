dtRates.bammdata.diversification <- function(ephy, segmat, tol, ism) {
    dtrates <- .Call("dtrates", ephy, segmat, tol, ism, 0L, PACKAGE = "BAMMtools")
    for (i in 1:2) {
        # recover original ordering of segmat. see dtRates source code
        names(dtrates[[i]]) <- rownames(segmat)
        dtrates[[i]] <- dtrates[[i]][as.character(1:nrow(segmat))]
        names(dtrates[[i]]) <- NULL
    }
    if (sum(is.na(dtrates[[1]]))) {
        warning(sprintf("Found %d NA speciation rates. Coercing to zero.", sum(is.na(dtrates[[1]]))))
        dtrates[[1]][is.na(dtrates[[1]])] <- 0
    }
    if (sum(is.na(dtrates[[2]]))) {
        warning(sprintf("Found %d NA extinction rates. Coercing to zero.", sum(is.na(dtrates[[2]]))))
        dtrates[[2]][is.na(dtrates[[2]])] <- 0
    }
    return(dtrates)
}

dtRates.bammdata.trait <- function(ephy, segmat, tol, ism) {
    dtrates <- .Call("dtrates", ephy, segmat, tol, ism, 1L, PACKAGE = "BAMMtools")
    # recover original ordering of segmat. see dtRates source code 
    names(dtrates) <- rownames(segmat)
    dtrates <- dtrates[as.character(1:nrow(segmat))]
    names(dtrates) <- NULL
    if (sum(is.na(dtrates))) {
        warning(sprintf("Found %d NA phenotypic rates. Coercing to zero.", sum(is.na(dtrates))))
        dtrates[is.na(dtrates)] <- 0
    }
}

dtRates.bammdata.binarystate <- function(ephy, segmat, tol, ism) {
    dtrates <- .Call("dtrates_binary", ephy, segmat, tol, ism, PACKAGE = "BAMMtools")
    for (i in 1:2) {
        # recover original ordering of segmat. see dtRates source code
        names(dtrates[[i]]) <- rownames(segmat)
        dtrates[[i]] <- dtrates[[i]][as.character(1:nrow(segmat))]
        names(dtrates[[i]]) <- NULL
    }
    if (sum(is.na(dtrates[[1]]))) {
        warning(sprintf("Found %d NA forward transition rates. Coercing to zero.", sum(is.na(dtrates[[1]]))))
        dtrates[[1]][is.na(dtrates[[1]])] <- 0
    }
    if (sum(is.na(dtrates[[2]]))) {
        warning(sprintf("Found %d NA reverse transition rates. Coercing to zero.", sum(is.na(dtrates[[2]]))))
        dtrates[[2]][is.na(dtrates[[2]])] <- 0
    }
    return(dtrates)
}

dtRates <- function (ephy, tau, ism = NULL, tmat = FALSE) {
    if (inherits(ephy, "bammdata")) {
        if (attributes(ephy)$order != "cladewise") {
            stop("Function requires tree in 'cladewise' order")
        }
        if (any(sapply(ephy$eventBranchSegs, function(x) is.unsorted(x[,1])))) {
            stop("Why are eventBranchSegs unordered?")
        }
        tH <- max(ephy$end)
        segmat <- segMap(ephy, tau)
        tol <- 0.00001
        if (is.null(ism)) {
            ism <- as.integer(1:length(ephy$eventBranchSegs))
        } else {
            ism <- as.integer(ism)
        }
        if (ism[length(ism)] > length(ephy$eventBranchSegs)) {
            warning("Sample index out of range")
            ism <- as.integer(1:length(ephy$eventBranchSegs))
        }
        # name the rows so we can recover original ordering
        # which corresponds to preorder traversal of tree (minus root)
        rownames(segmat) <- 1:nrow(segmat)
        # reorder by node number to match eventBranchSegs.
        # note that order of rows indexed by same node number
        # does not change, which preserves temporal sequence of
        # segments within each branch.
        segmat <- segmat[order(segmat[, 1]), ]
        if (inherits(ephy, "bammdata-diversification") || ephy$type == 'diversification') {
            dtrates <- dtRates.bammdata.diversification(ephy, segmat, tol, ism)
        } else if (inherits(ephy, "bammdata-trait") || ephy$type == 'trait') {
            dtrates <- dtRates.bammdata.trait(ephy, segmat, tol, ism)
        } else if (inherits(ephy, "bammdata-binarystate") || ephy$type == 'binarystate') {
            dtrates <- dtRates.bammdata.binarystate(ephy, segmat, tol, ism)
        } else {
            stop("Unrecognized model type.")
        }
        if (tmat) {
            segmat <- segmat[as.character(1:nrow(segmat)),]  # recover original ordering
            ephy$dtrates <- list(tau = tau, rates = dtrates, tmat = segmat)
            return(ephy)
        }
        ephy$dtrates <- list(tau = tau, rates = dtrates)
        return(ephy)
    }
    return(NULL)
}
