getMeanBranchLengthTree = function(ed, tau=0.01, rate='speciation') {
    ed = dtRates(ed, tau=tau, tmat=TRUE)
    obj = list()
    obj$phy = as.phylo(ed)

    r8ts = switch(ed$type,
        trait = tapply(ed$dtrates$rates, ed$dtrates$tmat[,1], mean)[as.character(ed$edge[,2])],
        diversification =
        if (rate == 'speciation') {
            tapply(ed$dtrates$rates[[1]], ed$dtrates$tmat[,1], mean)[as.character(ed$edge[,2])]
        } else if (rate == 'extinction') {
            tapply(ed$dtrates$rates[[2]], ed$dtrates$tmat[,1], mean)[as.character(ed$edge[,2])]
        } else if (rate == 'diversification') {
            tapply(ed$dtrates$rates[[1]] - ed$dtrates$rates[[2]], ed$dtrates$tmat[,1], mean)[as.character(ed$edge[,2])]
        } else {
            stop("Unrecognized rate type")
        },
        binarystate = tapply(ed$dtrates$rates[[1]], ed$dtrates$tmat[,1], mean)[as.character(ed$edge[,2])] * ed$edge.length
    )

    obj$phy$edge.length = r8ts
    obj$mean = mean(obj$phy$edge.length)
    obj$median = median(obj$phy$edge.length)
    return(obj)
}
