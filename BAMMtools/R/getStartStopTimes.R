#############################################################
#
#	getStartStopTimes(....)
#
#	adds begin and end times (absolute time) to each edge of 
#	phylogenetic tree

getStartStopTimes <- function(phy){
	res <- .Call('get_startstoptimes', phy$edge, phy$edge.length, Ntip(phy))
    phy$begin <- res[[1]][phy$edge[,2]]
    phy$end <- res[[2]][phy$edge[,2]]
    return(phy)
}
