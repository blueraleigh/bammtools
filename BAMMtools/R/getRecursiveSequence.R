#############################################################
#
#	getRecursiveSequence(....)
#
#	Private function, called by getEventDataDiversification

getRecursiveSequence = function(phy)
{
    L <- .C('get_recursivesequence', 
        phy$edge[,1], 
        phy$edge[,2], 
        Ntip(phy)+1L, 
        0, 
        nrow(phy$edge), 
        integer(2*Ntip(phy)-1), 
        integer(2*Ntip(phy)-1)
    )
    phy$downseq <- L[[6]]
    phy$lastvisit <- L[[7]]
	return(phy)
}