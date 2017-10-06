##' @title Create \code{bammdata} object from MCMC output
##'
##' @description \code{getEventData} Reads shift configuration data (the
##'     "event data" output) from a \code{BAMM} analysis and creates a
##'     \code{bammdata} object. The \code{bammdata} object is fundamental
##'     for extracting information about macroevolutionary rate variation
##'     through time and among lineages.
##'
##' @param phy An object of class \code{phylo} - specifically, the
##'     time-calibrated tree that was analyzed with \code{BAMM}.
##'     Alternatively, a character string specifying the path to a
##'     newick-formatted tree.
##' @param eventdata A character string specifying the path to a \code{BAMM}
##'     event-data file. Alternatively, an object of class \code{data.frame}
##'     that includes the event data from a \code{BAMM} run.
##' @param burnin A numeric indicating the fraction of posterior samples to
##'     discard as burn-in.
##' @param nsamples An integer indicating the number of posterior samples to
##'     include in the \code{bammdata} object. May be \code{NULL}.
##' @param verbose A logical. If \code{TRUE} progess is outputted to the
##'     console. Defaults to \code{FALSE}.
##' @param type A character string. Either "diversification" or "trait"
##'     depending on your \code{BAMM} analysis.
##'
##' @details In the \code{BAMM} framework, an "event" defines a
##'     macroevolutionary process of diversification or trait evolution. Every
##'     sample from the posterior includes at least one process, defined by
##'     such an "event". If a given sample includes just a single event, then
##'     the dynamics of diversification or trait evolution can be described
##'     entirely by a single time-constant or time-varying process that begins
##'     at the root of the tree. Any sample from the posterior distribution
##'     may include a complex mixture of distinct processes. To represent
##'     temporal heterogeneity in macroevolutionary rates, \code{BAMM} models
##'     a rate \eqn{R}, e.g. speciation, as a function that changes
##'     exponentially with time:
##'
##'     \eqn{R(t) = R(0)*exp(b*t)}.
##'
##'     Here \eqn{R(0)} is the initial rate and \eqn{b} is a parameter
##'     determining how quickly that rate grows or decays with time. 
##'
##'     The \code{eventdata} file (or data frame) is a record of events and
##'     associated parameters that were sampled with \code{BAMM} during
##'     simulation of the posterior with reversible jump MCMC. This complex,
##'     information-rich file is processed into a \code{bammdata} object,
##'     which serves as the core data object for numerous downstream analyses.
##'     From a \code{bammdata} object, you can summarize rate variation
##'     through time, among clades, extract locations of rate shifts,
##'     summarize clade-specific rates of speciation and extinction, and more.
##'
##'     In general, the user does not need to be concerned with the details of
##'     a \code{bammdata} object. The object is used as input by a number of
##'     \code{BAMMtools} functions. 
##'
##'     The parameter \code{nsamples} can be used to reduce the total amount
##'     of data included in the raw eventdata output from a \code{BAMM} run.
##'     The final \code{bammdata} object will consist of all data for
##'     \code{nsamples} from the posterior. These \code{nsamples} are equally
##'     spaced after discarding some \code{burnin} fraction as "burn-in". If
##'     \code{nsamples} is set to \code{NULL}, the \code{bammdata} object will
##'     include all samples in the posterior after discarding the
##'     \code{burnin} fraction.
##'
##' @return A list with many components:
##' \itemize{
##'     \item{edge} {See documentation for class \code{phylo} in package ape.}
##'     \item{Nnode} {See documentation for class \code{phylo} in package
##'         ape.}
##'     \item{tip.label} {See documentation for class \code{phylo} in package
##'         ape.}
##'     \item{edge.length} {See documentation for class \code{phylo} in
##'         package ape.}
##'     \item{begin} {The beginning time of each branch in absolute time (the
##'         root is set to time zero)}
##'     \item{end} {The ending time of each branch in absolute time.}
##'     \item{numberEvents} {An integer vector with the number of events
##'         contained in \code{phy} for each posterior sample. The length of
##'         this vector is equal to the number of posterior samples in the
##'         \code{bammdata} object.}
##'     \item{eventData} {A list of dataframes. Each element is a single
##'         posterior sample. Each row in a dataframe holds the data for a
##'         single event. Data associated with an event are: \code{node} - a
##'         node number. This identifies the branch where the event
##'         originates. \code{time} - this is the absolute time on that branch
##'         where the event originates (with the root at time 0). \code{lam1}
##'         - an initial rate of speciation or trait evolution. \code{lam2} -
##'         a decay/growth parameter. \code{mu1} - an initial rate of
##'         extinction. \code{mu2} - a decay/growth parameter. \code{index} -
##'         a unique integer associated with the event. See 'Details'.}
##'     \item{eventVectors} {A list of integer vectors. Each element is a
##'         single posterior sample. For each branch in \code{phy} the index
##'         of the event that occurs along that branch. Branches are ordered
##'         increasing here and elsewhere.}
##'     \item{eventBranchSegs} {A list of matrices. Each element is a single
##'         posterior sample. Each matrix has four columns: \code{Column 1}
##'         identifies a node in \code{phy}. \code{Column 2} identifies the
##'         beginning time of the branch or segment of the branch that
##'         subtends the node in \code{Column 1}. \code{Column 3} identifies
##'         the ending time of the branch or segment of the branch that
##'         subtends the node in \code{Column 1}. \code{Column 4} identifies
##'         the index of the event that occurs along the branch or segment of
##'         the branch that subtends the node in \code{Column 1}.}
##'     \item{tipStates} {A list of integer vectors. Each element is a single
##'         posterior sample. For each tip the index of the event that occurs
##'         along the branch subtending the tip. Tips are ordered increasing
##'         here and elsewhere.}
##'     \item{tipLambda} {A list of numeric vectors. Each element is a single
##'         posterior sample. For each tip the rate of speciation or trait
##'         evolution at the end of the terminal branch subtending that tip.}
##'     \item{tipMu} {A list of numeric vectors. Each element is a single
##'         posterior sample. For each tip the rate of extinction at the end
##'         of the terminal branch subtending that tip. Meaningless if working
##'         with \code{BAMM} trait results.}
##'     \item{meanTipLambda} {For each tip the mean of the marginal posterior
##'         density of the rate of speciation or trait evolution at the end of
##'         the terminal branch subtending that tip.}
##'     \item{meanTipMu} {For each tip the mean of the marginal posterior
##'         density of the rate of extinction at the end of the terminal
##'         branch subtending that tip. Meaningless if working with
##'         \code{BAMM} trait results.}
##'     \item{type} {A character string. Either "diversification" or "trait"
##'         depending on your \code{BAMM} analysis.}
##'     \item{downseq} {An integer vector holding the nodes of \code{phy}. The
##'         order corresponds to the order in which nodes are visited by a
##'         pre-order tree traversal.}
##'     \item{lastvisit} {An integer vector giving the index of the last node
##'         visited by the node in the corresponding position in
##'         \code{downseq}. \code{downseq} and \code{lastvisit} can be used to
##'         quickly retrieve the descendants of any node. e.g. the descendants
##'         of node 89 can be found by
##'         \code{downseq[which(downseq==89):which(downseq==lastvisit[89])}.}
##' }
##'
##' @note Currently the function does not check for duplicate tip labels in
##'     \code{phy}, which may cause the function to choke.
##'
##' @author Dan Rabosky, Mike Grundler
##'
##' @seealso \code{\link{summary.bammdata}}, \code{\link{plot.bammdata}},
##'     \code{\link{dtRates}}.
##' 
##' @references \url{http://bamm-project.org/}
##'
##' @examples
##' data(primates, events.primates)
##' xx <- getEventData(primates, events.primates, burnin=0.25, nsamples=500,
##'                    type = 'trait')
##' 
##' # compute mean phenotypic rate for primate body size evolution:
##' brates <- getCladeRates(xx)
##' mean(brates$beta)
##' 
##' # Plot rates:
##' plot(xx)
##' @keywords models
##' @export
getEventData <- function(phy, eventdata, burnin=0, nsamples=NULL, type=c("diversification","trait","binarystate")) {

    type <- match.arg(type)
    phy <- getRecursiveSequence(phy)
    phy <- getStartStopTimes(phy)

    if (class(eventdata) == 'data.frame') {
        eventspersample <- table(eventdata[,1])
        uniquegens <- as.integer(names(eventspersample))
        eventspersample <- unname(eventspersample)
    }
    else if (class(eventdata) == 'character') {
        eventdata <- read.csv(eventdata, header=TRUE, stringsAsFactors=FALSE)
        eventspersample <- table(eventdata[,1])
        uniquegens <- as.integer(names(eventspersample))
        eventspersample <- unname(eventspersample)
    } 
    else {
        err.string <- c('eventdata arg invalid\n\nType is ', class(eventdata), '\n', sep='')
        stop(err.string)
    }

    # check storage modes. sometimes a parameter vec that is all 0s
    # is stored as an integer, which causes problems in C API function
    # calls that assume it is a double
    for (i in 4:ncol(eventdata)) {
        storage.mode(eventdata[,i]) <- "double"
    }

    eventdata$eventnode <- getmrca(phy, as.integer(match(eventdata$leftchild, phy$tip.label)), as.integer(match(eventdata$rightchild, phy$tip.label, nomatch = 0L)))

    if (is.null(nsamples)) {
        nsamples <- length(uniquegens)
    }

    nsamples <- min(nsamples, length(uniquegens))
    drop <- max(nsamples * burnin, 1)
    from <- uniquegens[drop]
    to <- uniquegens[nsamples]
    firstgoodrow <- match(from, eventdata[,1], nomatch=0)
    lastgoodrow <- match(to, rev(eventdata[,1]), nomatch=0) + nrow(eventdata) - 1
    keeprows <- firstgoodrow:lastgoodrow

    # row indices
    ix <- 1:(length(firstgoodrow:lastgoodrow))

    # loop over the posterior samples and process each one
    res <- tapply(ix, eventdata[keeprows,1], .process.events, eventdata[keeprows,], phy, type, simplify=FALSE)

    # combine back into canonical form
    phy$eventData <- do.call(list, lapply(res, '[[', 'eventData'))
    phy$eventBranchSegs <- do.call(list, lapply(res, "[[", "eventBranchSegs"))
    phy$eventVectors <- do.call(list, lapply(res, "[[", "eventVectors"))
    phy$tipStates <- do.call(list, lapply(res, "[[", "tipStates"))
    phy$numberEvents <- eventspersample[drop:nsamples]

    # tip rates
    if (type == "diversification" || type == "trait") {
        stoptime <- max(phy$end)
        phy$tipLambda <- do.call(list, lapply(res, function(x) {
            tipstates <- x$tipStates
            ed <- x$eventData
            tiplam <- exponentialRate(stoptime - ed$time[tipstates], ed$lam1[tipstates], ed$lam2[tipstates])
            tiplam
        }))

        phy$tipMu <- do.call(list, lapply(res, function(x) {
            tipstates <- x$tipStates
            ed <- x$eventData
            tipmu <- ed$mu1[tipstates]
            tipmu
        }))

        phy$meanTipLambda <- colMeans(do.call(rbind, phy$tipLambda))
        phy$meanTipMu <- colMeans(do.call(rbind, phy$tipMu))
    }

    # finalize
    phy$type <- type
    class(phy) <- c('bammdata', sprintf('bammdata-%s', type))
    return(phy)
}




# `ix` indexes the rows of `ed` (the BAMM event 
# data output) corresponding to particular posterior
# sample. this is called on each posterior
# sample to build the bammdata object.
.process.events <- function(ix, ed, phy, type) {

    eventData <- ed[ix, c(ncol(ed), 4, 5:(ncol(ed)-1))]

    ## for backwards compatibility keep original column format when type is 'trait' or 'diversification'
    if (type == "diversification" || type == "trait") {
        colnames(eventData) <- c("node", "time", "lam1", "lam2", "mu1", "mu2")
    } else {
        colnames(eventData) <- c('node', 'time', colnames(ed)[5:(ncol(ed)-1)])
    }

    eventData <- eventData[order(eventData[,2]),]
    eventData$index <- 1:nrow(eventData)
    rownames(eventData) <- NULL
    
    statevec <- rep(1, nrow(phy$edge))

    if (nrow(eventData) > 1) {
        for (k in 2:nrow(eventData)) {
            s1 <- which(phy$downseq == eventData[k,1])
            s2 <- which(phy$downseq == phy$lastvisit[eventData[k,1]])
            descSet <- phy$downseq[s1:s2]
            isDescendantNode <- phy$edge[,2] %in% descSet
            statevec[isDescendantNode] <- k
        }
    }

    tipstates <- numeric(length(phy$tip.label))
    tipstates <- statevec[phy$edge[,2] <= phy$Nnode + 1]
    tipstates <- tipstates[order(phy$edge[phy$edge[,2] <= phy$Nnode + 1, 2])]

    eventBranchSegs <- matrix(0, nrow=(max(phy$edge[,1]) + nrow(eventData) - 2), ncol=4)
    
    pos <- 1
    
    is_noEventBranch <- !(phy$edge[,2] %in% eventData[,1])
    
    if (sum(is_noEventBranch) > 0) {
        
        eventBranchSegs[1:sum(is_noEventBranch), 1] <- phy$edge[,2][is_noEventBranch]
        eventBranchSegs[1:sum(is_noEventBranch), 2] <- phy$begin[is_noEventBranch]
        eventBranchSegs[1:sum(is_noEventBranch), 3] <- phy$end[is_noEventBranch]
        eventBranchSegs[1:sum(is_noEventBranch), 4] <- statevec[is_noEventBranch]
            
    } else {
        eventBranchSegs <- cbind(phy$edge[, 2], phy$begin, phy$end, statevec)
    }
    
    eventnodeset <- unique(eventData[-1,1])
    pos <- 1 + sum(is_noEventBranch)
    for (k in eventnodeset) {
        events.on.branch.index <- which(eventData[,1] == k)
        events.on.branch <- eventData[events.on.branch.index, ]
        events.on.branch <- events.on.branch[order(events.on.branch[,2]), ]
        
        fBranch <- match(k, phy$edge[,2])
        start.times <- c(phy$begin[fBranch], events.on.branch[,2])
        stop.times <- c(events.on.branch[,2], phy$end[fBranch])
        parent <- phy$edge[fBranch,1]
        if (parent == (Ntip(phy) + 1)) {
            # Parent is root:
            proc.set <- c(1, events.on.branch.index)   
        } else {
            proc.set <- c(statevec[match(parent, phy$edge[,2])], events.on.branch.index)
            #proc.set <- c(statevec[phy$edge[,2] == parent], events.on.branch.index)
        }
            
        zzindex = pos:(pos + nrow(events.on.branch))   
            
        eventBranchSegs[zzindex, 1] <- rep(k, length(zzindex))
        eventBranchSegs[zzindex, 2] <- start.times
        eventBranchSegs[zzindex, 3] <- stop.times
        eventBranchSegs[zzindex, 4] <- proc.set     
        pos <- pos + 1 + nrow(events.on.branch)
    }
    
    eventBranchSegs <- eventBranchSegs[order(eventBranchSegs[,1]),]

    return(list(eventData=eventData, eventBranchSegs=eventBranchSegs, eventVectors=statevec, tipStates=tipstates))
}
