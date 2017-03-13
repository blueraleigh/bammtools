######################################
#	Internal function called by dtRates(...)
#
#
segMap <- function(ephy, tau) {
	# tolerance for determining if segments fail to tesselate branch
	tol <- 0.0001
	tH <- max(ephy$end)
	# size to discretize branches into
	tau <- tH * tau
	# vector to track remainders when segments don't tesselate branch
	remainder <- numeric(max(ephy$edge[,1]))
	dtsegs <- vector("list", nrow(ephy$edge))
	for (i in 1:nrow(ephy$edge)) {
		if (remainder[ephy$edge[i,1]] > 0) {
			if (ephy$begin[i] + remainder[ephy$edge[i,1]] > ephy$end[i]) {
				remainder[ephy$edge[i,2]] <- (ephy$begin[i] + remainder[ephy$edge[i,1]]) - ephy$end[i]
				segs <- ephy$begin[i]
			}
			else {
				segs <- seq(ephy$begin[i] + remainder[ephy$edge[i,1]], ephy$end[i], tau)
				segs <- c(ephy$begin[i], segs)
			}
		}
		else {
			segs <- seq(ephy$begin[i], ephy$end[i], tau)
		}
		if (length(segs) > 1) {
			if (ephy$end[i] - tail(segs,1) > tol) {
				remainder[ephy$edge[i,2]] <- tau - (ephy$end[i] - tail(segs,1))
				segs <- c(segs, ephy$end[i])
			}
			segs <- cbind(segs[-length(segs)], segs[-1])
			segs <- cbind(rep(ephy$edge[i,2], nrow(segs)), segs)
		}
		else {
			if (remainder[ephy$edge[i,1]] == 0) {
				remainder[ephy$edge[i,2]] <- tau - (ephy$end[i] - tail(segs,1))
			}
			segs <- matrix(c(ephy$edge[i,2], ephy$begin[i], ephy$end[i]), nrow=1, ncol=3) 
		}
		dtsegs[[i]] <- segs
	}
	dtsegs <- do.call(rbind,dtsegs)
	return(dtsegs)
}



# segMap = function(nodes,begin,end,tau)
# {
	# foo = function(x,tau)
	# {
		# len = (x[3] - x[2])/tau if (len%%1 == 0) len = len+1
		# ret = seq(x[2],x[3],length.out=len)
		# if(length(ret) == 1) return(matrix(x,nrow=1))
		# #ret = seq(x[2],x[3],length.out=length(ret))
		# ret = rep(ret,each=2) ret=ret[-c(1,length(ret))]
		# ret = matrix(ret,ncol=2,byrow=TRUE)
		# return(cbind(matrix(rep(as.integer(x[1]),nrow(ret)),ncol=1), ret))
	# }
	# times = cbind(nodes,begin,end)
	# ret = apply(times,1,foo,tau)
	# return(do.call(rbind,ret))	
# }
