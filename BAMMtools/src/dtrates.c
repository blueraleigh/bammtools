#include <R.h>
#include <Rinternals.h>

//forward function declarations
SEXP getListElement(SEXP list, char *str);
SEXP dtrates(SEXP ephy, SEXP segmat, SEXP tol, SEXP sample, SEXP type);
SEXP dtrates_binary(SEXP ephy, SEXP segmat, SEXP tol, SEXP sample);
double getDblMatrixELT(SEXP matrix, int row, int col);
double getMeanRateExponential(double t1, double t2, double p1, double p2);
double getTimeIntegratedBranchRate(double t1, double t2, double p1, double p2);

// convenience function for accessing list element by name
SEXP getListElement(SEXP list, char *str) {
	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	int i;
	for (i = 0; i < LENGTH(list); i++) {
		if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
			elmt = VECTOR_ELT(list, i);
			break;
		}
	}
	return elmt;
}

// convenience function for getting an element from a
// flattened matrix using row column notation 
double getDblMatrixELT(SEXP matrix, int row, int col) {
	int nrow = INTEGER(getAttrib(matrix, R_DimSymbol))[0];
	return REAL(matrix)[row + nrow*col];
}

// exponential rate function used by bamm for speciation/extinction
// and continuous trait analysis
double getMeanRateExponential(double t1, double t2, double p1, double p2) {	
	if (p2 == 0.) {
		return p1;
	} else if (p2 < 0) {
		return (p1/p2)*(exp(t2*p2) - exp(t1*p2))/(t2 - t1);
	} else {
		return (p1/p2)*(2.*p2*(t2-t1) + exp(-t2*p2) - exp(-t1*p2))/(t2 - t1);
	}
}

// time integrated value of the above function
double getTimeIntegratedBranchRate(double t1, double t2, double p1, double p2) {
	if (p2 == 0.) {
		return (t2 - t1) * p1;
	} else if (p2 < 0) {
		return (p1/p2)*(exp(t2*p2) - exp(t1*p2));
	} else {
		return (p1/p2)*(2.*p2*(t2-t1) + exp(-t2*p2) - exp(-t1*p2));
	}
}

/*
** function to calculate rates through time along branches of phylogeny
**
** `ephy` 
** is the bammdata object, a list holding all the relevant data.
**
** `segmat` 
** is a matrix where the rows represent branches of a phylogeny
** broken up into many small segments.  each row indexes one such segment.
** columns 2 and 3 give the starting and ending times of that segment and
** column 1 is the node of the phylogeny to which that segement belongs.
**
** `tol` 
** is a precision parameter used for comparing starting and ending 
** times of approximating segments and starting and ending times of branches
** or branch segments on the phylogeny.
**
** `sample` 
** is a vector of posterior indices to calculate rates on
** 
** `type` 
** is a flag indicating the type of bamm analysis.
** this determines how the rates are calculated/interpreted
**   0 => speciation/extinction
**   1 => trait
**   2 => binary state
*/

SEXP dtrates(SEXP ephy, SEXP segmat, SEXP tol, SEXP sample, SEXP type) {
	double eps = REAL(tol)[0];
	
	int k, j, l, nprotect = 0;
	int nsamples = LENGTH(sample);	
	int nsegs = INTEGER(getAttrib(segmat, R_DimSymbol))[0];
	
	SEXP rates, erates;
	PROTECT(rates = allocVector(REALSXP, nsegs)); nprotect++;
	for(k = 0; k < nsegs; k++) {
		REAL(rates)[k] = 0.;
	}
	if (INTEGER(type)[0] == 0) {
		PROTECT(erates = allocVector(REALSXP, nsegs)); nprotect++;
		for(k = 0; k < nsegs; k++) {
			REAL(erates)[k] = 0.;
		}
	}
	
	int nrow, node, nnode, event, nxtevent, isGoodStart, isGoodEnd, place_holder;
	double begin, end, Start, lam1, lam2, mu1, mu2, relStart, relEnd, rightshift, leftshift, erightshift, eleftshift, ret;		
	
	for (k = INTEGER(sample)[0] - 1; k < INTEGER(sample)[nsamples - 1]; k++) {
		SEXP eventSegs, eventData;
		
		eventSegs = PROTECT(VECTOR_ELT(getListElement(ephy, "eventBranchSegs"), k)); nprotect++;
		eventData = PROTECT(VECTOR_ELT(getListElement(ephy, "eventData"), k)); nprotect++;
				
		nrow = INTEGER(getAttrib(eventSegs, R_DimSymbol))[0];
		place_holder = 0;
		//move down the rows of eventSegs. eventSegs is ordered by node number in first column
		for(j = 0; j < nrow; j++) {
			//eventSegs is 4 column matrix, node is in first column stored as double
			node = (int) REAL(eventSegs)[j + nrow * 0];
			
			//event index is in fourth column stored as double
			event = (int) REAL(eventSegs)[j  + nrow * 3];
			
			//begin and end of current branch segment are in second and third columns stored as doubles
			begin = REAL(eventSegs)[j + nrow * 1];
			end = REAL(eventSegs)[j + nrow * 2];
			
			//find next node to later check for shift point on branch
			if (j < (nrow-1)) {
				nnode = (int) REAL(eventSegs)[(j+1) + nrow * 0];
				nxtevent = (int) REAL(eventSegs)[(j+1) + nrow * 3];
				//Rprintf("%d\n", nxtevent);
			}
			
			//eventData is dataframe holding event parameters for the current tree
			//need to find the row that corresponds to the event index. in eventData
			//the rows are strictly ordered such that row 0 = event1, row 1 = event2, etc.
			Start = REAL(getListElement(eventData, "time"))[event-1];
			lam1 = REAL(getListElement(eventData, "lam1"))[event-1];
			lam2 = REAL(getListElement(eventData, "lam2"))[event-1];
			if (INTEGER(type)[0] == 0) {
			    mu1 = REAL(getListElement(eventData, "mu1"))[event-1];
			    mu2 = REAL(getListElement(eventData, "mu2"))[event-1];
			}
						
			//need to find which approximating dt segments match this branch segment
			//these are passed in strict order by node number so we only need to search top to bottom
			//and can ignore everything we've been over already
			for (l = place_holder; l < nsegs; l++) {
				if ( (int) getDblMatrixELT(segmat, l, 0) == node) {
					//isGoodStart = REAL(segbegin)[l] >= begin;
					isGoodStart = ( (getDblMatrixELT(segmat, l, 1) - begin) >= 0. || ( (getDblMatrixELT(segmat, l, 1) - begin) < 0. && (getDblMatrixELT(segmat, l, 1) - begin) >= -1.*eps) );
					//isGoodEnd = REAL(segend)[l] <= end;
					isGoodEnd =  ( (getDblMatrixELT(segmat, l, 2) - end) <= 0. || ( (getDblMatrixELT(segmat, l, 2) - end) > 0. && (getDblMatrixELT(segmat, l, 2) - end) <= eps) );

					if (isGoodStart && isGoodEnd) {					
						relStart = getDblMatrixELT(segmat, l, 1) - Start;
						relEnd = getDblMatrixELT(segmat, l, 2) - Start;
						ret = getMeanRateExponential(relStart,relEnd,lam1,lam2);
						REAL(rates)[l] += ret/((double) nsamples);
						if (INTEGER(type)[0] == 0) {
							ret = getMeanRateExponential(relStart,relEnd,mu1,mu2);
							REAL(erates)[l] += ret/((double) nsamples);
						}
					}
					//check for shift straddle
					if (node == nnode) {
						// there is a shift on the branch, need to see if the dtseg straddles it
						isGoodStart = getDblMatrixELT(segmat, l, 1) < end;
						isGoodEnd = getDblMatrixELT(segmat, l, 2) > end;
						if (isGoodStart && isGoodEnd) {	
							relStart = getDblMatrixELT(segmat, l, 1) - Start;
							relEnd = end - Start;
							leftshift = getTimeIntegratedBranchRate(relStart,relEnd,lam1,lam2);
							if (INTEGER(type)[0] == 0) {
								eleftshift = getTimeIntegratedBranchRate(relStart,relEnd,mu1,mu2);
							}
							relStart = 0.;
							relEnd = getDblMatrixELT(segmat, l, 2) - end;
							lam1 = REAL(getListElement(eventData, "lam1"))[nxtevent-1];
							lam2 = REAL(getListElement(eventData, "lam2"))[nxtevent-1];
							
							rightshift = getTimeIntegratedBranchRate(relStart,relEnd,lam1,lam2);
							ret = (leftshift+rightshift)/(getDblMatrixELT(segmat, l, 2) - getDblMatrixELT(segmat, l, 1));
							REAL(rates)[l] += ret/((double) nsamples);
							
							if (INTEGER(type)[0] == 0) {
								mu1 = REAL(getListElement(eventData, "mu1"))[nxtevent-1];
								mu2 = REAL(getListElement(eventData, "mu2"))[nxtevent-1];
								erightshift = getTimeIntegratedBranchRate(relStart,relEnd,mu1,mu2);
								ret = (eleftshift+erightshift)/(getDblMatrixELT(segmat, l, 2) - getDblMatrixELT(segmat, l, 1));
								REAL(erates)[l] += ret/((double) nsamples);
							}
							place_holder = l; place_holder++;
							break;
						}
					}
				}
				else {
					place_holder = l;
					break;
				}
			}			
		}
		UNPROTECT(2); nprotect -= 2; //protected eventSegs and eventData, which we no longer need
	}
	if (INTEGER(type)[0] == 0) {
		SEXP retlist;
		PROTECT(retlist = allocVector(VECSXP, 2)); nprotect++;
		SET_VECTOR_ELT(retlist, 0, rates);
		SET_VECTOR_ELT(retlist, 1, erates);
		UNPROTECT(nprotect);
		return retlist;
	}
	UNPROTECT(nprotect);
	return rates;
}




/*
** Same function as the above lazily copied
** to work with binary data
*/
SEXP dtrates_binary(SEXP ephy, SEXP segmat, SEXP tol, SEXP sample) {
	double eps = REAL(tol)[0];
	
	int k, j, l, nprotect = 0;
	int nsamples = LENGTH(sample);	
	int nsegs = INTEGER(getAttrib(segmat, R_DimSymbol))[0];
	
	SEXP q01, q10;
	PROTECT(q01 = allocVector(REALSXP, nsegs)); nprotect++;
	PROTECT(q10 = allocVector(REALSXP, nsegs)); nprotect++;
	memset(REAL(q01), 0.0, nsegs * sizeof(double));
	memset(REAL(q10), 0.0, nsegs * sizeof(double));
	
	int nrow, node, nnode, event, nxtevent, isGoodStart, isGoodEnd, place_holder;
	double begin, end, forward_rate, reverse_rate, rightshift_r, leftshift_r, rightshift_f, leftshift_f, ret;		
	
	for (k = INTEGER(sample)[0] - 1; k < INTEGER(sample)[nsamples - 1]; k++) {
		SEXP eventSegs, eventData;
		
		eventSegs = PROTECT(VECTOR_ELT(getListElement(ephy, "eventBranchSegs"), k)); nprotect++;
		eventData = PROTECT(VECTOR_ELT(getListElement(ephy, "eventData"), k)); nprotect++;
				
		nrow = INTEGER(getAttrib(eventSegs, R_DimSymbol))[0];
		place_holder = 0;
		//move down the rows of eventSegs. eventSegs is ordered by node number in first column
		for(j = 0; j < nrow; j++) {
			//eventSegs is 4 column matrix, node is in first column stored as double
			node = (int) REAL(eventSegs)[j + nrow * 0];
			
			//event index is in fourth column stored as double
			event = (int) REAL(eventSegs)[j  + nrow * 3];
			
			//begin and end of current branch segment are in second and third columns stored as doubles
			begin = REAL(eventSegs)[j + nrow * 1];
			end = REAL(eventSegs)[j + nrow * 2];
			
			//find next node to later check for shift point on branch
			if (j < (nrow-1)) {
				nnode = (int) REAL(eventSegs)[(j+1) + nrow * 0];
				nxtevent = (int) REAL(eventSegs)[(j+1) + nrow * 3];
				//Rprintf("%d\n", nxtevent);
			}
			
			//eventData is dataframe holding event parameters for the current tree
			//need to find the row that corresponds to the event index. in eventData
			//the rows are strictly ordered such that row 0 = event1, row 1 = event2, etc.
			forward_rate = REAL(getListElement(eventData, "q01"))[event-1];
			reverse_rate = REAL(getListElement(eventData, "q10"))[event-1];
						
			//need to find which approximating dt segments match this branch segment
			//these are passed in strict order by node number so we only need to search top to bottom
			//and can ignore everything we've been over already
			for (l = place_holder; l < nsegs; l++) {
				if ( (int) getDblMatrixELT(segmat, l, 0) == node) {
					isGoodStart = ( (getDblMatrixELT(segmat, l, 1) - begin) >= 0. || ( (getDblMatrixELT(segmat, l, 1) - begin) < 0. && (getDblMatrixELT(segmat, l, 1) - begin) >= -1.*eps) );
					isGoodEnd =  ( (getDblMatrixELT(segmat, l, 2) - end) <= 0. || ( (getDblMatrixELT(segmat, l, 2) - end) > 0. && (getDblMatrixELT(segmat, l, 2) - end) <= eps) );

					if (isGoodStart && isGoodEnd) {
						REAL(q01)[l] += forward_rate/((double) nsamples);
						REAL(q10)[l] += reverse_rate/((double) nsamples);
					}
					//check for shift straddle
					if (node == nnode) {
						// there is a shift on the branch, need to see if the dtseg straddles it
						isGoodStart = getDblMatrixELT(segmat, l, 1) < end;
						isGoodEnd = getDblMatrixELT(segmat, l, 2) > end;
						if (isGoodStart && isGoodEnd) {

							// calculate rates to the left of the shift
							leftshift_f = forward_rate * (end - getDblMatrixELT(segmat, l, 1));
							leftshift_r = reverse_rate * (end - getDblMatrixELT(segmat, l, 1)); 
							
							// calculate rates to the right of the shift
							forward_rate = REAL(getListElement(eventData, "q01"))[nxtevent-1];
							reverse_rate = REAL(getListElement(eventData, "q10"))[nxtevent-1];
							rightshift_f = forward_rate * (getDblMatrixELT(segmat, l, 2) - end);
							rightshift_r = reverse_rate * (getDblMatrixELT(segmat, l, 2) - end); 
							
							// weighted average of forward
							ret = (leftshift_f+rightshift_f)/(getDblMatrixELT(segmat, l, 2) - getDblMatrixELT(segmat, l, 1));
							REAL(q01)[l] += ret/((double) nsamples);

							// weighted average of reverse
							ret = (leftshift_r+rightshift_r)/(getDblMatrixELT(segmat, l, 2) - getDblMatrixELT(segmat, l, 1));
							REAL(q10)[l] += ret/((double) nsamples);
							
							place_holder = l; place_holder++;
							break;
						}
					}
				}
				else {
					place_holder = l;
					break;
				}
			}			
		}
		UNPROTECT(2); nprotect -= 2; //protected eventSegs and eventData, which we no longer need
	}
	SEXP res;
	PROTECT(res = allocVector(VECSXP, 2)); nprotect++;
	SET_VECTOR_ELT(res, 0, q01);
	SET_VECTOR_ELT(res, 1, q10);
	UNPROTECT(nprotect);
	return res;
}
