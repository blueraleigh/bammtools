#include <R.h>
#include <Rinternals.h>


SEXP get_startstoptimes(SEXP edge, SEXP edgelen, SEXP ntip);

/*
** Return beginning and ending times for each branch of phylogeny.
**
** `edge` is the edge matrix of ape's phylo class
** `edgelen` is the edge.length vector of ape's phylo class
** `ntip` is the number of terminal nodes in the phylogeny
*/
SEXP get_startstoptimes(SEXP edge, SEXP edgelen, SEXP ntip) {
    
    int nbranches = 2 * INTEGER(ntip)[0] - 2;  // root branch not included
    int nprotect = 0;

    // beginning and ending times for branches in phylogeny
    SEXP start_times;
    SEXP stop_times;
    PROTECT(start_times = allocVector(REALSXP, nbranches+1));
    nprotect++;
    PROTECT(stop_times = allocVector(REALSXP, nbranches+1));
    nprotect++;

    // initialize root. it has ape-index ntip + 1, c-index ntip
    REAL(start_times)[INTEGER(ntip)[0]] = 0.0;
    REAL(stop_times)[INTEGER(ntip)[0]] = 0.0;
    for (int i = 0; i < nbranches; i++) {
        int nodeix = INTEGER(edge)[i + 1*nbranches] - 1;  // get the i-th row, 2nd column
        int parentnodeix = INTEGER(edge)[i + 0*nbranches] - 1;  // get the i-th row, 1st column
        double el = REAL(edgelen)[i];
        REAL(start_times)[nodeix] = REAL(stop_times)[parentnodeix];
        REAL(stop_times)[nodeix] = REAL(stop_times)[parentnodeix] + el;
    }

    // return list
    SEXP res;
    PROTECT(res = allocVector(VECSXP, 2));
    nprotect++;
    SET_VECTOR_ELT(res, 0, start_times);
    SET_VECTOR_ELT(res, 1, stop_times);
    UNPROTECT(nprotect);
    return res;
}
