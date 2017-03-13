#include <R.h>
#include <Rinternals.h>


SEXP get_eventdata(SEXP x, SEXP nrowx, SEXP ncolx, SEXP edge, SEXP edgelen, SEXP tiplabels, SEXP ntip, 
    SEXP begin, SEXP end, SEXP downseq, SEXP lastvisit);

/*
** `x` is a dataframe holding one posterior sample
** `edge` is the edge matrix from ape's phylo class
** `edgelen` is the edge length vector from ape's phylo class
** `tiplabels` is the tip label vector from ape's phylo class
** `ntip` is the number of tips in the phylogeny
** `begin` is the vector of branch start times
** `end` is the vector of branch end times
** `downseq` is the vector of pre-order traversal node visitation indices
** `lastvisit` is the vector that stores the node-indices of last visits
*/

SEXP get_eventdata(SEXP x, SEXP nrowx, SEXP ncolx, SEXP edge, SEXP edgelen, SEXP tiplabels, SEXP ntip, 
    SEXP begin, SEXP end, SEXP downseq, SEXP lastvisit) 
{    
    
}



