#include <R.h>
#include <Rinternals.h>

void get_recursivesequence(int* anc, int* desc, int* node, int*pos, int* ne, int* downseq, int* lastvisit);


/*
** `anc` is the first column of the edge matrix in ape's phylo class
** `desc` is the second column of the edge matrix in ape's phylo class
** `node` is the focal node
** `pos` is the `node`'s visitation order
** `ne` is the number of edges in the phylogeny
**
** from R perspective:
** `downseq` stores visitation order. e.g. downseq[1] is ape-index of first node-visited
** `lastvisit` stores last node visited. e.g. lastvisit[40] returns ape-index of node last visited by node with ape-index 40
*/
void get_recursivesequence(int* anc, int* desc, int* node, int* pos, int* ne, int* downseq, int* lastvisit) {

    int i, d = 0;
    int* children;
    children = Calloc(2, int);

    // add node at its visitation position
    downseq[*pos] = *node;

    // find `node` in the edge matrix
    for (i = 0; i < *ne; i++) {
        if (anc[i] == *node) {
            children[d] = desc[i];
            d++;
        }
        if (d == 2) {
            // break out if we have both descendants
            break;
        }
    }
    if (children[0] != 0 && children[1] != 0) {
        int* child;
        child = Calloc(1, int);
        for (i = 0; i < 2; i++) {
            *child = children[i];
            (*pos)++;
            get_recursivesequence(anc, desc, child, pos, ne, downseq, lastvisit);
        }
        Free(child);
    }
    lastvisit[*node-1] = downseq[*pos];  
    Free(children); 
}
