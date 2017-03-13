#include <R.h>

void get_mrca(int* anc, int* desc, int* root, int* ne, int* npair, int* t1, int* t2, int* ret);

/*
** anc is phy$edge[,1]
** desc is phy$edge[,2]
** root is Ntip(phy) + 1
** ne is 2*Ntip(phy) - 2
** npair is length(t1)
** t1 and t2 are node pairs to compute mrcas on
** ret is the result
*/
void get_mrca(int* anc, int* desc, int* root, int* ne, int* npair, int* t1, int* t2, int* ret) {   
    int i,j,k;

    int cnt, node, mrca;
    int* path;
    
    for (k = 0; k < *npair; k++) {
        if (t2[k] == 0) {  // a 0 means no other node was passed as part of the pair
            ret[k] = t1[k];
            continue;
        }
        path = Calloc(*ne, int);
        cnt = 0; 
        mrca = 0;
        node = t1[k];
        // path to root for first node
        while (node != *root) {
            for (i = 0; i < *ne; i++) {
                if (desc[i] == node) {
                    node = anc[i];
                    path[cnt] = node;
                    cnt++;
                    break;
                }
            }
        }
        
        // path to root for second node, interrupted if intersects first path
        node = t2[k];
        while (node != *root) {
            for (i = 0; i < *ne; i++) {
                if (desc[i] == node) {
                    node = anc[i];
                    // check if ancestor in path of first node
                    for (j = 0; j < *ne; j++) {
                        if (node == path[j]) {
                            mrca = 1;
                            break;
                        }
                    }
                }
                if (mrca == 1) 
                    break;
            }
            if (mrca == 1) 
                break;
        }
        if (mrca == 1) {
            ret[k] = node;
        } else {
            ret[k] = *root;
        }
        Free(path);
    }
}
