#include <errno.h>
#include "pastml.h"
#include "output_states.h"
#include "output_tree.h"

void calc_correct(Node *nd, Node *root, int nbanno) {
    int i, j, tmpnum;
    double tmp_correct = 10.0, local_cor;

    /* if it's a tip, we have nothing to do */
    if (nd->nneigh == 1) {
        return;
    }

    // process children
    for (i = (nd == root) ? 0: 1; i < nd->nneigh; i++) {
        calc_correct(nd->neigh[i], root, nbanno);
    }

    /* local_fraction_optimize */
    for (i = 0; i < nbanno; i++) {
        local_cor = 0.;
        for (j = 0; j < nbanno; j++) {
            if (j <= i) {
                local_cor += pow(nd->marginal[j] - 1.0 / ((double) i + 1.0), 2);
            } else {
                local_cor += pow(nd->marginal[j], 2);
            }
        }
        if (tmp_correct > local_cor) {
            tmp_correct = local_cor;
            tmpnum = i + 1;
        }
    }
    for (i = 0; i < nbanno; i++) {
        nd->local_flag[i] = (i < tmpnum) ? 1: 0;
    }
}



void
order_marginal(Node *nd, Node *root, int num_annotations)
{
    /**
     * Recursively sets node's marginal field to marginal state probabilities in decreasing order,
     * and node's tmp_best field to the corresponding state indices.
     */
    int i;

    /* if it's a tip, we have nothing to do */
    if (nd->nneigh == 1) {
        return;
    }

    double* marginal = (nd == root) ? nd->mar_prob: nd->condlike_mar;

    // sort index array by marginal probabilities in marginal array, in decreasing order
    int* index = calloc(num_annotations, sizeof(int));
    for(i = 0; i < num_annotations; i++){
        index[i] = i;
    }
    int cmp(const void *a, const void *b){
        int ia = *(int *)a;
        int ib = *(int *)b;
        return -(marginal[ia] < marginal[ib] ? -1 : marginal[ia] > marginal[ib]);
    }
    qsort(index, num_annotations, sizeof(int), cmp);
    for (i = 0; i < num_annotations; i++) {
        nd->marginal[i] = marginal[index[i]];
        nd->tmp_best[i] = index[i];
    }
    free(index);

    // process children
    for (i = (nd == root) ? 0: 1; i < nd->nneigh; i++) {
        order_marginal(nd->neigh[i], root, num_annotations);
    }
}
