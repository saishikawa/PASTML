
#include "pastml.h"

void calc_correct(Node *nd, Node *root, int n) {
    /**
     * Chooses an optimal number of non-zero probabilities to keep, and sets all of them to be equal.
     */
    int i, j, best_num_states = n;
    double smallest_correction = INFINITY, correction_i, equal_p_i;

    /* if it's a tip, we have nothing to do */
    if (nd->nneigh == 1) {
        return;
    }

    // process children
    for (i = (nd == root) ? 0: 1; i < nd->nneigh; i++) {
        calc_correct(nd->neigh[i], root, n);
    }

    /* local fraction optimisation */
    for (i = 0; i < n; i++) {
        /* Calculate the squared difference between our probabilities 
         * and choosing (i + 1) states with probability 1 / (i + 1) each:  
         * correction_i = sum_{j <= i} (p_j - 1 / (i + 1))^2 + sum_{j > i} p_j^2 */
        correction_i = 0.;
        equal_p_i = 1.0 / ((double) i + 1.0);
        for (j = 0; j < n; j++) {
            if (j <= i) {
                correction_i += pow(nd->marginal[j] - equal_p_i, 2);
            } else {
                correction_i += pow(nd->marginal[j], 2);
            }
        }
        if (smallest_correction > correction_i) {
            smallest_correction = correction_i;
            best_num_states = i + 1;
        }
    }
    equal_p_i = 1.0 / ((double) best_num_states);
    for (i = 0; i < n; i++) {
        nd->marginal[i] = (i < best_num_states) ? equal_p_i: 0.0;
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

    /* put index array in the nd->best_states
     * and sort it by marginal probabilities in marginal array, in decreasing order */
    for(i = 0; i < num_annotations; i++){
        nd->best_states[i] = i;
    }
    int cmp(const void *a, const void *b){
        int ia = *(int *)a;
        int ib = *(int *)b;
        return -(nd->condlike_mar[ia] < nd->condlike_mar[ib] ? -1 : nd->condlike_mar[ia] > nd->condlike_mar[ib]);
    }
    qsort(nd->best_states, num_annotations, sizeof(int), cmp);
    for (i = 0; i < num_annotations; i++) {
        nd->marginal[i] = nd->condlike_mar[nd->best_states[i]];
    }

    // process children
    for (i = (nd == root) ? 0: 1; i < nd->nneigh; i++) {
        order_marginal(nd->neigh[i], root, num_annotations);
    }
}
