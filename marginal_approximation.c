
#include "pastml.h"



void
order_marginal(Tree* tree, int num_annotations)
{
    /**
     * Recursively sets node's marginal field to marginal state probabilities in decreasing order,
     * and node's tmp_best field to the corresponding state indices.
     */
    int i;
    Node* nd;
    int* indices = malloc(num_annotations * sizeof(int));
    for (i = 0; i < num_annotations; i++) {
        indices[i] = i;
    }

    for (int k = 0; k < tree->nb_nodes; k++) {
        nd = tree->nodes[k];
        /* put index array in the nd->best_states
         * and sort it by marginal probabilities in marginal array, in decreasing order */
        memcpy(nd->best_states, indices, num_annotations * sizeof(int));
        int cmp_index(const void *a, const void *b) {
            int ia = *(int *) a;
            int ib = *(int *) b;
            return -(nd->marginal[ia] < nd->marginal[ib] ? -1 : nd->marginal[ia] > nd->marginal[ib]);
        }
        qsort(nd->best_states, num_annotations, sizeof(int), cmp_index);
        int cmp_value(const void *a, const void *b) {
            double da = *(double *) a;
            double db = *(double *) b;
            return -(da < db ? -1 : da > db);
        }
        qsort(nd->marginal, num_annotations, sizeof(double), cmp_value);
    }
    free(indices);
}

void calc_correct(Tree *tree, int n) {
    /**
     * Chooses an optimal number of non-zero probabilities to keep, and sets all of them to be equal.
     */

    int i, j, best_num_states;
    double smallest_correction, correction_i, equal_p_i;
    Node* nd;

    for (int k = 0; k < tree->nb_nodes; k++) {
        nd = tree->nodes[k];
        smallest_correction = INFINITY;
        best_num_states = n;
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
            nd->marginal[i] = (i < best_num_states) ? equal_p_i : 0.0;
        }
    }
}

void choose_likely_states(Tree *tree, int n) {
    /**
     * Chooses an optimal number of non-zero probabilities to keep, and sets all of them to be equal.
     */

    // order marginal probabilities
    order_marginal(tree, n);

    //choose the most likely states to keep among those with non-zero probabilities
    calc_correct(tree, n);
}
