
#include "pastml.h"
#include "likelihood.h"

Node* nd;

int cmp_index(const void *a, const void *b) {
    size_t ia = *(size_t *) a;
    size_t ib = *(size_t *) b;
    return -(nd->result_probs[ia] < nd->result_probs[ib] ? -1 : nd->result_probs[ia] > nd->result_probs[ib]);
}

void
order_marginal(Tree* tree, size_t num_annotations)
{
    /**
     * Recursively sets node's marginal field to marginal state probabilities in decreasing order,
     * and node's tmp_best field to the corresponding state indices.
     */
    size_t i;

    for (i = 0; i < tree->nb_nodes; i++) {
        nd = tree->nodes[i];
        qsort(nd->best_states, num_annotations, sizeof(size_t), cmp_index);
    }
}

void calc_correct(Tree *tree, size_t n) {
    /**
     * Chooses an optimal number of non-zero probabilities to keep, and sets all of them to be equal.
     */

    size_t i, j, k, best_num_states;
    double smallest_correction, correction_i, equal_p_i;
    Node* nd;

    for (k = 0; k < tree->nb_nodes; k++) {
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
                    correction_i += pow(nd->result_probs[nd->best_states[j]] - equal_p_i, 2);
                } else {
                    correction_i += pow(nd->result_probs[nd->best_states[j]], 2);
                }
            }
            if (smallest_correction > correction_i) {
                smallest_correction = correction_i;
                best_num_states = i + 1;
            }
        }

        equal_p_i = 1.0 / ((double) best_num_states);
        for (i = 0; i < n; i++) {
            nd->result_probs[nd->best_states[i]] = (i < best_num_states) ? equal_p_i : 0.0;
        }
    }
}

void normalize_result_probabilities(Tree *tree, size_t n) {
    size_t i;
    for (i = 0; i < tree->nb_nodes; i++) {
        normalize(tree->nodes[i]->result_probs, n);
    }
}

void set_id_best_states(Tree *tree, size_t n) {
    size_t * indices = malloc(n * sizeof(size_t));
    size_t i;
    for (i = 0; i < n; i++) {
        indices[i] = i;
    }
    for (i = 0; i < tree->nb_nodes; i++) {
        memcpy(tree->nodes[i]->best_states, indices, n * sizeof(size_t));
    }
    free(indices);
}

void choose_likely_states(Tree *tree, size_t n) {
    /**
     * Chooses an optimal number of non-zero probabilities to keep, and sets all of them to be equal.
     */
    // order marginal probabilities
    order_marginal(tree, n);

    //choose the most likely states to keep among those with non-zero probabilities
    calc_correct(tree, n);
}


void choose_best_marginal_states(Tree *tree, size_t n) {
    /**
     * Chooses the state with the highest marginal probability for each node,
     * sets its result_probs to 1, and the others to 0.
     */
    size_t i, j;
    Node* nd;
    size_t best_state_i;
    double best_p, cur_p;

    for (i = 0; i < tree->nb_nodes; i++) {
        nd = tree->nodes[i];
        best_state_i = 0;
        best_p = 0.0;
        for (j = 0; j < n; j++) {
            cur_p = nd->result_probs[j];
            // set current prob to 0 (if it is the best we'll reset it to 1 below)
            nd->result_probs[j] = 0.0;
            if (best_p < cur_p) {
                best_p = cur_p;
                // set previous best to 0
                nd->result_probs[best_state_i] = 0.0;
                best_state_i = j;
                // set new best to 1
                nd->result_probs[best_state_i] = 1.0;
            }
        }
    }
}
