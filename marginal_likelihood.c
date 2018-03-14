#include "pastml.h"
#include "scaling.h"
#include "likelihood.h"

extern SIMULATION;

int *calculate_top_down_likelihoods(const Node *nd, const Node *root, size_t num_annotations, double *frequencies) {
    /**
     * The up-likelihood of a given node is computed based on the information
     * coming from all the tips that are not descending from the studied node
     * and contained in the “up-subtree” of that node.
     *
     *
     * For example, in a rooted tree ((T1, (T2, T3)N2)N1, (T4, T5)N3)ROOT;
     * to compute the up-likelihood of N1 = x,
     * we consider substitutions from its parent node N3
     * (note that we ignore the root and consider the tree as unrooted).
     * We then multiply each down likelihoods of N3 = y
     * by the probability of substitution from y to x during the time dist(N1) + dist(N3), and sum them.
     *
     * Lup(N1=x|D, params) = \sum_y( P(y->x, dist(N1) + dist(N3)) Ldown(N3=y) )
     *
     * The up-likelihood of N2 having a certain state x can then be calculated
     * by summing up-likelihoods of N1 having state y,
     * multiplying each of them by the probability of substitutions among N1 and T1 (y → tip_state(T1), in branch(T3)),
     * and N1 and N2 (y → x, in branch(N2)).
     *
     * We calculate up-likelihoods for N3 in same manner.
     */
    Node *father = nd->neigh[0];
    Node *other_child;
    int my_id = -1;
    int father_scaling_factors[num_annotations];
    int *scaling_factors = calloc(num_annotations, sizeof(int));
    double prob_father[num_annotations];
    double mu = get_mu(frequencies, num_annotations);
    int child_id, j, i, k;

    for (child_id = 0; child_id < father->nb_neigh; child_id++) {
        if (father->neigh[child_id] == nd) {
            my_id = child_id;
        }
    }
    // current node has annotation i
    for (i = 0; i < num_annotations; i++) {
        scaling_factors[i] = 0;

        if (father == root) {
            nd->top_down_likelihood[i] = 1.0;
            /* as our tree is rooted, there will be just one other child */
            for (child_id = 0; child_id < father->nb_neigh; child_id++) {
                if (child_id != my_id) {
                    other_child = father->neigh[child_id];
                    /* we ignore the root and consider the tree as unrooted,
                     * therefore the other_child becomes the parent of our nd:
                     * L_up(nd=i|D, params) = \sum_j( L_down(other_child=j) P(j->i, dist(nd) + dist(other_child)) )
                     */
                    double prob_up_i = 0.0;
                    for (j = 0; j < num_annotations; j++) {
                        prob_up_i += other_child->bottom_up_likelihood[j]
                                     * get_pij(frequencies, mu, nd->branch_len + other_child->branch_len, j, i);
                    }
                    nd->top_down_likelihood[i] *= prob_up_i;
                }
            }
        } else {
            nd->top_down_likelihood[i] = 0.0;
            /* we need to combine the up probability of our parent being in a state j
             * with the probabilities of all its children (including nd) evolving from j.
             * L_up(nd=i|D, params) =
             * \sum_j( L_up(father=j) P(j->i, dist(nd)) P(j->state(other_child_1), dist(other_child_1) ...) )
             */
            for (j = 0; j < num_annotations; j++) {
                prob_father[j] = nd->pij[j][i] * father->top_down_likelihood[j];
                father_scaling_factors[j] = 0;
                // as our father is not root, its first nb_neigh is our grandfather,
                // and we should iterate over children staring from 1
                for (child_id = 1; child_id < father->nb_neigh; child_id++) {
                    if (child_id != my_id) {
                        other_child = father->neigh[child_id];
                        double other_child_prob = 0.0;
                        for (k = 0; k < num_annotations; k++) {
                            other_child_prob += other_child->pij[j][k] * other_child->bottom_up_likelihood[k];
                        }
                        prob_father[j] *= other_child_prob;
                    }
                    if (prob_father[j] < LIM_P) {
                        int curr_scaler_pow1 = get_scaling_pow(prob_father[j]);
                        father_scaling_factors[j] += curr_scaler_pow1;
                        rescale(prob_father, j, curr_scaler_pow1);
                    }
                }
            }
            int max_father_factor = get_max(father_scaling_factors, num_annotations);
            for (j = 0; j < num_annotations; j++) {
                int curr_scaler_pow = max_father_factor - father_scaling_factors[j];
                if (curr_scaler_pow != 0) {
                    rescale(prob_father, j, curr_scaler_pow);
                }
                nd->top_down_likelihood[i] += prob_father[j];
            }
            scaling_factors[i] = max_father_factor;
        }
    }
    return scaling_factors;
}

void calculate_node_marginal_probabilities(Node *nd, Node *root, size_t num_annotations, double *frequency) {
    /**
     * The up-likelihood of a given node is computed based on the information
     * coming from all the tips that are not descending from the studied node
     * and contained in the “up-subtree” of that node.
     *
     *
     * For example, in a rooted tree ((T1, (T2, T3)N2)N1, (T4, T5)N3)ROOT;
     * to compute the up-likelihood of N1 = x,
     * we consider substitutions from its parent node N3
     * (note that we ignore the root and consider the tree as unrooted).
     * We then multiply each down likelihoods of N3 = y
     * by the probability of substitution from y to x during the time dist(N1) + dist(N3), and sum them.
     *
     * The up-likelihood of N2 having a certain state x can then be calculated
     * by summing up-likelihoods of N1 having state y,
     * multiplying each of them by the probability of substitutions among N1 and T1 (y → tip_state(T1), in branch(T3)),
     * and N1 and N2 (y → x, in branch(N2)).
     *
     * We calculate up-likelihoods for N3 in same manner.
     *
     * Finally, the marginal likelihood of a certain state can be computed by multiplying its up- and down-likelihoods,
     * and its prior probability (frequency).
     */
    int i;
    int curr_scaler_pow, max_factor;

    if (nd == root) {
        memcpy((void *) nd->marginal, (void *) nd->bottom_up_likelihood, num_annotations * sizeof(double));
    } else {
        int *tmp_factor = calculate_top_down_likelihoods(nd, root, num_annotations, frequency);
        max_factor = get_max(tmp_factor, num_annotations);
        for (i = 0; i < num_annotations; i++) {
            curr_scaler_pow = max_factor - tmp_factor[i];
            if (curr_scaler_pow != 0) {
                rescale(nd->top_down_likelihood, i, curr_scaler_pow);
            }
        }
        free(tmp_factor);

        // Finally, the marginal likelihood of a certain state can be computed
        // by multiplying its up-, down-likelihoods, and its frequency.
        for (i = 0; i < num_annotations; i++) {
            nd->marginal[i] = nd->top_down_likelihood[i] * nd->bottom_up_likelihood[i] * frequency[i];
        }
    }
    normalize(nd->marginal, num_annotations);
    if(SIMULATION == TRUE) {
        for (i = 0; i < num_annotations; i++) {
            nd->sim_marginal_prob[i] = nd->marginal[i];
        }
    }

    // recursively calculate marginal probabilities for the children
    for (i = (nd == root) ? 0 : 1; i < nd->nb_neigh; i++) {
        calculate_node_marginal_probabilities(nd->neigh[i], root, num_annotations, frequency);
    }
}

void calculate_marginal_probabilities(Tree *s_tree, size_t num_annotations, double *frequency) {
    /**
     * Calculates marginal probabilities of tree nodes.
     */
    calculate_node_marginal_probabilities(s_tree->root, s_tree->root, num_annotations, frequency);
}


