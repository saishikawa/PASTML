#include <assert.h>
#include "pastml.h"
#include "scaling.h"
#include "tree.h"
#include "logger.h"
#include "likelihood.h"


void _calculate_top_down_likelihood(const Node *nd, const Node *root, size_t num_frequencies) {
    /**
     * The up-likelihood of a given node is computed based on the information
     * coming from all the tips that are not descending from the studied node
     * and contained in the “up-subtree” of that node.
     *
     * To calculate it we assume that the tree is rooted in this node and combine the likelihoods of the “up-subtrees”,
     * e.g. to calculate the top-down likelihood of a node N1 being in a state i,
     * given that its parent node is P and its brother node is N2, we imagine that the tree is re-rooted in N1,
     * therefore P becoming the child of N1, and N2 its grandchild,
     * and then calculate the bottom-up likelihood from the P subtree:
     * L_top_down(N1, i) = \sum_j P(i -> j, dist(N1, P)) * L_top_down(P) * \sum_k P(j -> k, dist(N2, P)) * L_bottom_up (N2).
     *
     * For the root node we assume its top-down likelihood to be 1 for all the states.
     */
    Node *parent;
    Node *brother;
    int *parent_scaling_factors = calloc(num_frequencies, sizeof(int));
    int *scaling_factors = calloc(num_frequencies, sizeof(int));
    double prob_parent[num_frequencies];
    size_t child_id, j, i, k;

    if (nd == root) {
        for (i = 0; i < num_frequencies; i++) {
            root->top_down_likelihood[i] = 1.0;
            root->scaling_factor_up[0] = 0;
        }
    } else {
        parent = nd->neigh[0];
        nd->scaling_factor_up[0] = parent->scaling_factor_up[0];
        // current node has annotation i
        for (i = 0; i < num_frequencies; i++) {
            scaling_factors[i] = 0;
            nd->top_down_likelihood[i] = 0.0;
            /* we need to combine
             * the probability of the of a branch between our node and its parent to evolve from i to j,
             * with the up probability of our parent being in a state j,
             * with the probabilities of all its children (excluding nd) to evolve from j:
             * L_top_down(nd, i) = \sum_j P(i -> j, dist(nd, parent)) * L_top_down(parent, j) *
             * * \sum_k P(j -> k, dist(brother_1, parent)) * L_bottom_up (brother_1, k) ...
             */

            for (j = 0; j < num_frequencies; j++) {
                parent_scaling_factors[j] = 0;
                prob_parent[j] = nd->pij[i][j] * parent->top_down_likelihood[j];
                // if our parent is not root, its first nb_neigh is our grandfather,
                // and we should iterate over children staring from 1
                for (child_id = parent == root ? 0 : 1; child_id < parent->nb_neigh; child_id++) {
                    brother = parent->neigh[child_id];
                    if (brother == nd) {
                        continue;
                    }
                    double brother_prob = 0.0;
                    for (k = 0; k < num_frequencies; k++) {
                        brother_prob += brother->pij[j][k] * brother->bottom_up_likelihood[k];
                    }
                    prob_parent[j] *= brother_prob;
                    parent_scaling_factors[j] += rescale_if_needed(prob_parent, j);
                }
            }
            scaling_factors[i] = harmonise_scaling(prob_parent, parent_scaling_factors, num_frequencies);
            for (j = 0; j < num_frequencies; j++) {
                nd->top_down_likelihood[i] += prob_parent[j];
            }
        }
        nd->scaling_factor_up[0] += harmonise_scaling(nd->top_down_likelihood, scaling_factors, num_frequencies);
        for (child_id = parent == root ? 0 : 1; child_id < parent->nb_neigh; child_id++) {
            brother = parent->neigh[child_id];
            if (brother == nd) {
                continue;
            }
            nd->scaling_factor_up[0] += brother->scaling_factor_down[0];
        }
    }
    free(parent_scaling_factors);
    free(scaling_factors);

    // recursively calculate top down probabilities for the children
    for (i = (nd == root) ? 0 : 1; i < nd->nb_neigh; i++) {
        _calculate_top_down_likelihood(nd->neigh[i], root, num_frequencies);
    }
}

void calculate_top_down_likelihood(Tree *s_tree, size_t num_annotations) {
    _calculate_top_down_likelihood(s_tree->root, s_tree->root, num_annotations);
}

void _calculate_node_marginal_probabilities(Node *nd, Node *root, size_t num_annotations, const double *frequency) {
    /**
     * The marginal likelihood of a certain state can be computed by multiplying its up- and down-likelihoods,
     * and its prior probability (frequency).
     */
    size_t i;

    // Finally, the marginal likelihood of a certain state can be computed
    // by multiplying its up-, down-likelihoods, and its frequency.
    for (i = 0; i < num_annotations; i++) {
        nd->result_probs[i] = nd->top_down_likelihood[i] * nd->bottom_up_likelihood[i] * frequency[i];
    }
    // recursively calculate marginal probabilities for the children
    for (i = (nd == root) ? 0 : 1; i < nd->nb_neigh; i++) {
        _calculate_node_marginal_probabilities(nd->neigh[i], root, num_annotations, frequency);
    }
}

double get_marginal_likelihood(Node* nd, size_t num_annotations) {
    double lk = 0.0;
    size_t i;

    for (i = 0; i < num_annotations; i++) {
        lk += nd->result_probs[i];
    }
    lk = remove_upscaling_factors(log(lk), nd->scaling_factor_down[0] + nd->scaling_factor_up[0]);
    return round(lk * 10.0) / 10.0;
}

void check_marginal_probabilities(Tree *s_tree, size_t num_annotations) {
    /**
     * Sanity check: marginal likelihoods of all the nodes should be the same
     */
    size_t i;
    double node_lk, parent_lk;
    Node *node, *parent;

    for (i = 0; i < s_tree->nb_nodes; i++) {
        node = s_tree->nodes[i];
        if (node != s_tree->root && !isProblematic(node)) {
            parent = node->neigh[0];
            node_lk = get_marginal_likelihood(node, num_annotations);
            parent_lk = get_marginal_likelihood(parent, num_annotations);
            if (node_lk != parent_lk) {
                log_info("LH: %.2f\tvs\t%.2f (%s\tvs\t%s).\n", node_lk, parent_lk, node->name, parent->name);
            }
            assert(node_lk == parent_lk);
        }        
    }
}

void calculate_marginal_probabilities(Tree *s_tree, size_t num_annotations, double *frequencies) {
    /**
     * Calculates marginal probabilities of tree nodes.
     */
    _calculate_node_marginal_probabilities(s_tree->root, s_tree->root, num_annotations, frequencies);
//    check_marginal_probabilities(s_tree, num_annotations);
}


