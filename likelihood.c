#include "pastml.h"
#include "scaling.h"

int get_max(const int *array, size_t n) {
    /**
     * Finds the maximum in an array of positive integers.
     */
    int max_value = 0;
    size_t i;
    for (i = 0; i < n; i++) {
        if (max_value < array[i]) {
            max_value = array[i];
        }
    }
    return max_value;
}

void normalize(double *array, size_t n) {
    /**
     * Divides array members by their sum.
     */
    double sum = 0.0;
    size_t i;
    for (i = 0; i < n; i++) {
        sum += array[i];
    }
    for (i = 0; i < n; i++) {
        array[i] /= sum;
    }
}


double get_mu(const double *frequencies, size_t n) {
    /**
     * Calculates the mutation rate for F81 (and JC that is a simplification of it),
     * as \mu = 1 / (1 - sum_i \pi_i^2). This way the overall rate of mutation -\mu trace(\Pi Q) is 1.
     * See [Gascuel "Mathematics of Evolution and Phylogeny" 2005] for further details.
     */
    double sum = 0.0;
    size_t i;
    for (i = 0; i < n; i++) {
        sum += pow(frequencies[i], 2);
    }
    return 1.0 / (1.0 - sum);
}

double get_pij(const double *frequencies, double mu, double t, int i, int j) {
    /**
     * Calculate the probability of substitution i->j over time t, given the mutation rate mu:
     *
     * For K81 (and JC which is a simpler version of it)
     * Pxy(t) = \pi_y (1 - exp(-mu t)) + exp(-mu t), if x == y, \pi_y (1 - exp(-mu t)), otherwise
     * [Gascuel "Mathematics of Evolution and Phylogeny" 2005].
     */
    // if mu == inf (e.g. just one state) and t == 0, we should prioritise mu
    double exp_mu_t = (mu == INFINITY) ? 0.0: exp(-mu * t);

    double p_ij = frequencies[j] * (1.0 - exp_mu_t);
    if (i == j) {
        p_ij += exp_mu_t;
    }
    return p_ij;
}

double get_rescaled_branch_len(const Node *nd, double avg_br_len, double scaling_factor, double epsilon) {
    /**
     * Returns the node branch length multiplied by the scaling factor.
     *
     * For tips, before multiplying, the branch length is adjusted with epsilon:
     * bl = (bl + eps) / (avg_tip_len / (avg_tip_len + eps)).
     * This removes zero tip branches, while keeping the average branch length intact.
     */
    double bl = nd->branch_len;
    if (nd->nb_neigh == 1) { /*a tip*/
        bl = (bl + epsilon) * (avg_br_len / (avg_br_len + epsilon));
    }
    return bl * scaling_factor;
}

void set_p_ij(const Node *nd, double avg_br_len, size_t num_frequencies, const double *parameters, char* model) {
    /**
     * Sets node probabilities of substitution: p[i][j]:
     *
     * For K81 (and JC which is a simpler version of it)
     * Pxy(t) = \pi_y (1 - exp(-mu t)) + exp(-mu t), if x ==y, \pi_y (1 - exp(-mu t)), otherwise
     * [Gascuel "Mathematics of Evolution and Phylogeny" 2005]
     *
     * parameters = [frequency_1, .., frequency_n, scaling_factor, epsilon]
     */
    int i, j;

    double scaling_factor = parameters[num_frequencies];
    double epsilon = parameters[num_frequencies + 1];
    double t = get_rescaled_branch_len(nd, avg_br_len, scaling_factor, epsilon);
    double mu = get_mu(parameters, num_frequencies);

    if ((strcmp(model, JC) == 0) || (strcmp(model, F81) == 0)) {
        for (i = 0; i < num_frequencies; i++) {
            for (j = 0; j < num_frequencies; j++) {
                nd->pij[i][j] = get_pij(parameters, mu, t, i, j);
            }
        }
    }
}

int calculate_node_probabilities(const Node *nd, size_t num_annotations, size_t first_child_index, int is_marginal) {
    int factors = 0;
    double p_child_branch_from_i, cur_child_branch_from_i;
    size_t i, j, k;
    for (k = first_child_index; k < nd->nb_neigh; k++) {
        Node *child = nd->neigh[k];
        for (i = 0; i < num_annotations; i++) {
            /* Calculate the probability of having a branch from the node to its child node,
             * given that the node is in state i: p_child_branch_from_i = sum_j(p_ij * p_child_j)
             */
            p_child_branch_from_i = 0.;
            if (is_marginal) {
                // for marginal sum over all possible child states
                for (j = 0; j < num_annotations; j++) {
                    p_child_branch_from_i += child->pij[i][j] * child->bottom_up_likelihood[j];
                }
            } else {
                // for joint pick the best child state
                for (j = 0; j < num_annotations; j++) {
                    cur_child_branch_from_i = child->pij[i][j] * child->bottom_up_likelihood[j];
                    if (cur_child_branch_from_i > p_child_branch_from_i) {
                        p_child_branch_from_i = cur_child_branch_from_i;
                        child->best_states[i] = j;
                    }
                }
            }

            /* The probability of having the node in state i is a multiplication of
             * the probabilities of p_child_branch_from_i for all child branches.
             */
            if (k == first_child_index) {
                nd->bottom_up_likelihood[i] = p_child_branch_from_i;
            } else {
                nd->bottom_up_likelihood[i] *= p_child_branch_from_i;
            }
        }
        int add_factors = upscale_node_probs(nd->bottom_up_likelihood, num_annotations);
        // if all the probabilities are zero (shown by add_factors == -1), stop here
        if (add_factors == -1) {
            return -1;
        }
        factors += add_factors;
    }
    return factors;
}

int process_node(Node *nd, Tree *s_tree, size_t num_annotations, double *parameters, int is_marginal, char* model) {
    /**
     * Calculates node probabilities.
     * parameters = [frequency_char_1, .., frequency_char_n, scaling_factor, epsilon].
     */
    int factors = 0, add_factors;
    size_t i;
    Node* child;
    size_t first_child_index = (nd == s_tree->root) ? 0 : 1;
    
    /* recursively calculate probabilities for children if there are any children */
    for (i = first_child_index; i < nd->nb_neigh; i++) {
        child = nd->neigh[i];
        
        /* set probabilities of substitution */
        set_p_ij(child, s_tree->avg_tip_branch_len, num_annotations, parameters, model);

        add_factors = process_node(child, s_tree, num_annotations, parameters, is_marginal, model);
        // if all the probabilities are zero (shown by add_factors == -1), stop here
        if (add_factors == -1) {
            return -1;
        }
        factors += add_factors;
    }
    /* calculate own probabilities */
    add_factors = calculate_node_probabilities(nd, num_annotations, first_child_index, is_marginal);
    // if all the probabilities are zero (shown by add_factors == -1), stop here
    if (add_factors == -1) {
        return -1;
    }
    factors += add_factors;
    nd->scaling_factor_down[0] = factors;
    return factors;
}

double calculate_bottom_up_likelihood(Tree *s_tree, size_t num_annotations, double *parameters, int is_marginal,
                                      char* model) {
    /**
     * Calculates tree log likelihood.
     * parameters = [frequency_char_1, .., frequency_char_n, scaling_factor, epsilon].
     */
    double scaled_lk = 0;
    size_t i;

    int factors = process_node(s_tree->root, s_tree, num_annotations, parameters, is_marginal, model);

    /* if factors == -1, it means that the bottom_up_likelihood is 0 */
    if (factors != -1) {
        for (i = 0; i < num_annotations; i++) {
            /* multiply the probability by character frequency */
            scaled_lk += s_tree->root->bottom_up_likelihood[i] * parameters[i];
        }
    }
    return remove_upscaling_factors(log(scaled_lk), factors);
}


void
initialise_tip_probabilities(Tree *s_tree, char *const *tip_names, const int *states, size_t num_tips,
                             size_t num_annotations) {
    /**
     * Sets the state and likelihoods for a tip
     * by setting the likelihood of its real state (given in the metadata file) to 1
     * and the other to 0.
     */
    Node *nd;
    size_t j, i, k;
    double* all_ones_array = calloc(num_annotations, sizeof(double));
    for (j = 0; j < num_annotations; j++) {
       all_ones_array[j] = 1.0;
    }

    for (k = 0; k < s_tree->nb_nodes; k++) {
        nd = s_tree->nodes[k];
        /* if a tip, process it */
        if (nd->nb_neigh == 1) {
            for (i = 0; i < num_tips; i++) {
                if (strcmp(nd->name, tip_names[i]) == 0) {
                    // states[i] == -1 means that the annotation is missing
                    if (states[i] == -1) {
                        // and therefore any state is possible
                        memcpy(nd->bottom_up_likelihood, all_ones_array, num_annotations * sizeof(double));
                    } else {
                        nd->bottom_up_likelihood[states[i]] = 1.0;
                    }
                    break;
                }
            }
        }
    }

    free(all_ones_array);
}

void rescale_branch_lengths(Tree *s_tree, double scaling_factor, double epsilon) {
    /**
     * Rescales all the branches in the tree.
     */
    Node *nd;
    size_t i;
    for (i = 0; i < s_tree->nb_nodes; i++) {
        nd = s_tree->nodes[i];
        nd->branch_len = get_rescaled_branch_len(nd, s_tree->avg_tip_branch_len, scaling_factor, epsilon);
    }
}

void _choose_joint_states(Node *nd, Node *root, size_t best_state_i) {
    /**
     * Chooses best joint states and sets the corresponding result_probs to one.
     */
    size_t i;
    Node* child;

    nd->result_probs[best_state_i] = 1.0;


    /* recursively calculate best states for children if there are any children */
    for (i = (nd == root) ? 0: 1; i < nd->nb_neigh; i++) {
        child = nd->neigh[i];
        _choose_joint_states(child, root, child->best_states[best_state_i]);
    }
}


void choose_joint_states(Tree *s_tree, size_t num_annotations, const double* frequencies) {
    size_t best_root_state_i = 0;
    size_t i;
    double best_p = 0, cur_p;
    for (i = 0; i < num_annotations; i++) {
        cur_p = s_tree->root->bottom_up_likelihood[i] * frequencies[i];
        if (best_p < cur_p) {
            best_p = cur_p;
            best_root_state_i = i;
        }
    }
    _choose_joint_states(s_tree->root, s_tree->root, best_root_state_i);
}


