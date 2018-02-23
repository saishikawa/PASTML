#include "pastml.h"
#include "scaling.h"

extern Tree *s_tree;

int get_max(int *array, int n) {
    /**
     * Finds the maximum in an array of positive integers.
     */
    int max_value = 0;
    for (int i = 0; i < n; i++) {
        if (max_value < array[i]) {
            max_value = array[i];
        }
    }
    return max_value;
}

void normalize(double *array, int n) {
    /**
     * Divides array members by their sum.
     */
    double sum = 0.0;;
    for (int i = 0; i < n; i++) {
        sum += array[i];
    }
    for (int i = 0; i < n; i++) {
        array[i] /= sum;
    }
}

void set_p_ij(const Node *nd, int num_frequencies, const double *parameters);

int calculate_node_probabilities(const Node *nd, int num_annotations, int first_child_index);

double remove_upscaling_factors(double log_likelihood, int factors);

int process_node(Node *nd, Node* root, int num_annotations, double *parameters) {
    /**
     * Calculates node probabilities.
     * parameters = [frequency_char_1, .., frequency_char_n, scaling_factor, epsilon].
     */
    int factors = 0, add_factors;

    /* set probabilities of substitution */
    if (nd != root) {
        set_p_ij(nd, num_annotations, parameters);
    }

    /* not a tip */
    if (nd->nneigh != 1) {
        int first_child_index = (nd == root) ? 0 : 1;
        /* recursively calculate probabilities for children */
        for (int i = first_child_index; i < nd->nneigh; i++) {
            add_factors = process_node(nd->neigh[i], root, num_annotations, parameters);
            /* if all the probabilities are zero (shown by add_factors == -1),
             * there is no point to go any further
             */
            if (add_factors == -1) {
                return -1;
            }
            factors += add_factors;
        }
        /* calculate own probabilities */
        add_factors = calculate_node_probabilities(nd, num_annotations, first_child_index);
        /* if all the probabilities are zero (shown by add_factors == -1),
         * there is no point to go any further
         */
        if (add_factors == -1) {
            return -1;
        }
        factors += add_factors;
    }

    if (nd == root) {
        for (int i = 0; i < num_annotations; i++) {
            /* multiply the probability by character frequency */
            nd->bottom_up_likelihood[i] = nd->bottom_up_likelihood[i] * parameters[i];
        }
    }

    return factors;
}


double get_mu(const double* frequencies, int n) {
    /**
     * Calculates the mutation rate for F81 (and JC that is a simplification of it),
     * as \mu = 1 / (1 - sum_i \pi_i^2). This way the overall rate of mutation -\mu trace(\Pi Q) is 1.
     * See [Gascuel "Mathematics of Evolution and Phylogeny" 2005] for further details.
     */
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += pow(frequencies[i], 2);
    }
    return 1.0 / (1.0 - sum);
}

double calc_lik_bfgs(Node *root, int num_annotations, double *parameters) {
    /**
     * Calculates tree log likelihood.
     * parameters = [frequency_char_1, .., frequency_char_n, scaling_factor, epsilon].
     * mu = 1 / (1 - (frequency_1^2 + ... + frequency_n^2)).
     */
    double scaled_lk = 0;

    int factors = process_node(root, root, num_annotations, parameters);
    /* if factors == -1, it means that the bottom_up_likelihood is 0 */
    if (factors != -1) {
        for (int i = 0; i < num_annotations; i++) {
            scaled_lk += root->bottom_up_likelihood[i];
        }
    }
    return remove_upscaling_factors(log(scaled_lk), factors);
}


double remove_upscaling_factors(double log_likelihood, int factors) {
    /**
     * While calculating the node probabilities, the upscaling was done to avoid underflow problems:
     * if a certain node probability was too small, we multiplied this node probabilities by a scaling factor,
     * and kept the factor in mind to remove it from the final likelihood.
     *
     * Now its time to remove all the factors from the final likelihood.
     */
    int piecewise_scaler_pow;
    do {
        piecewise_scaler_pow = MIN(factors, 63);
        log_likelihood -= LOG2 * piecewise_scaler_pow;
        factors -= piecewise_scaler_pow;
    } while (factors != 0);
    return log_likelihood;
}

int calculate_node_probabilities(const Node *nd, int num_annotations, int first_child_index) {
    int factors = 0;
    for (int ii = first_child_index; ii < nd->nneigh; ii++) {
        Node *child = nd->neigh[ii];
        for (int i = 0; i < num_annotations; i++) {
            /* Calculate the probability of having a branch from the node to its child node,
             * given that the node is in state i: p_child_branch_from_i = sum_j(p_ij * p_child_j)
             */
            double p_child_branch_from_i = 0.;
            for (int j = 0; j < num_annotations; j++) {
                p_child_branch_from_i += child->pij[i][j] * child->bottom_up_likelihood[j];
            }

            /* The probability of having the node in state i is a multiplication of
             * the probabilities of p_child_branch_from_i for all child branches:
             * condlike_i = mult_ii(p_child_ii_branch_from_i)
             */
            if (ii == first_child_index) {
                nd->bottom_up_likelihood[i] = p_child_branch_from_i;
            } else {
                nd->bottom_up_likelihood[i] *= p_child_branch_from_i;
            }
        }
        int add_factors = upscale_node_probs(nd->bottom_up_likelihood, num_annotations);

        /* if all the probabilities are zero (shown by add_factors == -1),
         * there is no point to go any further
         */
        if (add_factors == -1) {
            fprintf(stderr, "Likelihood of one of the nodes is already 0, no need to go further.\n");
            return -1;
        }
        factors += add_factors;
    }
    return factors;
}

void
initialise_tip_probabilities(Node *nd, Node *root,
                             char *const *tipnames, const int *states, int num_tips, int num_annotations) {
    /**
     * Sets the state and likelihoods for a tip
     * by setting the likelihood of its real state (given in the metadata file) to 1
     * and the other to 0.
     */
    int i, j;

    /* if a tip, process it */
    if (nd->nneigh == 1) {
        for (i = 0; i < num_tips; i++) {
            if (strcmp(nd->name, tipnames[i]) == 0) {
                nd->fixed_state = states[i];

                // states[i] == num_annotations means that the annotation is missing
                if (states[i] == num_annotations) {
                    // and therefore any state is possible
                    for (j = 0; j < num_annotations; j++) {
                        nd->bottom_up_likelihood[j] = 1.0;
                    }
                } else {
                    nd->bottom_up_likelihood[states[i]] = 1.0;
                }
                break;
            }
        }
    /* if not a tip, call recursively on children */
    } else {
        if (nd->nneigh != 1) {
            for (i = (nd == root) ? 0 : 1; i < nd->nneigh; i++) {
                initialise_tip_probabilities(nd->neigh[i], root, tipnames, states, num_tips, num_annotations);
            }
        }
    }
}

double get_rescaled_branch_len(const Node *nd, double scaling_factor, double epsilon) {
    /**
     * Returns the node branch length multiplied by the scaling factor.
     *
     * For tips, before multiplying, the branch length is adjusted with epsilon:
     * bl = (bl + eps) / (avg_tip_len / (avg_tip_len + eps)).
     * This removes zero tip branches, while keeping the average branch length intact.
     */
    double bl = nd->brlen;
    if (nd->nneigh == 1) { /*a tip*/
        bl = (bl + epsilon) * (s_tree->tip_avg_branch_len / (s_tree->tip_avg_branch_len + epsilon));
    }
    return bl * scaling_factor;
}

void rescale_branch_lengths(Node *nd, Node *root, double scaling_factor, double epsilon) {
    /**
     * Rescales all the branches in the tree.
     */
    nd->brlen = get_rescaled_branch_len(nd, scaling_factor, epsilon);
    /* call recursively on children */
    for (int i = (nd == root) ? 0 : 1; i < nd->nneigh; i++) {
        rescale_branch_lengths(nd->neigh[i], root, scaling_factor, epsilon);
    }
}


double get_pij(const double *frequencies, double mu, double t, int i, int j) {
    /**
     * Calculate the probability of substitution i->j over time t, given the mutation rate mu:
     *
     * For K81 (and JC which is a simpler version of it)
     * Pxy(t) = \pi_y (1 - exp(-mu t)) + exp(-mu t), if x ==y, \pi_y (1 - exp(-mu t)), otherwise
     * [Gascuel "Mathematics of Evolution and Phylogeny" 2005].
     */
    double exp_mu_t = exp(-mu * t);

    double p_ij = frequencies[j] * (1.0 - exp_mu_t);
    if (i == j) {
        p_ij += exp_mu_t;
    }
    return p_ij;
}

void set_p_ij(const Node *nd, int num_frequencies, const double *parameters) {
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
    double t = get_rescaled_branch_len(nd, scaling_factor, epsilon);
    // TODO: do we need to recalculate mu each time
    // or we fix it in the beginning and play with the scaling factor instead?
    double mu = get_mu(parameters, num_frequencies);

    for (i = 0; i < num_frequencies; i++) {
        for (j = 0; j < num_frequencies; j++) {
            nd->pij[i][j] = get_pij(parameters, mu, t, i, j);
        }
    }
}

