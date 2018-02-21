#include "pastml.h"

extern Tree *s_tree;

int upscale_node_probs(const Node *nd, int num_annotations);

void set_p_ij(const Node *nd, int num_annotations, const double *parameters, double mu);

double get_rescaled_branch_len(const Node *nd, double scaling_factor, double epsilon);

void
initialise_tip_probabilities(Node *nd, char *const *tipnames, const int *states, int num_tips, int num_annotations);

int calculate_node_probabilities(const Node *nd, int num_annotations, int first_child_index);

double remove_upscaling_factors(double log_likelihood, int factors);

int process_node(Node *nd, Node* root, char *const *tipnames, const int *states, 
                 int num_tips, int num_annotations, double mu, double *parameters) {
    /**
     * Calculates node probabilities.
     * parameters = [frequency_char_1, .., frequency_char_n, scaling_factor, epsilon].
     * mu = 1 / (1 - (frequency_1^2 + ... + frequency_n^2)).
     */
    int factors = 0, add_factors;

    /* set probabilities of substitution */
    if (nd != root) {
        set_p_ij(nd, num_annotations, parameters, mu);
    }

    /* a tip */
    if (nd->nneigh == 1) {
        initialise_tip_probabilities(nd, tipnames, states, num_tips, num_annotations);
    } else {
        int first_child_index = (nd == root) ? 0 : 1;
        /* recursively calculate probabilities for children */
        for (int i = first_child_index; i < nd->nneigh; i++) {
            add_factors = process_node(nd->neigh[i], root, tipnames, states, num_tips, num_annotations, mu, parameters);
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
            nd->condlike[i] = nd->condlike[i] * parameters[i];
        }
    }

    if (nd->nneigh != 1) { // if not a tip
        memcpy( (void*)nd->up_like, (void*)nd->condlike, num_annotations * sizeof(double));
    }

    return factors;
}

double calc_lik_bfgs(Node *root, char *const *tipnames,const int *states, int num_tips, int num_annotations, double mu,
                     double *parameters) {
    /**
     * Calculates tree log likelihood.
     * parameters = [frequency_char_1, .., frequency_char_n, scaling_factor, epsilon].
     * mu = 1 / (1 - (frequency_1^2 + ... + frequency_n^2)).
     */
    double scaled_lk = 0;

    int factors = process_node(root, root, tipnames, states, num_tips, num_annotations, mu, parameters);
    /* if factors == -1, it means that the likelihood is 0 */
    if (factors != -1) {
        for (int i = 0; i < num_annotations; i++) {
            scaled_lk += root->condlike[i];
        }
        for (int i = 0; i < num_annotations; i++) {
            root->mar_prob[i] = root->condlike[i] / scaled_lk;
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
                p_child_branch_from_i += child->pij[i][j] * child->condlike[j];
            }

            /* The probability of having the node in state i is a multiplication of
             * the probabilities of p_child_branch_from_i for all child branches:
             * condlike_i = mult_ii(p_child_ii_branch_from_i)
             */
            if (ii == first_child_index) {
                nd->condlike[i] = p_child_branch_from_i;
            } else {
                nd->condlike[i] *= p_child_branch_from_i;
            }
        }
        int add_factors = upscale_node_probs(nd, num_annotations);

        /* if all the probabilities are zero (shown by add_factors == -1),
         * there is no point to go any further
         */
        if (add_factors == -1) {
            return -1;
        }
        factors += add_factors;
    }
    return factors;
}

void
initialise_tip_probabilities(Node *nd, char *const *tipnames, const int *states, int num_tips, int num_annotations) {
    /**
     * Sets the state and likelihoods for a tip
     * by setting the likelihood of its real state (given in the metadata file) to 1
     * and the other to 0.
     */
    int i, j;

    for (i = 0; i < num_tips; i++) {
        if (strcmp(nd->name, tipnames[i]) == 0) {
            nd->pupko_state = states[i];

            // states[i] == num_annotations means that the annotation is missing
            if (states[i] == num_annotations) {
                // and therefore any state is possible
                for (j = 0; j < num_annotations; j++) {
                    nd->condlike[j] = 1.0;
                    nd->up_like[j] = 1.0;
                }
            } else {
                nd->condlike[states[i]] = 1.0;
                nd->up_like[states[i]] = 1.0;
            }
            break;
        }
    }
}

void set_p_ij(const Node *nd, int num_annotations, const double *parameters, double mu) {
    /**
     * Sets node probabilities of substitution: p[i][j].
     */
    int i, j;

    double scaling_factor = parameters[num_annotations];
    double epsilon = parameters[num_annotations + 1];
    double bl = get_rescaled_branch_len(nd, scaling_factor, epsilon);
    double expmul = exp(-mu * bl);

    for (i = 0; i < num_annotations; i++) {
        for (j = 0; j < num_annotations; j++) {
            if (i == j) {
                nd->pij[i][j] = expmul + ((1.0 - expmul) * parameters[i]);
            } else {
                nd->pij[i][j] = parameters[j] * (1.0 - expmul);
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
    double bl = nd->br[0]->brlen;
    if (nd->nneigh == 1) { /*a tip*/
        bl = (bl + epsilon) * (s_tree->tip_avg_branch_len / (s_tree->tip_avg_branch_len + epsilon));
    }
    return bl * scaling_factor;
}

int upscale_node_probs(const Node *nd, int num_annotations) {
    /**
     * The rescaling is done to avoid underflow problems:
     * if a certain node probability is too small, we multiply this node probabilities by a scaling factor,
     * and keep the factor in mind to remove it from the final likelihood.
     */
    int j;
    int factors = 0;

    /* find the smallest non-zero probability */
    double smallest = 1.1;
    for (j = 0; j < num_annotations; j++) {
        if (nd->condlike[j] > 0.0 && nd->condlike[j] < smallest) {
            smallest = nd->condlike[j];
        }
    }

    if (smallest == 1.1) {
        fprintf(stderr, "Likelihood of the node %s is already 0, no need to go further.\n", nd->name);
        return -1;
    }

    if (smallest < LIM_P) {
        int curr_scaler_pow = (int) (POW * LOG2 - log(smallest)) / LOG2;
        int piecewise_scaler_pow;
        double curr_scaler;

        factors = curr_scaler_pow;
        do {
            piecewise_scaler_pow = MIN(curr_scaler_pow, 63);
            curr_scaler = ((unsigned long long) (1) << piecewise_scaler_pow);
            for (j = 0; j < num_annotations; j++) {
                nd->condlike[j] *= curr_scaler;
            }
            curr_scaler_pow -= piecewise_scaler_pow;
        } while (curr_scaler_pow != 0);
    }
    return factors;
}
