#include "pastml.h"
#include "scaling.h"
#include "logger.h"
#include "models.h"

extern char *global_model;

void *CAllocMem(long n, char *name, char *func, int showInfo)
{
	void *P;
	
	if ( (P=calloc(n, 1))==NULL ) {
		fprintf(stderr, "Out of memory allocating '%s': %s()\n", name, func);
		exit(0);
	}
	
	return P;
}

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

void set_p_ij(const Node *nd, double avg_br_len, size_t num_frequencies, const double *parameters) {
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
    // TODO: do we need to recalculate mu each time
    // or we fix it in the beginning and play with the scaling factor instead?
    double mu = get_mu(parameters, num_frequencies);
    double *P, *matrix[1];

    if ((strcmp(global_model, "JC") == 0) || (strcmp(global_model, "F81") == 0)) {
      for (i = 0; i < num_frequencies; i++) {
        for (j = 0; j < num_frequencies; j++) {
            nd->pij[i][j] = get_pij(parameters, mu, t, i, j);
        }
      }
    }
   
    if (strcmp(global_model, "HKY") == 0) {
      get_pij_hky(nd, num_frequencies, parameters, t);
    }

    if (strcmp(global_model, "JTT") == 0) {
      matrix[0] = CAllocMem(num_frequencies*num_frequencies*sizeof(double), "matrix", "CreateRates", 0);
      SetJTTMatrix(matrix[0], t);
      P=matrix[0];
      for(i=0;i<num_frequencies;i++){
        for(j=0;j<num_frequencies;j++){
          nd->pij[i][j] = (*P);
          P++;
        }
      }
    }
}

int calculate_node_probabilities(const Node *nd, size_t num_annotations, size_t first_child_index) {
    int factors = 0;
    size_t i, j, k;
    for (k = first_child_index; k < nd->nb_neigh; k++) {
        Node *child = nd->neigh[k];
        for (i = 0; i < num_annotations; i++) {
            /* Calculate the probability of having a branch from the node to its child node,
             * given that the node is in state i: p_child_branch_from_i = sum_j(p_ij * p_child_j)
             */
            double p_child_branch_from_i = 0.;
            for (j = 0; j < num_annotations; j++) {
                p_child_branch_from_i += child->pij[i][j] * child->bottom_up_likelihood[j];
            }

            /* The probability of having the node in state i is a multiplication of
             * the probabilities of p_child_branch_from_i for all child branches:
             * condlike_i = mult_ii(p_child_ii_branch_from_i)
             */
            if (k == first_child_index) {
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
            return -1;
        }
        factors += add_factors;
    }
    return factors;
}

int process_node(Node *nd, Tree *s_tree, size_t num_annotations, double *parameters) {
    /**
     * Calculates node probabilities.
     * parameters = [frequency_char_1, .., frequency_char_n, scaling_factor, epsilon].
     */
    int factors = 0, add_factors;
    size_t i, first_child_index;

    /* set probabilities of substitution */
    if (nd != s_tree->root) {
      set_p_ij(nd, s_tree->avg_tip_branch_len, num_annotations, parameters);
    }

    /* not a tip */
    if (nd->nb_neigh != 1) {
        first_child_index = (nd == s_tree->root) ? 0 : 1;
        /* recursively calculate probabilities for children */
        for (i = first_child_index; i < nd->nb_neigh; i++) {
	  add_factors = process_node(nd->neigh[i], s_tree, num_annotations, parameters);
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
    return factors;
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

double calculate_bottom_up_likelihood(Tree *s_tree, size_t num_annotations, double *parameters) {
    /**
     * Calculates tree log likelihood.
     * parameters = [frequency_char_1, .., frequency_char_n, scaling_factor, epsilon].
     */
    double scaled_lk = 0;
    size_t i;

    int factors = process_node(s_tree->root, s_tree, num_annotations, parameters);

    /* if factors == -1, it means that the bottom_up_likelihood is 0 */
    if (factors != -1) {
        for (i = 0; i < num_annotations; i++) {
            /* multiply the probability by character frequency */
            s_tree->root->bottom_up_likelihood[i] = s_tree->root->bottom_up_likelihood[i] * parameters[i];
            scaled_lk += s_tree->root->bottom_up_likelihood[i];
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

    for (k = 0; k < s_tree->nb_nodes; k++) {
        nd = s_tree->nodes[k];
        /* if a tip, process it */
        if (nd->nb_neigh == 1) {
            for (i = 0; i < num_tips; i++) {
                if (strcmp(nd->name, tip_names[i]) == 0) {
                    // states[i] == num_annotations means that the annotation is missing
                    if (states[i] == num_annotations) {
                        // and therefore any state is possible
                        for (j = 0; j < num_annotations; j++) {
                            nd->bottom_up_likelihood[j] = 1.0;
                            nd->joint_likelihood[j] = 1.0;
                        }
                    } else {
                        nd->bottom_up_likelihood[states[i]] = 1.0;
                        nd->joint_likelihood[states[i]] = 1.0;
			nd->best_joint_state = states[i];
                    }
                    break;
                }
            }
        }
    }
}

void rescale_branch_lengths(Tree *s_tree, double scaling_factor, double epsilon) {
    /**
     * Rescales all the branches in the tree.
     */
    Node *nd;
    size_t i;
    for (i = 0; i < s_tree->nb_nodes; i++) {
        nd = s_tree->nodes[i];
        nd->original_len = nd->branch_len;
        nd->branch_len = get_rescaled_branch_len(nd, s_tree->avg_tip_branch_len, scaling_factor, epsilon);
    }
}


