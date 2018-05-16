#include "pastml.h"
#include "marginal_likelihood.h"
#include "likelihood.h"
#include "marginal_approximation.h"
#include "param_minimization.h"
#include "states.h"
#include "logger.h"
#include "tree.h"
#include "parsimony.h"
#include <time.h>
#include <errno.h>

extern int* QUIET;

int calculate_frequencies(size_t num_annotations, size_t num_tips, int *states, char **character, char *model,
                          double *parameters) {
    /* we would need an additional spot in the count array for the missing data,
     * therefore num_annotations + 1*/
    int *count_array = calloc(num_annotations + 1, sizeof(int));
    size_t i;

    if (count_array == NULL) {
        fprintf(stderr, "Memory problems: %s\n", strerror(errno));
        fprintf(stderr, "Value of errno: %d\n", errno);
        return ENOMEM;
    }

    for (i = 0; i < num_tips; i++) {
        if (states[i] == -1) {
            states[i] = (int) num_annotations;
            sprintf(character[num_annotations], "?");
        }
        count_array[states[i]]++;
    }

    int sum_freq = 0;
    for (i = 0; i <= num_annotations; i++) {
        sum_freq += count_array[i];
    }
    log_info("INITIAL FREQUENCIES:\n\n");
    for (i = 0; i < num_annotations; i++) {
        if (strcmp(model, JC) == 0) {
            parameters[i] = ((double) 1) / num_annotations;
        } else if (strcmp(model, F81) == 0) {
            parameters[i] = ((double) count_array[i]) / sum_freq;
        }
        log_info("\t%s:\t%.10f\n", character[i], parameters[i]);
    }
    if (count_array[num_annotations] > 0.0) {
        log_info("\n\tMissing data:\t%.10f\n", (double) count_array[num_annotations] / (double) sum_freq);
    }
    log_info("\n");
    free(count_array);
    return EXIT_SUCCESS;
}

int is_valid_model(char* model) {
    return (strcmp(model, JC) == 0) || (strcmp(model, F81) == 0);
}

int is_marginal_method(char *prob_method) {
    return (strcmp(prob_method, MARGINAL) == 0) || (strcmp(prob_method, MARGINAL_APPROXIMATION) == 0)
           || (strcmp(prob_method, MAX_POSTERIORI) == 0);
}

int is_parsimonious_method(char *prob_method) {
    return (strcmp(prob_method, DOWNPASS) == 0) || (strcmp(prob_method, ACCTRAN) == 0)
           || (strcmp(prob_method, DELTRAN) == 0);
}

int is_ml_method(char *prob_method) {
    return (strcmp(prob_method, JOINT) == 0) || is_marginal_method(prob_method);
}

int is_valid_prediction_method(char *prob_method) {
    return is_parsimonious_method(prob_method) || is_ml_method(prob_method);
}

int runpastml(char *annotation_name, char *tree_name, char *out_annotation_name, char *out_tree_name,
              char *out_parameter_name, char *model,
              char* prob_method) {
    int i;
    int *states;
    double log_likelihood, sec;
    int minutes;
    double *parameters;
    char **character, **tips;
    size_t num_annotations, num_tips = 0;
    struct timespec time_start, time_end;
    int exit_val;
    Tree *s_tree;
    int is_marginal, is_parsimonious;

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);
    srand((unsigned) time(NULL));

    if (!is_valid_prediction_method(prob_method)) {
        fprintf(stderr, "Ancestral state prediction method \"%s\" is not valid", prob_method);
        return EINVAL;
    }

    is_parsimonious = is_parsimonious_method(prob_method);
    is_marginal = is_marginal_method(prob_method);

    if (!is_parsimonious && !is_valid_model(model)) {
        fprintf(stderr, "Model value \"%s\" is not valid", model);
        return EINVAL;
    }

    if (out_annotation_name == NULL) {
        out_annotation_name = calloc(256, sizeof(char));
        if (is_parsimonious) {
            sprintf(out_annotation_name, "%s.%s.pastml.out.csv", annotation_name, prob_method);
        } else {
            sprintf(out_annotation_name, "%s.%s.%s.pastml.out.csv", annotation_name, model, prob_method);
        }
    }
    if (out_tree_name == NULL) {
        out_tree_name = calloc(256, sizeof(char));
        if (is_parsimonious) {
            sprintf(out_tree_name, "%s.parsimonious.pastml.tree.nwk", annotation_name);
        } else {
            sprintf(out_tree_name, "%s.maxlikelihood.pastml.tree.nwk", annotation_name);
        }
    }
    if (out_parameter_name == NULL) {
        if (!is_parsimonious) {
            out_parameter_name = calloc(256, sizeof(char));
            sprintf(out_parameter_name, "%s.maxlikelihood.pastml.parameters.csv", annotation_name);
        }
    }

    log_info("ANCESTRAL STATE PREDICTION METHOD:\t%s\n\n", prob_method);
    if (!is_parsimonious) {
        log_info("MODEL:\t%s\n\n", model);
    }

    /* Allocate memory */
    states = calloc(MAXNSP, sizeof(int));
    tips = calloc(MAXNSP, sizeof(char *));
    for (i = 0; i < MAXNSP; i++) {
        tips[i] = calloc(MAXLNAME, sizeof(char));
    }
    states[0] = 0;

    size_t *num_anno_arr = calloc(1, sizeof(int));
    size_t *num_tips_arr = calloc(1, sizeof(int));
    character = read_annotations(annotation_name, tips, states, num_anno_arr, num_tips_arr);
    if (character == NULL) {
        return EXIT_FAILURE;
    }
    num_annotations = *num_anno_arr;
    num_tips = *num_tips_arr;
    free(num_anno_arr);
    free(num_tips_arr);

    s_tree = read_tree(tree_name, num_annotations);
    if (s_tree == NULL) {
        return EXIT_FAILURE;
    }

    if (s_tree->nb_taxa != num_tips) {
        fprintf(stderr, "Number of annotations (even empty ones) specified in the annotation file (%zd)"
                " and the number of tips (%zd) do not match", num_tips, s_tree->nb_taxa);
    }
    num_tips = s_tree->nb_taxa;

    initialise_tip_probabilities(s_tree, tips, states, num_tips, num_annotations);
    free(tips);

    if (is_parsimonious) {
        log_info("PREDICTING PARSIMONIOUS ANCESTRAL STATES...\n\n");
        parsimony(s_tree, num_annotations, prob_method);
        select_parsimonious_states(s_tree, num_annotations);
    } else {
        /* we would need two additional spots in the parameters array: for the scaling factor, and for the epsilon,
         * therefore num_annotations + 2*/
        parameters = calloc(num_annotations + 2, sizeof(double));
        if (parameters == NULL) {
            fprintf(stderr, "Memory problems: %s\n", strerror(errno));
            fprintf(stderr, "Value of errno: %d\n", errno);
            return ENOMEM;
        }

        exit_val = calculate_frequencies(num_annotations, num_tips, states, character, model, parameters);
        if (EXIT_SUCCESS != exit_val) {
            return exit_val;
        }
        parameters[num_annotations] = 1.0 / s_tree->avg_branch_len;
        parameters[num_annotations + 1] = s_tree->avg_tip_branch_len / 50.0;

        log_likelihood = calculate_bottom_up_likelihood(s_tree, num_annotations, parameters, is_marginal, model);
        if (log_likelihood == log(0)) {
            fprintf(stderr, "A problem occurred while calculating the bottom up likelihood: "
                    "Is your tree ok and has at least 2 children per every inner node?\n");
            return EXIT_FAILURE;
        }
        log_info("INITIAL LOG LIKELIHOOD:\t%.10f\n\n", log_likelihood);
        if (log_likelihood == log(1)) {
            log_info("INITIAL LIKELIHOOD IS PERFECT, CANNOT DO BETTER THAN THAT.\n\n");
        } else {
            log_info("OPTIMISING PARAMETERS...\n\n");
            log_likelihood = minimize_params(s_tree, num_annotations, parameters, character, model,
                                             0.001 / s_tree->avg_branch_len, 10.0 / s_tree->avg_branch_len,
                                             s_tree->avg_tip_branch_len / 100.0, s_tree->avg_tip_branch_len / 10.0);
            log_info("\n");

            log_info("OPTIMISED PARAMETERS:\n\n");
            if (0 == strcmp(F81, model)) {
                for (i = 0; i < num_annotations; i++) {
                    log_info("\tFrequency of %s:\t%.10f\n", character[i], parameters[i]);
                }
                log_info("\n");
            }
            log_info("\tScaling factor:\t%.10f \n", parameters[num_annotations]);
            log_info("\tEpsilon:\t%e\n", parameters[num_annotations + 1]);
            log_info("\n");
            log_info("OPTIMISED LOG LIKELIHOOD:\t%.10f\n", log_likelihood);
            log_info("\n");
        }

        exit_val = output_parameters(parameters, num_annotations, character, log_likelihood, model, out_parameter_name);
        if (EXIT_SUCCESS != exit_val) {
            return exit_val;
        }
        log_info("\tOptimised parameters are written to %s in csv format.\n", out_parameter_name);
        log_info("\n");

        rescale_branch_lengths(s_tree, parameters[num_annotations], parameters[num_annotations + 1]);

        if (is_marginal) {
            log_info("CALCULATING TOP-DOWN LIKELIHOOD...\n\n");
            calculate_top_down_likelihood(s_tree, num_annotations);

            log_info("CALCULATING MARGINAL PROBABILITIES...\n\n");
            calculate_marginal_probabilities(s_tree, num_annotations, parameters);

            normalize_result_probabilities(s_tree, num_annotations);
            set_id_best_states(s_tree, num_annotations);

            if (strcmp(prob_method, MARGINAL_APPROXIMATION) == 0) {
                log_info("PREDICTING MOST LIKELY ANCESTRAL STATES...\n\n");
                choose_likely_states(s_tree, num_annotations);
            } else if (strcmp(prob_method, MAX_POSTERIORI) == 0) {
                log_info("PREDICTING MOST LIKELY ANCESTRAL STATES...\n\n");
                choose_best_marginal_states(s_tree, num_annotations);
            }
        } else {
            // the branch lengths are already rescaled, so let's put scaling factor to 1, and epsilon to 0.
            parameters[num_annotations] = 1.0;
            parameters[num_annotations + 1] = 0.0;

            // calculate joint likelihood
            calculate_bottom_up_likelihood(s_tree, num_annotations, parameters, FALSE, model);

            log_info("PREDICTING MOST LIKELY ANCESTRAL STATES...\n\n");
            choose_joint_states(s_tree, num_annotations, parameters);
            set_id_best_states(s_tree, num_annotations);
        }
        free(parameters);
    }

    exit_val = write_nh_tree(s_tree, out_tree_name);
    if (EXIT_SUCCESS != exit_val) {
        return exit_val;
    }
    log_info("SAVING THE RESULTS...\n\n");
    log_info("\t%sree with internal node ids is written to %s.\n", (is_parsimonious) ? "T": "Scaled t", out_tree_name);

    exit_val = output_ancestral_states(s_tree, num_annotations, character, out_annotation_name);
    if (EXIT_SUCCESS != exit_val) {
        return exit_val;
    }
    log_info("\tState predictions are written to %s in csv format.\n", out_annotation_name);
    log_info("\n");

    //free all
    free(character);
    free(states);
    free_tree(s_tree, num_annotations);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_end);
    sec = (double) (time_end.tv_sec - time_start.tv_sec)
          + (time_end.tv_nsec - time_start.tv_nsec) / 1000.0 / 1000.0 / 1000.0;

    minutes = (int) (sec / 60.0);
    log_info("TOTAL EXECUTION TIME:\t%d minute%s %.2f seconds\n\n", minutes, (minutes != 1) ? "s" : "",
             sec - (60.0 * minutes));

    return EXIT_SUCCESS;
}
