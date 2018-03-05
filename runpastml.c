#include "pastml.h"
#include "marginal_likelihood.h"
#include "likelihood.h"
#include "marginal_approximation.h"
#include "param_minimization.h"
#include "output_tree.h"
#include "output_states.h"
#include "logger.h"
#include "make_tree.h"
#include <time.h>
#include <errno.h>

extern QUIET;

size_t tell_size_of_one_tree(char *filename) {
    /* the only purpose of this is to know about the size of a treefile (NH format)
     * in order to save memspace in allocating the string later on */
    size_t mysize = 0;
    int u;
    FILE *myfile = fopen(filename, "r");
    if (myfile) {
        while ((u = fgetc(myfile)) != ';') { /* termination character of the tree */
            if (feof(myfile)) {
                break;
            } /* shouldn't happen anyway */
            if (!isspace(u)) {
                mysize++;
            }
        }
        fclose(myfile);
    } /* end if(myfile) */
    return mysize + 1;
}

int copy_nh_stream_into_str(FILE *nh_stream, char *big_string) {
    int index_in_string = 0;
    int u;
    /* rewind(nh_stream);
     * DO NOT go to the beginning of the stream
     * if we want to make this flexible enough to read several trees per file */
    while ((u = fgetc(nh_stream)) != ';') { /* termination character of the tree */
        if (feof(nh_stream)) {
            big_string[index_in_string] = '\0';
            return EXIT_FAILURE;
        } /* error code telling that no tree has been read properly */
        if (index_in_string == MAX_TREELENGTH - 1) {
            fprintf(stderr, "Fatal error: tree file seems too big, are you sure it is a newick tree file?\n");
            return EXIT_FAILURE;
        }
        if (!isspace(u)) {
            big_string[index_in_string++] = (char) u;
        }
    }
    big_string[index_in_string++] = ';';
    big_string[index_in_string] = '\0';
    return EXIT_SUCCESS; /* leaves the stream right after the terminal ';' */
} /*end copy_nh_stream_into_str */

void free_node(Node *node, int count, size_t num_anno) {
    int j;

    if (node == NULL) return;
    if (node->name && count != 0) {
        free(node->name);
    }
    free(node->neigh);
    free(node->bottom_up_likelihood);
    free(node->marginal);
    free(node->best_states);
    free(node->top_down_likelihood);
    for (j = 0; j < num_anno; j++) {
        free(node->pij[j]);
    }
    free(node->pij);
    free(node);
}

void free_tree(Tree *tree, size_t num_anno) {
    int i;
    if (tree == NULL) return;
    for (i = 0; i < tree->nb_nodes; i++) {
        free_node(tree->nodes[i], i, num_anno);
    }
    free(tree->nodes);
    free(tree);
}


char **read_annotations(char *annotation_file_path, char **tips, int *states,
                        size_t *num_annotations, size_t *num_tips) {
    char annotation_line[MAXLNAME];
    int found_new_annotation;
    size_t i;
    size_t max_characters = 50;
    char **character = calloc(max_characters, sizeof(char *));
    for (i = 0; i < max_characters; i++) {
        character[i] = calloc(MAXLNAME, sizeof(char));
    }
    *num_annotations = 0;
    *num_tips = 0;

    char **annotations = calloc(MAXNSP, sizeof(char *));
    for (i = 0; i < MAXNSP; i++) {
        annotations[i] = calloc(MAXLNAME, sizeof(char));
    }

    /*Read annotation from file*/
    FILE *annotation_file = fopen(annotation_file_path, "r");
    if (!annotation_file) {
        fprintf(stderr, "Annotation file %s is not found or is impossible to access.", annotation_file_path);
        fprintf(stderr, "Value of errno: %d\n", errno);
        fprintf(stderr, "Error opening the file: %s\n", strerror(errno));
        return NULL;
    }

    while (fgets(annotation_line, MAXLNAME, annotation_file)) {
        sscanf(annotation_line, "%[^\n,],%[^\n\r]", tips[*num_tips], annotations[*num_tips]);
        char *annotation_value = annotations[*num_tips];
        if (strcmp(annotation_value, "") == 0) sprintf(annotation_value, "?");
        if (strcmp(annotation_value, "?") == 0) {
            states[*num_tips] = -1;
        } else {
            found_new_annotation = TRUE;
            for (i = 0; i < *num_tips; i++) {
                if (strcmp(annotations[i], "?") != 0 && strcmp(annotation_value, annotations[i]) == 0) {
                    states[*num_tips] = states[i];
                    strcpy(character[states[*num_tips]], annotation_value);
                    found_new_annotation = FALSE;
                    break;
                }
            }
            if (found_new_annotation == TRUE) {
                states[*num_tips] = (int) *num_annotations;
                if (*num_annotations >= max_characters) {
                    /* Annotations do not fit in the character array (of size max_characters) anymore,
                     * so we gonna double reallocate the memory for the array (of double size) and copy data there */
                    max_characters *= 2;
                    character = realloc(character, max_characters * sizeof(char *));
                    if (character == NULL) {
                        fprintf(stderr, "Problems with allocating memory: %s\n", strerror(errno));
                        fprintf(stderr, "Value of errno: %d\n", errno);
                        return NULL;
                    }
                    for (i = *num_annotations; i < max_characters; i++) {
                        character[i] = calloc(MAXLNAME, sizeof(char));
                    }
                }
                strcpy(character[*num_annotations], annotation_value);
                *num_annotations = *num_annotations + 1;
            }
        }
        *num_tips = *num_tips + 1;
    }
    fclose(annotation_file);
    free(annotations);
    return character;
}

int calculate_frequencies(size_t num_annotations, size_t num_tips, int *states, char **character, char *model,
                          double *parameters) {
    /* we would need an additional spot in the count array for the missing data,
     * therefore num_annotations + 1*/
    int *count_array = calloc(num_annotations + 1, sizeof(int));
    if (count_array == NULL) {
        fprintf(stderr, "Memory problems: %s\n", strerror(errno));
        fprintf(stderr, "Value of errno: %d\n", errno);
        return ENOMEM;
    }

    for (int i = 0; i < num_tips; i++) {
        if (states[i] == -1) {
            states[i] = (int) num_annotations;
            sprintf(character[num_annotations], "?");
        }
        count_array[states[i]]++;
    }

    int sum_freq = 0;
    for (int i = 0; i <= num_annotations; i++) {
        sum_freq += count_array[i];
    }
    log_info("MODEL:\t%s\n\n", model);
    log_info("INITIAL FREQUENCIES:\n\n");
    for (int i = 0; i < num_annotations; i++) {
        if (strcmp(model, "JC") == 0) {
            parameters[i] = ((double) 1) / num_annotations;
        } else if (strcmp(model, "F81") == 0) {
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

Tree *read_tree(char *nwk, size_t num_anno) {
    /**
     * Read a tree from newick file
     */

    Tree *s_tree;

    FILE *tree_file = fopen(nwk, "r");
    if (tree_file == NULL) {
        fprintf(stderr, "Tree file %s is not found or is impossible to access.\n", nwk);
        fprintf(stderr, "Value of errno: %d\n", errno);
        fprintf(stderr, "Error opening the file: %s\n", strerror(errno));
        return NULL;
    }

    size_t tree_file_size = 3 * tell_size_of_one_tree(nwk);
    if (tree_file_size > MAX_TREELENGTH) {
        fprintf(stderr, "Tree filesize for %s is more than %d bytes: are you sure it's a valid newick tree?\n",
                nwk, MAX_TREELENGTH / 3);
        return NULL;
    }

    void *retval;
    if ((retval = calloc(tree_file_size + 1, sizeof(char))) == NULL) {
        fprintf(stderr, "Not enough memory\n");
        return NULL;
    }
    char *c_tree = (char *) retval;

    if (EXIT_SUCCESS != copy_nh_stream_into_str(tree_file, c_tree)) {
        fprintf(stderr, "A problem occurred while parsing the reference tree.\n");
        return NULL;
    }
    fclose(tree_file);

    /*Make Tree structure*/
    s_tree = complete_parse_nh(c_tree, num_anno);
    if (NULL == s_tree) {
        fprintf(stderr, "A problem occurred while parsing the reference tree.\n");
        return NULL;
    }
    return s_tree;
}

int runpastml(char *annotation_name, char *tree_name, char *out_annotation_name, char *out_tree_name, char *model) {
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

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);
    srand((unsigned) time(NULL));

    if ((strcmp(model, "JC") != 0) && (strcmp(model, "F81") != 0)) {
        sprintf(stderr, "Model must be either JC or F81, not %s", model);
        return EINVAL;
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

    s_tree = read_tree(tree_name, num_annotations);
    if (s_tree == NULL) {
        return EXIT_FAILURE;
    }

    if (s_tree->nb_taxa != num_tips) {
        fprintf(stderr, "Number of annotations (even empty ones) specified in the annotation file (%zd)"
                " and the number of tips (%zd) do not match", num_tips, s_tree->nb_taxa);
    }
    num_tips = s_tree->nb_taxa;
    parameters[num_annotations] = 1.0;
    parameters[num_annotations + 1] = s_tree->min_branch_len;

    initialise_tip_probabilities(s_tree, tips, states, num_tips, num_annotations);
    free(tips);

    log_likelihood = calculate_bottom_up_likelihood(s_tree, num_annotations, parameters);
    if (log_likelihood == log(0)) {
        fprintf(stderr, "A problem occurred while calculating the bottom up likelihood: "
                "Is your tree ok and has at least 2 children per every inner node?\n");
        return EXIT_FAILURE;
    }
    log_info("INITIAL LOG LIKELIHOOD:\t%.10f\n\n", log_likelihood);


    log_info("OPTIMISING PARAMETERS...\n\n");
    log_likelihood = minimize_params(s_tree, num_annotations, parameters, character, model, 1.0 / 10000, 10000.0,
                                     MIN(s_tree->min_branch_len / 10.0, s_tree->avg_branch_len / 100.0),
                                     s_tree->avg_branch_len / 10.0);
    log_info("\n");

    log_info("OPTIMISED PARAMETERS:\n\n");
    if (0 == strcmp("F81", model)) {
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

    rescale_branch_lengths(s_tree, parameters[num_annotations], parameters[num_annotations + 1]);

    //Marginal bottom_up_likelihood calculation
    log_info("CALCULATING MARGINAL PROBABILITIES...\n\n");
    calculate_marginal_probabilities(s_tree, num_annotations, parameters);
    log_info("PREDICTING MOST LIKELY ANCESTRAL STATES...\n\n");
    choose_likely_states(s_tree, num_annotations);

    exit_val = write_nh_tree(s_tree, out_tree_name, parameters[num_annotations], parameters[num_annotations + 1]);
    if (EXIT_SUCCESS != exit_val) {
        return exit_val;
    }
    log_info("SAVING THE RESULTS...\n\n");
    log_info("\tScaled tree with internal node ids is written to %s.\n", out_tree_name);

    exit_val = output_state_ancestral_states(s_tree, num_annotations, character, out_annotation_name);
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
