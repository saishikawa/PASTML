#include "pastml.h"
#include "marginal_lik.h"
#include "lik.h"
#include "marginal_approxi.h"
#include "make_tree.h"
#include "param_minimization.h"
#include "output_tree.h"
#include "output_states.h"
#include <time.h>
#include <errno.h>

Tree *s_tree;
Node *root;

unsigned int tell_size_of_one_tree(char *filename) {
    /* the only purpose of this is to know about the size of a treefile (NH format) in order to save memspace in allocating the string later on */
    /* wew open and close this file independently of any other fopen */
    unsigned int mysize = 0;
    char u;
    FILE *myfile = fopen(filename, "r");
    if (myfile) {
        while ((u = fgetc(myfile)) != ';') { /* termination character of the tree */
            if (u == EOF) break; /* shouldn't happen anyway */
            if (isspace(u)) continue; else mysize++;
        }
        fclose(myfile);
    } /* end if(myfile) */
    return (mysize + 1);
}

int copy_nh_stream_into_str(FILE *nh_stream, char *big_string) {
    int index_in_string = 0;
    char u;
    /* rewind(nh_stream); DO NOT go to the beginning of the stream if we want to make this flexible enough to read several trees per file */
    while ((u = fgetc(nh_stream)) != ';') { /* termination character of the tree */
        if (u == EOF) {
            big_string[index_in_string] = '\0';
            return EXIT_FAILURE;
        } /* error code telling that no tree has been read properly */
        if (index_in_string == MAX_TREELENGTH - 1) {
            fprintf(stderr, "Fatal error: tree file seems too big, are you sure it is a newick tree file?\n");
            return EXIT_FAILURE;
        }
        if (isspace(u)) continue;
        big_string[index_in_string++] = u;
    }
    big_string[index_in_string++] = ';';
    big_string[index_in_string] = '\0';
    return EXIT_SUCCESS; /* leaves the stream right after the terminal ';' */
} /*end copy_nh_stream_into_str */

void free_node(Node *node, int count, int num_anno) {
    int j;

    if (node == NULL) return;
    if (node->name && count != 0) {
        free(node->name);
    }
    free(node->neigh);
    free(node->bottom_up_likelihood);
    free(node->condlike_mar);
    free(node->marginal);
    free(node->best_states);
    free(node->top_down_likelihood);
    for (j = 0; j < num_anno; j++) {
        free(node->pij[j]);
    }
    free(node->pij);
    free(node);
}

void free_tree(Tree *tree, int num_anno) {
    int i;
    if (tree == NULL) return;
    for (i = 0; i < tree->nb_nodes; i++) {
        free_node(tree->a_nodes[i], i, num_anno);
    }
    free(tree->a_nodes);
    free(tree);

}


char** read_annotations(char* annotation_name, char **tips, int *states, int* num_anno, int* num_tips) {
    char anno_line[MAXLNAME];
    int found_new_annotation, i;
    int max_characters = MAXCHAR;
    char **character = calloc(max_characters, sizeof(char *));
    for (i = 0; i < max_characters; i++) {
        character[i] = calloc(MAXLNAME, sizeof(char));
    }
    *num_anno = -1;
    *num_tips = 0;

    char **annotations = calloc(MAXNSP, sizeof(char *));
    for (i = 0; i < MAXNSP; i++) {
        annotations[i] = calloc(MAXLNAME, sizeof(char));
    }

    /*Read annotation from file*/
    FILE* annotationfile = fopen(annotation_name, "r");
    if (!annotationfile) {
        fprintf(stderr, "Annotation file %s is not found or is impossible to access.", annotation_name);
        fprintf(stderr, "Value of errno: %d\n", errno);
        fprintf(stderr, "Error opening the file: %s\n", strerror(errno));
        return NULL;
    }

    while (fgets(anno_line, MAXLNAME, annotationfile)) {
        sscanf(anno_line, "%[^\n,],%[^\n\r]", tips[*num_tips], annotations[*num_tips]);
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
                *num_anno = *num_anno + 1;
                states[*num_tips] = *num_anno;
                if (*num_anno >= max_characters) {
                    /* Annotations do not fit in the character array (of size max_characters) anymore,
                     * so we gonna double reallocate the memory for the array (of double size) and copy data there */
                    max_characters *= 2;
                    character = realloc(character, max_characters * sizeof(char *));
                    if (character == NULL) {
                        fprintf(stderr, "Problems with allocating memory: %s\n", strerror(errno));
                        fprintf(stderr, "Value of errno: %d\n", errno);
                        return NULL;
                    }
                    for (i = *num_anno; i < max_characters; i++) {
                        character[i] = calloc(MAXLNAME, sizeof(char));
                    }
                }
                strcpy(character[*num_anno], annotation_value);
            }
        }
        *num_tips = *num_tips + 1;
    }
    fclose(annotationfile);
    free(annotations);
    return character;
}

int calculate_frequencies(int num_anno, int num_tips, int* states, char** character, char* model, double* parameters) {
    /* we would need an additional spot in the count array for the missing data,
     * therefore num_anno + 1*/
    int* count_array = calloc(num_anno + 1, sizeof(int));
    if (count_array == NULL) {
        fprintf(stderr, "Memory problems: %s\n", strerror(errno));
        fprintf(stderr, "Value of errno: %d\n", errno);
        return ENOMEM;
    }

    for (int i = 0; i < num_tips; i++) {
        if (states[i] == -1) {
            states[i] = num_anno;
            sprintf(character[num_anno], "?");
        }
        count_array[states[i]]++;
    }

    int sum_freq = 0;
    for (int i = 0; i <= num_anno; i++) {
        sum_freq += count_array[i];
    }
    printf("\n*** Frequency of %d characters in the MODEL %s ***\n\n", num_anno, model);
    for (int i = 0; i < num_anno; i++) {
        if (strcmp(model, "JC") == 0) {
            parameters[i] = ((double) 1) / num_anno;
        } else if (strcmp(model, "F81") == 0) {
            parameters[i] = ((double) count_array[i]) / sum_freq;
        }
        printf("%s = %lf\n", character[i], parameters[i]);
    }

    printf("Frequency of missing data = %lf\n", (double) count_array[num_anno] / (double) sum_freq);
    free(count_array);
    return EXIT_SUCCESS;
}

Tree* read_tree(char* tree_name, int num_anno) {
    /*Read tree from file*/

    FILE* treefile = fopen(tree_name, "r");
    if (treefile == NULL) {
        fprintf(stderr, "Tree file %s is not found or is impossible to access.\n", tree_name);
        fprintf(stderr, "Value of errno: %d\n", errno);
        fprintf(stderr, "Error opening the file: %s\n", strerror(errno));
        return NULL;
    }

    unsigned int treefilesize = 3 * tell_size_of_one_tree(tree_name);
    if (treefilesize > MAX_TREELENGTH) {
        fprintf(stderr, "Tree filesize for %s is more than %d bytes: are you sure it's a valid newick tree?\n",
                tree_name, MAX_TREELENGTH / 3);
        return NULL;
    }

    void *retval;
    if( (retval=calloc(treefilesize + 1, sizeof(char))) == NULL ) {
        printf("Not enough memory\n");
        return NULL;
    }
    char* c_tree = (char *) retval;

    if (EXIT_SUCCESS != copy_nh_stream_into_str(treefile, c_tree)) {
        fprintf(stderr, "A problem occurred while parsing the reference tree.\n");
        return NULL;
    }
    fclose(treefile);

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
    double log_likelihood, nano, sec;
    double *parameters;
    char **character, **tips;
    int num_annotations, num_tips = 0;
    struct timespec samp_ini, samp_fin;
    int exit_val;

    if ((strcmp(model, "JC") != 0) && (strcmp(model, "F81") != 0)) {
        sprintf(stderr, "Model must be either JC or F81, not %s", model);
        return EINVAL;
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &samp_ini);
    srand((unsigned) time(NULL));

    /* Allocate memory */
    states = calloc(MAXNSP, sizeof(int));
    tips = calloc(MAXNSP, sizeof(char *));
    for (i = 0; i < MAXNSP; i++) {
        tips[i] = calloc(MAXLNAME, sizeof(char));
    }
    states[0] = 0;

    int* num_anno_arr = calloc(1, sizeof(int));
    int* num_tips_arr = calloc(1, sizeof(int));
    character = read_annotations(annotation_name, tips, states, num_anno_arr, num_tips_arr);
    if (character == NULL) {
        return EXIT_FAILURE;
    }
    num_annotations = *num_anno_arr;
    num_tips = *num_tips_arr;
    free(num_anno_arr);
    free(num_tips_arr);

    num_annotations++;

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

    root = s_tree->a_nodes[0];
    if (s_tree->nb_taxa != num_tips) {
        fprintf(stderr, "Number of annotations (even empty ones) specified in the annotation file (%d)"
                " and the number of tips (%d) do not match", num_tips, s_tree->nb_taxa);
    }
    num_tips = s_tree->nb_taxa;
    parameters[num_annotations] = 1.0;
    parameters[num_annotations + 1] = s_tree->tip_avg_branch_len / 100.0;

    initialise_tip_probabilities(root, root, tips, states, num_tips, num_annotations);

    log_likelihood = calc_lik_bfgs(root, num_annotations, parameters);
    if (log_likelihood == log(0)) {
        fprintf(stderr, "A problem occurred while calculating the bottom_up_likelihood: "
                "Is your tree ok and has at least 2 children per every inner node?\n");
        return EXIT_FAILURE;
    }
    printf("\n*** Initial log likelihood of the tree ***\n\n %lf\n\n", log_likelihood);

    log_likelihood = minimize_params(root, num_annotations, parameters, character, model, 1.0 / 10000, 10000.0,
                          s_tree->tip_avg_branch_len / 10000.0, s_tree->tip_avg_branch_len / 10.0);

    if (0 == strcmp("F81", model)) {
        printf("\n*** Optimized frequencies ***\n\n");
        for (i = 0; i < num_annotations; i++) {
            printf("%s = %.5f\n", character[i], parameters[i]);
        }
    }
    printf("\n*** Tree scaling factor ***\n\n %.5f \n\n*** Epsilon ***\n\n %.5e",
           parameters[num_annotations], parameters[num_annotations + 1]);
    printf("\n\n*** Optimised log likelihood ***\n\n %lf\n", log_likelihood);

    rescale_branch_lengths(root, root, parameters[num_annotations], parameters[num_annotations + 1]);

    //Marginal bottom_up_likelihood calculation
    printf("\n*** Calculating Marginal Likelihoods...\n");
    calculate_marginal_probabilities(root, root, num_annotations, parameters);

    printf("\n*** Predicting likely ancestral states by the Marginal Approximation method...\n\n");

    order_marginal(root, root, num_annotations);
    calc_correct(root, root, num_annotations);

    exit_val = write_nh_tree(root, out_tree_name, parameters[num_annotations], parameters[num_annotations + 1]);
    if (EXIT_SUCCESS != exit_val) {
        return exit_val;
    }
    printf("Scaled tree with internal node ids is written to %s\n", out_tree_name);

    exit_val = output_state_ancestral_states(root, num_annotations, character, out_annotation_name);
    if (EXIT_SUCCESS != exit_val) {
        return exit_val;
    }
    printf("Predictions for all internal nodes and the root are written to %s in csv format\n", out_annotation_name);

    //free all
    free(tips);
    free(character);
    free(states);
    free_tree(s_tree, num_annotations);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &samp_fin);
    sec = (double) (samp_fin.tv_sec - samp_ini.tv_sec);
    nano = (double) (samp_fin.tv_nsec - samp_ini.tv_nsec) / 1000.0 / 1000.0 / 1000.0;
    printf("\nTotal execution time = %3.10lf seconds\n\n", (sec + nano));

    return EXIT_SUCCESS;
}
