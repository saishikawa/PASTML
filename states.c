
#include <errno.h>
#include "pastml.h"


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

int output_ancestral_states(Tree *tree, size_t num_annotations, char **character, char *output_file_path) {
    FILE* outfile = fopen(output_file_path, "w");
    if (!outfile) {
        fprintf(stderr, "Output annotation file %s is impossible to access.", output_file_path);
        fprintf(stderr, "Value of errno: %d\n", errno);
        fprintf(stderr, "Error opening the file: %s\n", strerror(errno));
        return ENOENT;
    }
    size_t i, j, k;
    Node* nd;

    // print the header
    fprintf(outfile, "node ID");
    for (i = 0; i < num_annotations; i++) {
        fprintf(outfile, ",\"%s\"", character[i]);
    }
    fprintf(outfile, "\n");


    for (k = 0; k < tree->nb_nodes; k++) {
        nd = tree->nodes[k];

        fprintf(outfile, "%s", nd->name);

        for (i = 0; i < num_annotations; i++) {
            for (j = 0; j < num_annotations; j++) {
                if (strcmp(character[i], character[nd->best_states[j]]) == 0) {
                    fprintf(outfile, ",%.e", nd->result_probs[j]);
                }
            }
        }
        fprintf(outfile, "\n");
    }

    fclose(outfile);
    return EXIT_SUCCESS;
}


int output_parameters(double *parameters, size_t num_annotations, char **character, double log_lh, char *model,
                      char *output_file_path) {
    FILE* outfile = fopen(output_file_path, "w");
    if (!outfile) {
        fprintf(stderr, "Output parameter file %s is impossible to access.", output_file_path);
        fprintf(stderr, "Value of errno: %d\n", errno);
        fprintf(stderr, "Error opening the file: %s\n", strerror(errno));
        return ENOENT;
    }
    size_t i;


    fprintf(outfile, "parameter,value\n");
    fprintf(outfile, "model,%s\n", model);
    fprintf(outfile, "log likelihood,%.8f\n", log_lh);
    fprintf(outfile, "scaling factor,%.8f\n", parameters[num_annotations]);
    fprintf(outfile, "epsilon,%.e\n", parameters[num_annotations + 1]);
    for (i = 0; i < num_annotations; i++) {
        fprintf(outfile, "\"%s\",%.8f\n", character[i], parameters[i]);
    }
    fclose(outfile);
    return EXIT_SUCCESS;
}
