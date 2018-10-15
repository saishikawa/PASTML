
#include <errno.h>
#include "pastml.h"
#include "states.h"

size_t read_parameters(char* parameter_file_path, char **character, size_t num_annotations, double *parameters,
        bool read_frequencies) {
    char param_line[MAXLNAME];
    size_t i;
    double double_value = 0;
    FILE *param_file = fopen(parameter_file_path, "r");
    char *quoted_state = (char*)malloc(255 * sizeof(char));
    // we'll set its 0's bite to 1 if the frequencies are set, and its 1st bite to 1 if scaling is set
    size_t result = 0;
    size_t num_frequency_values = 0;

    if (!param_file) {
        fprintf(stderr, "Parameter file %s is not found or is impossible to access.", parameter_file_path);
        fprintf(stderr, "Value of errno: %d\n", errno);
        fprintf(stderr, "Error opening the file: %s\n", strerror(errno));
        return result;
    }

    while (fgets(param_line, MAXLNAME, param_file)) {
        char* name = strtok(param_line, ",");
        if (strcmp(name, "parameter") == 0 || strcmp(name, "model") == 0 || strcmp(name, "log likelihood") == 0) {
            continue;
        }
        char* value = (name != NULL) ? strtok(NULL, ","): NULL;
        if (strcmp(name, "scaling factor") == 0) {
            sscanf(value, "%lf", &double_value);
            parameters[num_annotations] = double_value;
            result |= SF_SET;
        } else if (read_frequencies) {
            for (i = 0; i < num_annotations; i++) {
                sprintf(quoted_state, "\"%s\"", character[i]);
                if (strcmp(character[i], name) == 0 || strcmp(quoted_state, name) == 0) {
                    break;
                }
            }
            if (i < num_annotations) {
                sscanf(value, "%lf", &double_value);
                parameters[i] = double_value;
                num_frequency_values++;
            }
        }
    }
    if (num_frequency_values == num_annotations) {
        result |= FREQUENCIES_SET;
    }
    fclose(param_file);
    return result;
}


char **read_annotations(char *annotation_file_path, char **tips, int *states,
                        size_t *num_annotations, size_t *num_tips) {
    char annotation_line[MAXLNAME];
    bool found_new_annotation;
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
            found_new_annotation = true;
            for (i = 0; i < *num_tips; i++) {
                if (strcmp(annotations[i], "?") != 0 && strcmp(annotation_value, annotations[i]) == 0) {
                    states[*num_tips] = states[i];
                    strcpy(character[states[*num_tips]], annotation_value);
                    found_new_annotation = false;
                    break;
                }
            }
            if (found_new_annotation) {
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

int output_ancestral_states(Tree *tree, size_t num_annotations, char **character, char *output_file_path, char* format) {
    FILE* outfile = fopen(output_file_path, "w");
    if (!outfile) {
        fprintf(stderr, "Output annotation file %s is impossible to access.", output_file_path);
        fprintf(stderr, "Value of errno: %d\n", errno);
        fprintf(stderr, "Error opening the file: %s\n", strerror(errno));
        return ENOENT;
    }
    size_t i, k;
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
            fprintf(outfile, ",");
            fprintf(outfile, format, nd->result_probs[i]);
        }
        fprintf(outfile, "\n");
    }

    fclose(outfile);
    return EXIT_SUCCESS;
}

int output_parameters(const double *parameters, size_t num_annotations, char **character, double log_lh,
        const char *model, size_t set_values, const char *output_file_path, Tree *tree) {
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
    fprintf(outfile, "scaling factor was fixed,%s\n", (set_values & SF_SET) == 0 ? "No": "Yes");
    fprintf(outfile, "scaling factor,%.8f\n", parameters[num_annotations]);
    fprintf(outfile, "state changes per avg branch length,%.8f\n",
            parameters[num_annotations] * tree->avg_branch_len);
    fprintf(outfile, "frequencies were fixed,%s\n", (set_values & FREQUENCIES_SET) == 0 ? "No": "Yes");
    for (i = 0; i < num_annotations; i++) {
        fprintf(outfile, "\"%s\",%.8f\n", character[i], parameters[i]);
    }
    fclose(outfile);
    return EXIT_SUCCESS;
}
