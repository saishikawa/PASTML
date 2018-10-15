/**
* Created by azhukova on 1/24/18.
**/

#ifndef PASTML_OUTPUT_STATES_H
#define PASTML_OUTPUT_STATES_H

#include "pastml.h"
#include <stdbool.h>

#define FREQUENCIES_SET 1
#define SF_SET 2

int output_ancestral_states(Tree *tree, size_t num_annotations, char **character, char *output_file_path, char* format);
int output_parameters(const double *parameters, size_t num_annotations, char **character, double log_lh,
                      const char *model, size_t set_values, const char *output_file_path, Tree* tree);
char **read_annotations(char *annotation_file_path, char **tips, int *states,
                        size_t *num_annotations, size_t *num_tips);

size_t read_parameters(char* parameter_file_path, char **character, size_t num_annotations, double *parameters,
        bool read_frequencies);

#endif //PASTML_OUTPUT_STATES_H
