//
// Created by azhukova on 1/24/18.
//

#ifndef PASTML_OUTPUT_STATES_H
#define PASTML_OUTPUT_STATES_H

#include "pastml.h"

int output_ancestral_states(Tree *tree, size_t num_annotations, char **character, char *output_file_path);
int output_parameters(double *parameters, size_t num_annotations, char **character, double log_lh, char *model,
                      char *output_file_path);
char **read_annotations(char *annotation_file_path, char **tips, int *states,
                        size_t *num_annotations, size_t *num_tips);

#endif //PASTML_OUTPUT_STATES_H
