//
// Created by azhukova on 2/20/18.
//

#ifndef PASTML_PARAM_MINIMIZATION_H
#define PASTML_PARAM_MINIMIZATION_H

#include <stdbool.h>

double minimize_params(Tree *s_tree, size_t num_annotations, double *parameters, char **character, size_t set_values,
                       double scale_low, double scale_up);
#endif //PASTML_PARAM_MINIMIZATION_H
