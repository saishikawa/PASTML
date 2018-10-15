//
// Created by azhukova on 2/20/18.
//

#ifndef PASTML_PARAM_MINIMIZATION_H
#define PASTML_PARAM_MINIMIZATION_H

static const double EPSILON_APPROXIMATING_ZERO = 1.e-8;

#include <stdbool.h>
#include <math.h>

void minimize_params(Tree *s_tree, size_t num_annotations, double *parameters, char **character, size_t set_values);
#endif //PASTML_PARAM_MINIMIZATION_H
