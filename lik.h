//
// Created by azhukova on 1/24/18.
//
#include "pastml.h"

#ifndef PASTML_LIK_H
#define PASTML_LIK_H

double
calc_lik_bfgs(Node *root, char *const *tipnames, const int *states, int num_tips, int num_annotations, double mu,
              double *parameters);

double get_rescaled_branch_len(const Node *nd, double scaling_factor, double epsilon);

#endif //PASTML_LIK_H
