//
// Created by azhukova on 1/24/18.
//
#include "pastml.h"

#ifndef PASTML_LIK_H
#define PASTML_LIK_H

double calc_lik_bfgs(Node *nd, char **tipnames, int *states, int num_tips, int num_annotations, double mu, char *model, double *parameters);

#endif //PASTML_LIK_H
