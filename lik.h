//
// Created by azhukova on 1/24/18.
//
#include "pastml.h"

#ifndef PASTML_LIK_H
#define PASTML_LIK_H

void calc_lik_bfgs(Node *nd, char **tipnames, int *states, int nb, int nbanno, double mu, char *model, double *p,
                   double *likelihood);

#endif //PASTML_LIK_H
