//
// Created by azhukova on 1/24/18.
//

#ifndef PASTML_FLETCHER_H
#define PASTML_FLETCHER_H

#include "pastml.h"

void frprmn(Node *nd, char **tipnames, int *states, int nb, int nbanno, double mu, char *model, double *p, int n,
            double ftol, int *iter, double *fret, char **character);

#endif //PASTML_FLETCHER_H
