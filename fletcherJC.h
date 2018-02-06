//
// Created by azhukova on 1/24/18.
//

#ifndef PASTML_FLETCHERJC_H
#define PASTML_FLETCHERJC_H

#include "pastml.h"

void frprmnJC(Node *nd, char **tipnames, int *states, int nb, int nbanno, double mu, char *model, double *p, int n,
              double ftol, int *iter, double *fret, char **character, double *frequency);

#endif //PASTML_FLETCHERJC_H
