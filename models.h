//
// Created by saishikawa on 3/12/18.
//

#ifndef PASTML_MODEL_H
#define PASTML_MODEL_H

#include "pastml.h"

void exchange_params(size_t num_tips, int *states, char **character, char *model);
void SetJTTMatrix(double *matrix, double len);
void get_pij_hky(const Node *nd, size_t num_frequencies, const double *frequencies, double bl);
void setJTTFrequencies(double* parameters, size_t n);
void setHKYFrequencies(double* parameters);

#endif //PASTML_MODEL_H
