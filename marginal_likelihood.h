/**
* Created by azhukova on 1/24/18.
**/

#ifndef PASTML_MARGINAL_LIK_H_H
#define PASTML_MARGINAL_LIK_H_H

#include "pastml.h"

void calculate_marginal_probabilities(Tree *s_tree, size_t num_annotations, double *frequencies);
void calculate_top_down_likelihood(Tree *s_tree, size_t num_annotations);

#endif //PASTML_MARGINAL_LIK_H_H
