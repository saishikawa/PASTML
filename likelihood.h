/**
* Created by azhukova on 1/24/18.
**/
#include "pastml.h"

#ifndef PASTML_LIK_H
#define PASTML_LIK_H

double
calculate_bottom_up_likelihood(Tree *s_tree, size_t num_annotations, double *parameters, int is_marginal,
                               char* model);

void rescale_branch_lengths(Tree *s_tree, double scaling_factor, double epsilon);
double get_mu(const double* frequencies, size_t n);
void
initialise_tip_probabilities(Tree *s_tree, char *const *tip_names, const int *states,
                             size_t num_tips, size_t num_annotations);
double get_pij(const double *frequencies, double mu, double t, int i, int j);
void normalize(double *array, size_t n);
int get_max(const int *array, size_t n);
void choose_joint_states(Tree *s_tree, size_t num_annotations, const double* frequencies);

#endif //PASTML_LIK_H
