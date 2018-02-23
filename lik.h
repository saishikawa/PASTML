//
// Created by azhukova on 1/24/18.
//
#include "pastml.h"

#ifndef PASTML_LIK_H
#define PASTML_LIK_H

double
calc_lik_bfgs(Node *root, int num_annotations, double *parameters);

double get_rescaled_branch_len(const Node *nd, double scaling_factor, double epsilon);
void rescale_branch_lengths(Node *nd, Node *root, double scaling_factor, double epsilon);
double get_mu(const double* frequencies, int n);
void
initialise_tip_probabilities(Node *nd, Node *root,
                             char *const *tipnames, const int *states, int num_tips, int num_annotations);
double get_pij(const double *frequencies, double mu, double t, int i, int j);
void normalize(double *array, int n);
int get_max(int *array, int n);

#endif //PASTML_LIK_H
