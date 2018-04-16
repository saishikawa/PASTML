//
// Created by azhukova on 2/22/18.
//

#ifndef PASTML_SCALING_H
#define PASTML_SCALING_H
int upscale_node_probs(double* array, size_t n);
void rescale(double *array, size_t i, int curr_scaler_pow);
int get_scaling_pow(double value);
double remove_upscaling_factors(double log_likelihood, int factors);
int harmonise_scaling(double* array, int* scaling_factors, size_t n);
int rescale_if_needed(double* array, size_t i);
#endif //PASTML_SCALING_H
