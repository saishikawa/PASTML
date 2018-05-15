//
// Created by azhukova on 2/22/18.
//



#include <stdio.h>
#include "pastml.h"
#include "likelihood.h"

int get_scaling_pow(double value) {
    return 128;
}

void rescale_array(double *array, size_t n, int curr_scaler_pow) {
    size_t i;

    for (i = 0; i < n; i++) {
      array[i] *= pow(2,curr_scaler_pow);
    }
}

int upscale_node_probs(double* array, size_t n) {
    /**
     * The rescaling is done to avoid underflow problems:
     * if a certain node probability is too small, we multiply this node probabilities by a scaling factor,
     * and keep the factor in mind to remove it from the final likelihood.
     */

  double smallest = 1.1, largest = 0.0;
    size_t i;
    int factors = 0;

    /* find the smallest non-zero probability */

    for (i = 0; i < n; i++) {
        if (array[i] > 0.0 && array[i] > largest) {
            largest = array[i];
        }
    }
    /* the whole array is zero */
    if (largest == 0.0) {
        return -1;
    }
    if (largest < LIM_P) {
        factors = get_scaling_pow(largest);
        rescale_array(array, n, factors);
    }

    return factors;
}

void rescale(double *array, size_t i, int curr_scaler_pow) {

    array[i] *= pow(2,curr_scaler_pow);
}

int rescale_if_needed(double* array, size_t i) {
    int curr_scaler_pow = 0;
    if (array[i] < LIM_P) {
        curr_scaler_pow = get_scaling_pow(array[i]);
        rescale(array, i, curr_scaler_pow);
    }
    return curr_scaler_pow;
}

int harmonise_scaling(double* array, int* scaling_factors, size_t n) {
    size_t i;
    int curr_scaler_pow;
    int max_factor = get_max(scaling_factors, n);
    for (i = 0; i < n; i++) {
        curr_scaler_pow = max_factor - scaling_factors[i];
        if (curr_scaler_pow != 0) {
            rescale(array, i, curr_scaler_pow);
        }
    }
    return max_factor;
}

double remove_upscaling_factors(double log_likelihood, int factors) {
    /**
     * While calculating the node probabilities, the upscaling was done to avoid underflow problems:
     * if a certain node probability was too small, we multiplied this node probabilities by a scaling factor,
     * and kept the factor in mind to remove it from the final likelihood.
     *
     * Now its time to remove all the factors from the final likelihood.
     */

     log_likelihood -= LOG2 * factors;

    return log_likelihood;
}
