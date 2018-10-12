//
// Created by azhukova on 2/22/18.
//



#include <stdio.h>
#include "pastml.h"
#include "likelihood.h"

int get_scaling_pow(double value) {
    return (int) ceil(POW - log(value) / LOG2);
}

void rescale_array(double *array, size_t n, int curr_scaler_pow) {
    int piecewise_scaler_pow;
    size_t i;
    unsigned long long curr_scaler;
    do {
        piecewise_scaler_pow = MIN(curr_scaler_pow, 63);
        curr_scaler = ((unsigned long long) (1) << piecewise_scaler_pow);
        for (i = 0; i < n; i++) {
            array[i] *= curr_scaler;
        }
        curr_scaler_pow -= piecewise_scaler_pow;
    } while (curr_scaler_pow != 0);
}

int upscale_node_probs(double* array, size_t n) {
    /**
     * The rescaling is done to avoid underflow problems:
     * if a certain node probability is too small, we multiply this node probabilities by a scaling factor,
     * and keep the factor in mind to remove it from the final likelihood.
     */
    bool smallest_set = false;
    double smallest;
    size_t i;
    int factors = 0;

    /* find the smallest non-zero probability */
    for (i = 0; i < n; i++) {
        if (array[i] > 0.0 && (!smallest_set || (array[i] < smallest))) {
            smallest_set = true;
            smallest = array[i];
        }
    }
    /* the whole array is zero */
    if (!smallest_set) {
        return -1;
    }
    if (smallest < LIM_P) {
        factors = get_scaling_pow(smallest);
        rescale_array(array, n, factors);
    }
    return factors;
}

void rescale(double *array, size_t i, int curr_scaler_pow) {
    int piecewise_scaler_pow;
    unsigned long long curr_scaler;
    do {
        piecewise_scaler_pow = MIN(curr_scaler_pow, 63);
        curr_scaler = ((unsigned long long) (1) << piecewise_scaler_pow);
        array[i] *= curr_scaler;
        curr_scaler_pow -= piecewise_scaler_pow;
    } while (curr_scaler_pow != 0);
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
    int piecewise_scaler_pow;
    do {
        piecewise_scaler_pow = MIN(factors, 63);
        log_likelihood -= LOG2 * piecewise_scaler_pow;
        factors -= piecewise_scaler_pow;
    } while (factors != 0);
    return log_likelihood;
}