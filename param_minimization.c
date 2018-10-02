#include <stdio.h>
#include <gsl/gsl_multimin.h>
#include "likelihood.h"
#include "logger.h"
#include "param_minimization.h"

#define GRADIENT_STEP 1.0e-5


double random_double(double min_bound, double max_bound) {
    /**
     * Generates a random double within given bounds.
     */
    return min_bound + (max_bound - min_bound) * rand() / RAND_MAX;
}

void softmax(double* xs, size_t n) {
    /**
     * transforms an array of n arbitrary values x in such a way that all of them become between 0 and 1 and sum to 1,
     * using the softmax function.
     */
    double exp_sum = 0.0;
    size_t i;
    for (i = 0; i < n; i++) {
        xs[i] = exp(xs[i]);
        exp_sum += xs[i];
    }
    for (i = 0; i < n; i++) {
        xs[i] /= exp_sum;
    }
}

double sigmoid(double x, const double lower_bound, const double upper_bound) {
    /**
     * transforms an arbitrary value x to be in between lower and upper bound, using a sigmoid function.
     */
    return lower_bound + (upper_bound - lower_bound) / (1. + exp(-x));
}

void optimised_parameters2real_parameters(const gsl_vector *v, size_t num_annotations, double scale_low,
                                          double scale_up,
                                          double epsilon_low, double epsilon_up, double *cur_parameters,
                                          bool optimize_frequencies) {
    size_t i;
    if (optimize_frequencies) {
        /* 1. Frequencies */
        for (i = 0; i < num_annotations - 1; i++) {
            cur_parameters[i] = gsl_vector_get(v, i);
        }
        cur_parameters[num_annotations - 1] = 0.0;
        softmax(cur_parameters, num_annotations);
    }
    size_t scaling_factor_index = optimize_frequencies ? (num_annotations - 1): 0;
    /* 2. Scaling factor */
    cur_parameters[num_annotations] = sigmoid(gsl_vector_get(v, scaling_factor_index), scale_low, scale_up);

    /* 3. Epsilon */
    cur_parameters[num_annotations + 1] = sigmoid(gsl_vector_get(v, scaling_factor_index + 1), epsilon_low, epsilon_up);
}

void
log_cur_parameter_values(size_t num_annotations, double *parameters, bool log_frequencies, size_t iter,
                         const gsl_multimin_fdfminimizer *s, char **character) {
    size_t i;

    if (character != NULL) {
        log_info ("\tstep\tlog-lh\t\t");
        if (log_frequencies) {
            for (i = 0; i < num_annotations; i++) {
                log_info("%s\t", character[i]);
            }
        }
        log_info ("scaling factor\tepsilon\n");
    }

    log_info("\t%3zd\t%5.10f\t\t", iter, -s->f);
    if (log_frequencies) {
        for (i = 0; i < num_annotations; i++) {
            log_info("%.10f\t", parameters[i]);
        }
    }
    log_info("%.10f\t%.e\n", parameters[num_annotations], parameters[num_annotations + 1]);
}


double anti_sigmoid(double x, const double lower_bound, const double upper_bound) {
    /**
     * undoes the transformation of an arbitrary value x in between lower and upper bound done with a sigmoid function.
     */
    return (upper_bound == lower_bound) ? -log(0): -log((upper_bound - lower_bound) / (x - lower_bound) - 1.0);
}

gsl_vector *real_parameters2optimised_parameters(size_t num_annotations, double scale_low, double scale_up,
                                                 double epsilon_low, double epsilon_up,
                                                 double *parameters, bool optimize_frequencies) {
    size_t i;
    size_t n = optimize_frequencies ? (num_annotations + 1) : 2;
    gsl_vector* x = gsl_vector_alloc(n);
    if (optimize_frequencies) {
        double exp_sum = 1.0 / parameters[num_annotations - 1];
        for (i = 0; i < n - 2; i++) {
            gsl_vector_set(x, i, log(parameters[i] * exp_sum));
        }
    }
    gsl_vector_set(x, n - 2, anti_sigmoid(parameters[num_annotations], scale_low, scale_up));
    gsl_vector_set(x, n - 1, anti_sigmoid(parameters[num_annotations + 1], epsilon_low, epsilon_up));
    return x;
}

double
minus_loglikelihood (const gsl_vector *v, void *params, double* cur_parameters, bool optimize_frequencies, Tree* s_tree)
{
    /**
     * Calculates the -log likelihood value.
     */
    double *p = (double *)params;
    size_t num_annotations = (size_t) p[0];
    double scale_low = p[1], scale_up = p[2], epsilon_low = p[3], epsilon_up = p[4];
    optimised_parameters2real_parameters(v, num_annotations, scale_low, scale_up, epsilon_low, epsilon_up,
                                         cur_parameters, optimize_frequencies);

    return -calculate_bottom_up_likelihood(s_tree, num_annotations, cur_parameters, TRUE);
}
void
d_minus_loglikelihood (const gsl_vector *v, void *params, gsl_vector *df, double* cur_parameters,
                       double cur_minus_log_likelihood, bool optimize_frequencies, Tree* s_tree)
{
    /** Fills in the gradient vector for each of the parameters.
     * parameters = [frequency_char_1, .., frequency_char_n, scaling_factor, epsilon].
     * cur_minus_log_likelihood in the given point can be pre-specified,
     * otherwise should be put to a negative value to show that it needs recalculation.
     */
    double *p = (double *)params;
    size_t num_annotations = (size_t) p[0];
    double scale_low = p[1], scale_up = p[2], epsilon_low = p[3], epsilon_up = p[4];
    double diff_log_likelihood;
    size_t i;

    // if the cur_minus_log_likelihood is already given, let's not recalculate it
    // otherwise it is negative to show that we need to recalculate it.
    if (cur_minus_log_likelihood < 0) {
        optimised_parameters2real_parameters(v, num_annotations, scale_low, scale_up, epsilon_low, epsilon_up,
                                             cur_parameters, optimize_frequencies);
        cur_minus_log_likelihood = -calculate_bottom_up_likelihood(s_tree, num_annotations, cur_parameters, TRUE);
    }

    size_t n = optimize_frequencies ? (num_annotations + 1): 2;
    for (i = 0; i < n; i++) {
        /* create a v + delta v array, where all the values but the i-th are the same as in the initial v,
         * and the i-th value is increased by the corresponding step.*/
        gsl_vector_set(v, i, gsl_vector_get(v, i) + GRADIENT_STEP);
        optimised_parameters2real_parameters(v, num_annotations, scale_low, scale_up, epsilon_low, epsilon_up,
                                             cur_parameters, optimize_frequencies);
        diff_log_likelihood = -calculate_bottom_up_likelihood(s_tree, num_annotations, cur_parameters, TRUE)
                              - cur_minus_log_likelihood;
        /* set the corresponding gradient value*/
        gsl_vector_set(df, i, diff_log_likelihood / GRADIENT_STEP);
        /* put the i-th value back to the one from the initial v.*/
        gsl_vector_set(v, i, gsl_vector_get(v, i) - GRADIENT_STEP);
    }
}

double my_f(const gsl_vector *v, void *params);
void my_df(const gsl_vector *v, void *params, gsl_vector *df);
void my_fdf(const gsl_vector *v, void *params, double *f, gsl_vector *df);



double _minimize_params(Tree* s_tree, size_t num_annotations, double *parameters, char **character,
        bool optimize_frequencies, double scale_low, double scale_up, double epsilon_low, double epsilon_up) {
    /**
     * Optimises the following parameters:
     * parameters = [frequency_char_1, .., frequency_char_n, scaling_factor, epsilon],
     * using BFGS algorithm.
     * The parameters variable is updated to contain the optimal parameters found.
     * The optimal value of the likelihood is returned.
     */

    size_t iter = 0;
    int status;

    // set initial values to random ones within bounds
    parameters[num_annotations] = random_double(scale_low, scale_up);
    parameters[num_annotations + 1] = random_double(epsilon_low, epsilon_up);

    if (optimize_frequencies) {
        double upper_bound = 1.0;
        for (int i = 0; i < num_annotations; i++) {
            parameters[i] = random_double(0.0, upper_bound);
            upper_bound -= parameters[i];
        }
    }

    log_info("Scaling factor can vary between %.10f and %.10f, starting at %.10f.\n", scale_low, scale_up,
             parameters[num_annotations]);
    log_info("Epsilon can vary between %.e and %.e, starting at %.e.\n", epsilon_low, epsilon_up,
             parameters[num_annotations + 1]);

    size_t n = (size_t) (optimize_frequencies ? (num_annotations + 1) : 2);

    gsl_multimin_fdfminimizer *s;
    /* Parameters: num_annotations, scale_low, scale_up, epsilon_low, epsilon_up, mu. */
    double par[5] = {(double) num_annotations, scale_low, scale_up, epsilon_low, epsilon_up};

    double my_f(const gsl_vector *v, void *params) {
        return minus_loglikelihood(v, params, parameters, optimize_frequencies, s_tree);
    }

    void my_df(const gsl_vector *v, void *params, gsl_vector *df) {
        return d_minus_loglikelihood(v, params, df, parameters, -1, optimize_frequencies, s_tree);
    }

    void my_fdf(const gsl_vector *v, void *params, double *f, gsl_vector *df) {
        *f = my_f(v, params);
        d_minus_loglikelihood(v, params, df, parameters, *f, optimize_frequencies, s_tree);
    }

    gsl_vector *x;
    gsl_multimin_function_fdf my_func;
    my_func.n = n;
    my_func.f = my_f;
    my_func.df = my_df;
    my_func.fdf = my_fdf;
    my_func.params = par;

    /* Starting point */
    x = real_parameters2optimised_parameters(num_annotations, scale_low, scale_up, epsilon_low,  epsilon_up,
                                             parameters,  optimize_frequencies);
    s = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs2, n);

    double step_size = .1;
    double tol = .1;
    gsl_multimin_fdfminimizer_set(s, &my_func, x, step_size, tol);

    log_cur_parameter_values(num_annotations, parameters, optimize_frequencies, iter, s, character);

    double epsabs = MIN(1e-3, epsilon_low / 2);
    do
    {
        iter++;
        status = gsl_multimin_fdfminimizer_iterate(s);

        if (status) {
            // if the iteration is not making progress towards solution let's try to reduce the step size
            if (GSL_ENOPROG == status && step_size > GRADIENT_STEP) {
                step_size /= 10.0;
                iter--;
                status = GSL_CONTINUE;
                gsl_multimin_fdfminimizer_set(s, &my_func, gsl_multimin_fdfminimizer_x(s), step_size, tol);
                log_info("\t\t(decreased the step size to %.1e)\n", step_size);
                continue;
            }
            log_info("\t\t(stopping minimization as %s)\n", gsl_strerror(status));
            break;
        }

        status = gsl_multimin_test_gradient(s->gradient, epsabs);
        log_cur_parameter_values(num_annotations, parameters, optimize_frequencies, iter, s, NULL);

        if (status == GSL_SUCCESS) {
            // let's adjust the tolerance to make sure we are at the minimum
            if (iter < 100 && epsabs > 1e-5) {
                epsabs /= 10.0;
                status = GSL_CONTINUE;
                log_info("\t\t(found an optimum candidate, but to be sure decreased the gradient tolerance to %.1e)\n",
                         epsabs);
            } else {
                log_info("\t\t(optimum found!)\n");
            }
        }
    }
    while (status == GSL_CONTINUE && iter < 500);

    /* Make sure that the parameters contain the best value */
    optimised_parameters2real_parameters(gsl_multimin_fdfminimizer_x(s), num_annotations, scale_low, scale_up,
                                         epsilon_low, epsilon_up, parameters, optimize_frequencies);
    double optimum = -gsl_multimin_fdfminimizer_minimum(s);

    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);
    return optimum;
}

double minimize_params(Tree* s_tree, size_t num_annotations, double *parameters, char **character,
        bool optimize_frequencies, double scale_low, double scale_up, double epsilon_low, double epsilon_up) {
    /**
     * Optimises the following parameters:
     * parameters = [frequency_char_1, .., frequency_char_n, scaling_factor, epsilon],
     * using BFGS algorithm.
     * The parameters variable is updated to contain the optimal parameters found.
     * The optimal value of the likelihood is returned.
     */
    double* best_parameters = calloc(num_annotations + 2, sizeof(double));
    double optimum;
    bool relaxed_bounds = FALSE;

    for (int i = 0; i < 5; i++) {
        log_info("\nOptimising parameters, iteration %d out of 5\n", i + 1);
        double res = _minimize_params(s_tree, num_annotations, parameters, character, optimize_frequencies,
                                          scale_low, scale_up, epsilon_low, epsilon_up);

        if (relaxed_bounds == FALSE) {
            // if we hit the bound with SF or Epsilon, let's relax the bound a bit and reoptimise.
            if ((fabs(scale_up - scale_low) > 1e-5)
            && ((fabs(parameters[num_annotations] - scale_up) < 1e-5)
            || (fabs(parameters[num_annotations] - scale_low) < 1e-5))) {
                    relaxed_bounds = TRUE;
                    scale_up *= 2;
                    scale_low /= 2;
            }
            if ((fabs(epsilon_up - epsilon_low) > 1e-5)
            && ((fabs(parameters[num_annotations + 1] - epsilon_low) < 1e-5)
            || (fabs(parameters[num_annotations + 1] - epsilon_up) < 1e-5))) {
                    relaxed_bounds = TRUE;
                    epsilon_low /= 2;
            }
        }

        if ((i == 0) || (res > optimum)) {
            memcpy(best_parameters, parameters, (num_annotations + 2) * sizeof(double));
            optimum = res;
        }
    }

    memcpy(parameters, best_parameters, (num_annotations + 2) * sizeof(double));
    return optimum;
}

