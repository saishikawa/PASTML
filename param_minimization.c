#include <stdio.h>
#include <gsl/gsl_multimin.h>
#include "likelihood.h"
#include "logger.h"
#include "param_minimization.h"
#include "states.h"

#define GRADIENT_STEP 1.0e-6


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

void optimised_parameters2real_parameters(const gsl_vector *v, size_t num_annotations, double scale_low, double scale_up,
                                          double *cur_parameters, size_t set_values) {
    size_t i = 0;
    if ((set_values & FREQUENCIES_SET) == 0) {
        /* 1. Frequencies */
        for (; i < num_annotations - 1; i++) {
            cur_parameters[i] = gsl_vector_get(v, i);
        }
        cur_parameters[num_annotations - 1] = 0.0;
        softmax(cur_parameters, num_annotations);
    }
    if ((set_values & SF_SET) == 0) {
        /* 2. Scaling factor */
        cur_parameters[num_annotations] = MIN(MAX(gsl_vector_get(v, i), scale_low), scale_up);
    }
}

void
log_cur_parameter_values(size_t num_annotations, const double *parameters, size_t set_values, size_t iter,
                         double res, char **character) {
    size_t i;

    if (character != NULL) {
        log_info ("\tstep\tlog-lh\t");
        if ((set_values & FREQUENCIES_SET) == 0) {
            for (i = 0; i < num_annotations; i++) {
                log_info("\t%s", character[i]);
            }
        }
        if ((set_values & SF_SET) == 0) {
            log_info("\tscaling factor (state changes per avg branch len)");
        }
        log_info ("\n");
    }

    log_info("\t%3zd\t%5.10f\t", iter, res);
    if ((set_values & FREQUENCIES_SET) == 0) {
        for (i = 0; i < num_annotations; i++) {
            log_info("\t%.10f", parameters[i]);
        }
    }
    if ((set_values & SF_SET) == 0) {
        log_info("\t%.10f", parameters[num_annotations]);
    }
    log_info ("\n");
}

size_t get_num_parameters(size_t num_annotations, size_t set_values) {
    size_t n = num_annotations;
    if ((set_values & FREQUENCIES_SET) != 0) {
        n -= (num_annotations - 1);
    }
    if ((set_values & SF_SET) != 0) {
        n--;
    }
    return n;
}

gsl_vector *real_parameters2optimised_parameters(size_t num_annotations, double scale_low, double scale_up,
                                                 const double *parameters, size_t set_values) {
    size_t i = 0;
    size_t n = get_num_parameters(num_annotations, set_values);
    gsl_vector* x = gsl_vector_alloc(n);
    if ((set_values & FREQUENCIES_SET) == 0) {
        double exp_sum = 1.0 / parameters[num_annotations - 1];
        for (; i < num_annotations - 1; i++) {
            gsl_vector_set(x, i, log(parameters[i] * exp_sum));
        }
    }
    if ((set_values & SF_SET) == 0) {
        gsl_vector_set(x, i, parameters[num_annotations]);
    }
    return x;
}

double
minus_loglikelihood (const gsl_vector *v, void *params, double* cur_parameters, size_t set_values, Tree* s_tree)
{
    /**
     * Calculates the -log likelihood value.
     */
    double *p = (double *)params;
    size_t num_annotations = (size_t) p[0];
    double scale_low = p[1], scale_up = p[2];
    optimised_parameters2real_parameters(v, num_annotations, scale_low, scale_up, cur_parameters, set_values);

    return -calculate_bottom_up_likelihood(s_tree, num_annotations, cur_parameters, true);
}

void
d_minus_loglikelihood (const gsl_vector *v, void *params, gsl_vector *df, const double* cur_parameters,
                       double cur_minus_log_likelihood, size_t set_values, Tree* s_tree)
{
    /** Fills in the gradient vector for each of the parameters.
     * parameters = [frequency_char_1, .., frequency_char_n, scaling_factor].
     * cur_minus_log_likelihood in the given point can be pre-specified,
     * otherwise should be put to a negative value to show that it needs recalculation.
     */
    double *p = (double *)params;
    size_t num_annotations = (size_t) p[0];
    double scale_low = p[1], scale_up = p[2];
    double diff_log_likelihood;
    size_t i;

    double *changed_params = malloc((num_annotations + 1) * sizeof(double));
    memcpy(changed_params, cur_parameters, (num_annotations + 1) * sizeof(double));

    // if the cur_minus_log_likelihood is already given, let's not recalculate it
    // otherwise it is negative to show that we need to recalculate it.
    if (cur_minus_log_likelihood < 0) {
        optimised_parameters2real_parameters(v, num_annotations, scale_low, scale_up, changed_params, set_values);
        cur_minus_log_likelihood = -calculate_bottom_up_likelihood(s_tree, num_annotations, changed_params, true);
    }

    size_t n = get_num_parameters(num_annotations, set_values);
    gsl_vector* v_copy = gsl_vector_alloc(n);
    gsl_vector_memcpy(v_copy, v);
    bool same_params;

    for (i = 0; i < n; i++) {
        /* create a v + delta v array, where all the values but the i-th are the same as in the initial v,
         * and the i-th value is increased by the corresponding step.*/
        gsl_vector_set(v_copy, i, gsl_vector_get(v, i) + GRADIENT_STEP);
        optimised_parameters2real_parameters(v_copy, num_annotations, scale_low, scale_up, changed_params, set_values);
        same_params = true;
        for (size_t j = 0; j < num_annotations + 1; j++) {
            if (fabs(cur_parameters[j] - changed_params[j]) > EPSILON_APPROXIMATING_ZERO) {
                same_params = false;
                break;
            }
        }
        diff_log_likelihood = (same_params) ? 0 :
                (-calculate_bottom_up_likelihood(s_tree, num_annotations, changed_params, true) - cur_minus_log_likelihood);
        /* set the corresponding gradient value*/
        gsl_vector_set(df, i, diff_log_likelihood / GRADIENT_STEP);
        /* put the i-th value back to the one from the initial v.*/
        gsl_vector_set(v_copy, i, gsl_vector_get(v, i));
    }
    gsl_vector_free(v_copy);
    free(changed_params);
}

double _minimize_params(double *parameters, char **character, size_t set_values, gsl_multimin_function_fdf my_func) {
    /**
     * Optimises the following parameters:
     * parameters = [frequency_char_1, .., frequency_char_n, scaling_factor],
     * using BFGS algorithm.
     * The parameters variable is updated to contain the optimal parameters found.
     * The optimal value of the likelihood is returned.
     */
    double* params = (double *) my_func.params;
    size_t num_annotations = (size_t) params[0];
    double scale_low = params[1];
    double scale_up = params[2];
    gsl_vector* x = real_parameters2optimised_parameters(num_annotations,
                                                         scale_low, scale_up, parameters, set_values);
    gsl_multimin_fdfminimizer* s = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs2,
                                                                   get_num_parameters(num_annotations, set_values));

    size_t iter = 0;
    int status;
    double step_size = .1;
    double tol = .1;
    gsl_multimin_fdfminimizer_set(s, &my_func, x, step_size, tol);

    log_cur_parameter_values(num_annotations, parameters, set_values, iter, -s->f, character);

    double epsabs = 1e-3;
    do
    {
        iter++;
        status = gsl_multimin_fdfminimizer_iterate(s);

        optimised_parameters2real_parameters(gsl_multimin_fdfminimizer_x(s), num_annotations, scale_low, scale_up,
                                             parameters, set_values);
        log_cur_parameter_values(num_annotations, parameters, set_values, iter, -s->f, NULL);

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

        if (status == GSL_SUCCESS) {
            // let's adjust the tolerance to make sure we are at the minimum
            if (iter < 100 && epsabs > EPSILON_APPROXIMATING_ZERO) {
                epsabs /= 10.0;
                status = GSL_CONTINUE;
                log_info("\t\t(found an optimum candidate, but to be sure decreased the gradient tolerance to %.1e)\n",
                         epsabs);
                continue;
            }
            log_info("\t\t(optimum found!)\n");
        }
    }
    while (status == GSL_CONTINUE && iter < 500);

    /* Make sure that the parameters contain the best value */
    optimised_parameters2real_parameters(gsl_multimin_fdfminimizer_x(s), num_annotations, scale_low, scale_up,
                                         parameters, set_values);
    double optimum = -gsl_multimin_fdfminimizer_minimum(s);

    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);
    return optimum;
}

double my_f(const gsl_vector *v, void *params);
void my_df(const gsl_vector *v, void *params, gsl_vector *df);
void my_fdf(const gsl_vector *v, void *params, double *f, gsl_vector *df);

void minimize_params(Tree *s_tree, size_t num_annotations, double *parameters, char **character, size_t set_values) {
    /**
     * Optimises the following parameters:
     * parameters = [frequency_char_1, .., frequency_char_n, scaling_factor],
     * using BFGS algorithm.
     * The parameters variable is updated to contain the optimal parameters found.
     * The optimal value of the likelihood is returned.
     */
    double* best_parameters = calloc(num_annotations + 1, sizeof(double));
    double optimum;

    double par[3] = {(double) num_annotations, 0.001, 10};

    double my_f(const gsl_vector *v, void *params) {
        return minus_loglikelihood(v, params, parameters, set_values, s_tree);
    }

    void my_df(const gsl_vector *v, void *params, gsl_vector *df) {
        return d_minus_loglikelihood(v, params, df, parameters, -1, set_values, s_tree);
    }

    void my_fdf(const gsl_vector *v, void *params, double *f, gsl_vector *df) {
        *f = my_f(v, params);
        d_minus_loglikelihood(v, params, df, parameters, *f, set_values, s_tree);
    }

    gsl_multimin_function_fdf my_func;
    my_func.n = get_num_parameters(num_annotations, set_values);
    my_func.f = my_f;
    my_func.df = my_df;
    my_func.fdf = my_fdf;
    my_func.params = par;

    bool some_optimum_found = false;

    if ((set_values & SF_SET) == 0) {
        log_info("Scaling factor can vary between %.3f and %.1f.\n", par[1], par[2]);
    }

    size_t max_rounds = 10;
    /* we will perform at least 10 iterations,
     * if after that our optimal scaling factor will not be within the proposed bounds,
     * we'll perform more iterations (up to 25) trying to find a better one within those bounds. */
    for (size_t i = 1; i <= MIN(max_rounds, 25); i++) {
        log_info("\nOptimising parameters, iteration %d out of %d\n", i, max_rounds);
        if (i == 1) {
            if ((set_values & FREQUENCIES_SET) == 0) {
                log_info("Starting from the observed frequencies...\n");
            }

            if ((set_values & SF_SET) == 0) {
                parameters[num_annotations] = 1.;
            }
        } else if (i == 2 && (set_values & FREQUENCIES_SET) == 0) {
            log_info("Starting from equal frequencies...\n");
            for (size_t j = 0; j < num_annotations; j++) {
                parameters[j] = 1. / num_annotations;
            }

            if ((set_values & SF_SET) == 0) {
                parameters[num_annotations] = 1.;
            }
        } else {
            if ((set_values & FREQUENCIES_SET) == 0) {
                log_info("Starting from random frequencies...\n");
            }
            if ((set_values & FREQUENCIES_SET) == 0) {
                for (size_t j = 0; j < num_annotations; j++) {
                    parameters[j] = random_double(0., 1.);
                }
                normalize(parameters, num_annotations);
            }

            if ((set_values & SF_SET) == 0) {
                parameters[num_annotations] = pow(10., random_double(-3., 1.));
            }
        }

        double res = _minimize_params(parameters, character, set_values, my_func);

        if (!some_optimum_found || (res > optimum)) {
            memcpy(best_parameters, parameters, (num_annotations + 1) * sizeof(double));
            optimum = res;
            log_info("...%s: %.10f...\n", some_optimum_found ? "improved the optimum" : "our first optimum candidate",
                     optimum);
            some_optimum_found = true;
        }
        // in case the optimal scaling factor that we found is not within the initial bounds, let's try for a bit longer.
        if (i >= 10 && i <= 25 && (set_values & SF_SET) == 0 &&
            (best_parameters[num_annotations] <= .001 || best_parameters[num_annotations] >= 10.)) {
            log_info("...the best scaling factor value found so far (%.5f) hits the bound, "
                     "so let's optimise our parameters for a bit longer...\n", parameters[num_annotations]);
            max_rounds += 1;
        }
    }

    memcpy(parameters, best_parameters, (num_annotations + 1) * sizeof(double));
    free(best_parameters);
}