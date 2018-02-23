#include <stdio.h>
#include <gsl/gsl_multimin.h>
#include "lik.h"

#define STEP_EPSILON 1.0e-6

double softmax(double* xs, int n) {
    /**
     * transforms an array of n arbitrary values x in such a way that all of them become between 0 and 1 and sum to 1,
     * using the softmax function.
     */
    double exp_sum = 0.0;
    for (int i = 0; i < n; i++) {
        xs[i] = exp(xs[i]);
        exp_sum += xs[i];
    }
    for (int i = 0; i < n; i++) {
        xs[i] /= exp_sum;
    }
}


double sigmoid(double x, const double lower_bound, const double upper_bound) {
    /**
     * transforms an arbitrary value x to be in between lower and upper bound, using a sigmoid function.
     */
    return lower_bound + (upper_bound - lower_bound) / (1.0 + exp(-x));
}

double anti_sigmoid(double x, const double lower_bound, const double upper_bound) {
    /**
     * undoes the transformation of an arbitrary value x in between lower and upper bound done with a sigmoid function.
     */
    return -log((upper_bound - lower_bound) / (x - lower_bound) - 1.0);
}

void *get_likelihood_parameters(const gsl_vector *v, int num_annotations, double scale_low, double scale_up,
                                double epsilon_low, double epsilon_up, double* cur_parameters, char* model) {
    size_t j;

    if (strcmp("F81", model) == 0) {
        /* 1. Frequencies */
        for (j = 0; j < num_annotations; j++) {
            cur_parameters[j] = gsl_vector_get(v, j);
        }
        softmax(cur_parameters, num_annotations);
    }
    int scaling_factor_index = (strcmp("F81", model) == 0) ? num_annotations: 0;
    /* 2. Scaling factor */
    cur_parameters[num_annotations] = sigmoid(gsl_vector_get(v, scaling_factor_index), scale_low, scale_up);

    /* 3. Epsilon */
    cur_parameters[num_annotations + 1] = sigmoid(gsl_vector_get(v, scaling_factor_index + 1), epsilon_low, epsilon_up);
}

double
minus_loglikelihood (const gsl_vector *v, void *params, double* cur_parameters, char* model, Node* root)
{
    /**
     * Calculates the -log likelihood value.
     */
    double *p = (double *)params;
    int num_annotations = (int) p[0];
    double scale_low = p[1], scale_up = p[2], epsilon_low = p[3], epsilon_up = p[4];
    get_likelihood_parameters(v, num_annotations, scale_low, scale_up, epsilon_low, epsilon_up, cur_parameters, model);

    return -calc_lik_bfgs(root, num_annotations, cur_parameters);
}
void
d_minus_loglikelihood (const gsl_vector *v, void *params, gsl_vector *df, double* cur_parameters,
                       double cur_minus_log_likelihood, char* model, Node* root)
{
    /** Fills in the gradient vector for each of the parameters.
     * parameters = [frequency_char_1, .., frequency_char_n, scaling_factor, epsilon].
     * cur_minus_log_likelihood in the given point can be pre-specified,
     * otherwise should be put to a negative value to show that it needs recalculation.
     */
    double *p = (double *)params;
    int num_annotations = (int) p[0];
    double scale_low = p[1], scale_up = p[2], epsilon_low = p[3], epsilon_up = p[4];
    double diff_log_likelihood;

    // if the cur_minus_log_likelihood is already given, let's not recalculate it
    // otherwise it is negative to show that we need to recalculate it.
    if (cur_minus_log_likelihood < 0) {
        get_likelihood_parameters(v, num_annotations, scale_low, scale_up, epsilon_low, epsilon_up,
                                  cur_parameters, model);
        cur_minus_log_likelihood = -calc_lik_bfgs(root, num_annotations, cur_parameters);
    }

    int n = (strcmp("F81", model) == 0) ? (num_annotations + 2): 2;
    for (size_t i = 0; i < n; i++) {
        /* create a next_step_parameters array, where all the values but the i-th are the same as in parameters,
         * and the i-th value is increased by the corresponding step.*/
        if (i > 0) {
            gsl_vector_set(v, i - 1, gsl_vector_get(v, i - 1) - STEP_EPSILON);
        }
        gsl_vector_set(v, i, gsl_vector_get(v, i) + STEP_EPSILON);
        get_likelihood_parameters(v, num_annotations, scale_low, scale_up, epsilon_low, epsilon_up,
                                  cur_parameters, model);

        diff_log_likelihood = -calc_lik_bfgs(root, num_annotations, cur_parameters) - cur_minus_log_likelihood;

        /* calculate the gradients*/
        gsl_vector_set(df, i, diff_log_likelihood / STEP_EPSILON);
    }
    gsl_vector_set(v, n - 1, gsl_vector_get(v, n - 1) - STEP_EPSILON);
}

double minimize_params(Node* root, int num_annotations, double *parameters, char **character, char *model, double scale_low,
                       double scale_up, double epsilon_low, double epsilon_up) {
    /**
     * Optimises the following parameters:
     * parameters = [frequency_char_1, .., frequency_char_n, scaling_factor, epsilon],
     * using BFGS algorithm.
     * If model is JC, the frequences are not optimised.
     * The parameters variable is updated to contain the optimal parameters found.
     * The optimal value of the likelihood is returned.
     */

    size_t iter = 0;
    int status;

    size_t n = (size_t) ((strcmp("JC", model) == 0) ? 2 : (num_annotations + 2));

    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;
    /* Parameters: num_annotations, scale_low, scale_up, epsilon_low, epsilon_up, mu. */
    double par[5] = {(double) num_annotations, scale_low, scale_up, epsilon_low, epsilon_up};

    double my_f(const gsl_vector *v, void *params) {
        return minus_loglikelihood(v, params, parameters, model, root);
    }

    void my_df(const gsl_vector *v, void *params, gsl_vector *df) {
        return d_minus_loglikelihood (v, params, df, parameters, -1, model, root);
    }

    void my_fdf (const gsl_vector *v, void *params, double *f, gsl_vector *df) {
        *f = minus_loglikelihood(v, params, parameters, model, root);
        d_minus_loglikelihood (v, params, df, parameters, *f, model, root);
    }

    gsl_vector *x;
    gsl_multimin_function_fdf my_func;
    my_func.n = n;
    my_func.f = my_f;
    my_func.df = my_df;
    my_func.fdf = my_fdf;
    my_func.params = par;

    /* Starting point */
    x = gsl_vector_alloc(n);
    if (strcmp("F81", model) == 0) {
        for (size_t j = 0; j < num_annotations; j++) {
            gsl_vector_set(x, j, log(parameters[j]));
        }
    }
    gsl_vector_set(x, n - 2, anti_sigmoid(parameters[num_annotations], scale_low, scale_up));
    gsl_vector_set(x, n - 1, anti_sigmoid(parameters[num_annotations + 1], epsilon_low, epsilon_up));


    T = gsl_multimin_fdfminimizer_vector_bfgs2;
    s = gsl_multimin_fdfminimizer_alloc(T, n);

    gsl_multimin_fdfminimizer_set(s, &my_func, x, 1, .1);

    printf ("step\tlog-lh\t\t");
    if (strcmp("F81", model) == 0) {
        for (size_t j = 0; j < num_annotations; j++) {
            printf("%s\t", character[j]);
        }
    }
    printf ("scaling\tepsilon\n");
    do
    {
        iter++;
        status = gsl_multimin_fdfminimizer_iterate(s);

        if (status) {
            printf("Stopping minimization as %s.\n", gsl_strerror(status));
            break;
        }

        status = gsl_multimin_test_gradient(s->gradient, .1);

        if (status == GSL_SUCCESS) {
            printf("Minimum found!\n");
        }
        get_likelihood_parameters(s->x, num_annotations, scale_low, scale_up, epsilon_low, epsilon_up, parameters,
                                  model);

        printf ("%5d\t%10.5f\t\t", iter, -s->f);
        if (strcmp("F81", model) == 0) {
            for (size_t j = 0; j < num_annotations; j++) {
                printf("%.5f\t", parameters[j]);
            }
        }
        printf ("%.5f\t%.5f\n", parameters[num_annotations], parameters[num_annotations + 1]);
    }
    while (status == GSL_CONTINUE && iter < 150);
    double optimum = -s->f;

    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);
    return optimum;
}
