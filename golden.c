#include "pastml.h"
#include "lik.h"

#define SIGMA 1e-5
#define gratio 1.6180339887498948482045868343656

extern Tree *s_tree;
extern Node *root;

double f(double x, char *const *tipnames, const int *states, int num_tips, int num_annotations, double mu,
         double *parameters);

void golden(char **tipnames, int *states, int num_tips, int num_annotations, double mu, double *parameters,
            double lower_bound, double upper_bound) {
    /**
     * Maximises the scaling factor with the golden-section search,
     * a technique for finding the extremum of a strictly unimodal function.
     *
     * parameters = [frequency_char_1, .., frequency_char_n, scaling_factor, epsilon].
     * mu = 1 / (1 - (frequency_1^2 + ... + frequency_n^2)).
     */
    double x_start, x_end, x_mid, x_next, fx_start, fx_end, fx_mid, fx_next;

    x_start = lower_bound;
    x_end = upper_bound;
    x_mid = lower_bound + (x_end - x_start) / gratio;
    fx_start = f(x_start, tipnames, states, num_tips, num_annotations, mu, parameters);
    fx_end = f(x_end, tipnames, states, num_tips, num_annotations, mu, parameters);
    fx_mid = f(x_mid, tipnames, states, num_tips, num_annotations, mu, parameters);
    
    /* if current max is at x_end, return it */
    if (fx_mid > fx_start && fx_end > fx_mid) {
        parameters[num_annotations] = x_end;
        return;
    }
    /* if current max is at x_start, return it */
    if (fx_mid < fx_start && fx_end < fx_mid) {
        parameters[num_annotations] = x_start;
        return;
    }

    /* the max is inside our interval so let's search for it with the golden section */
    while (TRUE) {
        x_next = x_mid + (x_end - x_mid) / gratio;
        fx_next = f(x_next, tipnames, states, num_tips, num_annotations, mu, parameters);
        if (fabs(fx_next - fx_mid) < SIGMA) {
            parameters[num_annotations] = x_mid;
            return;
        }
        if (fx_next > fx_mid) {            
            x_mid = x_next;
            fx_mid = fx_next;
        } else {
            x_end = x_next;
        }
    }
}

double f(double x, char *const *tipnames, const int *states, int num_tips, int num_annotations, double mu,
         double *parameters) {
    double fx_start;
    parameters[num_annotations] = x;
    fx_start = -calc_lik_bfgs(root, tipnames, states, num_tips, num_annotations, mu, parameters);
    return fx_start;
}
