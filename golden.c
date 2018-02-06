#include "pastml.h"
#include "lik.h"

#define SIGMA 1e-5
#define sqrt5 2.236067977499789696
#define ALPHA 1e-5
#define gratio 1.6180339887498948482045868343656

extern Tree *s_tree;
extern Node *root;

int Stopping_Rule(double x0, double x1) {
    double xm = 0.5 * fabs(x0 + x1);

    if (xm <= 1.0) {
        if (fabs(x1 - x0) < ALPHA) {
            return 0;
        } else {
            return 1;
        }
    } else {
        if (fabs(x1 - x0) < ALPHA * xm) {
            return 0;
        } else {
            return 1;
        }
    }
}

void golden(Node *nd, char **tipnames, int *states, int nb, int nbanno, double mu, char *model, double *p, double ub) {
    double lb = 0.05 / s_tree->avgbl, x1, x2, fx1 = 0., fx2 = 0., diff;
    int count = 0;

    x1 = ub + (lb - ub) / gratio;
    p[nbanno] = x1;
    calc_lik_bfgs(root, tipnames, states, nb, nbanno, mu, model, p, &fx1);
    fx1 = fx1 * (-1.0);
    x2 = lb + (ub - lb) / gratio;
    p[nbanno] = x2;
    calc_lik_bfgs(root, tipnames, states, nb, nbanno, mu, model, p, &fx2);
    fx2 = fx2 * (-1.0);
    diff = fx1 - fx2;

    while (Stopping_Rule(lb, ub) == 1) {
        if (fabs(diff) < SIGMA) break;
        if (diff < 0.0) {
            lb = x1;
            x1 = x2;
            fx1 = fx2;
            x2 = lb + (ub - lb) / gratio;
            p[nbanno] = x2;
            calc_lik_bfgs(root, tipnames, states, nb, nbanno, mu, model, p, &fx2);
            fx2 = fx2 * (-1.0);
            diff = fx1 - fx2;
        } else if (diff > 0.0) {
            ub = x2;
            x2 = x1;
            fx2 = fx1;
            x1 = ub + (lb - ub) / gratio;
            p[nbanno] = x1;
            calc_lik_bfgs(root, tipnames, states, nb, nbanno, mu, model, p, &fx1);
            fx1 = fx1 * (-1.0);
            diff = fx1 - fx2;
        }
        count++;
    }
    p[nbanno] = (x1 + x2) / 2.0;
    //printf("fx1 = %lf, fx2 = %lf\n",fx1, fx2);
}
