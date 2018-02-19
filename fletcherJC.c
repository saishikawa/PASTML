#include "pastml.h"
#include "nrutil.h"
#include "lik.h"

int ncomjc;
double *pcomjc, *xicomjc, scale_upjc, scale_lowjc, *store_pjc, *best_pjc;
extern Tree *s_tree;
extern Node *root;

void SHFT2(double *a, double *b, double *c, const double *d) {

    *a = (*b);
    *b = (*c);
    *c = (*d);

}

void gradientJC(char **tipnames, int *states, int nb, int nbanno, double mu, char *model, double *p, double *g, int n,
                double *frequency) {
    int i, j;
    double pluslnl, oldlnl, *likp;

    likp = vector(0, nbanno + 1);
    for (i = 0; i < nbanno; i++) {
        likp[i] = frequency[i];
    }
    likp[nbanno] = p[0];
    likp[nbanno + 1] = p[1];
    oldlnl = calc_lik_bfgs(root, tipnames, states, nb, nbanno, mu, likp);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (j == i) {
                if (j == 0) {
                    likp[nbanno] = p[j] + STEP2;
                } else if (j == 1) {
                    likp[nbanno + 1] = p[j] + STEP3;
                }
            }
        }
        pluslnl = calc_lik_bfgs(root, tipnames, states, nb, nbanno, mu, likp);
        if (i == 0) {
            g[i] = (fabs(pluslnl) - fabs(oldlnl)) / (STEP2);
        } else if (i == 1) {
            g[i] = (fabs(pluslnl) - fabs(oldlnl)) / (STEP3);
        }
    }
    free_vector(likp, 0, nbanno + 1);
}


double f1dimJC(char **tipnames, int *states, int nb, int nbanno, double mu, char *model, double x, const double *frequency) {
    int i, j;
    double f, *xt, e1, e2, exp1, exp2, *likp;

    xt = vector(0, ncomjc - 1);
    likp = vector(0, nbanno + 1);

    if (xicomjc[0] == 0.0) {
        exp1 = 1.0;
    } else {
        e1 = log(fabs(xicomjc[0]));
        exp1 = pow(10.0, e1);
    }
    if (xicomjc[1] == 0.0) {
        exp2 = 1.0;
    } else {
        e2 = log(fabs(xicomjc[1]));
        exp2 = pow(10.0, e2);
    }
    xt[0] = (pcomjc[0] * exp1 + x * xicomjc[0]) / exp1;
    xt[1] = (pcomjc[1] * exp2 + x * xicomjc[1]) / exp2;
    if (xt[0] > scale_upjc || xt[0] < scale_lowjc) xt[0] = pcomjc[0];
    if (fabs(xt[1]) > 1.0 || fabs(xt[1]) < SCAL_MIN || xt[1] < 0.0) xt[1] = pcomjc[1];

    for (i = 0; i < nbanno; i++) {
        likp[i] = frequency[i];
    }
    likp[nbanno] = xt[0];
    likp[nbanno + 1] = xt[1];
    f = calc_lik_bfgs(root, tipnames, states, nb, nbanno, mu, likp);
    for (j = 0; j < ncomjc; j++) {
        store_pjc[j] = xt[j];
    }
    free_vector(xt, 0, ncomjc - 1);
    free_vector(likp, 0, nbanno + 1);
    return fabs(f);
}

void
mnbrakJC(char **tipnames, int *states, int nb, int nbanno, double mu, char *model, double *ax, double *bx, double *cx,
         double *fa, double *fb, double *fc, double *frequency) {
/*Given a function func , and given distinct initial points ax and bx , this routine searches in
the downhill direction (defined by the function as evaluated at the initial points) and returns
new points ax , bx , cx that bracket a minimum of the function. Also returned are the function
values at the three points, fa , fb , and fc .*/
    double ulim, u, r, q, fu, dum = 0.0, *adum, *bdum, *cdum, tmp;

    adum = &dum;
    *fa = f1dimJC(tipnames, states, nb, nbanno, mu, model, *ax, frequency);
    *fb = f1dimJC(tipnames, states, nb, nbanno, mu, model, *bx, frequency);
    if (*fb > *fa) {
        SHFT2(adum, ax, bx, adum);
        SHFT2(adum, fb, fa, adum);
    }
    *cx = (*bx) + GOLD * (*bx - *ax);
    *fc = f1dimJC(tipnames, states, nb, nbanno, mu, model, *cx, frequency);
    while (*fb > *fc) {
        r = ((*bx) - (*ax)) * ((*fb) - (*fc));
        q = ((*bx) - (*cx)) * ((*fb) - (*fa));
        u = (*bx) - (((*bx) - (*cx)) * q - ((*bx) - (*ax)) * r) / (2.0 * SIGN(FMAX(fabs(q - r), TINY), q - r));
        ulim = (*bx) + GLIMIT * (*cx - *bx);
        if ((*bx - u) * (u - *cx) > 0.0) {
            fu = f1dimJC(tipnames, states, nb, nbanno, mu, model, u, frequency);
            if (fu < *fc) {
                *ax = (*bx);
                *bx = u;
                *fa = (*fb);
                *fb = fu;
                return;
            } else if (fu > *fb) {
                *cx = u;
                *fc = fu;
                return;
            }
            u = (*cx) + GOLD * ((*cx) - (*bx));
            fu = f1dimJC(tipnames, states, nb, nbanno, mu, model, u, frequency);
        } else if ((*cx - u) * (u - ulim) > 0.0) {
            fu = f1dimJC(tipnames, states, nb, nbanno, mu, model, u, frequency);
            if (fu < *fc) {
                tmp = (*cx) + GOLD * ((*cx) - (*bx));
                bdum = &tmp;
                SHFT2(bx, cx, &u, &bdum);
                tmp = f1dimJC(tipnames, states, nb, nbanno, mu, model, u, frequency);
                cdum = &tmp;
                SHFT2(fb, fc, &fu, &cdum);
            }
        } else if ((u - ulim) * (ulim - *cx) >= 0.0) {
            u = ulim;
            fu = f1dimJC(tipnames, states, nb, nbanno, mu, model, u, frequency);
        } else {
            u = (*cx) + GOLD * ((*cx) - (*bx));
            fu = f1dimJC(tipnames, states, nb, nbanno, mu, model, u, frequency);
        }
        SHFT2(ax, bx, cx, &u);
        SHFT2(fa, fb, fc, &fu);
    }
}

double
brentJC(char **tipnames, int *states, int nb, int nbanno, double mu, char *model, double ax, double bx, double cx,
        double tol, double *xmin, double *frequency) {
/*Given a function f , and given a bracketing triplet of abscissas ax , bx , cx (such that bx is
between ax and cx , and f(bx) is less than both f(ax) and f(cx) ), this routine isolates
the minimum to a fractional precision of about tol using Brent’s method. The abscissa of
the minimum is returned as xmin , and the minimum function value is returned as brent , the
returned function value.*/

    int iter, i;
    double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
    double e = 0.0;

    a = (ax < cx ? ax : cx);
    b = (ax > cx ? ax : cx);
    x = w = v = bx;
    fw = fv = fx = f1dimJC(tipnames, states, nb, nbanno, mu, model, x, frequency);
    for (iter = 0; iter < ITMAX; iter++) {
        xm = 0.5 * (a + b);
        tol2 = 2.0 * (tol1 = tol * fabs(x) + SIGMA);
        if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) {
            *xmin = x;
            return (-1.0) * fx;
        }
        if (fabs(e) > tol1) {
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0.0) p = -p;
            q = fabs(q);
            etemp = e;
            e = d;
            if (fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)) {
                d = CGOLD * (e = (x >= xm ? a - x : b - x));
            } else {
                d = p / q;
                u = x + d;
                if (u - a < tol2 || b - u < tol2) d = SIGN(tol1, xm - x);
            }
        } else {
            d = CGOLD * (e = (x >= xm ? a - x : b - x));
        }
        u = (fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
        fu = f1dimJC(tipnames, states, nb, nbanno, mu, model, u, frequency);
        if (fu <= fx) {
            if (u >= x) a = x; else b = x;
            SHFT2(&v, &w, &x, &u);
            SHFT2(&fv, &fw, &fx, &fu);
            for (i = 0; i < ncomjc; i++) {
                best_pjc[i] = store_pjc[i];
            }
        } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v = u;
                fv = fu;
            }
        }
    }
    nrerror("Warnng: Too many iterations >200 in Brent's method");
    *xmin = x;
    return (-1.0) * fx;
}

void linminJC(char **tipnames, int *states, int nb, int nbanno, double mu, char *model, double *p, double *xi, int n,
              double *fret, double *frequency) {
/*Given an n -dimensional point p[1..n] and an n -dimensional direction xi[1..n] , moves and
resets p to where the function func(p) takes on a minimum along the direction xi from p ,
and replaces xi by the actual vector displacement that p was moved. Also returns as fret
the value of func at the returned location p . This is actually all accomplished by calling the
routines mnbrak and brent .*/
    int i, j;
    double xx, xmin, fx, fb, fa, bx, ax, *pold, f_check, *likp;

    ncomjc = n;
    pcomjc = vector(0, n - 1);
    pold = vector(0, n - 1);
    xicomjc = vector(0, n - 1);
    store_pjc = vector(0, n - 1);
    best_pjc = vector(0, n - 1);
    likp = vector(0, nbanno + 1);
    for (j = 0; j < n; j++) {
        pcomjc[j] = p[j];
        xicomjc[j] = xi[j];
        pold[j] = p[j];
        best_pjc[j] = p[j];
    }
    ax = 0.0;
    xx = 1.0;
    mnbrakJC(tipnames, states, nb, nbanno, mu, model, &ax, &xx, &bx, &fa, &fx, &fb, frequency);
    *fret = brentJC(tipnames, states, nb, nbanno, mu, model, ax, xx, bx, TOL, &xmin, frequency);

    for (j = 0; j < n; j++) {
        p[j] = best_pjc[j];
    }
    if (p[0] > scale_upjc || p[0] < scale_lowjc) p[0] = pold[0];
    if (p[1] > 1.0 || p[1] < SCAL_MIN) p[1] = pold[1];
    for (i = 0; i < nbanno; i++) {
        likp[i] = frequency[i];
    }
    likp[nbanno] = p[0];
    likp[nbanno + 1] = p[1];
    f_check = calc_lik_bfgs(root, tipnames, states, nb, nbanno, mu, likp);
    free_vector(pcomjc, 0, n - 1);
    free_vector(xicomjc, 0, n - 1);
    free_vector(pold, 0, n - 1);
    free_vector(store_pjc, 0, n - 1);
    free_vector(best_pjc, 0, n - 1);
    free_vector(likp, 0, nbanno + 1);
}

void frprmnJC(char **tipnames, int *states, int nb, int nbanno, double mu, char *model, double *p, int n, double ftol,
              int *iter, double *fret, double *frequency) {
/*Given a starting point p[1..n] , Fletcher-Reeves-Polak-Ribiere minimization is performed on a
function func , using its gradient as calculated by a routine dfunc . The convergence tolerance
on the function value is input as ftol . Returned quantities are p (the location of the minimum),
iter (the number of iterations that were performed), and fret (the minimum value of the
function). The routine linmin is called to perform line minimizations.*/

    int i, j, its;
    double gg, gam, fp, dgg, fl;
    double *g, *h, *xi, *likp;

    scale_upjc = 5.0 / s_tree->avg_branch_len;
    scale_lowjc = 0.05 / s_tree->avg_branch_len;
    printf("Scalingfactor upbound = %lf, start point = %.5e\n\n***Fletcher-Reeves-Polak-Ribiere minimization ...\n\n",
           scale_upjc, p[0]);
    g = vector(0, n - 1);
    h = vector(0, n - 1);
    xi = vector(0, n - 1);
    likp = vector(0, nbanno + 1);
    for (i = 0; i < nbanno; i++) {
        likp[i] = frequency[i];
    }
    likp[nbanno] = p[0];
    likp[nbanno + 1] = p[1];
    fp = calc_lik_bfgs(root, tipnames, states, nb, nbanno, mu, likp);
    gradientJC(tipnames, states, nb, nbanno, mu, model, p, xi, n, frequency);
    for (j = 0; j < n; j++) {
        g[j] = -xi[j];
        xi[j] = h[j] = g[j];
    }
    for (its = 0; its < ITMAX_O; its++) {
        *iter = its;
        printf("Step%d, ", its);
        linminJC(tipnames, states, nb, nbanno, mu, model, p, xi, n, fret, frequency);
        for (i = 0; i < n; i++) {
            if (i == 0) printf("Scaling=%.6e, ", p[i]);
            if (i == 1) printf("Epsilon=%.6e, ", p[i]);
        }
        for (i = 0; i < nbanno; i++) {
            likp[i] = frequency[i];
        }
        likp[nbanno] = p[0];
        likp[nbanno + 1] = p[1];
        fl = calc_lik_bfgs(root, tipnames, states, nb, nbanno, mu, likp);
        printf("lnL = %lf\n", fl);
        if (fabs(fl - fp) <= ftol) {
            free_vector(g, 0, n - 1);
            free_vector(h, 0, n - 1);
            free_vector(xi, 0, n - 1);
            free_vector(likp, 0, nbanno + 1);
            return;
        }
        fp = fl;
        gradientJC(tipnames, states, nb, nbanno, mu, model, p, xi, n, frequency);
        dgg = gg = 0.0;
        for (j = 0; j < n; j++) {
            gg += g[j] * g[j];
            dgg += (xi[j] + g[j]) * xi[j]; /* Polak-Ribiere */
        }
        if (g == 0) {
            free_vector(g, 0, n - 1);
            free_vector(h, 0, n - 1);
            free_vector(xi, 0, n - 1);
            free_vector(likp, 0, nbanno + 1);
            return;
        }
        gam = dgg / gg;
        for (j = 0; j < n; j++) {
            g[j] = -xi[j];
            xi[j] = h[j] = g[j] + gam * h[j];
        }
    }
    nrerror("parameter or likelihood is NOT converged within 1000 steps\n*** PASTML returns a possible best solution found during optimization ***\n");
    free_vector(g, 0, n - 1);
    free_vector(h, 0, n - 1);
    free_vector(xi, 0, n - 1);
    free_vector(likp, 0, nbanno + 1);
}

