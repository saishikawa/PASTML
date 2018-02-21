#include "pastml.h"
#include "nrutil.h"

#define NR_END 1
#define FREE_ARG char*

void nrerror(char error_text[]) {
    /**
     * Numerical Recipes standard error handler
     */
    fprintf(stderr, "%s\n", error_text);
}

double *vector(long nl, long nh) {
    /**
     * allocate a double vector with subscript range v[nl..nh]
     */
    double *v;
    v = (double *) malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(double)));
    if (!v) nrerror("allocation failure in vector()");
    return v - nl + NR_END;
}

void free_vector(double *v, long nl) {
    /**
     * free a double vector allocated with vector()
     */
    free((FREE_ARG) (v + nl - NR_END));
}

