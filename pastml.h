#ifndef PASTML_H
#define PASTML_H

/* this is needed for time.h */
#define _POSIX_C_SOURCE 199309L

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdbool.h>


#define MAXLNAME 255
#define MAXNSP 50000
#define MAX_TREELENGTH    10000000 /* more or less 10MB for a tree file in NH format */
#define MAX_NAMELENGTH        255    /* max length of a taxon name */
#define POW (-500)
#define LIM_P pow(2, POW)
#define LOG2 0.69314718055994528623
#define MIN(a, b) ((a)<(b)?(a):(b))
#define MAX(a, b) ((a)>(b)?(a):(b))
#define JOINT "joint"
#define MARGINAL "marginal"
#define MARGINAL_APPROXIMATION "marginal_approx"
#define MAX_POSTERIORI "max_posteriori"
#define DOWNPASS "downpass"
#define ACCTRAN "acctran"
#define DELTRAN "deltran"
#define JC "JC"
#define F81 "F81"
#define EFT "EFT"

typedef struct __Node {
    char *name;
    int id;            /* unique id attributed to the node */
    size_t nb_neigh;    /* number of neighbours */
    struct __Node **neigh;    /* neighbour nodes */

    double **pij;           /* probability of substitution from i to j */
    double *bottom_up_likelihood;       /* conditional likelihoods at the node*/
    long *parsimony_states;
    long *down_parsimony_states;
    long *up_parsimony_states;
    double *result_probs;
    size_t *best_states;
    double *top_down_likelihood;
    int *scaling_factor_down;
    int *scaling_factor_up;
    double branch_len;
    bool unknown_state;
} Node;


typedef struct __Tree {
    Node **nodes;            /* array of node pointers */
    Node *root;            /* the root or pseudo-root node */
    size_t nb_nodes;
    size_t nb_edges;
    size_t nb_taxa;
    int next_avail_node_id;
    double avg_branch_len;
    double max_branch_len;
    double min_branch_len;
    double avg_tip_branch_len;
    int num_zero_tip_branches;
    double min_tip_branch_len;
} Tree;

#endif // PASTML_H
