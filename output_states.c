#include "pastml.h"

extern double global_like;
extern double global_factor;
extern Tree *s_tree;
extern Node *root;


void output_state_anc_PP(Node *nd, int nb, int nbanno, char **character, FILE *outfile) {
    int i, j, count = 0, node_start, num = 0;
    double sum = 0.0;

    if (nd == root) {
        fprintf(outfile, "Internal NodeID");
        for (i = 0; i < nbanno; i++) {
            fprintf(outfile, ", %s", character[i]);
        }
        fprintf(outfile, "\n");
    }
    if (nd->nneigh == 1) {
        return;
    }

    if (nd == root) {
        node_start = 0;
    } else {
        node_start = 1;
    }

    for (i = 0; i < nbanno; i++) {
        if (nd->local_flag[i] == 1) {
            num++;
        } else {
        }
    }

    for (i = node_start; i < nd->nneigh; i++) {
        output_state_anc_PP(nd->neigh[i], nb, nbanno, character, outfile);
    }

    //printf("%s", nd->name);
    for (i = 0; i < nbanno; i++) {
      //printf(", %.5f", nd->marginal[j]);
        if (nd->local_flag[i] == 1) {
            nd->marginal[i] = (double) 1.0 / num;
        } else {
        }
    }
    //printf("\n");

    fprintf(outfile, "%s", nd->name);
        for (i = 0; i < nbanno; i++) {
            for (j = 0; j < nbanno; j++) {

                if (strcmp(character[i], character[nd->tmp_best[j]]) == 0) {
                    if (nd->local_flag[j] == 1) {
                        fprintf(outfile, ", %.5f", nd->marginal[j]);
                    } else {
                        fprintf(outfile, ", 0.0");
                    }
                }
            }
        }
    fprintf(outfile, "\n");

    return;
}

void output_state_tip_PP(Node *nd, int nb, int nbanno, char **character, FILE *outfile) {
    int i, node_start;

    if (nd->nneigh == 1) {
        fprintf(outfile, "%s", nd->name);
        for (i = 0; i < nbanno; i++) {
            if (i == nd->pupko_state) {
                fprintf(outfile, ", 1.0");
            } else {
                fprintf(outfile, ", 0.0");
            }
        }
        fprintf(outfile, "\n");
        return;
    }

    if (nd == root) {
        node_start = 0;
    } else {
        node_start = 1;
    }

    for (i = node_start; i < nd->nneigh; i++) {
        output_state_tip_PP(nd->neigh[i], nb, nbanno, character, outfile);
    }

    return;
}

