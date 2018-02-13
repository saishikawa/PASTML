
#include <errno.h>
#include "pastml.h"
#include "output_states.h"
#include "output_tree.h"

extern Tree *s_tree;
extern Node *root;


void calc_correct(Node *nd, int nb, int nbanno, char **character) {
    int i, j, tmpnum, node_start;
    double tmp_correct = 10.0, local_cor;

    if (nd->nneigh == 1) {
        return;
    }

    if (nd == root) {
        node_start = 0;
    } else {
        node_start = 1;
    }

    for (i = node_start; i < nd->nneigh; i++) {
        calc_correct(nd->neigh[i], nb, nbanno, character);
    }

    /*local_fraction_optimize*/
    for (i = 0; i < nbanno; i++) {
        local_cor = 0.;
        for (j = 0; j < nbanno; j++) {
            if (j <= i) {
                local_cor +=
                        (nd->marginal[j] - 1.0 / ((double) i + 1.0)) * (nd->marginal[j] - 1.0 / ((double) i + 1.0));
            } else {
                local_cor += (nd->marginal[j] - 0.0) * (nd->marginal[j] - 0.0);
            }
        }
        if (tmp_correct > local_cor) {
            tmp_correct = local_cor;
            tmpnum = i + 1;
        }
    }
    for (i = 0; i < nbanno; i++) {
        if (i < tmpnum) {
            nd->local_flag[i] = 1;
        } else {
            nd->local_flag[i] = 0;
        }
    }
    if (nd == root) {
        return;
    }
}

void order_marginal(Node *nd, int nb, int nbanno) {
    int i, j, k, tmpnum, zerostate[nbanno], zero_count = 0, node_start;
    double tmpmarginal[nbanno], tmpbest;
    static int count = 0;

    zero_count = 0;
    if (nd->nneigh == 1) {
        count++;
        return;
    } else if (nd == root) {
        for (j = 0; j < nbanno; j++) {
            tmpmarginal[j] = root->mar_prob[j];
            if (tmpmarginal[j] == 0.0) {
                zerostate[zero_count] = j;
                zero_count++;
            }
        }
        zero_count = 0;
        for (j = 0; j < nbanno; j++) {
            tmpbest = 0.;
            tmpnum = -1;
            for (k = 0; k < nbanno; k++) {
                if (tmpmarginal[k] > tmpbest) {
                    tmpbest = tmpmarginal[k];
                    tmpnum = k;
                }
            }
            nd->marginal[j] = tmpbest;
            tmpmarginal[tmpnum] = 0.;
            nd->tmp_best[j] = tmpnum;
            if (tmpbest == 0.0) {

                nd->tmp_best[j] = zerostate[zero_count];
                zero_count++;
            }
        }
        count++;
    } else {
        for (j = 0; j < nbanno; j++) {
            tmpmarginal[j] = nd->condlike_mar[j];
            if (nd->condlike_mar[j] == 0.0) {
                zerostate[zero_count] = j;
                zero_count++;
            }
        }
        zero_count = 0;
        for (j = 0; j < nbanno; j++) {
            tmpbest = 0.;
            tmpnum = -1;
            for (k = 0; k < nbanno; k++) {
                if (tmpmarginal[k] > tmpbest) {
                    tmpbest = tmpmarginal[k];
                    tmpnum = k;
                }
            }
            nd->marginal[j] = tmpbest;
            tmpmarginal[tmpnum] = 0.;
            nd->tmp_best[j] = tmpnum;
            if (tmpbest == 0.0) {
                nd->tmp_best[j] = zerostate[zero_count];
                zero_count++;
            }
        }
        count++;
    }

    if (nd == root) {
        node_start = 0;
    } else {
        node_start = 1;
    }

    for (i = node_start; i < nd->nneigh; i++) {
        order_marginal(nd->neigh[i], nb, nbanno);
    }

    if (nd == root) count = 0;
}

int make_samples(int num_tips, int num_anno, char **character, double *parameter, char *out_annotation_file_name,
                 char *out_tree_name) {

    FILE *fptre, *fpanno;

    order_marginal(root, num_tips, num_anno);
    calc_correct(root, num_tips, num_anno, character);

    fptre = fopen(out_tree_name, "w");
    if (!fptre) {
        fprintf(stderr, "Output tree file %s is impossible to access.", out_tree_name);
        fprintf(stderr, "Value of errno: %d\n", errno);
        fprintf(stderr, "Error opening the file: %s\n", strerror(errno));
        return ENOENT;
    }
    /*write_result_tree*/
    if (EXIT_SUCCESS != write_nh_tree(s_tree, fptre, parameter[num_anno], parameter[num_anno + 1])) {
        return EXIT_FAILURE;
    }
    fclose(fptre);
    printf("Scaled tree with internal node ids is written to %s\n", out_tree_name);

    fpanno = fopen(out_annotation_file_name, "w");
    if (!fpanno) {
        fprintf(stderr, "Output annotation file %s is impossible to access.", out_annotation_file_name);
        fprintf(stderr, "Value of errno: %d\n", errno);
        fprintf(stderr, "Error opening the file: %s\n", strerror(errno));
        return ENOENT;
    }
    output_state_anc_PP(root, num_tips, num_anno, character, fpanno);
    output_state_tip_PP(root, num_tips, num_anno, character, fpanno);
    fclose(fpanno);
    printf("Predictions for all internal nodes and the root are written to %s in csv format\n",
           out_annotation_file_name);

    return EXIT_SUCCESS;

}
