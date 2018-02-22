#include "pastml.h"
#include "scaling.h"


int *calculate_sum_down(const Node *nd, const Node *root, int num_annotations);

int get_max(int *array, int n) {
    /**
     * Finds the maximum in an array of positive integers.
     */
    int max_value = 0;
    for (int i = 0; i < n; i++) {
        if (max_value < array[i]) {
            max_value = array[i];
        }
    }
    return max_value;
}

void normalize(double *array, int n) {
    /**
     * Divides array members by their sum.
     */
    double sum = 0.0;;
    for (int i = 0; i < n; i++) {
        sum += array[i];
    }
    for (int i = 0; i < n; i++) {
        array[i] /= sum;
    }
}

void down_like_marginal(Node *nd, Node *root, int num_tips, int num_annotations, double scale, double *frequency) {
    int child_i, cur_annot_i;
    int curr_scaler_pow, max_factor;

    /*a tip, therefore nothing to do*/
    if (nd->nneigh == 1) {
        return;
    }

    if (nd == root) {
        memcpy((void *) nd->sum_down, (void *) frequency, num_annotations * sizeof(double));
    } else {
        int *tmp_factor = calculate_sum_down(nd, root, num_annotations);
        max_factor = get_max(tmp_factor, num_annotations);
        for (cur_annot_i = 0; cur_annot_i < num_annotations; cur_annot_i++) {
            curr_scaler_pow = max_factor - tmp_factor[cur_annot_i];
            if (curr_scaler_pow != 0) {
                rescale(nd->sum_down, cur_annot_i, curr_scaler_pow);
            }
        }
        free(tmp_factor);

        for (cur_annot_i = 0; cur_annot_i < num_annotations; cur_annot_i++) {
            nd->condlike_mar[cur_annot_i] =
                    nd->sum_down[cur_annot_i] * nd->up_like[cur_annot_i] * frequency[cur_annot_i];
        }

        upscale_node_probs(nd->condlike_mar, num_annotations);
        normalize(nd->condlike_mar, num_annotations);
    }

    for (child_i = (nd == root) ? 0 : 1; child_i < nd->nneigh; child_i++) {
        down_like_marginal(nd->neigh[child_i], root, num_tips, num_annotations, scale, frequency);
    }
}

int *calculate_sum_down(const Node *nd, const Node *root, int num_annotations) {
    Node *father = nd->neigh[0];
    Node *other_child;
    int my_id = -1;
    int tmp_father_factor[num_annotations];
    int *tmp_factor = calloc(num_annotations, sizeof(int));
    double tmp_father_lik[num_annotations];

    for (int child_i = 0; child_i < father->nneigh; child_i++) {
        if (father->neigh[child_i] == nd) {
            my_id = child_i;
        }
    }
    //currentNODE has cur_annot_i
    for (int cur_annot_i = 0; cur_annot_i < num_annotations; cur_annot_i++) {
        nd->sum_down[cur_annot_i] = 1.0;
        tmp_factor[cur_annot_i] = 0;

        if (father == root) {
            for (int child_i = 0; child_i < father->nneigh; child_i++) {
                if (child_i != my_id) {
                    other_child = father->neigh[child_i];
                    double prob_down_son = 0.0;
                    //assume other sons have other_annot_i
                    for (int other_annot_i = 0; other_annot_i < num_annotations; other_annot_i++) {
                        /* Calculate the probability of having a branch from the root to a child node child_i,
                         * given that the root is in state cur_annot_i: p_child_branch_from_i = sum_j(p_ij * p_child_j)
                         */
                        nd->rootpij[cur_annot_i][other_annot_i] = 0.0;
                        for (int father_annotation_i = 0;
                             father_annotation_i < num_annotations; father_annotation_i++) {
                            nd->rootpij[cur_annot_i][other_annot_i] +=
                                    nd->pij[cur_annot_i][father_annotation_i] *
                                    other_child->pij[father_annotation_i][other_annot_i];
                        }
                        prob_down_son += nd->rootpij[cur_annot_i][other_annot_i] * other_child->up_like[other_annot_i];
                    }
                    nd->sum_down[cur_annot_i] *= prob_down_son;
                }
            }
        } else {
            //assume father has father_annotation_i
            for (int father_annotation_i = 0; father_annotation_i < num_annotations; father_annotation_i++) {
                tmp_father_lik[father_annotation_i] =
                        nd->pij[cur_annot_i][father_annotation_i] * father->sum_down[father_annotation_i];
                tmp_father_factor[father_annotation_i] = 0;
                // as our father is not root, its first nneigh is our grandfather,
                // and we should iterate over children staring at 1
                for (int child_i = 1; child_i < father->nneigh; child_i++) {
                    if (child_i != my_id) {
                        other_child = father->neigh[child_i];
                        double prob_down_son = 0.0;
                        //assume other sons have other_annot_i
                        for (int other_annot_i = 0; other_annot_i < num_annotations; other_annot_i++) {
                            prob_down_son += other_child->pij[father_annotation_i][other_annot_i]
                                             * other_child->up_like[other_annot_i];
                        }
                        tmp_father_lik[father_annotation_i] *= prob_down_son;
                    }
                    if (tmp_father_lik[father_annotation_i] < LIM_P) {
                        int curr_scaler_pow1 = get_scaling_pow(tmp_father_lik[father_annotation_i]);
                        tmp_father_factor[father_annotation_i] += curr_scaler_pow1;
                        rescale(tmp_father_lik, father_annotation_i, curr_scaler_pow1);
                    }
                }
            }
            int max_father_factor = get_max(tmp_father_factor, num_annotations);
            for (int father_annotation_i = 0; father_annotation_i < num_annotations; father_annotation_i++) {
                int curr_scaler_pow = max_father_factor - tmp_father_factor[father_annotation_i];
                if (curr_scaler_pow != 0) {
                    rescale(tmp_father_lik, father_annotation_i, curr_scaler_pow);
                }
                nd->sum_down[cur_annot_i] += tmp_father_lik[father_annotation_i];
            }
            tmp_factor[cur_annot_i] = max_father_factor;
        }
    }
    return tmp_factor;
}


