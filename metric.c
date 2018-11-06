#include <stdbool.h>
#include <errno.h>
#include "pastml.h"
#include "tree.h"

int node_check_tips(Node *nd, Tree *s_tree, size_t num_annotations) {
    /**
     * Check tips associated with the current node.
     */
    int my_tips=0, tmp=0;
    size_t i, j, k, l;
    Node* child;
    size_t first_child_index = (nd == s_tree->root) ? 0 : 1;
    double tmp_dist;
    int num_states;
    
    /* recursively calculate tips for children if there are any children */
    for (i = first_child_index; i < nd->nb_neigh; i++) {
        child = nd->neigh[i];
        if (isTip(child)) {
           my_tips ++;
        } else {
           my_tips += node_check_tips(child, s_tree, num_annotations);
        }
    }
    nd->color_tips = calloc(nd->nb_neigh, sizeof(int *));
    nd->local_dist_tips = calloc(nd->nb_neigh, sizeof(double *));
    for (i = first_child_index; i < nd->nb_neigh; i++){
        nd->color_tips[i] = calloc(my_tips, sizeof(int));
        nd->local_dist_tips[i] = calloc(my_tips, sizeof(double));
    }
    for (i = first_child_index; i < nd->nb_neigh; i++){
        for (j = 0; j < my_tips; j++){
            nd->color_tips[i][j] = -1;
            nd->local_dist_tips[i][j] = 0.0;
        }
    }
    tmp=0;
    for (i = first_child_index; i < nd->nb_neigh; i++){
        child = nd->neigh[i];
        tmp=0;
        if (isTip(child)){
            nd->color_tips[i][tmp] = child->tip_state;
            tmp ++;
        } else {
            for (k = 1; k < child->nb_neigh; k++){
                for (j = 0; j < child->num_tips; j++){
                    if (child->color_tips[k][j] > -1){
                        nd->color_tips[i][tmp] = child->color_tips[k][j];
                        tmp ++;
                    }
                }
            }
        }
    }
    tmp=0;
    for (i = first_child_index; i < nd->nb_neigh; i++){
        child = nd->neigh[i];
        tmp=0;
        if (isTip(child)){
            tmp_dist = 0;
            num_states = 0;
            for(j = 0; j < num_annotations; j++){
                if(nd->result_probs[j] > 0){
                    num_states ++;
                    if(j == child->tip_state){
                        tmp_dist += 0;
                    } else {
                        tmp_dist ++;
                    }
                }
            }
            tmp_dist /= num_states;
            nd->local_dist_tips[i][tmp] = tmp_dist;
            tmp ++;
        } else {
            tmp_dist = 0;
            num_states = 0;
            for(k = 0; k < num_annotations; k++){
                if(nd->result_probs[k] > 0){
                    for(l = 0; l < num_annotations; l++){
                        if(child->result_probs[l] > 0){
                            num_states ++;
                            if(l == k){
                                tmp_dist += 0;
                            } else {
                                tmp_dist ++;
                            }
                        }
                    }
                }
            }
            tmp_dist /= num_states;
            for (k = 1; k < child->nb_neigh; k++){
                for (j = 0; j < child->num_tips; j++){
                    if (child->color_tips[k][j] > -1){
                        nd->local_dist_tips[i][tmp] = child->local_dist_tips[k][j] + tmp_dist;
                        tmp ++;
                    }
                }
            }
        }
    }
    nd->num_tips = my_tips;
    /*printf("%s tips %d neigh = %d, tips = %d\n", nd->name, nd->num_tips, nd->nb_neigh);
    for(i = first_child_index; i < nd->nb_neigh; i++){
        for(j = 0; j < nd->num_tips; j++) {
            if(nd->color_tips[i][j] > -1) printf("tip%d side%d = %d, dist = %.3f\n", j, i, nd->color_tips[i][j], nd->local_dist_tips[i][j]);
        }
    }*/
    return my_tips;
}

void calculate_metric_color_pairs(Node *nd, Tree *s_tree, size_t num_annotations, double** metric_color_pairs, int** num_pairs_per_color){
    /**
     * calculate distance between two tips with distinct colors.
     */
    size_t i, j, k, l, m, n;
    Node* child, child_1, child_2;
    size_t first_child_index = (nd == s_tree->root) ? 0 : 1;
    
    /* recursively calculate tips for children if there are any children */
    for (i = first_child_index; i < nd->nb_neigh; i++) {
        child = nd->neigh[i];
        if (isTip(child)) {
        } else {
            calculate_metric_color_pairs(child, s_tree, num_annotations, metric_color_pairs, num_pairs_per_color);
        }
    }

    for (i = 0; i < num_annotations; i++) {
        for (j = i+1; j < num_annotations; j++) {
            for(k = first_child_index; k < nd->nb_neigh; k++){
                for(l = first_child_index; l < nd->nb_neigh; l++){
                    if(k != l){
                        for(m = 0; m < nd->num_tips; m++){
                            for(n = 0; n < nd->num_tips; n++){
                                if((nd->color_tips[k][m] > -1) && (nd->color_tips[l][n] > -1)){
                                    if((nd->color_tips[k][m]==i) && (nd->color_tips[l][n]==j)) {
                                        metric_color_pairs[i][j] += nd->local_dist_tips[k][m] + nd->local_dist_tips[l][n];
                                        num_pairs_per_color[i][j] ++;
                                    }                                    
                                } 
                            }
                        }
                    }
                }
            }
        }
    }

}

void compute_metric(Tree *s_tree, size_t num_annotations, char** character) {
    /*
     * Calculates tree metric vector.
     * by reference to Kendall et al., bioRxiv, DOI: https://doi.org/10.1101/251710
     * but the principle is a bit modified to evaluate the PastML output "prediction"
     */
    double scaled_lk = 0;
    int i, j;
    double **metric_color_pairs; 
    int **num_pairs_per_color;
    FILE *output_list, *output_pairs;
    char *filename_list, *filename_pairs;

    int my_tips = node_check_tips(s_tree->root, s_tree, num_annotations);
    //printf("all tips counted at the root = %d\n",my_tips);
    
    metric_color_pairs = calloc(num_annotations, sizeof(double *));
    for (i = 0; i < num_annotations; i++) {
        metric_color_pairs[i] = calloc(num_annotations, sizeof(double));
    }
    num_pairs_per_color = calloc(num_annotations, sizeof(int *));
    for (i = 0; i < num_annotations; i++) {
        num_pairs_per_color[i] = calloc(num_annotations, sizeof(int));
    }

    for (i = 0; i < num_annotations; i++) {
        for (j = 0; j < num_annotations; j++) {
            metric_color_pairs[i][j] = 0.0;
            num_pairs_per_color[i][j] = 0;
        }
    }

    calculate_metric_color_pairs(s_tree->root, s_tree, num_annotations, metric_color_pairs, num_pairs_per_color);

    filename_list = calloc(256, sizeof(char));
    sprintf(filename_list, "list_pairs.csv");
    output_list = fopen(filename_list, "w");
    if (!output_list) {
        fprintf(stderr, "Output file %s is impossible to access.", filename_list);
        fprintf(stderr, "Value of errno: %d\n", errno);
        fprintf(stderr, "Error opening the file: %s\n", strerror(errno));
        return ENOENT;
    }

    for (i = 0; i < num_annotations; i++) {
        for (j = i+1; j < num_annotations; j++) {
                fprintf(output_list, "%s - %s,", character[i], character[j]);
        }
    }
    fprintf(output_list, "\n");
    fclose(output_list);

    filename_pairs = calloc(256, sizeof(char));
    sprintf(filename_pairs, "dist_pairs.csv");
    output_pairs = fopen(filename_pairs, "w");
    if (!output_pairs) {
        fprintf(stderr, "Output file %s is impossible to access.", filename_pairs);
        fprintf(stderr, "Value of errno: %d\n", errno);
        fprintf(stderr, "Error opening the file: %s\n", strerror(errno));
        return ENOENT;
    }
    for (i = 0; i < num_annotations; i++) {
        for (j = i+1; j < num_annotations; j++) {
            if(num_pairs_per_color[i][j] > 0){
                metric_color_pairs[i][j] /= num_pairs_per_color[i][j];
                fprintf(output_pairs, "%.5f,", metric_color_pairs[i][j]);
            }
        }
    }
    fprintf(output_pairs, "\n");
    fclose(output_pairs);

    for(i = 0; i < num_annotations; i++){
        free(metric_color_pairs[i]);
        free(num_pairs_per_color[i]);
    }
    free(metric_color_pairs);
    free(num_pairs_per_color);
    free(filename_list); free(filename_pairs);

}