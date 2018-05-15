
#include <errno.h>
#include "pastml.h"

void name_simulation(Node *nd, Node *root, size_t first_child_index){
  size_t i;
  static int count = 0;

  if(nd->nb_neigh==1){
    return;
  }

  first_child_index = (nd == root) ? 0 : 1;
  
  for(i=first_child_index; i<nd->nb_neigh; i++){
    name_simulation(nd->neigh[i], root, first_child_index);
  }

  if(nd != root){
    count++;
    sprintf(nd->sim_name, "Node%d", count);
  } else {
    count=0;
  }

}

void output_sim_node_states(Node *nd, Node *root, size_t num_annotations, char **character, FILE *outfile, size_t method_num, size_t first_child_index, char **ID, char **CHAR, int num_nodes){
  size_t i, j, k, tmp_map;
  double tmp_prob;
  int current, parent;
  static double node_brier = 0.0;
  static double edge_brier = 0.0;
  static double num_pred = 0.0;

  if(nd->nb_neigh==1){
    return;
  }

  first_child_index = (nd == root) ? 0 : 1;
  
  for(i=first_child_index; i<nd->nb_neigh; i++){
    output_sim_node_states(nd->neigh[i], root, num_annotations, character, outfile, method_num, first_child_index, ID, CHAR, num_nodes);
  }

  //Node identification
  for(i=0;i<num_nodes;i++){
      if(strcmp(nd->sim_name,ID[i])==0){
        current = i;
      }
      if(nd != root) {
        if(strcmp(nd->neigh[0]->sim_name,ID[i])==0){
          parent = i;
        }
      }
  }
  //printf("Current = %s,%s,%s , Parent = %s,%s,%s\n", nd->sim_name, ID[current], CHAR[current], nd->neigh[0]->sim_name, ID[parent], CHAR[parent]);

  //Calculation of the node BS
  for(i=0;i<num_annotations;i++){
      if(strcmp(character[i],CHAR[current])==0){
        node_brier += (nd->result_probs[i] - 1.0)*(nd->result_probs[i] - 1.0);
      } else {
        node_brier += (nd->result_probs[i] - 0.0)*(nd->result_probs[i] - 0.0);
      }
  }
  //Calculation of the edge BS
  if(nd != root) {
    for(i=0;i<num_annotations;i++){
      for(j=0;j<num_annotations;j++){
        if(strcmp(character[i],CHAR[current])==0 && strcmp(character[j],CHAR[parent])==0){
          edge_brier += (nd->result_probs[i]*nd->neigh[0]->result_probs[j] - 1.0)*(nd->result_probs[i]*nd->neigh[0]->result_probs[j] - 1.0);
        } else {
          edge_brier += (nd->result_probs[i]*nd->neigh[0]->result_probs[j] - 0.0)*(nd->result_probs[i]*nd->neigh[0]->result_probs[j] - 0.0);
        }
      }
    }
  }
  for(i=0;i<num_annotations;i++){    
      if(nd->result_probs[i] > 0) num_pred += 1.0;
  }

  if(nd == root){
    node_brier /= (double)num_nodes;
    fprintf(outfile, "%.5lf\n", node_brier);
    node_brier = 0.0;
    edge_brier /= (double)(num_nodes - 1);
    fprintf(outfile, "%.5lf\n", edge_brier);
    edge_brier = 0.0;
    num_pred /= (double)num_nodes;
    fprintf(outfile, "%.5lf\n", num_pred);
  }

  return;

}



int output_simulation(Tree *tree, size_t num_annotations, char **character, char *output_file_path, size_t method_num, char **ID, char **CHAR, int num_nodes) {
    FILE* outfile = fopen(output_file_path, "w");
    if (!outfile) {
        fprintf(stderr, "Output file %s is impossible to access.", output_file_path);
        fprintf(stderr, "Value of errno: %d\n", errno);
        fprintf(stderr, "Error opening the file: %s\n", strerror(errno));
        return ENOENT;
    }
    output_sim_node_states(tree->root, tree->root, num_annotations, character, outfile, method_num, 0, ID, CHAR, num_nodes);

    fclose(outfile);
    return EXIT_SUCCESS;
}

