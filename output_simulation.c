
#include <errno.h>
#include "pastml.h"

void output_sim_node_states(Node *nd, Node *root, size_t num_annotations, char **character, FILE *outfile, size_t method_num, size_t first_child_index){
  size_t i, k, tmp_map;
  double tmp_prob;
  static int count = 0;

  if(nd->nb_neigh==1){
    return;
  }

  first_child_index = (nd == root) ? 0 : 1;
  
  for(i=first_child_index; i<nd->nb_neigh; i++){
    output_sim_node_states(nd->neigh[i], root, num_annotations, character, outfile, method_num, first_child_index);
  }

  if(nd != root){
    count++;
    sprintf(nd->sim_name, "Node%d", count);
  } else {
    count=0;
  }

  fprintf(outfile, "%s", nd->sim_name);
  if(nd->neigh[0] == root) fprintf(outfile, "_child", nd->sim_name);

  if(method_num == 0){
    //OUTPUT_JOINT
    for (i = 0; i < num_annotations; i++) {
      if(nd->result_probs[i] == 1.0) fprintf(outfile, ",%s\n", character[i]);
    }
  } else if (method_num == 1) {
    //OUTPUT_MARGINAL
    for (i = 0; i < num_annotations; i++) {
      fprintf(outfile, ",%s=%.10f", character[i], nd->result_probs[i]);
    }
    fprintf(outfile, "\n");
  } else if (method_num == 2) {
    //OUTPUT_MAP
    for (i = 0; i < num_annotations; i++) {
      if(nd->result_probs[i] == 1.0) fprintf(outfile, ",%s\n", character[i]);
    }
  } else if (method_num == 3) {
    //OUTPUT_MARGINAL_APPROXIMATION
    for (i = 0; i < num_annotations; i++) {
      if(nd->result_probs[i] > 0.0) fprintf(outfile, ",%s", character[i]);
    }
    fprintf(outfile, "\n"); 
  } else if (method_num == 4) {
    //OUTPUT_PARSIMONY
    for (i = 0; i < num_annotations; i++) {
      if(nd->result_probs[i] > 0.0) fprintf(outfile, ",%s", character[i]);
    }
    fprintf(outfile, "\n"); 
  }

  return;

}



int output_simulation(Tree *tree, size_t num_annotations, char **character, char *output_file_path, size_t method_num) {
    FILE* outfile = fopen(output_file_path, "w");
    if (!outfile) {
        fprintf(stderr, "Output file %s is impossible to access.", output_file_path);
        fprintf(stderr, "Value of errno: %d\n", errno);
        fprintf(stderr, "Error opening the file: %s\n", strerror(errno));
        return ENOENT;
    }
    output_sim_node_states(tree->root, tree->root, num_annotations, character, outfile, method_num, 0);

    fclose(outfile);
    return EXIT_SUCCESS;
}

