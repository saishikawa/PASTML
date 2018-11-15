
#include <errno.h>
#include "pastml.h"

void put_true_scenario(Tree *s_tree, Node *nd, Node *root, size_t first_child_index, size_t num_annotations, char **character, char **ID, char **CHAR){
  int i, j;
  static int count = 0;

  if(nd->nb_neigh==1){
    return;
  }

  first_child_index = (nd == root) ? 0 : 1;

  if(nd != root)
    count++;

  //printf("Node %d,", count);
    for(j=0;j<num_annotations;j++){
        if(strcmp(character[j],CHAR[count])==0) {
            nd->result_probs[j] = 1.0;
        } else {
            nd->result_probs[j] = 0.0;
        }
        //printf("%s %.1f,",character[j],nd->result_probs[j]);
    }
  //printf("\n");

  for(i=first_child_index; i<nd->nb_neigh; i++){
    put_true_scenario(s_tree, nd->neigh[i], root, first_child_index, num_annotations, character, ID, CHAR);
  }

}
