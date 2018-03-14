#include "pastml.h"
#include "likelihood.h"

void pick_best_joint(Node *nd, Node *root, int ancestor, size_t first_child_index){

    /**
     * The top-down tree traversal to pick up the joint estimation of each node
     *
     */

  int i;
  static int count = 0;
  
  if(nd->nb_neigh==1){
    return;
  }

  first_child_index = (nd == root) ? 0 : 1;

  if(nd == root){
    nd->best_joint_state = ancestor;
  } else {
    nd->best_joint_state = nd->joint_state[ancestor];
    ancestor=nd->best_joint_state;
  }
  
  for(i=first_child_index; i<nd->nb_neigh; i++){
    pick_best_joint(nd->neigh[i], root, ancestor, first_child_index);
  }
  
  if(nd != root){
    count++;
    sprintf(nd->sim_name, "Node%d", count);
  } else {
    count=0;
  }
  return;
}

void calculate_node_joint_probabilities(Node *nd, Node *root, size_t num_annotations, double *frequency, size_t first_child_index){
    /**
     * The joint-likelihood of a given node is computed based on the information
     * coming from all the tips descending from the studied node
     * using a dynamic programming proposed by Pupko et al 2000.
     *
     */
  size_t i,ii,j,best_root_state;
  double tmp_prob[num_annotations], curr_scaler, best_joint_lik, log_lik, smallest;
  static int factors=0;
  int curr_scaler_pow, piecewise_scaler_pow;

  //if tips
  if(nd->nb_neigh==1){
     for(i=0;i<num_annotations;i++){
       tmp_prob[i]=nd->joint_likelihood[i];
     }
     for(i=0;i<num_annotations;i++){
       nd->joint_likelihood[i]=0.0;
       /*assume state i at the ancestral node and state j at this tip*/
       for(j=0;j<num_annotations;j++){
         nd->joint_likelihood[i]+=nd->pij[i][j]*tmp_prob[j];
       }
     }
     return;
  }

  first_child_index = (nd == root) ? 0 : 1;
  
  for(i=first_child_index;i<nd->nb_neigh;i++){
    calculate_node_joint_probabilities(nd->neigh[i], root, num_annotations, frequency, first_child_index);
  }

  if(nd == root){ /* ROOT */
    for(i=0;i<num_annotations;i++){
      /*collect joint likelihoods from all descendant nodes assuming state i at the root*/
      for(ii=first_child_index;ii<nd->nb_neigh;ii++){      
        if(ii==first_child_index) {
          tmp_prob[i] = nd->neigh[ii]->joint_likelihood[i];
        } else {
          tmp_prob[i] *= nd->neigh[ii]->joint_likelihood[i];
        }
      }
      nd->joint_likelihood[i] = tmp_prob[i] * frequency[i];
    }

    best_joint_lik=0.0;
    for(i=0;i<num_annotations;i++){
      if(best_joint_lik < tmp_prob[i]){
        best_joint_lik = tmp_prob[i];
        best_root_state = i;
      }
    }

    /*rescale likelihood*/
    log_lik=log(best_joint_lik);
    do {
      piecewise_scaler_pow = MIN(factors,63);
      curr_scaler = ((unsigned long long)(1) << piecewise_scaler_pow);
      log_lik -= LOG2*piecewise_scaler_pow;
      factors -= piecewise_scaler_pow;
    } while(factors != 0);
    printf("Joint Likelihood = %.5f\n",log_lik);
    pick_best_joint(root, root, best_root_state, 0);
    factors=0;
    return;

  } else { /*internal nodes*/

    for(i=0;i<num_annotations;i++){
       /*assume state i at the ancestral node*/
       nd->joint_likelihood[i]=0.;
       for(j=0;j<num_annotations;j++){
         /*collect joint likelihoods from all descendant nodes assuming state j at this node*/ 
         for(ii=first_child_index;ii<nd->nb_neigh;ii++){      
           if(ii==first_child_index) {
             tmp_prob[j] = nd->neigh[ii]->joint_likelihood[j];
           } else {
             tmp_prob[j] *= nd->neigh[ii]->joint_likelihood[j];
           }
         }
         tmp_prob[j] *= nd->pij[i][j];
         /*find state j at this node giving the largest joint probability when its ancestor represents state i*/
         if(nd->joint_likelihood[i] < tmp_prob[j]) {
           nd->joint_likelihood[i] = tmp_prob[j];
           nd->joint_state[i] = j;
         }
       }
    }
    
    /*likelihood scaling*/
    smallest = 1.0;
    for(i=1;i<num_annotations;i++){
       if(nd->joint_likelihood[i] > 0.0 && nd->joint_likelihood[i] < smallest){
         smallest=nd->joint_likelihood[i];
       }
    }
    if(smallest < LIM_P){
         curr_scaler_pow = (int)(POW*LOG2-log(smallest))/LOG2;
         curr_scaler     = ((unsigned long long)(1) << curr_scaler_pow);
         factors+=curr_scaler_pow;
         do {
           piecewise_scaler_pow = MIN(curr_scaler_pow,63);
           curr_scaler = ((unsigned long long)(1) << piecewise_scaler_pow);
           for(i=0;i<num_annotations;i++){
             nd->joint_likelihood[i] *= curr_scaler;
           }
           curr_scaler_pow -= piecewise_scaler_pow;
         } while(curr_scaler_pow != 0);
    }
  } 
  return;
}

void calculate_joint_probabilities(Tree *s_tree, size_t num_annotations, double *frequency) {
    /**
     * Calculates joint probabilities of tree nodes.
     */
  calculate_node_joint_probabilities(s_tree->root, s_tree->root, num_annotations, frequency, 0);
}

