#include "asrml.h"

extern Tree* s_tree;
extern Node* root;


void down_like_marginal(Node* nd, int nb, int nbanno, double mu, double scale, double* frequency){
  int i,j,k, ii;
  double mul, expmul, sum=0., prob_up=0., prob_down=0., prob_down_son=0., smallest, scaled_lk=0., logroot, root_my_BL=0., bl;
  Node* father;
  int factor=0;
  double curr_scaler;
  int curr_scaler_pow, piecewise_scaler_pow, my_id, node_start;

  //printf("reached at %s, ", nd->name);
  if(nd->nneigh==1){
     return;
  }
  
  if(nd==root){
  } else {
    father = nd->neigh[0];
    //printf("%s, father=%s ", nd->name, father->name);

    for(i=0;i<father->nneigh;i++){
        if(father->neigh[i] == nd) {
          my_id=i;
          //printf("my_id = %d,",my_id);
        }
    }
    for(j=0;j<nbanno;j++){
        nd->sum_down[j]=0.;
        //assume this node is j
        for(k=0;k<nbanno;k++){
          //assume father is k
          prob_down=1.0;
          if(father==root) {
            node_start=0;
          } else {
            node_start=1;
          }
          nd->down_factor=0;
          for(i=node_start;i<father->nneigh;i++){
            if(i!=my_id){
              nd->down_factor+=father->neigh[i]->up_factor;
              prob_down_son=0.;
              for(ii=0;ii<nbanno;ii++){
                // assume another son is ii
                //printf("%s, pij = %lf, up = %.5e\n", father->neigh[i]->name,father->neigh[i]->pij[k][ii],father->neigh[i]->up_like[ii]);
                prob_down_son+=father->neigh[i]->pij[k][ii] * father->neigh[i]->up_like[ii];  
              }
              prob_down = prob_down * prob_down_son;
              //printf("prob_down = %.5e, ", prob_down);
            }
          }
          if(father != root) {
            prob_down = prob_down * father->sum_down[k];
            nd->down_factor+=father->down_factor;
          }
          nd->sum_down[j] += nd->pij[k][j] * prob_down;
          if(nd->calc_flag[j]==0) nd->sum_down[j]=0.;
          //printf("sum_down %d = %.5e, ", j, nd->sum_down[j]);
        }
      }

      smallest=1.0;
      for(j=0;j<nbanno;j++){
        if(nd->sum_down[j] > 0.0) {
          if(nd->sum_down[j] < smallest) smallest=nd->sum_down[j];
        }
      }
      if(smallest < LIM_P){
        curr_scaler_pow = (int)(POW*LOG2-log(smallest))/LOG2;
        curr_scaler     = ((unsigned long long)(1) << curr_scaler_pow);
        nd->down_factor+=curr_scaler_pow;
        do {
           piecewise_scaler_pow = MIN(curr_scaler_pow,63);
           curr_scaler = ((unsigned long long)(1) << piecewise_scaler_pow);
           for(j=0;j<nbanno;j++){
             nd->sum_down[j] *= curr_scaler;
           }
           curr_scaler_pow -= piecewise_scaler_pow;
        } while(curr_scaler_pow != 0);
      }

      for(j=0;j<nbanno;j++){
        nd->condlike_mar[j] = frequency[j] * nd->up_like[j] * nd->sum_down[j];
        sum+=nd->condlike_mar[j] * frequency[j];
      }
      for(j=0;j<nbanno;j++){
        nd->condlike_mar[j] = nd->condlike_mar[j] * frequency[j] / sum;
        //printf("marginal %d=%lf, ", j, nd->condlike_mar[j]);
      }
      //printf("\n");

  }

  if(nd==root) {
    node_start=0;
  } else {
    node_start=1;
  }
  for(i=node_start;i<nd->nneigh;i++){
    down_like_marginal(nd->neigh[i], nb, nbanno, mu, scale, frequency);
  }

  return;
}

