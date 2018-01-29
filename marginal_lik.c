#include "asrml.h"

extern Tree *s_tree;
extern Node *root;


void down_like_marginal(Node* nd, int nb, int nbanno, double mu, double scale, double* frequency){
  int i,j,k, ii;
  double mul, expmul, sum=0., prob_up=0., prob_down=0., prob_down_son=0., smallest, scaled_lk=0., logroot, root_my_BL=0., bl, tmp_br, sum_mu;
  Node* father;
  int factor=0, tmp_factor;
  double curr_scaler, sumlog;
  int curr_scaler_pow, piecewise_scaler_pow, my_id, node_start;

  sum_mu=0.0;
  for(i=0;i<nbanno;i++){
      sum_mu+=frequency[i]*frequency[i];
  }
  mu=1/(1-sum_mu);
  if(nd->nneigh==1){
     return;
  }
  
  if(nd==root){
    for(j=0;j<nbanno;j++){
      nd->sum_down[j]=frequency[j];
      nd->down_factor=0;
    }
  } else {
    father = nd->neigh[0];
    for(i=0;i<father->nneigh;i++){
        if(father->neigh[i] == nd) {
          my_id=i;
        }
    }

    //currentNODE has j
     for(j=0;j<nbanno;j++){
        if(father==root) {
          node_start=0;
          nd->sum_down[j] =1.0;
          prob_down = 1.0;
          for(i=node_start;i<father->nneigh;i++){
            if(i!=my_id){
              prob_down_son = 0.0;
              tmp_br = (nd->br[0]->brlen + father->neigh[i]->br[0]->brlen) * scale;
              mul=-1.*mu*tmp_br;
              expmul=exp(mul);
              //assume other sons have ii
              for(ii=0;ii<nbanno;ii++){
                 //assume root has k
                nd->rootpij[j][ii]=0.0;
                for(k=0;k<nbanno;k++){
                    nd->rootpij[j][ii]+=nd->pij[j][k] * father->neigh[i]->pij[k][ii];
                }
                prob_down_son += nd->rootpij[j][ii] * father->neigh[i]->up_like[ii];
              }
              nd->sum_down[j] *= prob_down_son;
	    }
          }
        } else {
          node_start=1;
          nd->sum_down[j] = 0.0;
          //assume father has k
          for(k=0;k<nbanno;k++){
            prob_down = nd->pij[j][k] * father->sum_down[k];
            for(i=node_start;i<father->nneigh;i++){
              if(i!=my_id){
                prob_down_son = 0.0;
                //assume other sons have ii
                for(ii=0;ii<nbanno;ii++){
                  prob_down_son += father->neigh[i]->pij[k][ii] * father->neigh[i]->up_like[ii];
                }
                prob_down *= prob_down_son;
	      }
            }
            nd->sum_down[j] += prob_down;
          }
        }
      }
      
      nd->down_factor=nd->up_factor;
      for(j=0;j<nbanno;j++){
        nd->condlike_mar[j] = nd->sum_down[j] * nd->up_like[j] * frequency[j];
      }

      smallest=1.0;
      for(j=0;j<nbanno;j++){
        if(nd->condlike_mar[j] > 0.0) {
          if(nd->condlike_mar[j] < smallest) smallest=nd->condlike_mar[j];
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
             nd->condlike_mar[j] *= curr_scaler;
           }
           curr_scaler_pow -= piecewise_scaler_pow;
        } while(curr_scaler_pow != 0);
      }
 
      sum=0.0;
      for(j=0;j<nbanno;j++){
        sum+=nd->condlike_mar[j];
      }
      //printf("%s, SUM=%lf, factor=%d\n",nd->name,log(sum),nd->down_factor);
      for(j=0;j<nbanno;j++){
        nd->condlike_mar[j] = (nd->condlike_mar[j]) / sum;
      }

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

