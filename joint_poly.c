#include "asrml.h"

extern Tree* s_tree;
extern Node* root;
extern have_miss;

void pick_scenario(Node* nd, int nb, int anc){
  static int count=0;
  int i, node_start;
  
  if(nd->nneigh==1){
    for(i=0;i<s_tree->nb_nodes;i++){
      if(strcmp(s_tree->a_nodes[i]->name,nd->name)==0){
        nd->pupko_state=nd->best_char[anc];
      }
    }
    
    count++;
    return;
  }
  if(nd==root){
    for(i=0;i<s_tree->nb_nodes;i++){
      if(strcmp(s_tree->a_nodes[i]->name,nd->name)==0){
        nd->pupko_flag=anc;
        count++;
      }
    }
  } else {
    for(i=0;i<s_tree->nb_nodes;i++){
      if(strcmp(s_tree->a_nodes[i]->name,nd->name)==0){
        nd->pupko_flag=nd->best_char[anc];
        count++;
      }
    }
    anc=nd->best_char[anc];
  }

  if(nd==root) {
    node_start=0;
  } else {
    node_start=1;
  }
  
  for(i=node_start;i<nd->nneigh;i++){
    pick_scenario(nd->neigh[i],nb,anc);
  }

  //printf("%s,%d\n",nd->nom,nd->pupko_flag);
  if(nd==root) count=0;
  return;
}

void joint(Node* nd, char** tipnames, int* states, int nb, int nbanno, double* frequency){
  int i,ii,j,k,root_char,node_start;
  double prob_i[nbanno], root_best, smallest, scaled_lk, logroot, prob_sons[MAXPOLY];
  static int factors=0;
  double curr_scaler;
  int curr_scaler_pow, piecewise_scaler_pow;

  for(i=0;i<nbanno;i++) prob_i[i]=0.;

  if(nd->nneigh==1){
     /*tip probability*/
     for(i=0;i<nb;i++){
       if(strcmp(nd->name, tipnames[i])==0){
         if(states[i]==have_miss){
           //printf("have a missing data at tips, give equal probabilities\n");
           for(j=0;j<nbanno;j++){
             prob_i[j]=(double)1.0/(double)nbanno;
           }
         } else {
           prob_i[states[i]]=1.0;
         }
         for(j=0;j<nbanno;j++) nd->best_char[j]=states[i];
         break;
       }
     }
     for(i=0;i<nbanno;i++){
       /*assume i at the ancestral node*/
       nd->prob[i]=0.;
       for(j=0;j<nbanno;j++){
         nd->prob[i]+=nd->pij[i][j]*prob_i[j];
       }
       nd->best_i[i]=nd->prob[i];
     }
     return;
  }

  if(nd==root) {
    node_start=0;
  } else {
    node_start=1;
  }
  
  for(i=node_start;i<nd->nneigh;i++){
    joint(nd->neigh[i], tipnames, states, nb, nbanno, frequency);
  }

  if(nd==root){ /* ROOT */
    root_best=0.;
    root_char=256;
    for(i=0;i<nbanno;i++){
      for(ii=node_start;ii<nd->nneigh;ii++){      
        if(ii==node_start) {
          nd->prob[i] = nd->neigh[ii]->best_i[i];
        } else {
          nd->prob[i] = nd->prob[i] * nd->neigh[ii]->best_i[i];
        }
        //printf("condlike%d = %lf\n", i, nd->condlike[i]);
      }
      nd->prob[i] = nd->prob[i] * frequency[i];
    }

    root_best=DBL_MIN;
    for(i=0;i<nbanno;i++){
      if(root_best < nd->prob[i]){
        root_best=nd->prob[i];
        root_char=i;
      }
    }
    scaled_lk=root_best;
    logroot=log(scaled_lk);
    s_tree->pupko_like=scaled_lk;
    s_tree->pupko_factor=factors;
    //printf("sum scaling = %d, scaled lk = %.12f\n",factors,logroot);
    do {
      piecewise_scaler_pow = MIN(factors,63);
      curr_scaler = ((unsigned long long)(1) << piecewise_scaler_pow);
      logroot = logroot - LOG2*piecewise_scaler_pow;
      //printf("scale down by %d to %.12f\n", piecewise_scaler_pow, logroot);
      factors -= piecewise_scaler_pow;
    } while(factors != 0);
    //printf("%s, Best state = %d, L = %.12f\n", nd->nom, root_char, logroot);
    s_tree->pupko_value=logroot;
    pick_scenario(root,nb,root_char);
    factors=0;
    return;

  } else { /*internal nodes*/

    for(i=0;i<nbanno;i++){
       /*assume i at the ancestral node*/
       nd->best_i[i]=DBL_MIN;
       for(j=0;j<nbanno;j++){
         for(ii=node_start;ii<nd->nneigh;ii++){      
           if(ii==node_start) {
             nd->prob[j] = nd->neigh[ii]->best_i[j];
           } else {
             nd->prob[j] = nd->prob[j] * nd->neigh[ii]->best_i[j];
           }
           //printf("condlike%d = %lf\n", i, nd->condlike[i]);
         }
         nd->prob[j]=nd->pij[i][j] * nd->prob[j];
       }

       for(k=0;k<nbanno;k++){
         if(nd->best_i[i] < nd->prob[k]){
           nd->best_i[i]=nd->prob[k];
           nd->best_char[i]=k;
         }
       }
    }
    
       /*scaling factor*/
    smallest=DBL_MAX;
    for(j=1;j<nbanno;j++){
       if(nd->best_i[j] > DBL_MIN){
           if(nd->best_i[j] < smallest) smallest=nd->best_i[j];
       }
    }
    if(smallest < LIM_P){
         curr_scaler_pow = (int)(POW*LOG2-log(smallest))/LOG2;
         curr_scaler     = ((unsigned long long)(1) << curr_scaler_pow);
         factors+=curr_scaler_pow;
         do {
           piecewise_scaler_pow = MIN(curr_scaler_pow,63);
           curr_scaler = ((unsigned long long)(1) << piecewise_scaler_pow);
           for(j=0;j<nbanno;j++){
             nd->best_i[j] *= curr_scaler;
           }
           curr_scaler_pow -= piecewise_scaler_pow;
         } while(curr_scaler_pow != 0);
    }
    //printf("%s", nd->nom);
    //for(i=0;i<nbanno;i++) printf(", L%d=%.12f, C%d=%d", i, log(nd->best_i[i]), i, nd->best_char[i]);
    //printf("\n");
  } 
  return;
}

