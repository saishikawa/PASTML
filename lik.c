#include "asrml.h"
 extern int have_miss;
extern Tree *s_tree;
extern Node *root;

double global_like;
double global_factor;

void calc_lik_bfgs(Node *nd, char** tipnames, int* states, int nb, int nbanno, double mu, char* model, double* p, double* likelihood){
  int i,j,ii;
  double mul, expmul, sum=0., prob_left=0., prob_right=0., bl, smallest, scaled_lk, logroot, prob_sons[MAXPOLY];
  static int factors=0;
  double curr_scaler;
  int curr_scaler_pow, piecewise_scaler_pow, node_start;
  double ts=8.0, tv=1.0, beta, freqTC, freqAG, mul2, rig, lef;

  for(i=0;i<nbanno;i++) {
    nd->prob[i]=0.;
  }
  sum=0.;
  if(nd->nneigh==1){ /*tips*/
     //printf("Calculating at %s:%lf\n",nd->name,nd->br[0]->brlen);
     bl=nd->br[0]->brlen;
     if(bl==0.0){
       //printf("zero branches\n");
       bl = (nd->br[0]->brlen+p[nbanno+1]) * (s_tree->avgbl / (s_tree->avgbl + p[nbanno+1]));
     } else {
       bl=nd->br[0]->brlen*p[nbanno];
       //printf("bl = %lf\n",bl);
     }
     mul=-1.*mu*bl;
     expmul=exp(mul);
     /*tip probability*/
     for(i=0;i<nb;i++){
       if(strcmp(nd->name, tipnames[i])==0){
         nd->pupko_state=states[i];
         if(states[i]==have_miss){
           //printf("have a missing data at tips, give equal probabilities\n");
           for(j=0;j<nbanno;j++){
             nd->condlike[j]=1.0;
             nd->up_like[j]=1.0;
           }
           //printf("probability = %lf\n",nd->condlike[0]);
         } else {
           nd->condlike[states[i]]=1.0;
           nd->up_like[states[i]]=1.0;
         }
         break;
       }
     }
     /*Pij*/
     for(i=0;i<nbanno;i++){
       for(j=0;j<nbanno;j++){
         if(i==j){
           nd->pij[i][j]=expmul + (( 1.0 - expmul ) * p[i] );
         } else {
           nd->pij[i][j]=p[j] * (1.0 - expmul);
         }
       }
     }
     for(i=0;i<nbanno;i++){
         //printf("condlike%d = %lf, ", i, nd->condlike[i]);
     }
     //printf("\n");
     return;
  }
  if(nd==root) {
    node_start=0;
  } else {
    node_start=1;
  }
  
  for(i=node_start;i<nd->nneigh;i++){
    calc_lik_bfgs(nd->neigh[i], tipnames, states, nb, nbanno, mu, model, p, likelihood);
  }

  if(nd==root){
    //printf("Calculating at %s\n",nd->name);
    for(i=0;i<nbanno;i++){
      for(ii=node_start;ii<nd->nneigh;ii++){      
        prob_sons[ii]=0.;
        for(j=0;j<nbanno;j++){
          prob_sons[ii] += nd->neigh[ii]->pij[i][j] * nd->neigh[ii]->condlike[j];
        }
        //printf("prob_sons%d = %lf, ",ii,prob_sons[ii]);
        if(ii==node_start) {
          nd->condlike[i] = prob_sons[ii];
        } else {
          nd->condlike[i] = nd->condlike[i] * prob_sons[ii];
        }
        //printf("condlike%d = %lf\n", i, nd->condlike[i]);
      }
      nd->condlike[i] = nd->condlike[i] * p[i];
      if(nd->calc_flag[i]==0) nd->condlike[i]=0.;
    }
    for(i=0;i<nbanno;i++){
         //printf("condlike%d = %.6e, ", i, nd->condlike[i]);
    }
    //printf("\n");
    smallest=1.0;
    for(j=0;j<nbanno;j++){
      if(nd->condlike[j] > 0.0) {
        if(nd->condlike[j] < smallest) smallest=nd->condlike[j];
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
             nd->condlike[j] *= curr_scaler;
           }
           curr_scaler_pow -= piecewise_scaler_pow;
      } while(curr_scaler_pow != 0);
    }
    scaled_lk=0.;
    for(i=0;i<nbanno;i++){
      sum += nd->condlike[i]*p[i];
      scaled_lk += nd->condlike[i];
    }
    for(i=0;i<nbanno;i++){
      nd->mar_prob[i] = nd->condlike[i]*p[i]/sum;
    }
    global_like=scaled_lk;
    global_factor=factors;
    nd->up_factor=factors;
    for(i=0;i<nbanno;i++){
      nd->up_like[i]=nd->condlike[i];
    }
    logroot=log(scaled_lk);
    //if(factors>0) printf("likelihood scaling is on\n");
    do {
      piecewise_scaler_pow = MIN(factors,63);
      curr_scaler = ((unsigned long long)(1) << piecewise_scaler_pow);
      logroot = logroot - LOG2*piecewise_scaler_pow;
      factors -= piecewise_scaler_pow;
    } while(factors != 0);
    factors=0;
    //printf("scaled lk = %.12e, log = %lf\n",scaled_lk,log(scaled_lk));
    *likelihood = logroot;
    //return logroot;
    return;
  } else {
    //printf("Calculating at %s:%lf\n",nd->name,nd->br[0]->brlen);
    bl=nd->br[0]->brlen;
    if(bl==0.0){
       //printf("zero branches\n");
       bl = (nd->br[0]->brlen+p[nbanno+1]) * (s_tree->avgbl / (s_tree->avgbl + p[nbanno+1]));
    } else {
       bl=nd->br[0]->brlen*p[nbanno];
       //printf("bl = %lf\n",bl);
    }
    mul=-1.*mu*bl;
    expmul=exp(mul);
    /*Pij*/
    for(i=0;i<nbanno;i++){
      for(j=0;j<nbanno;j++){
        if(i==j){
          nd->pij[i][j]=expmul + (( 1.0 - expmul ) * p[i] );
        } else {
          nd->pij[i][j]=p[j] * (1.0 - expmul );
        }
      }
    }

    for(i=0;i<nbanno;i++){
      for(ii=node_start;ii<nd->nneigh;ii++){      
        prob_sons[ii]=0.;
        for(j=0;j<nbanno;j++){
          prob_sons[ii] += nd->neigh[ii]->pij[i][j] * nd->neigh[ii]->condlike[j];
        }
        //printf("prob_sons%d = %lf, ",ii,prob_sons[ii]);
        if(ii==node_start) {
          nd->condlike[i] = prob_sons[ii];
        } else {
          nd->condlike[i] = nd->condlike[i] * prob_sons[ii];
        }
        //printf("condlike%d = %lf\n", i, nd->condlike[i]);
      }
      if(nd->calc_flag[i]==0) nd->condlike[i]=0.;
    }
    /*for(i=0;i<nbanno;i++){
        printf("%s condlike%d = %.5e, ", nd->name, i, nd->condlike[i]);
    }
    printf("\n");*/
    /*scaling*/
    smallest=1.0;
    for(j=0;j<nbanno;j++){
      if(nd->condlike[j] > 0.0) {
        if(nd->condlike[j] < smallest) smallest=nd->condlike[j];
      } else {
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
             nd->condlike[j] *= curr_scaler;
           }
           curr_scaler_pow -= piecewise_scaler_pow;
      } while(curr_scaler_pow != 0);
    }
    nd->up_factor=factors;
    for(i=0;i<nbanno;i++){
      nd->up_like[i]=nd->condlike[i];
      //printf("%s up_like%d = %.5e, ", nd->name, i, nd->up_like[i]);
    }
    //printf("\n");
  } 
  return;
}
