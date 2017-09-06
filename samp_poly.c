#include "asrml.h"

extern Tree* s_tree;
extern Node* root;

void output_samp(Node* nd, int nb, int nbanno, char **character, FILE* outfile){
  int i, node_start;
  static int count=0;

  if(nd->nneigh==1){
     count++;
     return;
  }
  count++;

  if(nd==root) {
    node_start=0;
  } else {
    node_start=1;
  }
  
  for(i=node_start;i<nd->nneigh;i++){
    output_samp(nd->neigh[i],nb,nbanno,character,outfile);
  }

  //fprintf(outfile,"%s",nd->name);
  fprintf(outfile,"%s\n",character[nd->prob_sampled]);

  return;
}

int prob_sampling(Node* nd, char** tipnames, int* states, int nb, int nbanno, int ancstate, double* frequency){
  int i,j,tempanc,tmpfactor,node_start;
  double sumlik=0., random, border, logroot;
  static double sampval=1.;
  static int factors=0;
  double curr_scaler;
  int curr_scaler_pow, piecewise_scaler_pow;
  static int count=0;
  
  if(nd->nneigh==1){
    /*sampling at the tips*/
    sampval=sampval*nd->pij[ancstate][nd->best_char[ancstate]];
    if(sampval < LIM_P){
         curr_scaler_pow = (int)(POW*LOG2-log(sampval))/LOG2;
         curr_scaler     = ((unsigned long long)(1) << curr_scaler_pow);
         factors+=curr_scaler_pow;
         do {
           piecewise_scaler_pow = MIN(curr_scaler_pow,63);
           curr_scaler = ((unsigned long long)(1) << piecewise_scaler_pow);
           sampval *= curr_scaler;
           curr_scaler_pow -= piecewise_scaler_pow;
         } while(curr_scaler_pow != 0);
    }
    count++;
    return 0;
  }
  if(nd==root){
      for(i=0;i<nbanno;i++){
        sumlik+=nd->condlike[i];
      }
      for(i=0;i<nbanno;i++){
        nd->sortedlike[i]=nd->condlike[i]/sumlik;
        nd->sortedstates[i]=i;
      }
      random=(double)rand()/RAND_MAX;
      border=0.;
      for(i=0;i<nbanno;i++){
        border+=nd->sortedlike[i];
        if(random <= border){
          ancstate=nd->sortedstates[i];
          nd->prob_sampled=nd->sortedstates[i];
          count++;
          sampval=sampval*frequency[ancstate];
          break;
        }
      }
      if(ancstate==-1){
        printf("error in sampling \n");
        exit(0);
      }
  } else {
      for(j=0;j<nbanno;j++){
        sumlik+=nd->condlike[j]*nd->pij[ancstate][j];
      }
      for(j=0;j<nbanno;j++){
        nd->sortedlike[j]=nd->condlike[j]*nd->pij[ancstate][j]/sumlik;
        nd->sortedstates[j]=j;
      }
      random=(double)rand()/RAND_MAX;
      border=0.;
      for(i=0;i<nbanno;i++){
        border+=nd->sortedlike[i];
        if(random <= border){
          sampval=sampval*nd->pij[ancstate][nd->sortedstates[i]];
          ancstate=nd->sortedstates[i];
          nd->prob_sampled=nd->sortedstates[i];
	      count++;
          break;
        }
      }
      if(sampval < LIM_P){
        curr_scaler_pow = (int)(POW*LOG2-log(sampval))/LOG2;
        curr_scaler     = ((unsigned long long)(1) << curr_scaler_pow);
        factors+=curr_scaler_pow;
        do {
          piecewise_scaler_pow = MIN(curr_scaler_pow,63);
          curr_scaler = ((unsigned long long)(1) << piecewise_scaler_pow);
          sampval *= curr_scaler;
          curr_scaler_pow -= piecewise_scaler_pow;
        } while(curr_scaler_pow != 0);
      }
  }
  tempanc=ancstate;

  if(nd==root) {
    node_start=0;
  } else {
    node_start=1;
  }
  
  for(i=node_start;i<nd->nneigh;i++){
    prob_sampling(nd->neigh[i], tipnames, states, nb, nbanno, tempanc, frequency);
  }

  tmpfactor=factors;
  if(nd==root){
    s_tree->sample_like=sampval;
    logroot=log(sampval);
    do {
      piecewise_scaler_pow = MIN(factors,63);
      curr_scaler = ((unsigned long long)(1) << piecewise_scaler_pow);
      logroot = logroot - LOG2*piecewise_scaler_pow;
      factors -= piecewise_scaler_pow;
    } while(factors != 0);
    s_tree->sample_value=logroot;
    sampval=1.;
    logroot=0.;
    count=0;
  }
  return tmpfactor;
}


void make_samp_prob (char** tipnames, int* states, int num_tips, int num_anno, char **character, int check_out, double* frequency){
  int tmpfactor;
  FILE *fpn;
  char tmpsamp[100];

    tmpfactor=prob_sampling(root, tipnames, states, num_tips, num_anno, -1, frequency);

    s_tree->scenario_values=s_tree->sample_value;
    s_tree->scenario_like=s_tree->sample_like;
    s_tree->scenario_factors=tmpfactor;
    
    if(check_out==1){
      sprintf(tmpsamp,"scenariosamp.txt");
      fpn=fopen(tmpsamp, "w");
      output_samp(root, num_tips, num_anno, character, fpn);
      fclose(fpn);
    }


  return;
}
