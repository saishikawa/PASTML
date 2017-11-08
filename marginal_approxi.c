#include "asrml.h"
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>

extern double global_like;
extern double global_factor;
double global_correct, local_corrects_opt=0.0, local_corrects_best=0.0, local_possibs=0.0;

extern Tree* s_tree;
extern Node* root;
int output_count=0;
double root_corrects;
double edge_global=0., edge_global_best=0., edge_global_PPscore=0.;


double calc_correct(Node* nd, int nb, int nbanno, char **character){
  static int count=0, edge_count=0;
  static double possib=0;
  int i, j, k, tmpnum, node_start;
  double tmppossib, tmpmax=0.0, num_pos;
  static double correctness=0.0;
  int check=0;
  double root_correct=0.;
  int num_node, local_count, tmp_max;
  double local_numpos, edge_numpos, tmp_max_lik=0., tmp_correct=100.0, tmp_numpos, local_cor, node_corr=0.;
  double edge_cor=0.;
  
  num_node = s_tree->nb_nodes - s_tree->nb_taxa;

  if(nd->nneigh==1){
    return 0.0;
  }

  if(nd==root) {
    node_start=0;
  } else {
    node_start=1;
  }

  for(i=node_start;i<nd->nneigh;i++){
    tmppossib = calc_correct(nd->neigh[i],nb,nbanno,character);
  }

  //while(nd->marginal[nd->relax_count]==0.0){
    //nd->relax_count--;
  //}
  //possib+=nd->relax_count+1;
  //num_pos = 1.0 / (double) (nd->relax_count + 1);

 /*local_fraction_optimize*/
  for(i=0;i<nbanno;i++){
    local_cor=0.;
    for(j=0;j<nbanno;j++){
      if(j<=i){
        local_cor+=(nd->marginal[j]-1.0/((double)i+1.0))*(nd->marginal[j]-1.0/((double)i+1.0));
      } else {
        local_cor+=(nd->marginal[j]-0.0)*(nd->marginal[j]-0.0);
      }
    }
    if(tmp_correct > local_cor){
      tmp_correct = local_cor;
      tmpnum = i+1;
    }
  }
  for(i=0;i<nbanno;i++){
    if(i<tmpnum) {
      nd->local_flag[i]=1;
    } else {
      nd->local_flag[i]=0;
    }
  }
  local_numpos = 1.0 / (double) tmpnum;
  nd->pppos=tmpnum;
  nd->count=count;

  if(nd==root){
    edge_count=0;
    correctness=0.0;
    possib=0.0;
    output_count=1;
    return tmppossib;
  }
  count++;
  return 0.0;
}

void order_marginal(Node* nd, int nb, int nbanno){
  int i,j,k,ii,jj,tmpnum,stop=0,zerostate[nbanno],zero_count=0,node_start;
  double tmpmarginal[nbanno], tmpbest;
  static int count=0;
  int pupko;

  zero_count=0;
  if(nd->nneigh==1){
    count++;
    return;
  } else if (nd==root){
      pupko=1;
      for(j=0;j<nbanno;j++) {
        tmpmarginal[j]=root->mar_prob[j];
        if(tmpmarginal[j]==0.0){
          zerostate[zero_count]=j;
          zero_count++;
        }
      }
      zero_count=0;
      for(j=0;j<nbanno;j++){
	tmpbest=0.;
	tmpnum=-1;
	for(k=0;k<nbanno;k++){
	  /*if(k==nd->pupko_flag&&pupko==1){
            pupko=0;
            tmpbest=tmpmarginal[k];
	    tmpnum=k;
            break;
          }*/
	  if(tmpmarginal[k]>tmpbest) {
	    tmpbest=tmpmarginal[k];
	    tmpnum=k;
	  }
        }
	nd->marginal[j]=tmpbest;
	tmpmarginal[tmpnum]=0.;
	nd->mar_state[j]=tmpnum;
        nd->tmp_best[j]=tmpnum;
        if(tmpbest==0.0){
	  nd->mar_state[j]=zerostate[zero_count];
          nd->tmp_best[j]=zerostate[zero_count];
          zero_count++;
        }
      }
      count++;
  } else {
      pupko=1;
      for(j=0;j<nbanno;j++){
        tmpmarginal[j]=nd->condlike_mar[j];
        if(nd->condlike_mar[j]==0.0){
          zerostate[zero_count]=j;
          zero_count++;
        }
      }
      zero_count=0;
      for(j=0;j<nbanno;j++){
        tmpbest=0.;
        tmpnum=-1;
        for(k=0;k<nbanno;k++){
          /*if(k==nd->pupko_flag&&pupko==1){
            pupko=0;
	        tmpbest=tmpmarginal[k];
	        tmpnum=k;
            break;
          }*/
          if(tmpmarginal[k]>tmpbest) {
            tmpbest=tmpmarginal[k];
            tmpnum=k;
          }
        }
        nd->marginal[j]=tmpbest;
        tmpmarginal[tmpnum]=0.;
        nd->mar_state[j]=tmpnum;
        nd->tmp_best[j]=tmpnum;
        if(tmpbest==0.0){
	  nd->mar_state[j]=zerostate[zero_count];
          nd->tmp_best[j]=zerostate[zero_count];
          zero_count++;
        }
      }
      count++;
  }

  if(nd==root) {
    node_start=0;
  } else {
    node_start=1;
  }
  
  for(i=node_start;i<nd->nneigh;i++){
    order_marginal(nd->neigh[i], nb, nbanno);
  }

  if(nd==root) count=0;
  return;
}

void make_samples (char** tipnames, int* states, int num_tips, int num_anno, char **character, double* parameter){
   int i,j,k,ii,tmpfactor=0,power,border,out025=1,out05=1,out075=1,out090=1,depth,endloop;
   double fraction, samp_error, lik_marginal, factor_marginal;
   FILE *fp, *tru;
  int  true_count=0;
  char fname[100], tmp[100];
  double sum_lk=0., maxlnl=0., possib, mem_root_correct;
  int max_step;

  order_marginal(root, num_tips, num_anno);  
  printf("###for develop ### Ordering end\n");
  possib=calc_correct(root, num_tips, num_anno, character);

  printf("###for develop ### Marginal approximation end\n");
  sprintf(fname,"Result_treeIDs.%d.taxa.%d.states.tre",num_tips,num_anno);
  fp=fopen(fname, "w");

  /*write_result_tree*/
  write_nh_tree(s_tree, fp, parameter[num_anno], parameter[num_anno + 1]);
  fclose(fp);
  printf("Scaled tree with internal NodeID is written in Result_treeIDs.%d.taxa.%d.states.tre\n",num_tips,num_anno);

  sprintf(fname,"Result_states_probs.FULL.%d.taxa.%d.states.txt",num_tips,num_anno);
  fp=fopen(fname, "w");
  //output_state_PP(root,num_tips,num_anno,character,fp);
  output_state_anc_PP(root,num_tips,num_anno,character,fp);
  output_state_tip_PP(root,num_tips,num_anno,character,fp);
  fclose(fp);
  printf("Predictions for all internal Nodes and the Root are written in Result_states_probs.FULL.%d.taxa.%d.states.txt as a table format\n",num_tips,num_anno);

  /*write_simplified_tree*/
  /*sprintf(fname,"Result_summaryTree.%d.taxa.%d.states.tre",num_tips,num_anno);
  fp=fopen(fname, "w");
  summary(num_tips,num_anno,character,fp,1);
  fclose(fp);    
  printf("Summarized tree and prediction are written in Result_summaryTree.%d.taxa.%d.states.tre\n",num_tips,num_anno);

  sprintf(fname,"Result_states_probs.COLLAPSED.%d.taxa.%d.states.txt",num_tips,num_anno);
  fp=fopen(fname, "w");
  output_state_PP(root,num_tips,num_anno,character,fp);
  fclose(fp);
  printf("Predictions for only internal Nodes of collapsed tree are written in Result_states_probs.COLLAPSED.%d.taxa.%d.states.txt as a table format\n",num_tips,num_anno);

  sprintf(fname,"Result_IDs.COLLAPSED.%d.taxa.%d.states.txt",num_tips,num_anno);
  fp=fopen(fname, "w");
  output_state_IDs(root,num_tips,num_anno,character,fp,root->name);
  fclose(fp);
  printf("Relationships between collapsed nodes are written in Result_IDs.COLLAPSED.%d.taxa.%d.states.txt as a table format\n",num_tips,num_anno);*/
  return;
 
}
