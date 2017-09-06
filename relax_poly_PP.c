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

  while(nd->marginal[nd->relax_count]==0.0){
    nd->relax_count--;
  }
  possib+=nd->relax_count+1;
  num_pos = 1.0 / (double) (nd->relax_count + 1);

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


int sampling(Node* nd, int nb, int nbanno, int ancstate, int tmpnode, int tmpstate, int rep, char** tipnames, int* states, char str[100], char **character) {
  int i,j,ii, tmpfactor=0, count, num_anno, node_start;
  double mul, expmul,  prob_left=0., prob_right=0., smallest, scaled_lk=0., logroot;
  static int factors=0;
  static int ndcount=0;
  double curr_scaler;
  int curr_scaler_pow, piecewise_scaler_pow;
  char str2[50];

  count=ndcount;
  if(ndcount==tmpnode){
      for(i=0;i<=nd->relax_count;i++){
        if(i==nd->relax_count){
          nd->marginal[i]=0.;
          nd->calc_flag[nd->mar_state[i]]=0;
          nd->mar_state[i]=-1;
        }
      }
      nd->relax_count--;
      sprintf(str, "%s:%s", nd->name, character[nd->tmp_best[0]]);
      for(i=1;i<=nd->relax_count;i++){
        sprintf(str2, "->%s", character[nd->tmp_best[i]]);
        strcat(str,str2);
      }
  }
  num_anno=nd->relax_count+1;
  ndcount++;

  if(nd->nneigh==1){
     return 0;
  }

  if(nd==root) {
    node_start=0;
  } else {
    node_start=1;
  }
  
  for(i=node_start;i<nd->nneigh;i++){
    tmpfactor=sampling(nd->neigh[i], nb, nbanno, ancstate, tmpnode, tmpstate, rep, tipnames, states, str, character);
  }

  if(nd==root){
    ndcount=0;
    factors=0;
    return tmpfactor;
  } else {
  } 
  return tmpfactor;
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
        if(root->mar_prob[j]==0.0){
          zerostate[zero_count]=j;
          zero_count++;
        }
      }
      zero_count=0;
      for(j=0;j<nbanno;j++){
	tmpbest=0.;
	tmpnum=-1;
	for(k=0;k<nbanno;k++){
	  if(k==nd->pupko_flag&&pupko==1){
            pupko=0;
            tmpbest=tmpmarginal[k];
	    tmpnum=k;
            break;
          }
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
          if(k==nd->pupko_flag&&pupko==1){
            pupko=0;
	        tmpbest=tmpmarginal[k];
	        tmpnum=k;
            break;
          }
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

void nd_count(Node* nd, int nbanno){
  int i, node_start;
  
  if(nd->nneigh==1){
    return;
  } else {
    nd->relax_count=nbanno-1;
  }

  if(nd==root) {
    node_start=0;
  } else {
    node_start=1;
  }
  
  for(i=node_start;i<nd->nneigh;i++){
    nd_count(nd->neigh[i], nbanno);
  }

  return;
}

void reset_flags(Node* nd, int nbanno){
  int i, node_start;
  
  if(nd->nneigh==1){
    return;
  } else {
    for(i=0;i<nbanno;i++) nd->calc_flag[i]=1;
  }

  if(nd==root) {
    node_start=0;
  } else {
    node_start=1;
  }
  
  for(i=node_start;i<nd->nneigh;i++){
    reset_flags(nd->neigh[i], nbanno);
  }

  return;
}

void marginal_reset(Node* nd, int nbanno){
  int i, node_start;
  
  if(nd->nneigh==1){
    return;
  } else {
    for(i=0;i<nbanno;i++){
      nd->marginal[i]=0.0;
    }
  }

  if(nd==root) {
    node_start=0;
  } else {
    node_start=1;
  }
  
  for(i=node_start;i<nd->nneigh;i++){
    marginal_reset(nd->neigh[i], nbanno);
  }

  return;
}

void reset_nd(Node* nd, int nb, int nbanno){
  int i, j, node_start;

  if(nd->nneigh==1){
    return;
  } else {
    nd->up_factor=0;
    nd->down_factor=0;
    for(i=0;i<nbanno;i++){
      nd->condlike[i]=0.;
      nd->condlike_mar[i]=0.;
      nd->up_like[i]=0.;
      nd->sum_down[i]=0.;
    }
  }

  if(nd==root) {
    node_start=0;
  } else {
    node_start=1;
  }
  
  for(i=node_start;i<nd->nneigh;i++){
    reset_nd(nd->neigh[i], nb, nbanno);
  }

  return;
}

int relaxing(Node* nd, int nb, int nbanno, int rep, char** tipnames, int* states, char **character ) {
  int tmpfactor=0, judge=1, i, node_start;
  static int count=0, tmpstate=-1, tmpnode=-1;
  static double tmpmar=10.;
  static double prob=1.0;
  char str[100];
  
  if(nd->nneigh==1){
    count++;
    return tmpfactor;
  } else {
    //printf("%s=%d", nd->name, nd->relax_count);
    //for(i=0;i<=nd->relax_count;i++) printf(", %d = %.12f", nd->mar_state[i], nd->marginal[i]);
    //printf("\n");
    if(nd->relax_count!=0&&nd->marginal[nd->relax_count]!=0.0){
      if(nd->marginal[nd->relax_count]<prob){
         prob=nd->marginal[nd->relax_count];
         tmpstate=nd->mar_state[nd->relax_count];
         tmpnode=count;
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
    tmpfactor=relaxing(nd->neigh[i], nb, nbanno, rep, tipnames, states, character);
  }

  if(nd==root){
    if(tmpstate!=-1){
      tmpfactor=sampling(root, nb, nbanno, -1, tmpnode, tmpstate, rep, tipnames, states, str, character);
      printf("\nStep %d, Relax at %s, ", rep+1, str);
    } else {
      tmpfactor=-10;
    }
    count=0;
    tmpstate=-1;
    tmpnode=-1;
    tmpmar=10.;
    prob=1.0;
  }
  return tmpfactor;
}

void make_samples (char** tipnames, int* states, int num_tips, int num_anno, char **character, double mu, double scale, double* frequency, char* model, double border_frac){
   int i,j,k,ii,tmpfactor=0,power,border,out025=1,out05=1,out075=1,out090=1,depth,endloop;
   double fraction, samp_error, lik_marginal, factor_marginal;
   FILE *fp, *tru;
  int  true_count=0;
  char fname[100], tmp[100];
  double sum_lk=0., maxlnl=0., possib, mem_root_correct;
  int max_step;

  nd_count(root,num_anno);
  reset_flags(root,num_anno);
  marginal_reset(root,num_anno);
  reset_nd(root,num_tips,num_anno);
  calc_lik(root, tipnames, states, num_tips, num_anno, mu, scale, model, frequency, &maxlnl);
  down_like_marginal(root, num_tips, num_anno, mu, scale, frequency);
  order_marginal(root, num_tips, num_anno);
  tmpfactor=sampling(root,num_tips,num_anno,-1,-1,-1,0,tipnames,states,tmp,character);
  fraction=1.0;
  lik_marginal=global_like;
  factor_marginal=global_factor;
  
  //sprintf(fname,"Fraction_Correctness.txt");
  //fp=fopen(fname, "a");

  if(border_frac==0.0) {
    possib=calc_correct(root, num_tips, num_anno, character);
    sprintf(fname,"Result_treeIDs.%d.taxa.%d.states.tre",num_tips,num_anno);
    fp=fopen(fname, "w");

    /*write_result_tree*/
    write_nh_tree(s_tree, fp, scale);
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
    sprintf(fname,"Result_summaryTree.%d.taxa.%d.states.tre",num_tips,num_anno);
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
    printf("Relationships between collapsed nodes are written in Result_IDs.COLLAPSED.%d.taxa.%d.states.txt as a table format\n",num_tips,num_anno);
    return;
  }

  for(i=0;i<((num_tips - 1)*(num_anno - 1));i++){
    printf("Likelihood fraction methods\n");
    tmpfactor=relaxing(root, num_tips, num_anno, i, tipnames, states, character);
    printf("Relaxing\n");
    reset_nd(root,num_tips,num_anno);
    calc_lik(root, tipnames, states, num_tips, num_anno, mu, scale, model, frequency, &maxlnl);
    down_like_marginal(root, num_tips, num_anno, mu, scale, frequency);
    order_marginal(root, num_tips, num_anno);

    power=factor_marginal - global_factor;
    sum_lk=global_like*pow(2.0,(double)power);
    fraction=sum_lk/lik_marginal;

    if(fraction <= border_frac){
      printf("Stopped at the likelihood fraction = %lf...\n",border_frac);
      sprintf(fname,"Result_treeIDs.%d.taxa.%d.states.tre",num_tips,num_anno);
      fp=fopen(fname, "w");
      /*write_result_tree*/
      write_nh_tree(s_tree, fp, scale);
      fclose(fp);
      printf("Scaled tree with internal NodeID is written in Result_treeIDs.%d.taxa.%d.states.tre\n",num_tips,num_anno);
      
      sprintf(fname,"Result_states_probs.FULL.%d.taxa.%d.states.txt",num_tips,num_anno);
      fp=fopen(fname, "w");
      output_state_anc(root,num_tips,num_anno,character,fp);
      //output_state_tip(root,num_tips,num_anno,character,fp);
      fclose(fp);
      printf("Predictions for all internal Nodes and the Root are written in Result_states_probs.FULL.%d.taxa.%d.states.txt as a table format\n",num_tips,num_anno);

      /*write_simplified_tree*/
      sprintf(fname,"Result_summaryTree.%d.taxa.%d.states.tre",num_tips,num_anno);
      fp=fopen(fname, "w");
      summary(num_tips,num_anno,character,fp,0);
      fclose(fp);    
      printf("Summarized tree and prediction are written in Result_summaryTree.%d.taxa.%d.states.tre\n",num_tips,num_anno);

      sprintf(fname,"Result_states_probs.COLLAPSED.%d.taxa.%d.states.txt",num_tips,num_anno);
      fp=fopen(fname, "w");
      output_state_anc(root,num_tips,num_anno,character,fp);
      //output_state_tip(root,num_tips,num_anno,character,fp);
      fclose(fp);
      printf("Predictions for only internal Nodes of collapsed tree are written in Result_states_probs.COLLAPSED.%d.taxa.%d.states.txt as a table format\n",num_tips,num_anno);
      break;
    }
  }
 
  return;
}
