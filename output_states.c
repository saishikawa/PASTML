#include "asrml.h"
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>

extern double global_like;
extern double global_factor;
extern Tree* s_tree;
extern Node* root;

void output_state_PP(Node* nd, int nb, int nbanno, char **character, FILE* outfile){
  int i, j, count=0, node_start, num=0;
  double sum=0.0;
  static int tip_counts[MAXCHAR], node_counts=0, num_collapsed_tips=0, root_node_counts;

  if(nd==root){
    node_counts=0;
    num_collapsed_tips=0;
    for(i=0;i<nbanno;i++){
      tip_counts[i]=0;
      root_node_counts=0;
    }
    fprintf(outfile,"NodeID");
    for(i=0;i<nbanno;i++){
      fprintf(outfile, ", %s", character[i]);
    }
    fprintf(outfile,", Size, Is_Tip");
    //for(i=0;i<nbanno;i++){
      //fprintf(outfile, ", Tips_%s", character[i]);
    //}
    fprintf(outfile,"\n");
  }
  if(nd->nneigh==1){
     for(i=0;i<nbanno;i++){
       if(i==nd->state_flag) tip_counts[i]++;
     }
     return;
  }

  if(nd==root) {
    node_start=0;
  } else {
    node_start=1;
  }

  for(i=0;i<nbanno;i++){
    if(nd->local_flag[i]==1) {
      num++;
      //nd->marginal[i]=nd->marginal[i]/sum;
    } else {
    }
  }
  for(i=0;i<nbanno;i++){
    if(nd->local_flag[i]==1) {
      nd->marginal[i]=(double)1.0/num;
    } else {
    }
  }
  
  for(i=node_start;i<nd->nneigh;i++){
    output_state_PP(nd->neigh[i],nb,nbanno,character,outfile);
  }

  if(nd->node_flag==1) {
    /*if(nd->neigh[0]==root){
      root_node_counts+=node_counts;
      //in full output root decendants are not represent, problem
    } else if (nd == root) {
      fprintf(outfile,"%s",nd->name);
      for(i=0;i<nbanno;i++){
        for(j=0;j<nbanno;j++){
          if(strcmp(character[i],character[nd->tmp_best[j]])==0){
            //printf("%s output PP of %s\n",nd->name,character[j]);
            if(nd->local_flag[j]==1) {
              fprintf(outfile, ", %.8f",nd->marginal[j]);
            } else {
              fprintf(outfile, ", 0.0");
            }
          }   
        }
      }
      fprintf(outfile,", %d, F",root_node_counts);
    } else {*/
      fprintf(outfile,"%s",nd->name);
      node_counts++;
      for(i=0;i<nbanno;i++){
        for(j=0;j<nbanno;j++){
          if(strcmp(character[i],character[nd->tmp_best[j]])==0){
            //printf("%s output PP of %s\n",nd->name,character[j]);
            if(nd->local_flag[j]==1) {
              fprintf(outfile, ", %.8f",nd->marginal[j]);
            } else {
              fprintf(outfile, ", 0.0");
            }
          }   
        }
      }
      fprintf(outfile,", %d, F",node_counts);
      for(i=0;i<nbanno;i++){
        if(tip_counts[i]>0){
          nd->tip_counts[i]=tip_counts[i];
          num_collapsed_tips++;
          fprintf(outfile,"\nTip%d",num_collapsed_tips);
          sprintf(nd->tip_names[i],"Tip%d",num_collapsed_tips);
          for(j=0;j<nbanno;j++){
            if(j==i) {
              fprintf(outfile, ", 1");
            } else {
              fprintf(outfile, ", 0");
            }
          }
          fprintf(outfile,", %d, T",tip_counts[i]);
        }
      }
      fprintf(outfile,"\n");
      node_counts=0;
      for(i=0;i<nbanno;i++){
        tip_counts[i]=0;
      }     
    //} 
  } else {
    node_counts++;
  }
  return;
}

void output_state_IDs(Node* nd, int nb, int nbanno, char **character, FILE* outfile, char* parent_name){

  int i, j, node_start;
  char* current_parent;
  //struct __Node* current;
  if(nd==root){
    fprintf(outfile,"NodeID, ParentID\n");
  }
  if(nd->nneigh==1){
     return;
  }
  if(nd->node_flag==1) {
    if(nd==root){
      current_parent=parent_name;
    } else {
      /*if(nd->neigh[0]==root){
        current_parent=parent_name;
      } else {*/

      fprintf(outfile,"%s, %s\n",nd->name,parent_name);
      for(i=0;i<nbanno;i++){
        if(nd->tip_counts[i]>0){
          fprintf(outfile,"%s, %s\n",nd->tip_names[i],nd->name);
        }
      }
      current_parent=nd->name;

      //}
    }
  } else {
    current_parent=parent_name;
  }

  if(nd==root) {
    node_start=0;
  } else {
    node_start=1;
  }

  for(i=node_start;i<nd->nneigh;i++){
    output_state_IDs(nd->neigh[i],nb,nbanno,character,outfile,current_parent);
  }

  return;
}


void output_state_anc(Node* nd, int nb, int nbanno, char **character, FILE* outfile){
  int i, j, count=0, node_start;
  double sum=0.0;

  if(nd==root){
    fprintf(outfile,"Internal NodeID");
    for(i=0;i<nbanno;i++){
      fprintf(outfile, ", %s", character[i]);
    }
    fprintf(outfile,"\n");
  }
  if(nd->nneigh==1){
     return;
  }

  if(nd==root) {
    node_start=0;
  } else {
    node_start=1;
  }
  
  for(i=node_start;i<nd->nneigh;i++){
    output_state_anc(nd->neigh[i],nb,nbanno,character,outfile);
  }
  if(nd->node_flag==1) {
    fprintf(outfile,"%s",nd->name);
    for(i=0;i<nbanno;i++){
     for(j=0;j<nbanno;j++){
       if(strcmp(character[i],character[nd->tmp_best[j]])==0){
        fprintf(outfile, ", %.8f",nd->marginal[j]);
       }
     }
    }
    fprintf(outfile,"\n");
  }
  return;
}

void output_state_tip(Node* nd, int nb, int nbanno, char **character, FILE* outfile){
  int i, node_start;

  if(nd->nneigh==1){
    fprintf(outfile,"%s",nd->name);
     for(i=0;i<nbanno;i++){
       if(i==nd->pupko_state){
         fprintf(outfile, ", 1.0");
       } else {
         fprintf(outfile, ", 0.0");
       }
     }
     fprintf(outfile,"\n");
     return;
  }

  if(nd==root) {
    node_start=0;
  } else {
    node_start=1;
  }
  
  for(i=node_start;i<nd->nneigh;i++){
    output_state_tip(nd->neigh[i],nb,nbanno,character,outfile);
  }

  return;
}

void output_state_anc_PP(Node* nd, int nb, int nbanno, char **character, FILE* outfile){
  int i, j, count=0, node_start, num=0;
  double sum=0.0;

  if(nd==root){
    fprintf(outfile,"Internal NodeID");
    for(i=0;i<nbanno;i++){
      fprintf(outfile, ", %s", character[i]);
    }
    fprintf(outfile,"\n");
  }
  if(nd->nneigh==1){
     return;
  }

  if(nd==root) {
    node_start=0;
  } else {
    node_start=1;
  }
  
  for(i=0;i<nbanno;i++){
    if(nd->local_flag[i]==1) {
      num++;
      //nd->marginal[i]=nd->marginal[i]/sum;
    } else {
    }
  }
  for(i=0;i<nbanno;i++){
    if(nd->local_flag[i]==1) {
      nd->marginal[i]=(double)1.0/num;
    } else {
    }
  }

  for(i=node_start;i<nd->nneigh;i++){
    output_state_anc_PP(nd->neigh[i],nb,nbanno,character,outfile);
  }

  if(nd->node_flag==1) {
    fprintf(outfile,"%s",nd->name);
    for(i=0;i<nbanno;i++){
     for(j=0;j<nbanno;j++){

       if(strcmp(character[i],character[nd->tmp_best[j]])==0){
            if(nd->local_flag[j]==1) {
              fprintf(outfile, ", %.8f",nd->marginal[j]);
            } else {
              fprintf(outfile, ", 0.0");
            }
       }
     }
    }
    fprintf(outfile,"\n");
  }
  return;
}

void output_state_tip_PP(Node* nd, int nb, int nbanno, char **character, FILE* outfile){
  int i, node_start;

  if(nd->nneigh==1){
    fprintf(outfile,"%s",nd->name);
     for(i=0;i<nbanno;i++){
       if(i==nd->pupko_state){
         fprintf(outfile, ", 1.0");
       } else {
         fprintf(outfile, ", 0.0");
       }
     }
     fprintf(outfile,"\n");
     return;
  }

  if(nd==root) {
    node_start=0;
  } else {
    node_start=1;
  }
  
  for(i=node_start;i<nd->nneigh;i++){
    output_state_tip_PP(nd->neigh[i],nb,nbanno,character,outfile);
  }

  return;
}

