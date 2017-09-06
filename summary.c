#include "asrml.h"
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>

extern Tree* s_tree;
extern Node* root;

void *check_alloc(int nbrelt, int sizelt);

int ctreewrite(Node* node, char* ctree, char **character, int num_tips, int nbanno, int check){
  int i,j,same=1,tmpsame,here=-1,come=-1;
  char str[100], str2[100];
  int num;
  node->node_flag=0;
  if (node->nneigh==1){
    node->state_flag=node->pupko_state;
    sprintf(str,"%s", node->name);
    strcat(ctree,str);
    return 0;
  } else if (node==root){
    same=1;
    strcat(ctree,"(");
    num=ctreewrite(node->neigh[0], ctree, character, num_tips, nbanno, check);
  } else{
    if(check==0){
      here=node->relax_count;
      come=node->neigh[0]->relax_count;
    } else {
      for(i=0;i<nbanno;i++){
        if(node->local_flag[i]==1)here++;
        if(node->neigh[0]->local_flag[i]==1)come++;
      }
    }
    if(here==come){
      for(i=0;i<=come;i++){
        tmpsame=1;
        for(j=0;j<=here;j++){
          if(node->neigh[0]->tmp_best[i]==node->tmp_best[j]) tmpsame=0;
        }
        if(tmpsame==0){
          same=0;
        } else {
          same=1;
        }
      }
    } else {
      same=1;
    }
   
    if(node->neigh[0]==root){
      strcat(ctree,"(");
    } else if(same!=0){
      strcat(ctree,"(");
    }
    num=ctreewrite(node->neigh[1], ctree, character, num_tips, nbanno, check);
  }
  
  strcat(ctree,",");
  if(node==root) {
    num=ctreewrite(node->neigh[1], ctree, character, num_tips, nbanno, check);
  } else {
    num=ctreewrite(node->neigh[2], ctree, character, num_tips, nbanno, check);
  }
  
  node->node_flag=0;
  //if(node!=root){
   if(same!=0 || node->neigh[0]==root){
     strcat(ctree,")");
     sprintf(str,"%s", node->name);
     if(check==0){
       node->node_flag=1;
       for(j=0;j<=node->relax_count;j++){
          sprintf(str2, "_%s", character[node->tmp_best[j]]);
          strcat(str,str2);
       }
     } else {
      node->node_flag=1;
      for(i=0;i<nbanno;i++){
        for(j=0;j<nbanno;j++){
          if(strcmp(character[i],character[node->tmp_best[j]])==0){
            if(node->local_flag[j]==1) {
              sprintf(str2, "_%s", character[node->tmp_best[j]]);
              strcat(str,str2);
            }
          }
        }
      }
     }
     strcat(ctree,str);
   }
   //}
  return 0;
}

void summary(int nb, int nbanno, char **character, FILE *outfile, int check){
  char *ctree;
  int num;

  ctree=(char*)check_alloc(200*nb, sizeof(char));
  num=ctreewrite(root, ctree, character, nb, nbanno, check);
  fprintf(outfile, "%s\n", ctree);
}

