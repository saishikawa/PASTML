#include "asrml.h"

#define SIGMA 1e-5
#define sqrt5 2.236067977499789696
#define ALPHA 1e-5
#define gratio 1.6180339887498948482045868343656

extern Tree* s_tree;
extern Node* root;

int Stopping_Rule(double x0, double x1) {
  double xm=0.5*fabs(x0+x1);
  
  if(xm<=1.0){
     if ( fabs( x1 - x0 ) < ALPHA ){
       return 0;
     } else {
       return 1;
     }
  } else {
     if ( fabs( x1 - x0 ) < ALPHA*xm){
       return 0;
     } else {
       return 1;
     }    
  }
}

/*double golden(char** tipnames, int* states, int nb, int nbanno, double mu, double ub, char* model, double* frequency, double* scale){
  double lb=0.0, x1, x2, fx1=0., fx2=0., diff;
  int count=0;

  x1=ub+(lb-ub)/gratio;
  x2=lb+(ub-lb)/gratio;
  calc_lik(root,tipnames,states,nb,nbanno,mu,x1, model, frequency, &fx1);
  calc_lik(root,tipnames,states,nb,nbanno,mu,x2, model, frequency, &fx2);
  *scale=(x1+x2)/2;
  diff=fx1-fx2;

  while(Stopping_Rule(lb,ub)==1){
    if(fabs(diff)<SIGMA) break;
    if(diff < 0.0){
      lb=x1;
      x1=x2;
      fx1=fx2;
      x2=lb+(ub-lb)/gratio;
      calc_lik(root,tipnames,states,nb,nbanno,mu,x2, model, frequency,&fx2);
      *scale=(x1+x2)/2;
      diff=fx1-fx2;
    } else if (diff > 0.0) {
      ub=x2;
      x2=x1;
      fx2=fx1;
      x1=ub+(lb-ub)/gratio;
      calc_lik(root,tipnames,states,nb,nbanno,mu,x1, model, frequency,&fx1);
      *scale=(x1+x2)/2;
      diff=fx1-fx2;
    }
    count++;
  }
  *scale=(x1+x2)/2.0;
  return (fx1+fx2)/2.0;
}*/

void golden(char** tipnames, int* states, int nb, int nbanno, double mu, double ub, char* model, double* frequency, double* scale){
  double lb=0.0, x1, x2, fx1=0., fx2=0., diff;
  int count=0;

  x1=ub+(lb-ub)/gratio;
  s_tree->scale_B=x1;
  calc_lik(root,tipnames,states,nb,nbanno,mu,x1, model, frequency, &fx1);
  x2=lb+(ub-lb)/gratio;
  s_tree->scale_B=x2;
  calc_lik(root,tipnames,states,nb,nbanno,mu,x2, model, frequency, &fx2);
  *scale=(x1+x2)/2;
  diff=fx1-fx2;

  while(Stopping_Rule(lb,ub)==1){
    if(fabs(diff)<SIGMA) break;
    if(diff < 0.0){
      lb=x1;
      x1=x2;
      fx1=fx2;
      x2=lb+(ub-lb)/gratio;
      s_tree->scale_B=x2;
      calc_lik(root,tipnames,states,nb,nbanno,mu,x2, model, frequency,&fx2);
      *scale=(x1+x2)/2;
      diff=fx1-fx2;
    } else if (diff > 0.0) {
      ub=x2;
      x2=x1;
      fx2=fx1;
      x1=ub+(lb-ub)/gratio;
      s_tree->scale_B=x1;
      calc_lik(root,tipnames,states,nb,nbanno,mu,x1, model, frequency,&fx1);
      *scale=(x1+x2)/2;
      diff=fx1-fx2;
    }
    count++;
  }
  *scale=(x1+x2)/2.0;
  //printf("fx1 = %lf, fx2 = %lf\n",fx1, fx2);
  s_tree->gold_1=fx1;
  s_tree->gold_2=fx2;
}
