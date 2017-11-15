#include "asrml.h"
#include "nrutil.h"
#define ALF 1.0e-4
#define ITMAX 200
#define EPS 3.0e-8
#define TOLX 1.0e-6
#define TOL 1.0e-7
#define SIGMA 1.0e-3
#define STPMX 100.0
#define STEP 1.0e-5
#define STEP2 1.0e-4
#define STEP3 1.0e-9
#define SIMPLEX 1.0e-3
#define SCAL_MAX 10.0
#define SCAL_MIN 1.0e-4

#define FREEALL free_vector(xi,1,n);free_vector(pnew,1,n);free_matrix(hessin,1,n,1,n);free_vector(hdg,1,n);free_vector(g,1,n);free_vector(dg,1,n);

extern Tree* s_tree;
extern Node* root;

void gradient(Node *nd, char** tipnames, int* states, int nb, int nbanno, double mu, char* model, double* p, double* g, int n) {
  int i,j;
  double pluslnl, oldlnl, minuslnl, tmp_p[n], sum=0.0;

  calc_lik_bfgs(root, tipnames, states, nb, nbanno, mu, model, p, &oldlnl);
  //printf("Present lnL = %.15f, scale = %lf, epsilon = %e\n",oldlnl,p[nbanno],p[nbanno+1]);
  for(i=0;i<n;i++){
    sum=0.0;
    for(j=0;j<n;j++){
      if(j==i){
        if(j < nbanno) {
          tmp_p[j] = p[j] + STEP;
        } else if (j==nbanno) {
          tmp_p[j] = p[j] + STEP2;
        } else if (j==nbanno+1) {
          tmp_p[j] = p[j] + STEP3;
        }
      } else {
        tmp_p[j] = p[j] ;
      }
      if(j < nbanno) sum+=tmp_p[j];
    }
    for(j=0;j<nbanno;j++) {
      tmp_p[j] = tmp_p[j]/sum;
    }
    calc_lik_bfgs(root, tipnames, states, nb, nbanno, mu, model, tmp_p, &pluslnl);
    //if(i==n-1) printf("plus = %.15f, param = %.10e\n",pluslnl,tmp_p[i]);
    if(i<nbanno) {
      g[i] = (fabs(pluslnl) - fabs(oldlnl)) / (STEP);
    } else if (i==nbanno) {
      g[i] = (fabs(pluslnl) - fabs(oldlnl)) / (STEP2*10);
    } else if (i==nbanno+1) {
      g[i] = (fabs(pluslnl) - fabs(oldlnl)) / (STEP3*1.0e+4);
    }
    //printf("diff = %e, %d th gradient = %.10f\n",fabs(pluslnl) - fabs(oldlnl),i,g[i]);
  }
  
}

void lnsrch(Node *nd, char** tipnames, int* states, int nb, int nbanno, double mu, char* model, int n, double* xold, double fold, double* g, double* p, double* x, double *f, double stpmax, int *check, double* p_x) {

/*Given an n-dimensional point xold[1..n] , the value of the function and gradient there, fold
and g[1..n] , and a direction p[1..n] , finds a new point x[1..n] along the direction p from
xold where the function func has decreased “sufficiently.” The new function value is returned
in f . stpmax is an input quantity that limits the length of the steps so that you do not try to
evaluate the function in regions where it is undefined or subject to overflow. p is usually the
Newton direction. The output quantity check is false (0) on a normal exit. It is true (1) when
x is too close to xold . In a minimization algorithm, this usually signals convergence and can
be ignored. However, in a zero-finding algorithm the calling program should check whether the
convergence is spurious. Some “difficult” problems may require double precision in this routine.*/

  int i,j;
  double a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,sum,temp,test,tmplam,sum_freq=0.0,exp_pi[nbanno],fold2;
  double* local_xold, sum_x, tmp_px, scale_up, scale_down;

  scale_up = 1.0 / s_tree->avgbl;
  scale_down = 1.0e-4 / s_tree->min_bl;  
  alam = alam2 = f2 = fold2 = tmplam = .0;
  local_xold = calloc(n,sizeof(double));
  for(i=0;i<n;i++) local_xold[i] = xold[i];
  *check = 0;
  for(sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
  sum = sqrt(sum);
  if (sum > stpmax) {
    for(i=0;i<n;i++) p[i] *= stpmax/sum;
  }
  for(slope=0.0,i=0;i<n;i++) {
    slope += g[i]*p[i];
  }
  if(slope >= 0.0) nrerror("Roundoff problem in lnsrch\n");
  test=0.0;
  for(i=0;i<n;i++) {
    temp=fabs(p[i])/FMAX(fabs(local_xold[i]),1.0);
    if (temp > test) test=temp;
  }
  alamin=TOL/test;
  alam=1.0;
  for(;;) {
    sum_freq=0.0;
    sum_x=0.0;
    for(i=0;i<nbanno;i++) {
      //p_x[i] = ((p_x[i]*100) + alam*p[i])/100;
      //p_x[i] = ((p_x[i]) + alam*p[i]);
      
      /*if(i==0) { 
        p_x[i] = 0.0;
      } else {
        p_x[i] = ((p_x[i])*100 + alam*p[i])/100;
      }
      tmp_px = (1.0 / (1.0 + exp(p_x[i])) )+ SIGMA;
      sum_x+=tmp_px;
      p_x[i] = tmp_px;*/

      x[i]=(xold[i]*100+alam*p[i])/100;
      //printf("next vector %d = %lf, old %lf, alam %e, p %e\n",i,x[i],xold[i],alam,p[i]);
      /*if(i==nbanno) {
        if(x[i] < SCAL_MIN) x[i] = xold[i];
        if(x[i] > SCAL_MAX) x[i] = xold[i];
      }
      if(i==nbanno+1) {
        if(x[i] < DBL_MIN*1.0e+4) x[i] = xold[i];
        if(x[i] > s_tree->min_bl*1.0e+4) x[i] = xold[i];
      }*/
      sum_freq+=fabs(x[i]);
      //xold[i]=x[i];
    }
    for(i=0;i<nbanno;i++) {
      //x[i] =  p_x[i] / sum_x;
      x[i] = fabs(x[i]) / sum_freq;
      //printf("next vector %d = %lf, old %lf, alam %e, p %e\n",i,x[i],xold[i],alam,p[i]);
    }
    x[nbanno] = ((xold[nbanno]*10) + alam*p[i])/10;
    //printf("next vector %d = %lf, old %lf, alam %e, p %e\n",nbanno,x[nbanno],xold[nbanno],alam,p[nbanno]);
    x[nbanno+1] =  ((xold[nbanno+1]*1.0e+4) + alam*p[i])*1.0e-4;
    //printf("next vector %d = %lf, old %lf, alam %e, p %e\n",nbanno+1,x[nbanno+1],xold[nbanno+1],alam,p[nbanno+1]);

    if(x[nbanno] < SCAL_MIN || x[nbanno] > scale_up) x[nbanno] = xold[nbanno];
    if(x[nbanno+1] < DBL_MIN || x[nbanno+1] > s_tree->min_bl) x[nbanno+1] = xold[nbanno+1];
    for(i=0;i<n;i++){
      xold[i] = x[i];
    }
    /*for(i=0;i<nbanno;i++) {
      x[i] = fabs(x[i]) / sum_freq;
    }
    x[nbanno] = x[nbanno] / 10;
    x[nbanno+1] = x[nbanno+1] * 1.0e-4;*/
    

    calc_lik_bfgs(root, tipnames, states, nb, nbanno, mu, model, x, f); /*likelihood calculation*/  

    /*for(i=0;i<nbanno;i++) {
      x[i] = x[i]*100;
    }
    x[nbanno] = x[nbanno] * 10;
    x[nbanno+1] = x[nbanno+1] * 1.0e+4;*/

    if (alam < alamin) {
      for(i=0;i<n;i++) xold[i]=local_xold[i];
      *check = 1;
      free(local_xold);
      return;
    } else if (*f <= fold + ALF * alam * slope) {
      for(i=0;i<n;i++) xold[i]=local_xold[i];
      free(local_xold);
      return;
    } else {
      if(alam==1.0) {
        tmplam = -slope/(2.0*(*f - fold - slope));
      } else {
        rhs1 = *f-fold-alam*slope;
        rhs2=f2-fold2-alam2*slope;
        a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
        b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
        if (a == 0.0) {
          tmplam = -slope/(2.0*b);
        } else {
          disc=b*b-3.0*a*slope;
        }
        if (disc < 0.0) {
          tmplam=0.5*alam;
        } else if (b <= 0.0) {
          tmplam=(-b+sqrt(disc))/(3.0*a);
        } else {
          tmplam=-slope/(b+sqrt(disc));
        }
        if (tmplam > 0.5*alam) tmplam=0.5*alam;
      }
    }
    alam2=alam;
    f2=*f;
    fold2=fold;
    alam=FMAX(tmplam,0.1*alam);
  }
}

void frprmn(Node *nd, char** tipnames, int* states, int nb, int nbanno, double mu, char* model, double *p, int n, int *iter, double *fret) {

/*Given a starting point p[1..n] that is a vector of length n , the Broyden-Fletcher-Goldfarb-
Shanno variant of Davidon-Fletcher-Powell minimization is performed on a function func , using
its gradient as calculated by a routine dfunc . The convergence requirement on zeroing the
gradient is input as gtol . Returned quantities are p[1..n] (the location of the minimum),
iter (the number of iterations that were performed), and fret (the minimum value of the
function). The routine lnsrch is called to perform approximate line minimizations.*/

  int check,i,its,j;
  double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test,fpold,tmp_fpbest,scale_up,scale_down;
  double *dg,*g,*hdg,**hessin,*pnew,*xi,*pold,diff,*init_p, *tmp_pbest, *p_x;

  scale_up = 1.0 / s_tree->avgbl;
  //scale_down = 1.0e-4 / s_tree->min_bl;
  p[nbanno] = (scale_up - SCAL_MIN) / 2.0;
  printf("upbound = %lf, start = %lf\n",scale_up,p[nbanno]);
  check=0;
  dg=vector(0,n-1);
  g=vector(0,n-1);
  hdg=vector(0,n-1);
  hessin=matrix(0,n-1,0,n-1);
  pnew=vector(0,n-1);
  pold=vector(0,n-1);
  init_p=vector(0,n-1);
  xi=vector(0,n-1);
  p_x=vector(0,n-1);
  tmp_pbest=vector(0,n-1);
  tmp_fpbest=DBL_MAX;

  calc_lik_bfgs(root, tipnames, states, nb, nbanno, mu, model, p, &fp); /*likelihood calculation*/  
  fpold = fp;
  for(i=0;i<n;i++) {
    p_x[i]=0.0;
    g[i]=0.0;
    init_p[i]=p[i];
  }
  gradient(root, tipnames, states, nb, nbanno, mu, model, p, g, n); /*gradient calculation*/

  for(i=0;i<n;i++){
    /*if(i<nbanno) {
      p[i] = p[i]*100;
    } else if(i==nbanno) {
      p[i] = p[i]*10;
    } else if(i==nbanno+1) {
      p[i] = p[i]*1.0e+4;
    } else {
    }*/
    for(j=0;j<n;j++) hessin[i][j]=0.0;
    hessin[i][i]=1.0;
    xi[i] = -g[i];
    sum += p[i]*p[i];
    //printf("\n. BFGS x[%d]: %f p[i]: %f",i,xi[i],p[i]); 
  }
  stpmax=STPMX*FMAX(sqrt(sum),(double)n);
  for(its=0;its<ITMAX;its++) {
    *iter=its; 
    for(i=0;i<n;i++) pold[i]=p[i];

    lnsrch(root, tipnames, states, nb, nbanno, mu, model, n, p, fp, g, xi, pnew, fret, stpmax, &check, p_x);
    /*The new function evaluation occurs in lnsrch; save the function value in fp for the next line search. It is usually safe to ignore the value of check.*/

    diff = (*fret) - fp;
    diff = fabs(diff);
    printf("lnL = %lf, ", *fret);
    fpold = fp;
    fp = *fret;
    for(i=0;i<n;i++){
      xi[i]=pnew[i]-pold[i];
      if(i==nbanno) {
        if(pnew[i] < SCAL_MIN) pnew[i] = pold[i];
        if(pnew[i] > scale_up) pnew[i] = pold[i];
      }
      if(i==nbanno+1) {
        if(pnew[i] < DBL_MIN) pnew[i] = pold[i];
        if(pnew[i] > s_tree->min_bl) pnew[i] = pold[i];
      }
      p[i]=pnew[i];
    }
    test=0.0;
    for(i=0;i<n;i++){
      temp=fabs(xi[i])/FMAX(fabs(p[i]),1.0);
      if(temp > test) test=temp;
    }
    if (fabs(fp) > fabs(fpold)) {
      //for(i=0;i<n;i++) p[i] = pold[i];
    } else {
      if (test < TOLX || (diff < SIMPLEX && its > 1)) {
        //FREEALL
        if(fp < fpold - SIGMA) {

              for(i=0;i<n;i++) p[i] = pold[i];
              //printf("Final likelihood is not converged, bigger than before ... \n");
        }
        if(fabs(fp) > fabs(tmp_fpbest)){
          for(i=0;i<n;i++) p[i] = tmp_pbest[i];
        /*for(i=0;i<nbanno;i++){
          p[i] = p[i]*100;
        }
        p[nbanno] = p[nbanno] * 10;
        p[nbanno+1] = p[nbanno+1] * 1.0e+4;*/
        }
        return;
      }
    }
    for(i=0;i<n;i++) dg[i]=g[i];
    /*for(i=0;i<nbanno;i++){
      p[i] = p[i]/100;
    }
    p[nbanno] = p[nbanno] / 10;
    p[nbanno+1] = p[nbanno+1] * 1.0e-4;*/
    for(i=0;i<n;i++) printf("param%d=%.5e, ",i,p[i]);
    if(fabs(tmp_fpbest) > fabs(fp)){
      tmp_fpbest = fp;
      for(i=0;i<n;i++) tmp_pbest[i] = p[i];
      //printf("tmplnl = %lf",tmp_fpbest);
    }
    printf("\n");
    gradient(root, tipnames, states, nb, nbanno, mu, model, p, g, n); /*gradient calculation*/
    /*for(i=0;i<nbanno;i++){
      p[i] = p[i]*100;
    }
    p[nbanno] = p[nbanno] * 10;
    p[nbanno+1] = p[nbanno+1] * 1.0e+4;*/

    test=0.0;
    den=FMAX(*fret,1.0);
    for(i=0;i<n;i++) {
      temp=fabs(g[i])*FMAX(fabs(p[i]),1.0)/den;
      if (temp > test) test=temp;
    }
    if (fabs(fp) < fabs(fpold)) {
      if (test < TOLX || (diff < SIMPLEX && its > 1)) {
        if(fabs(fp) > fabs(tmp_fpbest)){
          for(i=0;i<n;i++) p[i] = tmp_pbest[i];
          /*for(i=0;i<nbanno;i++){
            p[i] = p[i]*100;
          }
          p[nbanno] = p[nbanno] * 10;
          p[nbanno+1] = p[nbanno+1] * 1.0e+4;*/
        }
        //FREEALL
        //printf("finally likelihood optimised = %lf\n",*fret);
        return;
      }
    }
    for(i=0;i<n;i++) dg[i]=g[i]-dg[i];
    for(i=0;i<n;i++) {
      hdg[i]=0.0;
      for (j=0;j<n;j++) hdg[i] += hessin[i][j]*dg[j];
    }
    fac=fae=sumdg=sumxi=0.0;
    for(i=0;i<n;i++) {
      fac += dg[i]*xi[i];
      fae += dg[i]*hdg[i];
      sumdg += SQR(dg[i]);
      sumxi += SQR(xi[i]);
    }
    if (fac > sqrt(EPS*sumdg*sumxi)) {
      fac=1.0/fac;
      fad=1.0/fae;
      for(i=0;i<n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
      for(i=0;i<n;i++) {
        for (j=0;j<n;j++) {
          hessin[i][j] += fac*xi[i]*xi[j]-fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
        }
      }
    }
    for(i=0;i<n;i++) {
      xi[i]=0.0;
      for (j=0;j<n;j++) xi[i] -= hessin[i][j]*g[j];
    }
  }
  nrerror("parameter or likelihood is NOT converged untill 1000th step\n*** PASTML returns a possible best solution found during optimization ***\n");
  for(i=0;i<n;i++) p[i] = tmp_pbest[i];
  //exit(0);
  return;
  //FREEALL;
}
