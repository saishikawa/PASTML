#include "pastml.h"
#include "nrutil.h"
#include "lik.h"

int ncom;
double *pcom, *xicom, scale_up, scale_low, *store_p, *best_p;
extern Tree* s_tree;
extern Node* root;

void SHFT(double* a, double* b, double* c, double* d){

  *a=(*b);
  *b=(*c);
  *c=(*d);

}

void gradient(Node *nd, char** tipnames, int* states, int nb, int nbanno, double mu, char* model, double* p, double* g, int n) {
  int i,j;
  double pluslnl, oldlnl, minuslnl, tmp_p[n], sum=0.0;

  calc_lik_bfgs(root, tipnames, states, nb, nbanno, mu, model, p, &oldlnl);
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
    if(i<nbanno) {
      g[i] = (fabs(pluslnl) - fabs(oldlnl)) / (STEP);
    } else if (i==nbanno) {
      g[i] = (fabs(pluslnl) - fabs(oldlnl)) / (STEP2);
    } else if (i==nbanno+1) {
      g[i] = (fabs(pluslnl) - fabs(oldlnl)) / (STEP3);
    }
  }
  
}


double f1dim(Node *nd, char** tipnames, int* states, int nb, int nbanno, double mu, char* model, double x) {
  int j;
  double f,*xt,x_sum=0.0,tmp_x,sumfreq,e1,e2,exp1,exp2;
  
  xt=vector(0,ncom-1);
  sumfreq=0.0;
  x_sum=0.0;
  for(j=0;j<nbanno;j++) {
    xt[j]=(pcom[j]*10+x*xicom[j])/10;
    sumfreq+=fabs(xt[j]);
    //if(j==0) xt[j]=0.0;
    tmp_x = (1.0 / (1.0 + exp(xt[j])) )+ SIGMA;
    x_sum+=tmp_x;
    //xt[j]=tmp_x;
  }
  for(j=0;j<nbanno;j++) {
    //xt[j]=fabs(xt[j])/x_sum;
    //xt[j]=xt[j]/x_sum;
    xt[j] = fabs(xt[j]) / sumfreq;
  }
  if(xicom[nbanno]==0.0){
    exp1=1.0;
  } else {
    e1=log(fabs(xicom[nbanno]));
    exp1=pow(10.0,e1);
  }
  if(xicom[nbanno+1]==0.0){
    exp2=1.0;
  } else {
    e2=log(fabs(xicom[nbanno+1]));
    exp2=pow(10.0,e2);
  }
  xt[nbanno]=(pcom[nbanno]*exp1+x*xicom[nbanno])/exp1;
  xt[nbanno+1]=(pcom[nbanno+1]*exp2+x*xicom[nbanno+1])/exp2;
  if(xt[nbanno] >  scale_up || xt[nbanno] < scale_low) xt[nbanno] = pcom[nbanno];
  if(fabs(xt[nbanno+1]) >  1.0 || fabs(xt[nbanno+1]) < SCAL_MIN || xt[nbanno+1] < 0.0) xt[nbanno+1] = pcom[nbanno+1];
  
  calc_lik_bfgs(root, tipnames, states, nb, nbanno, mu, model, xt, &f);
  for(j=0;j<ncom;j++) {
    store_p[j] = xt[j];
  }
  free_vector(xt,0,ncom-1);
  return fabs(f);
}

void mnbrak(Node *nd, char** tipnames, int* states, int nb, int nbanno, double mu, char* model, double *ax, double *bx, double *cx, double *fa, double *fb, double *fc) {
/*Given a function func , and given distinct initial points ax and bx , this routine searches in
the downhill direction (defined by the function as evaluated at the initial points) and returns
new points ax , bx , cx that bracket a minimum of the function. Also returned are the function
values at the three points, fa , fb , and fc .*/
  double ulim,u,r,q,fu,dum=0.0,*adum,*bdum,*cdum,tmp;

  adum=&dum;
  *fa=f1dim(root, tipnames, states, nb, nbanno, mu, model, *ax);
  *fb=f1dim(root, tipnames, states, nb, nbanno, mu, model, *bx);
  if (*fb > *fa) {
    SHFT(adum,ax,bx,adum);
    SHFT(adum,fb,fa,adum);
  }
  *cx=(*bx)+GOLD*(*bx-*ax);
  *fc=f1dim(root, tipnames, states, nb, nbanno, mu, model, *cx);
  while (*fb > *fc) {
    r=((*bx)-(*ax))*((*fb)-(*fc));
    q=((*bx)-(*cx))*((*fb)-(*fa));
    u=(*bx)-(((*bx)-(*cx))*q-((*bx)-(*ax))*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
    ulim=(*bx)+GLIMIT*(*cx-*bx);
    if ((*bx-u)*(u-*cx) > 0.0) {
      fu=f1dim(root, tipnames, states, nb, nbanno, mu, model, u);
      if (fu < *fc) {
        *ax=(*bx);
        *bx=u;
        *fa=(*fb);
        *fb=fu;
        return;
      } else if (fu > *fb) {
        *cx=u;
        *fc=fu;
        return;
      }
      u=(*cx)+GOLD*((*cx)-(*bx));
      fu=f1dim(root, tipnames, states, nb, nbanno, mu, model, u);
    } else if ((*cx-u)*(u-ulim) > 0.0) {
        fu=f1dim(root, tipnames, states, nb, nbanno, mu, model, u);
        if (fu < *fc) {
          tmp=(*cx)+GOLD*((*cx)-(*bx));
           bdum=&tmp;
          SHFT(bx,cx,&u,&bdum);
          tmp=f1dim(root, tipnames, states, nb, nbanno, mu, model, u);
           cdum=&tmp;
          SHFT(fb,fc,&fu,&cdum);
        }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
        u=ulim;
        fu=f1dim(root, tipnames, states, nb, nbanno, mu, model, u);
    } else {
        u=(*cx)+GOLD*((*cx)-(*bx));
        fu=f1dim(root, tipnames, states, nb, nbanno, mu, model, u);
    }
    SHFT(ax,bx,cx,&u);
    SHFT(fa,fb,fc,&fu);
  }
}

double brent(Node *nd, char** tipnames, int* states, int nb, int nbanno, double mu, char* model, double ax, double bx, double cx, double tol, double *xmin) {
/*Given a function f , and given a bracketing triplet of abscissas ax , bx , cx (such that bx is
between ax and cx , and f(bx) is less than both f(ax) and f(cx) ), this routine isolates
the minimum to a fractional precision of about tol using Brent’s method. The abscissa of
the minimum is returned as xmin , and the minimum function value is returned as brent , the
returned function value.*/

  int iter, i;
  double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;

  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=f1dim(root, tipnames, states, nb, nbanno, mu, model, x);
  for (iter=0;iter<ITMAX;iter++) {
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+SIGMA);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin=x;
      return (-1.0)*fx;
    }
    if (fabs(e) > tol1) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabs(q);
      etemp=e;
      e=d;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)) {
        d=CGOLD*(e=(x >= xm ? a-x : b-x));
      } else {
        d=p/q;
        u=x+d;
        if (u-a < tol2 || b-u < tol2) d=SIGN(tol1,xm-x);
      }
    } else {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu=f1dim(root, tipnames, states, nb, nbanno, mu, model, u);
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      SHFT(&v,&w,&x,&u);
      SHFT(&fv,&fw,&fx,&fu);
      for(i=0; i<ncom;i++){
        best_p[i] = store_p[i];
      }
    } else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
        v=w;
        w=u;
        fv=fw;
        fw=fu;
      } else if (fu <= fv || v == x || v == w) {
        v=u;
        fv=fu;
      }
    }
  }
  nrerror("Warnng: Too many iterations >200 in Brent's method");
  *xmin=x;
  return (-1.0)*fx;
}

void linmin(Node *nd, char** tipnames, int* states, int nb, int nbanno, double mu, char* model, double* p, double* xi, int n, double* fret) {
/*Given an n -dimensional point p[1..n] and an n -dimensional direction xi[1..n] , moves and
resets p to where the function func(p) takes on a minimum along the direction xi from p ,
and replaces xi by the actual vector displacement that p was moved. Also returns as fret
the value of func at the returned location p . This is actually all accomplished by calling the
routines mnbrak and brent .*/
  int i,j;
  double xx,xmin,fx,fb,fa,bx,ax,*pold,tmp_x,x_sum,f_check;

  ncom=n;
  pcom=vector(0,n-1);
  pold=vector(0,n-1);
  xicom=vector(0,n-1);
  store_p=vector(0,n-1);
  best_p=vector(0,n-1);
  for (j=0;j<n;j++) {
    pcom[j]=p[j];
    xicom[j]=xi[j];
    pold[j]=p[j];
    best_p[j] = p[j];
  }
  ax=0.0;
  xx=1.0;
  mnbrak(root, tipnames, states, nb, nbanno, mu, model, &ax,&xx,&bx,&fa,&fx,&fb);
  *fret=brent(root, tipnames, states, nb, nbanno, mu, model, ax,xx,bx,TOL,&xmin);
  for(j=0;j<n;j++) {
    p[j]=best_p[j];
  }
  if(p[nbanno] >  scale_up || p[nbanno] < scale_low) p[nbanno] = pold[nbanno];
  if(p[nbanno+1] >  1.0 || p[nbanno+1] < SCAL_MIN) p[nbanno+1] = pold[nbanno+1];
  calc_lik_bfgs(root, tipnames, states, nb, nbanno, mu, model, p, &f_check);
  free_vector(pcom,0,n-1); free_vector(xicom,0,n-1); free_vector(pold,0,n-1); free_vector(store_p,0,n-1); free_vector(best_p,0,n-1);
}

void frprmn(Node *nd, char** tipnames, int* states, int nb, int nbanno, double mu, char* model, double* p, int n, double ftol, int* iter, double* fret, char **character) {
/*Given a starting point p[1..n] , Fletcher-Reeves-Polak-Ribiere minimization is performed on a
function func , using its gradient as calculated by a routine dfunc . The convergence tolerance
on the function value is input as ftol . Returned quantities are p (the location of the minimum),
iter (the number of iterations that were performed), and fret (the minimum value of the
function). The routine linmin is called to perform line minimizations.*/

  int i,j,its;
  double gg,gam,fp,dgg,fl;
  double *g, *h, *xi, sum;

  scale_up = 5.0 / s_tree->avgbl;
  scale_low = 0.05 / s_tree->avgbl;
  printf("Scalingfactor upbound = %lf, start point = %.5e\n\n***Fletcher-Reeves-Polak-Ribiere minimization ...\n\n",scale_up,p[nbanno]);
  g=vector(0,n-1);
  h=vector(0,n-1);
  xi=vector(0,n-1);
  calc_lik_bfgs(root, tipnames, states, nb, nbanno, mu, model, p, &fp);
  gradient(root, tipnames, states, nb, nbanno, mu, model, p, xi, n);
  for (j=0;j<n;j++) {
    g[j] = -xi[j];
    xi[j]=h[j]=g[j];
  }
  for (its=0;its<ITMAX_O;its++) {
    *iter=its;
    printf("Step%d, ",its);
    linmin(root, tipnames, states, nb, nbanno, mu, model, p,xi,n,fret);
    for (i=0;i<n;i++) {
      if(i < nbanno) printf("%s=%.6e, ",character[i],p[i]);
      if(i == nbanno) printf("Scaling=%.6e, ",p[i]);
      if(i == nbanno+1) printf("Epsilon=%.6e, ",p[i]);
    }
    calc_lik_bfgs(root, tipnames, states, nb, nbanno, mu, model, p, &fl);
    printf("lnL = %lf\n",fl);
    if (fabs(*fret-fp) <= ftol) {
      free_vector(g,0,n-1); free_vector(h,0,n-1); free_vector(xi,0,n-1);
      return;
    }
    fp= *fret;
    gradient(root, tipnames, states, nb, nbanno, mu, model, p, xi, n);
    dgg=gg=0.0;
    for (j=0;j<n;j++) {
      gg += g[j]*g[j];
      //dgg += xi[j]*xi[j]; /* Fletcher-Reeves */
      dgg += (xi[j]+g[j])*xi[j]; /* Polak-Ribiere */
    }
    if(g == 0) {
      free_vector(g,0,n-1); free_vector(h,0,n-1); free_vector(xi,0,n-1);
      return;
    }
    gam=dgg/gg;
    for(j=0; j<n; j++) {
      g[j] = -xi[j];
      xi[j] = h[j] = g[j]+gam*h[j];
    }
  }
  nrerror("parameter or likelihood is NOT converged within 1000 steps\n*** PASTML returns a possible best solution found during optimization ***\n");
  free_vector(g,0,n-1); free_vector(h,0,n-1); free_vector(xi,0,n-1);
}

