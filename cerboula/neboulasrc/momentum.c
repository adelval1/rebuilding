# define min(x,y)  (((x) <= (y)) ? (x) : (y))
# define max(x,y)  (((x) >= (y)) ? (x) : (y))
extern double *vector(int);
extern void free_vector(double *);
extern void tridiag_inv(int,double *,double,double,double,double *,
                        double *,double *,double *,double *,double *,double *);

/*********************************************************************
     subroutine momentum
     Integrating the momentum equation 
     Returns: f', dimensionless tangential velocity
*********************************************************************/

void momentum(double *u,double *up1,double *rho,double *l0,
              double *l0s,double *V,double *result,double xit2,double lam0)
{

   extern int neq;
   extern double beta,rhodelta,deta;
   extern double a_e,K_bl;

   int i;
   double f1,g1,h1,*a,*b,*c,*d,*F,*dF;
   double an,bn,cn,dn,fn,dfn,pt,deta22,K_bl2;

   a=vector(neq); b=vector(neq); c=vector(neq); d=vector(neq);
   F=vector(neq); dF=vector(neq);

   deta22 = deta*deta;
   K_bl2=K_bl*K_bl;


   for(i=0;i<neq;i++)
   {
/*    Here the equation is written in the form 
      an*dw/dcsi+bn*dw/deta=cn*d^2w/deta^2+dn+F   */
      an = xit2*u[i]*lam0;
      bn = V[i]-K_bl2*l0s[i];
      cn = K_bl2*l0[i];
      dn = beta*rhodelta/rho[i]+0.5*rhodelta/rho[i]*a_e;
      fn = -beta*u[i]*u[i];
      dfn = -2.0*beta*u[i];

/*    Here the equation is written in the form
      ai*d^2w/deta^2+bi*dw/deta+ci*w=d   */
      a[i] = -cn/deta22;
      b[i] = bn/deta;
      /*c[i] = an-dfn;*/
      c[i] = an;
      /*d[i] = dn-u[i]*up1[i] +fn-dfn*u[i];*/
      d[i] = dn-u[i]*up1[i];
      F[i] = fn;
      dF[i] = dfn;
   }
      
   f1 = 0.0;
   g1 = 1.0;
   h1 = 0.0;
      
   tridiag_inv(neq,u,f1,g1,h1,a,b,c,d,F,dF,result);

   free_vector(a); free_vector(b); free_vector(c); free_vector(d);
   free_vector(F); free_vector(dF);
}
