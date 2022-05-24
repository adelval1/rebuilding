# include <math.h>
# define min(x,y)  (((x) <= (y)) ? (x) : (y))
# define max(x,y)  (((x) >= (y)) ? (x) : (y))
#include "externalvariables.h"
extern double *vector(int);
extern void free_vector(double *);
extern void tridiag_inv(int,double *,double,double,double,double *,
                        double *,double *,double *,double *,double *,double *);


/*********************************************************
     subroutine energy
     Integrating energy equation
     Returns: h, dimensionless enthalpy
**********************************************************/

void energy(double *u,double *us,double *h_old,double *hp1,double **cs,double **css,double *V,double *rho,double *l0,double *l3,double *l3s,double **h,
            double **hs,double **J,double **Js,double *result,double xit2,double duedx,int jstart,double lam0,int total_therm_cond/*NANNI*/)
{
   extern int neq,nspec;
   extern double udelta,rdelta,dcsidelta,rhodelta,mudelta,hdelta,dhdelta,deta,beta;
   extern double hwall,K_bl;
   extern int axisymm,two_d,cone,flat_plate; /* Geometric informations */

   int i,j;
   double ud2,f1,g1,h1,*aa,*bb,*cc,*dd,*F,*dF,tmp1,tmp2,fact;
   double an,bn,cn,dn,pt,deta22,K_bl2;

   aa=vector(neq); bb=vector(neq); cc=vector(neq); dd=vector(neq);
   F=vector(neq); dF=vector(neq);
   
   ud2 = udelta*udelta;
   deta22 = deta*deta;
   K_bl2=K_bl*K_bl;
   
   if(jstart != 0)
   {
      fact=sqrt(xit2)*rdelta/dcsidelta;
   }
   else
   {
      if(axisymm)
         if(cone)
            fact=1.0;  /* It is included in the definition of Ji */
         else
            fact=1.0/sqrt(2.0*rhodelta*mudelta*duedx);

      if(two_d)
         if(flat_plate)
            fact=1.0;  /* It is included in the definition of Ji */
         else
            fact=1.0/sqrt(rhodelta*mudelta*duedx);
   }
   
   
   /*if (total_therm_cond==1)        /*NANNI*
      {
       for (i=0;i<neq;i++)        /*NANNI*  
        {
         for (j=0;j<nspec;j++)    /*NANNI* 
	  {
	  J[j][i]=0;              /*NANNI*
	  Js[j][i]=0;             /*NANNI*
	  }
	}  
      }*/
 
   for(i=0;i<neq;i++)
   {
/*    Here the equation is written in the form 
      an*dw/dcsi+bn*dw/deta=cn*d^2w/deta^2+dn+F   */
      an = xit2*u[i]*lam0;
      bn = V[i]-K_bl2*l3s[i];
      cn = K_bl2*l3[i];
      tmp1=0.0;
      tmp2=0.0;
      for(j=0;j<nspec;j++)
      {
         if (total_therm_cond==0) tmp1=tmp1+Js[j][i]*h[j][i]/hdelta+J[j][i]*hs[j][i]/hdelta;
         if (total_therm_cond==1) tmp1=0;    /*NANNI*/
	 tmp2=tmp2+l3s[i]*cs[j][i]*h[j][i]/hdelta
                  +l3[i]*css[j][i]*h[j][i]/hdelta
                  +l3[i]*cs[j][i]*hs[j][i]/hdelta;
      }
      tmp1=tmp1*K_bl;
      tmp2=tmp2*K_bl2;
      /*printf("\n fact*tmp1 and tmp2 in energy %le %le",fact*tmp1,tmp2);*/
      /*dn = -xit2*u[i]/hdelta*dhdelta*h_old[i]-beta*ud2/hdelta*rhodelta/rho[i]*u[i]
           +l0[i]*ud2/hdelta*us[i]*us[i]-fact*tmp1-tmp2;*/
      dn = -beta*ud2/hdelta*rhodelta/rho[i]*u[i]+l0[i]*ud2/hdelta*us[i]*us[i]*K_bl2-fact*tmp1-tmp2;
      /*dn = -beta*ud2/hdelta*rhodelta/rho[i]*u[i]
           +l0[i]*ud2/hdelta*us[i]*us[i]-tmp1-tmp2;*/ /* The diffusion fluxes are
                                                       already multiplied by fact */

/*    Here the equation is written in the form
      ai*d^2w/deta^2+bi*dw/deta+ci*w=d   */
      aa[i] = -cn/deta22; 
      bb[i] = bn/deta;
      /*cc[i] = an;*/
      cc[i] = an+xit2*u[i]/hdelta*dhdelta;
      dd[i] = dn-u[i]*hp1[i];
      F[i] = 0.0;
      dF[i] = 0.0;
   }
   f1 = 0.0;
   g1 = 1.0;
   h1 = hwall/hdelta;
 
   tridiag_inv(neq,h_old,f1,g1,h1,aa,bb,cc,dd,F,dF,result);

   free_vector(aa); free_vector(bb); free_vector(cc); free_vector(dd);
   free_vector(F); free_vector(dF);
}
