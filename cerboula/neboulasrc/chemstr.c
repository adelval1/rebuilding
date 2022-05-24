#include <math.h>
#include "peg_def.h"
# define max(x,y)  (((x) >= (y)) ? (x) : (y))
# define min(x,y)  (((x) <= (y)) ? (x) : (y))
#include "redefine.h"
#include "externalvariables.h"

void nerror(char *);
double *vector(int);
double **matrix(int,int);
void free_matrix(double **,int);
void free_vector(double *);
void f_chemistry_(double *,double *,double *,double *,double *,double *,double *,
                  int *,int *,int *,int *);


/******************************************************************
     Calculates the species mass production rate and its jacobian
*******************************************************************/

void prod_and_dprod_rate(double t,double *c,double rho,double M,double *Mw,
                         double *h_spec,double cp,double *mapr,double **dmapr)
{
   int i,j,k,l;
   double *rhosp,**drhospdc,*w,*dwdt,**dwdrhosp,*dtdc;
   double kfr,kbr,f_prod,b_prod,tmp,**temp_dmapr;
   double *dwdrp,*dwdrm,tneq[4],dummy;
   double delta(int,int);

   rhosp=vector(nspec); drhospdc=matrix(nspec,nspec); dtdc=vector(nspec);
   w=vector(nspec); dwdt=vector(nspec); dwdrhosp=matrix(nspec,nspec);
   dwdrp=vector(nspec*nspec); dwdrm=vector(nspec*nspec);
   temp_dmapr=matrix(nspec,nspec);
   
   for(i=0;i<4;i++) tneq[i]=t;
   for(i=0;i<nspec;i++) rhosp[i]=c[i]*rho;

   for(i=0;i<nspec;i++) dtdc[i] = -h_spec[i]/cp;
   for(l=0;l<nspec;l++)
   {
      for(k=0;k<nspec;k++)
      {
         /*drhospdc[l][k] = rho*delta(l,k)+rhosp[l]*(-M/Mw[k]+h_spec[k]/cp/t);*/

         /*Paolo*/
         drhospdc[l][k] = rho*delta(l,k)+rhosp[l]*(-M/Mw[k]);
      }
   }

/* Set to zero the mass production rate and its jacobian */
   for(i=0;i<nspec;i++)
   {
      w[i]=0.0;
      dwdt[i]=0.0;
      dwdrp[i]=0.0;     /*NANNI*/
      dwdrm[i]=0.0;     /*NANNI*/
      for(j=0;j<nspec;j++){dwdrhosp[i][j]=0.0;}
   }
   
/* Calls the Pegase routines for chemical reaction rates and jacobians */
   f_chemistry_(&t,tneq,rhosp,w,dwdrp,dwdrm,dwdt,&flg_anha,&flg_neq,&flg_stop,&flg_termo);

   for(i=0;i<nspec;i++)
   {
      for(j=0;j<nspec;j++)
      {
         dwdrhosp[i][j] = dwdrp[i*nspec+j]+dwdrm[i*nspec+j];
      }
   }
   for(i=0;i<nspec;i++)
   {
      mapr[i] = w[i];
      for(k=0;k<nspec;k++)
      {
         tmp = 0.0;
         for(l=0;l<nspec;l++)
         {
            tmp = tmp+dwdrhosp[i][l]*drhospdc[l][k];
         }
	 if (oper==2)
         {
            dmapr[i][k] = dwdt[i]*dtdc[k]+tmp;
         }
	 else
	 {
	    temp_dmapr[i][k]=0.0; /*NANNI*/ 
            dmapr=temp_dmapr;  /*NANNI*/
	 }
         
      }
   } 

   free_vector(rhosp); free_matrix(drhospdc,nspec); free_vector(dtdc);
   free_vector(w); free_vector(dwdt); free_matrix(dwdrhosp,nspec);
   free_vector(dwdrp); free_vector(dwdrm); free_matrix(temp_dmapr,nspec);
}








/*double mixt_enthalpy(double *t,double h,double *c)
{
   int i;
   double mixture_ent_(double *,double *,double *,int *,int *,int *,
                            int *,int *,double *);
   double c_mixture_enthalpy_(double *,double *,double *,int *,int *,int *,
                            int *,int *,double *);
   double mixh[8],tneq[4];

   for(i=0;i<8;i++) mixh[i] = 0.0;
   for(i=0;i<4;i++) tneq[i] = *t;

   c_mixture_enthalpy_(t,tneq,c,&flg_anha,&mode,&flg_neq,&flg_stop,&flg_termo,mixh);

   return (mixh[0]-h);
} */



double delta(int a,int b)
{
   if(a == b) return 1.0;
   else return 0.0;
}
