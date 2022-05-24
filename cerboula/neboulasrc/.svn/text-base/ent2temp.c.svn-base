#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "redefine.h"
#include "peg_def.h"
#include "externalvariables.h"

double *vector(int);
void free_vector(double *);

/**********************************************************************/
/***   Given the enthalpy computes the corresponding temperature    ***/
/**********************************************************************/

void ent2temp(double *t,double *h,double **c,double hdelta,double twall)
{
   int i,j;
   double t0,t_1,t_newt,htmp,*ctmp;
   void get_temp(double *,double,double *,double,double,double,int);

   ctmp = vector(nspec);
   for(i=neq-1;i>=1;i--)
   {
      for(j=0;j<nspec;j++) ctmp[j] = c[j][i];
      t0=t[i];
      t_1=t0+1.0;
      htmp=h[i]*hdelta;
      get_temp(ctmp,htmp,&t_newt,t0,t_1,macheps,newton_it);
      t[i]=t_newt;
   }
   t[0]=twall;
   free_vector(ctmp);
}



void get_temp(double *Y,double h,double *T_newt,double T0,double T_1,
              double macheps,int newton_it)
{
   int i,counter;
   double e_n[8],e_n_1[8],T_n1,T_n,T_n_1,Tneq[4];
   double C1,C2,f_n,f_n_1,steptol,fvectol;
   void c_mixture_energy_(double *,double *,double *,int *,int *,int *,
                          int *,int *,double *);
   double c_mixture_rspecific_(double *,int *,int *);

   steptol=pow(macheps,(2./3.));
   fvectol=pow(macheps,(1./3.));

   C2=c_mixture_rspecific_(Y,&mode,&flg_stop);
   C1=h/C2;

   for(i=0;i<4;i++) Tneq[i]=T0;
   T_n=T0;
   c_mixture_energy_(&T_n,Tneq,Y,&flg_anha,&mode,&flg_neq,&flg_stop,
                     &flg_termo,e_n);
   f_n = -T_n-e_n[0]/C2+C1;

   for(i=0;i<4;i++) Tneq[i]=T_1;
   T_n_1=T_1;
   c_mixture_energy_(&T_n_1,Tneq,Y,&flg_anha,&mode,&flg_neq,&flg_stop,
                     &flg_termo,e_n_1);
   f_n_1 = -T_n_1-e_n_1[0]/C2+C1;

   T_n1 = T_n - f_n*(T_n-T_n_1)/(f_n-f_n_1);

   if(fabs(f_n)<=(1.e-2*fvectol))
   {
      /*write(*,*) 'abs(fn)<=1.d-2*fvectol'*/
      *T_newt = T_n1;
      return;
   }

   counter=0.0;
   while(counter<=newton_it)
   {
      counter=counter+1;
      T_n_1 = T_n;
      T_n = T_n1;
      f_n_1 = f_n;
      c_mixture_energy_(&T_n,Tneq,Y,&flg_anha,&mode,&flg_neq,
                        &flg_stop,&flg_termo,e_n);
      f_n = -T_n-e_n[0]/C2+C1;

      T_n1 = T_n - f_n*(T_n-T_n_1)/(f_n-f_n_1);

      if(T_n1<=0.0)
      {
         printf("\n  Found negative temperature in get_temp %d %le \n",counter,T_n1);
         exit(1);
      }

      if(fabs(f_n)<=fvectol)
      {
         *T_newt = T_n1;
         /*write(*,*) 'abs(f)<=fvectol'
         write(*,*) counter*/
         return;
      }

      if((fabs(T_n1-T_n)/fabs(T_n1))<=steptol)
      {
         *T_newt = T_n1;
         /*write(*,*) 'abs(Tn1-Tn)/abs(Tn1)<=steptol'
         write(*,*) counter*/
         return;
      }

   }

   if(counter>newton_it)
   {
      printf("\n  Maximum number of iterations exceeded in get_temp \n");
      exit(1);
   }

}
