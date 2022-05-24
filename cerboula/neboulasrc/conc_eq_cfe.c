#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "redefine.h"
#include "externalvariables.h"

double *vector(int);
void free_vector(double *);



void eq_conc_cfe(double **c,double *t)
{
#include "peg_def.h"
   int i,j,flg_log,flg_guess,niter;
   double Tneq[4],*y,*y_guess,eps,t_act;

   void c_massfrac_(double *,double *,double *,double *,double *,int *,int *,int *,
                    int *,int *,int *,double *,int *);

   flg_log = 0;
   flg_guess = 0;
   niter=100;/*before there wannns 10000*/
   eps=1.e-10;

   y=vector(nspec); y_guess=vector(nspec);
   
   t_act=0.;
   
   for(j=0;j<neq;j++)
     {
       for(i=0;i<nspec;i++) y[i]=0.0;
       for(i=0;i<nspec;i++) y_guess[i]=0.0;  
       for(i=0;i<4;i++) Tneq[i]=t[j];
  
       t_act=t[j];
       c_massfrac_(y,&pdelta,&t_act,Tneq,&eps,&flg_log,&flg_anha,&flg_neq,&flg_stop,
		   &flg_termo,&flg_guess,y_guess,&niter);
    
       for(i=0;i<nspec;i++) c[i][j] = y[i];
     }
       free_vector(y); free_vector(y_guess);
   }

