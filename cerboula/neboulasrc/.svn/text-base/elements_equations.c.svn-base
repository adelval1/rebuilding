#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define min(x,y)  (((x) <= (y)) ? (x) : (y))
#define max(x,y)  (((x) >= (y)) ? (x) : (y))
#include "peg_def.h"
#include "redefine.h"
#include "externalvariables.h"

extern void ableit(double *,double *,double,int);
extern double *vector(int);
extern double **matrix(int,int);
extern double ***trimatrix(int,int,int);
extern void free_vector(double *);
extern void free_matrix(double **,int);
extern void free_trimatrix(double ***,int,int);
void block_tridiag(double **,double **,double **,double **,double **,double ***,double *,
                   double *,double *,double **,double **,double **,int,int);


/***********************************************************************************
     subroutine elements_eq
     Integrating concentration equations in terms of nuclei
     Returns number of moles of the elements per unit volume
************************************************************************************/

void elements_eq(double **c,double **cs,double **css,double **cp1,double *u,double *V,
		 double **J,double **Js,double ***Dij,double *l0,double *l0s,double **mapr,
		 double ***dmapr,double *M,double *Ms,double *rho,double *Mw,double **result,
		 double xit2,double duedx,int jstart,double lam0,double *Mw_scaled,
		 double **csbis, double *t/*NANNI*/,double *var_add,int it)/*, double *Jel)*/
{
   /*extern int n_elements;*/

   extern int neq,nspec;
   extern double deta,udelta,rdelta,dcsidelta,rhodelta,mudelta,K_bl,pdelta;


   void c_massfrac_eq_(double *,double *,double *,double *,double *,int *,int *,
                        int *,int *,int *,int *,double *,int *,double *);
   void wall_elements(double *,double *,double *,double **, double *,
		      double *,double **,double **,double *, double *,
		      double *,double ***,double **,double **,
		      double ,double ,double ,int,int);/*,double *);*/

   void stef_max_sutton(double ,double *,double *,double *,double **,
                       double ,double ,double *,double *,double);
 
   void stef_max(double **,double **,double *,double ***,double *,double *,
              double *,double **,double **,double,double,int,
              void (*)(double ,double *,double *,double *,double **,
                       double ,double ,double *,double *,double) );   
    
   int i,j,k,i1;
   int flg_log,flg_guess,niter;
   double **aa,**bb,**cc,**dd,**F,***dF,***RHS,**res,tneq[4];
   double *f1,*g1,*h1,**dh1dc,*rhoiw,tmp_pie;

   /*   double fact,fact1,fact2,deta22,K_bl2;*/
   double fact1,deta22,K_bl2;
   double *y,*y_guess;

   double t_now, **cxi,**cxi1,**cxi1s, **dcxi, **rho_spec, **J_fake;
   double **J_fake_s, **J1, **J1s, **sum_Js, **Xc, *Xc_now;
   double *rho_s, **cxi_s, var_add_now, **c1, **c1s, *M1, *M1s, *rho1;
   double maxdiff_cxi,diff,**cxi_old;
   double *Mel,*Mel_scaled,max, eps;
   double *Xc_dernow;
   
   const double Pr=0.72,Le=1.2,meps=1.e-5,eps1=1.e-5;



   
   FILE *fp;
   
   /*- Vectors Allocation -*/ 
   rho1=vector(neq);

   f1=vector(n_elements); 
   g1=vector(n_elements); 
   h1=vector(n_elements); 
   Xc_now=vector(n_elements);
   rhoiw=vector(nspec);
   y=vector(nspec);
   y_guess=vector(nspec);
   M1=vector(neq);
   M1s=vector(neq);
   Xc_dernow=vector(n_elements);
   rho_s=vector(neq);
   Mel=vector(n_elements);
   Mel_scaled=vector(n_elements);
   
   /*- Matrix Allocation -*/

   dF=trimatrix(neq,n_elements,n_elements);
   RHS=trimatrix(neq,n_elements,n_elements);

   aa=matrix(neq,n_elements); 
   bb=matrix(neq,n_elements); 
   cc=matrix(neq,n_elements); 
   dd=matrix(neq,n_elements);
   F=matrix(neq,n_elements); 
   res=matrix(neq,n_elements);
   Xc=matrix(n_elements,neq);
   dh1dc=matrix(n_elements,n_elements);
   J1=matrix(nspec,neq);
   J1s=matrix(nspec,neq);
   cxi=matrix(n_elements,neq);
   cxi1=matrix(n_elements,neq);
   cxi1s=matrix(n_elements,neq);
   dcxi=matrix(n_elements,neq);
   c1=matrix(nspec,neq);
   c1s=matrix(nspec,neq);
   cxi_s=matrix(n_elements,neq);
   cxi_old=matrix(n_elements,neq);
   sum_Js=matrix(n_elements,neq);
   J_fake=matrix(n_elements,neq);
   J_fake_s=matrix(n_elements,neq);
   rho_spec=matrix(nspec,neq);


   flg_log = 0;
   flg_guess = 0;
   niter=100;
   eps=1.e-8;
   max=0.0; 
   
   deta22 = deta*deta;
   K_bl2=K_bl*K_bl;

/*A[rows][column]*/
/*Calculates the species partial densities & the actual values of cxi*/

   for(i=0;i<neq;i++)
     {
       for(j=0;j<nspec;j++)rho_spec[j][i]=rho[i]*c[j][i];
       for(k=0;k<n_elements;k++)
	 {
	   cxi[k][i]=0.0;
	   for(j=0;j<nspec;j++)cxi[k][i]= cxi[k][i]+lambda_matrx[j][k]*rho_spec[j][i]/Mw[j];
	 }
     }

/*Calculates the derivatives of cxi*/

   for(i=0;i<n_elements;i++)
     {
       /*cxi[i][neq-1];/**/
       ableit(cxi[i],cxi_s[i],deta,neq);
     }
/*Defines the fake additional term.....*/	    

  for(i=0;i<neq;i++)for(j=0;j<n_elements;j++)J_fake[j][i]=(1.0/rho[i])*(Le/Pr)*cxi_s[j][i]; 

/*..... and its first derivative*/

  for(i=0;i<n_elements;i++)ableit(J_fake[i],J_fake_s[i],deta,neq);

/*Calculates the derivative of rho */

  ableit(rho,rho_s,deta,neq);


         /*for(i=0;i<neq;i++)
         { 
            double dummy1;
            dummy1=0.0;
            for(j=0;j<nspec;j++)
            {
               dummy1=dummy1+Js[j][i];
               printf(" %le",Js[j][i]);
            }
            printf(" %le\n",dummy1);
         }		
	 getchar();*/





/*Evaluates the sum of J derivatives Js*/  
  
  for(i=0;i<neq;i++)
    for(k=0;k<n_elements;k++)
      {
	sum_Js[k][i]=0.0;
	for(j=0;j<nspec;j++)sum_Js[k][i]=sum_Js[k][i]+lambda_matrx[j][k]*Js[j][i]/Mw[j];
      }
 
  /*fact = 1.0/(2.0*duedx);*/
            
  fact1 = 1.0/sqrt(2.0*rhodelta*mudelta*duedx);
            
  /*fact2 = 1.0/sqrt(2.0*duedx/rhodelta/mudelta);*/
        
  for(i=0;i<n_elements;i++)for(j=0;j<n_elements;j++)dh1dc[i][j] = 0.0;
   
/***************************************************************************/      
/* 	COMPUTES THE COUPLING TERMS FOR THE NUCLEI CONSERVATION EQS (BEGIN)*/      
/***************************************************************************/      

/* Calculating the dcxi : dcxi[j][i] infinitesimal variation of the number of 
   moles of element j at the position i */

  /*for(i=0;i<neq;i++)for(j=0;j<n_elements;j++)dcxi[j][i]= meps*cxi[j][i];*/

 for(k=0;k<n_elements;k++)for(i=0;i<neq;i++)cxi1[k][i]=cxi[k][i];
  /*  --------------------------------------------------------- */
  /*  |   |   |   |   |   |   |   |   |   |   |   |   |   |   | */
  /*  --------------------------------------------------------- */

  /*      for(k=0;k<n_elements;k++)
          { 
	    cxi1[k][i]=0.0;
	    if(k==j) 
	      { 
		cxi1[k][i]=cxi[k][i]+dcxi[k][i];
	      }   
	    else if(k!=j) 
	      {
	        cxi1[k][i]=cxi[k][i];/*+dcxi[k][i]; /*!!!!  Ok  !!!!!!!
	      }
            Xc_dernow[k]=cxi1[k][i]*R*t[i]/pdelta;
	   }/* Calculating the derivative of the NUMBER fluxes wrt cxi */ 

  for(j=0;j<n_elements;j++)
  { 
   for(i=0;i<neq;i++) 
     { 
       
      dcxi[j][i]= meps*fabs(cxi[j][i]);
      cxi1[j][i]=cxi[j][i]+dcxi[j][i];
      if(fabs(cxi[j][i])<=1.0)
      {
         dcxi[j][i] = eps1;
         cxi1[j][i] = cxi[j][i]+eps1;
      }
      dcxi[j][i] = cxi1[j][i]-cxi[j][i];
      
      for(k=0;k<n_elements;k++)Xc_dernow[k]=cxi1[k][i]*R*t[i]/pdelta;
      
      /*printf("\n pdelta %le",pdelta);/**/ 
      for(k=0;k<4;k++)tneq[k]=t[i];

      t_now=t[i];

      /*printf("\n OLA %p %le \n ",&t[i],t[i]);*/

      c_massfrac_eq_(y,&pdelta,&t_now,tneq,&eps,&flg_log,&flg_anha,&flg_neq,&flg_stop,
                     &flg_termo,&flg_guess,y_guess,&niter,Xc_dernow);

      /*void c_massfrac_eq_(double *,double *,double *,double *,double *,int *,int *,
	                    int *,int *,int *,int *,double *,int *,double *);*/
      for (k=0;k<nspec;k++)c1[k][i]=y[k];      
     }
     
     for (k=0;k<n_elements;k++)ableit(cxi1[k],cxi1s[k],deta,neq);
      
     for (k=0;k<nspec;k++)ableit(c1[k],c1s[k],deta,neq); 
     
     /* Calculating the new value of M */
      
      for (i=0;i<neq;i++) 
	{
	  M1[i]=0.0;
	  for (k=0;k<nspec;k++)M1[i]= M1[i]+c1[k][i]/Mw[k];
	  M1[i]=1.0/M1[i];
	  rho1[i]=pdelta*M1[i]/(R*t[i]);	  
	}

      ableit(M1,M1s,deta,neq);
 
      /*ableit(rho1,rho1s,deta,neq);*/
       
     /*New diffusion fluxes with Stefan-Maxwell*/ 
      /*printf("\n xit2 = %le", xit2);*/
      stef_max(c1,c1s,rho1,Dij,Mw,M1,M1s,J1,J1s,duedx,xit2,jstart,stef_max_sutton);
      
     /* Computes the new eta derivatives of the number fluxes for the elements */
		    
     for (i=0;i<neq;i++)/* Over the cells */ 
     { 
     for(k=0;k<n_elements;k++) /* Over the elements */
      { 
       tmp_pie=0.0;
       /* New dNF/dcxi for element k thanks to the changement of cxi[j][i]=cxi[j][i]+dcxi[i]*/
       for(i1=0;i1<nspec;i1++)tmp_pie=tmp_pie+lambda_matrx[i1][k]*J1s[i1][i]/Mw[i1];
       RHS[i][k][j]=0.0;/*(fact1*K_bl)*(-tmp_pie+sum_Js[k][i])/dcxi[j][i];/**/
      } 
     } 
     for(i=0;i<neq;i++)cxi1[j][i]=cxi[j][i];
    }  
      
/***************************************************************************/      
/* 	COMPUTES THE COUPLING TERMS FOR THE NUCLEI CONSERVATION EQS (END)  */      
/***************************************************************************/             

   for(i=0;i<neq;i++)
   {
      for(j=0;j<n_elements;j++)
      {
/*    Here the equation is written in the form
      ai*d^2w/deta^2+bi*dw/deta+ci*w=d   */
/*       aa[i][j] = -K_bl2*(1.0/rho[i])*Le/Pr/deta22; 
         bb[i][j] = (1.0/rho[i]*V[i]+K_bl2*Le/Pr*rho_s[i]/(rho[i]*rho[i]))/deta;
         cc[i][j] = -V[i]*rho_s[i]/(rho[i]*rho[i]);
         dd[i][j] = -K_bl2*J_fake_s[j][i]/cxi_delta[j];*/

         aa[i][j] = -(1.0/rho[i])*Le/Pr/deta22; 
         bb[i][j] = (V[i]/rho[i]+Le/Pr*rho_s[i]/(rho[i]*rho[i]))/deta;
         cc[i][j] = -V[i]*rho_s[i]/(rho[i]*rho[i]);
         dd[i][j] = -J_fake_s[j][i];

	 F[i][j] = -fact1*K_bl*sum_Js[j][i]; 
	
	 for(k=0;k<n_elements;k++)dF[i][j][k] =  RHS[i][j][k];/**/
      }
   }

/******* Boundary condition using a prediction-correction method *********/
 
 wall_elements(f1,g1,h1,dh1dc,t,rho,c,cs,M,Ms,Mw,Dij,cxi,cxi_s,
	       fact1,xit2,duedx,jstart,it);/*,Jel);*/

 /*wall_cxi2(f1,g1,h1,dh1dc,t,rho,c,cs,M,Ms,Mw,Dij,cxi,cxi_s,
   fact2/rho[0],xit2,duedx,jstart,it,Jel);*/ 
   

/*************************************************************************/   
   
 for(j=0;j<n_elements;j++)for(i=0;i<nspec;i++)Mel[j]=Mel[j]+lambda_matrx[i][j]*Mw[i];
	     
 for(j=0;j<n_elements;j++)
    {  
       for(k=0;k<n_elements;k++)
	{
	  /*printf("\n Mel[%d]=%le",j,Mel[k]);*/
	if(max<Mel[k])  max=Mel[k];
	}    
       Mel_scaled[j]=Mel[j]/max;/**/
       /*Mel_scaled[j]=1.0;/**/
       /*printf("\n max=%le",Mel_scaled[j]);*/
      } 

 for(i=0;i<neq;i++)
      for(j=0;j<n_elements;j++)
      {
         aa[i][j]=aa[i][j]/Mel_scaled[j];
         bb[i][j]=bb[i][j]/Mel_scaled[j];
         cc[i][j]=cc[i][j]/Mel_scaled[j];
         dd[i][j]=dd[i][j]/Mel_scaled[j];
         F[i][j]=F[i][j]/Mel_scaled[j];
         for(k=0;k<n_elements;k++)dF[i][j][k]=dF[i][j][k]/Mel_scaled[j];
      }

   for(i=0;i<n_elements;i++)
   {
      f1[i]=f1[i]/Mel_scaled[i];
      g1[i]=g1[i]/Mel_scaled[i];
      h1[i]=h1[i]/Mel_scaled[i];
      for(j=0;j<n_elements;j++)dh1dc[i][j]=dh1dc[i][j]/Mel_scaled[i];
   }/**/
 
   maxdiff_cxi=0.0;
   diff=0.0;
  
   for (i=0;i<neq;i++)for(j=0;j<n_elements;j++)cxi_old[j][i]=cxi[j][i];
       
  
   
   block_tridiag(aa,bb,cc,dd,F,dF,f1,g1,h1,dh1dc,cxi,res,neq,n_elements);

   /* for(i=0;i<neq-1;i++)for(j=0;j<n_elements;j++)cxi[j][i]=urf*res[i][j]*cxi_delta[j]+(1.0-urf)*cxi_old[j][i];*/

   for(i=0;i<neq-1;i++)
     {
       /*printf("\n");*/
       for(j=0;j<n_elements;j++)
	 {
	   cxi[j][i]=res[i][j];
	   /*if(i==98)
	     {
	       printf("\n res[%d][%d] = %le ",i,j,res[i][j]);
	       printf("\n cxi[%d][%d] = %le ",i,j,cxi[j][i]);
	       }*/
	   /*printf(" res[%d][%d] = %le ",i,j,res[i][j]);
	     printf("\n cxi[%d][%d] = %le ",i,j,cxi[j][i]);*/
	 }
     }
   /*
   for(j=0;j<n_elements;j++)printf("\n res[99][%d] = %le ",j,res[99][j]);
   for(j=0;j<n_elements;j++)printf("\n cxi[99][%d] = %le ",j,cxi[j][99]);  
    printf("\n res[99][0] = %le ",res[99][0]);*/
   /* printf("\n cxi[%d][99] = %le ",j,cxi[j][99]);	   
 /* Saves the new values of cxi in cxi_sl*/

   /* for(i=0;i<neq;i++)for(j=0;j<n_elements;j++)cxi_sl[j][i]=cxi[j][i];*/


/* Writes the values of wall nuclei fractions*/


   for (i=0;i<neq;i++)
    for (k=0;k<n_elements;k++)
      {
       Xc[k][i]=cxi[k][i]*R*t[i]/pdelta;
      }    
   
   for (j=0;j<nspec;j++)y_guess[j]=0.0;/*c[j][0]*M[i]/Mw[j];/* THIS IS x_guess!!!!!!*/
      
   
   for (i=0;i<neq;i++)
    {       
      var_add_now=0.0;

      t_now=t[i];

      for (j=0;j<n_elements;j++) Xc_now[j]=Xc[j][i];

      for (j=0;j<4;j++) tneq[j]=t_now;

      c_massfrac_eq_(y,&pdelta,&t_now,tneq,&eps,&flg_log,&flg_anha,&flg_neq,&flg_stop,
                   &flg_termo,&flg_guess,y_guess,&niter,Xc_now);/*,&var_add_now*/

      for(j=0;j<nspec;j++)result[j][i] = y[j];
      for (j=0;j<nspec;j++)y_guess[j]=0.0;/*y[j]*M[i]/Mw[j]; /* THIS IS x_guess!!!!!!*/
      var_add[i]=var_add_now;
    }   

   /*- Deallocation -*/

   free_trimatrix(dF,neq,n_elements);  
   free_trimatrix(RHS,neq,n_elements);
 
   free_matrix(aa,neq);
   free_matrix(bb,neq);
   free_matrix(cc,neq);
   free_matrix(dd,neq);
   free_matrix(F,neq);
   free_matrix(res,neq);
   free_matrix(dh1dc,n_elements);
   free_matrix(cxi_old,n_elements);
   free_matrix(Xc,n_elements);
   free_matrix(J1,nspec);
   free_matrix(J1s,nspec);
   free_matrix(cxi,n_elements);
   free_matrix(dcxi,n_elements);
   free_matrix(cxi_s,n_elements);
   free_matrix(cxi1,n_elements);
   free_matrix(cxi1s,n_elements);
   free_matrix(c1,nspec);
   free_matrix(c1s,nspec);
   free_matrix(sum_Js,n_elements); 
   free_matrix(J_fake,n_elements);
   free_matrix(J_fake_s,n_elements); 
   free_matrix(rho_spec,nspec);

   free_vector(f1); 
   free_vector(g1); 
   free_vector(h1);
   free_vector(rhoiw);
   free_vector(Mel);
   free_vector(Mel_scaled);
   free_vector(Xc_now);
   free_vector(y);
   free_vector(y_guess);
   free_vector(rho_s);
   free_vector(Xc_dernow);
   free_vector(M1);
   free_vector(M1s);
   free_vector(rho1);
   
}

