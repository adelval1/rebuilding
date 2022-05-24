#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "peg_def.h"
#include "redefine.h"
#include "externalvariables.h"

/* Boundary condition for the nuclei fraction based on a prediction-correction
    method */
    

void wall_elements(double *f1,double *g1,double *h1,double **dh1dc, double *t,
                  double *rho,double **c,double **cs,double *M, double *Ms,
		  double *Mw,double ***Dij,double **cxi,double **cxi_s,
		   double fact,double xit2,double duedx,int jstart, int it)/*,double *Jel)*/
{

 /*printf("\n Inside wall_cxi2 \n");*/
 
 void ableit(double *,double *,double,int);
 double *vector(int);
 double **matrix(int,int);
 void free_vector(double *);
 void free_matrix(double **,int);
 void stef_max_sutton(double ,double *,double *,double *,double **,
                       double ,double ,double *,double *,double);
 
 void stef_max(double **,double **,double *,double ***,double *,double *,
              double *,double **,double **,double,double,int,
              void (*)(double ,double *,double *,double *,double **,
                       double ,double ,double *,double *,double) );

 extern double K_bl;


 int i,j,k;/*,nit;*/
 int flg_log,flg_guess,niter;
 double tmp,t_now,tneq[4];
 double **x, *x1, *Di,*Di1, *Dj_tilde, *Dj_tilde1, **Jbc, **Jsbc, *Jw, **Dijw;
 double **dcxi, **cxi1, *rho1, **c1, **c1s, *M1, *M1s;
 double  *h1b, **rhoi, **cxi1s;
 double  *Jwb, *Xc_dernow;
 double	 *y, *y_guess, *Jelwall, ampli;

 const double Pr=0.72,Le=1.2;

 const double eps=1.e-8,meps=1.e-5,eps1=1.e-5;  /*eps is used to calculate the equilibrium
				       composition in c_massfrac_eq_!!!*/

 FILE *fp;

 x1=vector(nspec); 
 Di1=vector(nspec);
 Di=vector(nspec);
 Dj_tilde=vector(n_elements);
 Dj_tilde1=vector(n_elements);
 Jw=vector(nspec); 
 rho1=vector(neq); 
 M1=vector(neq);
 M1s=vector(neq);  
 h1b=vector(n_elements); 
 Jwb=vector(nspec);
 Xc_dernow=vector(n_elements); 
 y=vector(nspec); 
 y_guess=vector(nspec);
 Jelwall=vector(n_elements);


 Jbc=matrix(nspec,neq); 
 Jsbc=matrix(nspec,neq);
 Dijw=matrix(nspec,nspec); 
 c1=matrix(nspec,neq);
 c1s=matrix(nspec,neq); 
 cxi1=matrix(n_elements,neq); 
 rhoi=matrix(nspec,neq); 
 cxi1s=matrix(n_elements,neq);
 dcxi=matrix(n_elements,neq); 
 x=matrix(nspec,neq);


 flg_log = 0;
 flg_guess = 0;
 niter=100;
 ampli=1.;

 tmp=0.0;
/*Dj_tilde=0.0;*/

 /*nit=0;*/

 for(k=0;k<neq;k++)
   {
     for(i=0;i<nspec;i++)
       {
	 x[i][k]=c[i][k]*M[k]/Mw[i];
	 rhoi[i][k]=c[i][k]*rho[k];
       }
   }

 /*ableit(rho,rhos,deta,neq);*/


 /* Calculates the multicomponent diffusion coefficients at the wall */

 for(i=0;i<nspec;i++)
      {
	tmp=0.0;
	for(j=0;j<nspec;j++)
	  {
	    Dijw[i][j]=Dij[0][i][j];
	    if(j != i)tmp=tmp+x[j][0]/Dijw[i][j];	
	  }
	Di[i]=(1.0-x[i][0])/tmp;
      }

/*Calculating the Dj_tilde */ 
 
/* if(it<nit)
  for(j=0;j<n_elements;j++)for(i=0;i<nspec;i++)Dj_tilde[j]=Dj_tilde[j]+lambda_matrx[i][j]*Di[i]*10000.0; */ 
 
 /*if(it>=nit)*/
   for(j=0;j<n_elements;j++)
     {
       Dj_tilde[j]=0.0;
       for(i=0;i<nspec;i++)Dj_tilde[j]=Dj_tilde[j]+lambda_matrx[i][j]*Di[i]*ampli;/**/  
     }

/* Calculates f1 and g1 */

for(j=0;j<n_elements;j++)
         {
	   /*f1[j]=-Dj_tilde[j]/deta/rho[0]/Di[0];*/
	   f1[j]=Le/Pr/rho[0]/deta;
           g1[j]=0.0;/*K_bl*Dj_tilde[j]*rhos[0]/rho[0];*/
	 }
	 
/* Calculates the diffusion fluxes using the Stef-Max model */

stef_max(c,cs,rho,Dij,Mw,M,Ms,Jbc,Jsbc,duedx,xit2,jstart,stef_max_sutton);

for(i=0;i<nspec;i++)Jw[i]=Jbc[i][0];

/****************************************************************************/

/* Calculates h1 */
   
for(j=0;j<n_elements;j++)
    {
     tmp=0.0;
     for(i=0;i<nspec;i++)tmp=tmp+lambda_matrx[i][j]*Jw[i]/Mw[i];
      
      Jelwall[j]=tmp;

      /*      h1[j]=-tmp*fact*K_bl-Dj_tilde[j]*cxi_s[j][0]/rho[0]/Di[0];/*+K_bl*Dj_tilde[j]*rhos[0]/rho[0]*cxi[j][0];*/

      h1[j]=tmp*fact*K_bl+Le/Pr*cxi_s[j][0]/rho[0];

      /*h1[j]= -tmp*fact*K_bl-cxi_s[j][0]/rho[0];*/

      /* printf("tmp*fact : %le \n",tmp*fact);*/
    }


/* Writes the values of the wall diff. fluxes for the elements */

 if(fmod(it,50)==0.0)
   {
   if((fp=fopen("conv_Jwall_el.dat","a"))==NULL)
      {
       printf("Cannot open conv_Jwall_el.dat");
       exit(1);
      }
    else 
     { 
      fprintf(fp,"%d  %le   %le\n",it,Jelwall[0],Jelwall[1]);
     }
    fclose(fp);
   }

/* Calculating the derivative of h1[i] wrt cxi[j] */
/* for(i=0;i<nspec;i++) rhoiw1[i]=rhoiw[i];
   for(kder=0;kder<nspec;kder++)
    {
      dri[kder] = meps*fabs(rhoiw[kder]);
      rhoiw1[kder] = rhoiw[kder]+dri[kder];
      if(fabs(rhoiw[kder])<=1.0)
      {
         dri[kder] = eps;
         rhoiw1[kder] = rhoiw[kder]+eps;
      }
      dri[kder] = rhoiw1[kder]-rhoiw[kder];*/

/* Calculating the dcxi */

/*for(i=0;i<neq;i++)for(j=0;j<n_elements;j++)dcxi[j][i]= meps*fabs(cxi[j][i]);*/

/* Calculating the derivative of the  diff fluxes wrt cxi */ 

/*   for(k=0;k<n_elements;k++)
          {
            cxi1[k][i]=0.0;
            if(k==j) 
              {
        	cxi1[k][i]=cxi[k][i]+dcxi[k][i];/*+meps;/*
              }	   
            else   
              { 
		cxi1[k][i]=cxi[k][i];/*+dcxi[k][i]; /*!!!!  Ok  !!!!!!!
	      }
	  }*/

for(k=0;k<n_elements;k++)for(i=0;i<neq;i++)cxi1[k][i]=cxi[k][i];

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
      
      t_now=t[i];
     
      for(k=0;k<4;k++)tneq[k]=t_now;
      
      for(k=0;k<nspec;k++)
	{
	  y_guess[k]=0.0;
	  y[k]=0.0;
	}

      c_massfrac_eq_(y,&pdelta,&t_now,tneq,&eps,&flg_log,&flg_anha,&flg_neq,&flg_stop,
                   &flg_termo,&flg_guess,y_guess,&niter,Xc_dernow);
	
      for(k=0;k<nspec;k++)c1[k][i]=y[k];

      }
     
      for (k=0;k<n_elements;k++)ableit(cxi1[k],cxi1s[k],deta,neq);
      
      for (k=0;k<nspec;k++)ableit(c1[k],c1s[k],deta,neq);
     
     
     /* Calculating the new value of M */
      
      for (i=0;i<neq;i++) 
         {
          for (k=0;k<nspec;k++)M1[i]= M1[i]+c1[k][i]/Mw[k];
	  M1[i]=1.0/M1[i];
	  rho1[i]=pdelta*M1[i]/(R*t[i]);
	  }
      
      ableit(M1,M1s,deta,neq);
      
     /* Calculating the new value of rho  */
      
     /*ableit(rho1,rho1s,deta,neq);*/
      
     /*New diffusion fluxes with Stefan-Maxwell*/ 
      
      stef_max(c1,c1s,rho1,Dij,Mw,M1,M1s,Jbc,Jsbc,duedx,xit2,jstart,stef_max_sutton);
	
      for(k=0;k<nspec;k++)Jwb[k]=Jbc[k][0];
		    
      /*--------------------------------------------------------------------*/
      /*----- In this part you find the computation of the new Dj_tilde ----*/
      /*--------------------------------------------------------------------*/      
      
      /* --- Compute the mole fractions at the wall --- */
      for(k=0;k<nspec;k++) x1[k]=c1[k][0]*M1[0]/Mw[k];
     
      for(i=0;i<nspec;i++)
      {
       tmp=0.0;
       for(k=0;k<nspec;k++)if(k != i)tmp=tmp+x1[k]/Dijw[i][k] ;
       Di1[i]=(1.0-x1[i])/tmp ;
      }
 
      /*      if(it<nit)
      for(k=0;k<n_elements;k++)for(i=0;i<nspec;i++)Dj_tilde1[k] = Dj_tilde1[k]+lambda_matrx[i][k]*Di1[i];*/  
 
      /*if(it>=nit)*/
      for(k=0;k<n_elements;k++)
	{
	  Dj_tilde1[k]=0.0;
	  for(i=0;i<nspec;i++)Dj_tilde1[k] = Dj_tilde1[k]+lambda_matrx[i][k]*Di1[i]*ampli;/**/
	} 
      /*--------------------------------------------------------------------*/
     
      for(k=0;k<n_elements;k++)
      {
       tmp=0.0;
       for(i=0;i<nspec;i++)tmp=tmp+lambda_matrx[i][k]*Jwb[i]/Mw[i];      
       
       /* h1b[k]=-tmp*fact*K_bl-Dj_tilde1[k]*cxi1s[k][0]/rho1[0]/Di1[0];/*+K_bl*Dj_tilde[k]*rho1s[0]/rho1[0]*cxi1[k][0];*/
       
       h1b[k]=tmp*fact*K_bl+Le/Pr*cxi_s[k][0]/rho1[0];/*for CO2*/
       
       /* h1b[k]= -tmp*fact*K_bl-cxi1s[k][0]/rho1[0];/*FOR AIR cx1 maybe*/
       /*printf("\n  h1[%d], h1b[%d] : %le %le\n",k,k,h1[k],h1b[k]);*/
       dh1dc[k][j]=(h1b[k]-h1[k])/dcxi[j][0];
      }
 
      /* for (k=0;k<n_elements;k++)dh1dc[k][j]=(h1b[k]-h1[k])/dcxi[j][0]; /**/
      
      for(i=0;i<neq;i++)cxi1[j][i]=cxi[j][i];
  }

/*for(i=0;i<n_elements;i++)Jel[i]=Jelwall[i];*/

/* for(j=0;j<n_elements;j++)printf("\n f1[%d], g1[%d], h1[%d], h1b[%d] : %le %le %le %le\n",j,j,j,j,f1[j],g1[j],h1[j],h1b[j]);*/
 
/*for(j=0;j<n_elements;j++)printf("\n dcxi[%d][0] : %le \n",j,dcxi[j][0]);*/
 
/*for(j=0;j<n_elements;j++)for(i=0;i<n_elements;i++)printf("\n dh1dc[%d][%d] = %le : ",j,i,dh1dc[j][i]);*/
 /**/
 free_vector(x1);
 free_vector(Di);
 free_vector(Di1);
 free_vector(Jw);

 free_vector(rho1);
 free_vector(M1);free_vector(M1s);
 free_vector(h1b);


 free_vector(Jwb);
 
 free_vector(Dj_tilde);
 free_vector(Dj_tilde1);
 free_vector(y); 
 free_vector(y_guess); 
 free_vector(Jelwall); 
 free_vector(Xc_dernow);

 free_matrix(Jbc,nspec);
 free_matrix(Jsbc,nspec);
 free_matrix(Dijw,nspec);
 free_matrix(c1,nspec); 
 free_matrix(c1s,nspec);
 free_matrix(cxi1,n_elements); 

 free_matrix(rhoi,nspec);
 free_matrix(cxi1s,n_elements);


 free_matrix(dcxi,n_elements);
 free_matrix(x,nspec); 

}


