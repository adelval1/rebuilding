#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "peg_def.h"
#include "redefine.h"
#include "externalvariables.h"

double *vector(int);
double **matrix(int,int);
void free_vector(double *);
void free_matrix(double **,int);
void ableit(double *,double *,double,int);

void prod_and_dprod_rate(double,double *,double,double,double *,double *,double,double *,double **);
void diffusion_coeff(double **,double *,double *);
double number_density_(double *,double *);
void c_cp_frozen_(double *,double *,double *,int *,int *,int *,int *,int *,double *);
void species_enthalpy_(int *,double *,double *,int *,int *,int *,int *,int *,double *);
void species_cp_(int *,double *,double *,int *,int *,int *,int *,int *,double *);
void set_traco_(double *,double *,double *,double *,double *,double *,
                double *,int *);
double viscosity_();
double thermal_cond_frz_(double *);
/*double thermal_cond_tot_b_(double *);    /*NANNI*/
double c_mixture_enthalpy_(double *,double *,double *,int *,int *,int *,int *,int *,double *);
double c_mix_den_mass_(double *,double *,double *,int *);



void coefficients(double *l0,double ***Dij,double *l3,double **mapr,double ***dmapr,
                  double **h,double *rho,double *y,double *u,double *up1,double **c,
                  double **cs,double *t,double *M,double *Ms,double *Mw,
                  double xit2,double lam0,double twalldelta)
{
   int i,j,j1,l,nreact=0;
   double Pr,*mu,*k,cp,cp_m[8],cp_s[8],*cpr,*cpv,*cpe,*n_den,*rho_sp,*p_sp;
   double tmp1,tmp2,*ctmp,htmp[8],*h_spec,*k_tot/*NANNI*/,*k_frz/*NANNI*/;
   double *dum_mapr,**dum_dmapr,tneq[4],k_t/*NANNI*/;
   double *xtmp,*h_ord, h_nanni[nspec];
   double *Cpf;
 
   
 
   ctmp=vector(nspec);
   for (i=0;i<nspec;i++) h_nanni[i]=0.0;
   
   k_t=0.;
    
    h_spec=vector(nspec); 
   for (i=0;i<nspec;i++) h_spec[i]=h_nanni[i];

   h_ord=vector(nspec);
   dum_mapr=vector(nspec); dum_dmapr=matrix(nspec,nspec);
   cpr=vector(nspec); cpv=vector(nspec); cpe=vector(nspec);
   n_den=vector(nspec); rho_sp=vector(nspec); p_sp=vector(nspec);
   xtmp=vector(nspec);
   Cpf=vector(neq); k=vector(neq); mu=vector(neq);

   k_tot=vector(neq); k_frz=vector(neq);             /*NANNI*/
   
   /*printf("\n Starting coefficients \n"); /*NANNI*/

/* Computes the enthalpy at the wall */
   for(i=0;i<8;i++) htmp[i] = 0.0;
   for(i=0;i<4;i++) tneq[i] = twalldelta;
   for(i=0;i<nspec;i++){ctmp[i] = c[i][0];}
   
   /*printf("\n mixture enthalpy at the wall\n"); /*NANNI*/
   
     c_mixture_enthalpy_(&twalldelta,tneq,ctmp,&flg_anha,&mode,&flg_neq,&flg_stop,&flg_termo,htmp);
  
   hwall = htmp[0];


   /*printf("\n molecular weight\n"); /*NANNI*/ 

/* Computes the molecular weight */   
   for(i=0;i<neq;i++)
   {
      M[i]=0.0;
      for(j=0;j<nspec;j++)
      {
         M[i]=M[i]+c[j][i]/Mw[j];
      }
      M[i]=1.0/M[i];
   }

/* Computes the first derivative of molecular weight */
   for(i=0;i<neq;i++)
   {
      tmp1=0.0;
      for(j=0;j<nspec;j++)
      {
         tmp1=tmp1+cs[j][i]/Mw[j];
      }
      Ms[i]=-M[i]*M[i]*tmp1;
   }
   rhodelta = pdelta*M[neq-1]/t[neq-1]/R;
   
   
   for(i=0;i<neq;i++)
   {
      rho[i] = pdelta*M[i]/t[i]/R;
      for(j=0;j<nspec;j++)
      {
         ctmp[j]=c[j][i];
         xtmp[j]=ctmp[j]*M[i]/Mw[j];
         rho_sp[j] = rho[i]*ctmp[j];
         p_sp[j]=rho_sp[j]*R/Mw[j]*t[i];
         n_den[j]=number_density_(&p_sp[j],&t[i]);
      }
      for(j=0;j<4;j++) tneq[j]=t[i];
      for(j=0;j<8;j++) cp_m[j]=0.0;
      
      /*printf("\n call to cp_frozen Pegase routine mixture \n"); /*NANNI*/
      
      /* call to cp_frozen Pegase routine mixture */
      c_cp_frozen_(&t[i],tneq,ctmp,&flg_anha,&mode,&flg_neq,&flg_stop,&flg_termo,cp_m);
      cp=cp_m[0];
      Cpf[i]=cp;


      /*printf("\n before finding species cp \n"); /*NANNI*/
      
      /*printf("\n mixture cp in coeff %le \n",cp);*/
      for(j=0;j<nspec;j++)
      {
         int mode1;
         mode1 = 1;
         for(l=0;l<8;l++)
         {
            htmp[l]=0.0;
            cp_s[l]=0.0;
         }
         j1=j+1;
         /* call to species_enthalpy pegase routine */
         
	 /*printf("\n t[%d]= %le \n" ,i,t[i]); /*NANNI*/
	 species_enthalpy_(&j1,&t[i],tneq,&flg_anha,&mode,&flg_neq,&flg_stop,
                           &flg_termo,htmp);
         h[j][i] = htmp[0]+htmp[7];

         /*printf("\n h[][]= %le \n" , h[j][i]); /*NANNI*
	 getchar();/*NANNI*/




/***********************************************************************************
          call to routine used to compute species cp 
          Be careful all the cp' s are for unit mole and note unit mass as the 
            other quantities you compute 
*************************************************************************************/
         species_cp_(&j1,&t[i],tneq,&flg_anha,&mode1,&flg_neq,&flg_stop,
                     &flg_termo,cp_s);

         if(flg_anha==0)
         {
            cpr[j] = cp_s[2];
            cpv[j] = cp_s[3];
            cpe[j] = cp_s[4];
         }
         else
         {
            cpr[j] = cp_s[0];
            cpv[j] = 0.0;
            cpe[j] = 0.0;
         }
      }

      for(l=0;l<nspec;l++)
      {
         h_ord[l] = 0.0;
      }
      set_traco_(xtmp,n_den,cpr,cpv,cpe,h_ord,tneq,&oper_traco);


      /*printf("\n after set_traco subroutine \n"); /*NANNI*/
      
      mu[i] = viscosity_();
      
      /*if (total_therm_cond==0) 
         {
      	  k_frz[i] = thermal_cond_frz_(tneq);  /*NANNI*
          k_tot[i]= thermal_cond_tot_(tneq);  /*NANNI*
          k[i]=k_frz[i];
	 }
      
      if (total_therm_cond==1) 
         {
      	  k_tot[i] = thermal_cond_tot_(tneq);  /*NANNI*
          k_frz[i]= thermal_cond_frz_(tneq);  /*NANNI*
          k[i]=k_tot[i];
	 }
      
      /*------------- NEW IMPLEMENTATION FOR K tot -------------------*/
      
      
      if (total_therm_cond==0) 
         {
	k_frz[i] = thermal_cond_frz_(tneq);         /*NANNI*/
	k[i] = k_frz[i];                            /*NANNI*/ }
	 
	 
      if (total_therm_cond==1) 	  
        {
	 k_frz[i] = thermal_cond_frz_(tneq);         /*NANNI*/
	 printf(" change the flag ");
	 getchar();
/*	 k_tot[i] = thermal_cond_tot_b_(tneq);   /*NANNI*/
	 k[i] = k_tot[i];                            /*NANNI*/
        }
      
      /*---------------------------------------------------------------*/
      
      diffusion_coeff(Dij[i],tneq,&pdelta);
      
      Pr=mu[i]*cp/k[i];
      /*printf("\n Prandtl number in coeff %le ",Pr);*/
      l0[i]=rho[i]*mu[i]/(rhodelta*mudelta);
      /*printf(" l0 number in coeff %le ",l0[i]);*/
      l3[i]=l0[i]/Pr;
      for(l=0;l<nspec;l++) h_spec[l] = h[l][i];
      /*printf("\n h_spec = %le %le %le %le %le \n",h_spec[0],h_spec[1],h_spec[2],h_spec[3],h_spec[4]); /*NANNI*/ 
      
      /*Mass production rates are calculated  in the nonequilibrium 
      	case only, otherwise they are set to 0. NANNI */
      
      /*printf("\n before calculating the mass production rates \n"); /*NANNI*/
      
      /*if (oper==2) /*NANNI*
      {*/
       /*printf("\n prod_and_dprod_rate \n"); /*NANNI*/
       
       prod_and_dprod_rate(t[i],ctmp,rho[i],M[i],Mw,h_spec,cp,dum_mapr,dum_dmapr);
       
       for(j=0;j<nspec;j++)
        {
         /*printf("\n mapr \n"); /*NANNI*/
	 mapr[i][j]=dum_mapr[j];
            
         for(l=0;l<nspec;l++)
         {
               /*printf("\n dmapr \n"); /*NANNI*/
	       dmapr[i][j][l]=dum_dmapr[j][l];
         }
        }
     /*} 
       else /*NANNI*
        {
         for(j=0;j<nspec;j++)  /*NANNI*
          {
	  printf("\n before mapr \n"); /*NANNI* 
	  mapr[i][j]=1.0;  /*NANNI*
            
           for(l=0;l<nspec;l++)  /*NANNI*
           {
	       printf("\n before mapr \n"); /*NANNI*
               dmapr[i][j][l]=1.0; /*NANNI*
           }
          }
        }*/
     /*printf("\n after calculating the mass production rates \n"); /*NANNI*/ 
      
      
      y[i] = -(xit2*lam0+1.0)*u[i]-up1[i];
    
     
   }
   /*for(i=0;i<neq;i++) printf("\n k[%d]=%le   k_frz[%d]=%le   k_tot[%d]=%le",
                             i,k[i],i,k_frz[i],i,k_tot[i]);   /*NANNI*/
   
   
   cpwall=Cpf[0];
   kwall=k[0];
   muwall=mu[0];
   for(i=0;i<nspec;i++) hspwall[i] = h[i][0];
   
   
   free_vector(ctmp);    free_vector(xtmp); 
   free_vector(h_spec);    free_vector(h_ord);
   free_vector(dum_mapr); free_matrix(dum_dmapr,nspec);
   free_vector(cpr); free_vector(cpv);
   free_vector(cpe); free_vector(n_den); free_vector(rho_sp);
   free_vector(p_sp);
   free_vector(Cpf); free_vector(k); free_vector(mu);
   free_vector(k_tot); free_vector(k_frz);
}




void diffusion_coeff(double **Dij,double *tneq,double *p)
{
   int i,j;
   double *Dij_vec;
   void diff_coeff_f_(double *,double *,double *,int *);

   Dij_vec=vector(nspec*nspec);
   for(i=0;i<nspec*nspec;i++) Dij_vec[i] = 0.0;

   diff_coeff_f_(Dij_vec,tneq,p,&nspec);
   for(i=0;i<nspec;i++)
   {
      for(j=0;j<nspec;j++)
      {
         Dij[i][j]=Dij_vec[i*nspec+j];
      }
   } 
   free_vector(Dij_vec);
}


void der1(double *u,double *us,double **c,double **cs,double **css)
{
   int i;
   ableit(u,us,deta,neq);
   for(i=0;i<nspec;i++)
   {
      ableit(c[i],cs[i],deta,neq);
      ableit(cs[i],css[i],deta,neq);
   }
}


void der2(double *l0,double *l0s,double *l3,double *l3s,double **h,double **hs)
{
   int i;
   ableit(l0,l0s,deta,neq);
   ableit(l3,l3s,deta,neq);
   for(i=0;i<nspec;i++)
   {
      ableit(h[i],hs[i],deta,neq);
   }
}
