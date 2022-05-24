#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define WL 80
#include "redefine.h"
#include "externalvariables.h"

double *vector(int);
void free_vector(double *);
double **matrix(int,int);
void free_matrix(double **,int);
# define min(x,y)  (((x) <= (y)) ? (x) : (y))


void wall_flux(double *Jw,double *rhoiw,double *Mw)
{
   int i,k,r,iter;
   double *M_inc,tmp1,tmp2,tmp3;
   double ga_1_O,ga_2_O,ga_CO,ga_temp;/*,*JJJw; /*pietro*/
   double M_maximpinging;
   char dummystr[WL];
   FILE *fp;
   
   M_inc=vector(nspec);
   
   for(i=0;i<nspec;i++){
      M_inc[i]=-rhoiw[i]*sqrt(R*twalldelta/(2.0*pi*Mw[i]));
   }

   if(choice==1)
   {
   	/*FIRST PART OF THE IF CONDITION*/   
   	for(i=0;i<nspec;i++)
   	{
      	M_inc[i]=-(M_inc[i]/Mw[i]);
      	/*printf("\n %le %le %le",Mw[i],rhoiw[i],M_inc[i]);*/
   	}   
   	/****************************************************/
   	/*                  3 WALL REACTIONS                */
   	/****************************************************/
   
   	/****************************************************/
   	/* Recombination probability for each wall reactions*/
   	/****************************************************/
   
   	/*ga_CO=(2.*gam[1]/(2.+gam[1]));
   	ga_3_C=(2.*gam[1]/(2.+gam[1]));
   	ga_3_O=ga_3_C*0.5*M_inc[4]/M_inc[3];
   	ga_2_O=ga_CO*M_inc[2]/M_inc[3];
   	ga_1_O=(2.*gam[1]/(2.+gam[1]))-ga_2_O-ga_3_O;
   	/****************************************************/
   	/* Recombination probability for each wall reactions*/
   	/****************************************************/ 
   	/* (0) CO2 (1) O2 (2) CO (3) O  (6) C */
   	/*Jw[0]=(gam[1]/ga_CO)*Mw[0]*M_inc[3]*(ga_2_O+2*ga_3_O);
   	Jw[1]=(gam[1]/ga_CO)*ga_1_O*M_inc[3]*Mw[1]*0.5;
   	Jw[2]=-ga_2_O*(gam[1]/ga_CO)*M_inc[3]*Mw[2];
   	Jw[3]=-(2.*ga_1_O+ga_2_O+2.*2.*ga_3_O)*(gam[1]/ga_CO)*M_inc[3]*Mw[3];
   	Jw[4]=-2.*ga_3_O*M_inc[3];*/   
   	/****************************************************/
   	/*                  2 WALL REACTIONS                */
   	/****************************************************/
   	/*if(nrwall==2)*/
   	ga_CO=(2.*gamdelta[1]/(2.+gamdelta[1])); /*pietro*/   
   	ga_2_O=(2.*gamdelta[1]/(2.+gamdelta[1]))*(M_inc[2]/M_inc[3]); /*pietro*/
   	ga_1_O=(2.*gamdelta[1]/(2.+gamdelta[1]))-ga_2_O; /*pietro*/
   
   	/*if(ga_1_O<0)
   	{
    	Deltagamma=1.001;
    	ga_temp=gam[1]; 
    	iter=0;
    	while(ga_1_O<0. && ga_temp<1.)
     	{
      	ga_temp=ga_temp*Deltagamma;
      	ga_CO=(2.*ga_temp/(2.+ga_temp));    
      	ga_2_O=(2.*ga_temp/(2.+ga_temp))*(M_inc[2]/M_inc[3]); 
      	ga_1_O=(2.*ga_temp/(2.+ga_temp))-ga_2_O;    
      	iter++;
     	}
    	printf("\n gamma_temp %le ga_1_O %le ga_2_0 %le ",ga_temp,ga_1_O,ga_2_O);
   	if(ga_1_O>=0)gam[1]=ga_temp;  
   	if(ga_1_O>=0)gam[2]=ga_temp;
   	}

   	if(ga_temp>=1.)
   	{
   	ga_CO=(2.*gam[1]/(2.+gam[1]));    
   	ga_2_O=(2.*gam[1]/(2.+gam[1]))*(M_inc[2]/M_inc[3]); 
   	ga_1_O=(2.*gam[1]/(2.+gam[1]))-ga_2_O; 
   	}
    
   	/*pietro*/ 
   	/*printf("\n ga_1_O, ga_2_O, ga_CO \n");
   	printf(" %le %le %le ",ga_1_O*1.,ga_2_O*1.,ga_CO*1.);
   	/*printf("\n ga_1_O, ga_2_O, ga_CO \n");*/
   	if(ga_1_O<0) {printf("\n CIAO BRODO");
//	getchar();
	}/*pietro*/ 
   	/*for(i=0;i<nrwall;i++)
   	{
   	/*printf("\n gamma(i), i, nrwall \n");
   	printf(" %le %le %le ",gam[i],i*1.0,nrwall*1.0);/*pietro
   	}*/
   	/*ga_choice=1;*/
   
   	if(ga_choice==1)
   	{
    		if(ga_2_O+ga_1_O!=0.)
    		{
     		Jw[0]=ga_2_O*(gamdelta[1]/(ga_2_O+ga_1_O))*M_inc[3]*Mw[0];
     		Jw[1]=ga_1_O*(gamdelta[1]/(ga_2_O+ga_1_O))*M_inc[3]*Mw[1]*0.5;
     		Jw[2]=-ga_2_O*(gamdelta[1]/(ga_2_O+ga_1_O))*M_inc[3]*Mw[2];
     		Jw[3]=-(ga_2_O+ga_1_O)*(gamdelta[1]/(ga_2_O+ga_1_O))*M_inc[3]*Mw[3];
     		Jw[4]=0.0;   
    		}
    		else
    		{for(i=1;i<nspec;i++)Jw[i]=0.0;}
   	}
   	else
   	{
   	Jw[0]=ga_2_O*M_inc[3]*Mw[0];
   	Jw[1]=ga_1_O*M_inc[3]*Mw[1]*0.5;
   	Jw[2]=-ga_2_O*M_inc[3]*Mw[2];
   	Jw[3]=-(ga_2_O+ga_1_O)*M_inc[3]*Mw[3];
   	Jw[4]=0.0;   
   	}
   
//  	printf("\n These are the fluxes .... and the sum\n");
/*  	for(i=0;i<nspec;i++)
  	{
  	printf(" %e",Jw[i]/Mw[i]);
  	}
  */
  	/*
  	printf("\n Sum NEW %le",Jw[0]+Jw[1]+Jw[2]+Jw[3]+Jw[4]);/**/
     	/*printf(" \n Mw %le %le %le %le %le  \n",Mw[0],Mw[1],Mw[2],Mw[3],Mw[4]);*/
     	/*printf("\n Jw wall %le %le %le %le %le \n",Jw[0],Jw[1],Jw[2],Jw[3],Jw[4]);/**/ 
   	/*printf("\n C nuclei NEW %le ", Jw[0]/Mw[0]+Jw[2]/Mw[2]+Jw[4]/Mw[4]);
   	printf("\n O nuclei NEW %le ", 2*Jw[0]/Mw[0]+Jw[2]/Mw[2]+2*Jw[1]/Mw[1]+Jw[3]/Mw[3]);
   	/*-----------------------------------------------------------*/
   	/*      END MODIFICATIONS BY PIETRO FOR NEW CATALYCITY       */
   	/*-----------------------------------------------------------*/  
   }
   if(choice==2)
   {
//    for(r=0;r<nrwall;r++){printf("\n gamdelta %f ",gamdelta[r]);};
//    printf("\n");
	
      for(i=0;i<nspec;i++)
      {
         tmp1=0.0;
         for(r=0;r<nrwall;r++){tmp1=tmp1+nu[r][i]*gamdelta[r];}
         tmp2=0.0;
         for(r=0;r<nrwall;r++)
         {
            tmp3=0.0;
            for(k=0;k<nspec;k++){tmp3=tmp3+mu[r][i][k]*gamdelta[r]*M_inc[k];}
            tmp2=tmp2+tmp3;
         }
         Jw[i] = tmp1*M_inc[i] - tmp2;
//		printf("\n jw %d %e",i,Jw[i]);
      }
//		printf("\n");

   }
   /* END SECOND PART OF THE IF CONDITION */



if(choice==100) /*this is a generic bc for customization in the code right here*/{

	double gasconcentration[nspec],adsorb[nspec+1],wallmoleproduction[5],catalycities[nspec],cata;
	int jj,cattype=7;
	for(jj=0;jj<nspec;jj++){
		gasconcentration[jj]=rhoiw[jj]/Mw[jj];	
	
        }
        cata=gamdelta[0];

        catc_zgbplusl_(&nspec, catalycities,adsorb, &twalldelta,gasconcentration,Mw, wallmoleproduction);
//	catc_catalycity_(&cattype,&nspec,&cata, gasconcentration,Mw, &twalldelta, wallmoleproduction);	
	for(jj=0;jj<nspec;jj++){
		Jw[jj]=wallmoleproduction[jj]*Mw[jj];	
        }
/*
	for(jj=0;jj<nspec;jj++){
			printf("wmp%e \n",wallmoleproduction[jj]);
	}
	*/
//	getchar();
//	gamdelta[0]=-Jw[0]/Mw[0]/M_inc[2]*Mw[2];
//	gamdelta[1]=-(Jw[0]/Mw[0]+Jw[1]/Mw[1])/(M_inc[1]/Mw[1]+M_inc[3]/Mw[3]);


	
	if((fp = fopen("neboulaoutput/surfacecoverage.dat","w"))==NULL){
			printf("Cannot open file %s\n","neboulaoutput/surfacecoverage.dat");exit(1);
	}
	else{
		fprintf(fp,"surface species\n");
		fprintf(fp,"%e %e %e %e %e %e\n",adsorb[0],adsorb[1],adsorb[2],adsorb[3],adsorb[4],adsorb[5]);
		for(jj=0;jj<nspec;jj++){
			fprintf(fp,"gamma %d \n",jj);
			fprintf(fp,"%e\n",catalycities[jj]);
		}
	}
	fclose(fp);

   }
   free_vector(M_inc);
}


/* Paolo. It is the new boundary condition */

void wall_bc(double *f1,double *g1,double *h1,double **dh1dc,double rhow,double *cw,
             double M,double Msw,double *Mw,double *Jiw,double *Jbw,
             double l0,double fact)
{
   void wall_flux(double *,double *,double *);
   double delta(int,int);
   extern double K_bl,rhodelta,mudelta;

   int i,j,k,l,kk,kder;
   double *Jw,*rhoiw,*rhoiw1,*h1b,*Di,*x,tmp,tmp1,tmp2;
   double **drdc,*dri,**dh1dr,rhow1,M1,*cw1;
   const double Pr=0.72,Le=1.2;
   const double eps=1.e-5,meps=1.e-5;

   Jw=vector(nspec); rhoiw=vector(nspec); rhoiw1=vector(nspec); h1b=vector(nspec);
   Di=vector(nspec); x=vector(nspec); drdc=matrix(nspec,nspec); dri=vector(nspec);
   dh1dr=matrix(nspec,nspec); cw1=vector(nspec);

   for(i=0;i<nspec;i++)
   {
      rhoiw[i] = rhow*cw[i];
      x[i] = cw[i]*M/Mw[i];
   }
   tmp1=0.0;
   tmp2=0.0;
   for(i=0;i<nspec;i++)
   {
      tmp1=tmp1+x[i];
      tmp2=tmp2+cw[i];
   }

   wall_flux(Jw,rhoiw,Mw);

   for(i=0;i<nspec;i++) h1[i]=Jw[i]*fact;


   for(i=0;i<nspec;i++)
   {
      f1[i] = -l0*Le/Pr/deta*K_bl;
      g1[i] = 0.0;
   }


   /* Computes derivatives of partial densities with respec to mass fractions */

   for(i=0;i<nspec;i++)
   {
      for(j=0;j<nspec;j++)
      {
         drdc[i][j] = rhow*delta(i,j)+rhoiw[i]*(-M/Mw[j]);
      }
   }

   /* Computes derivatives of h1 with respect to partial densities */

   for(i=0;i<nspec;i++) rhoiw1[i]=rhoiw[i];

   for(kder=0;kder<nspec;kder++)
   {
      dri[kder] = meps*fabs(rhoiw[kder]);
      rhoiw1[kder] = rhoiw[kder]+dri[kder];
      if(fabs(rhoiw[kder])<=1.0)
      {
         dri[kder] = eps;
         rhoiw1[kder] = rhoiw[kder]+eps;
      }
      dri[kder] = rhoiw1[kder]-rhoiw[kder];

   /* Computes new density  */
      rhow1=0.0;
      for(i=0;i<nspec;i++) rhow1 = rhow1+rhoiw1[i];
      for(i=0;i<nspec;i++) cw1[i] = rhoiw1[i]/rhow1;
   /* Computes the new wall flux */   
      wall_flux(h1b,rhoiw1,Mw);
      for(i=0;i<nspec;i++) h1b[i]=h1b[i]*fact;

      for(kk=0;kk<nspec;kk++) dh1dr[kk][kder] = (h1b[kk]-h1[kk])/dri[kder];
      rhoiw1[kder] = rhoiw[kder];
   }

   for(i=0;i<nspec;i++)
   {
      for(k=0;k<nspec;k++)
      {
         tmp = 0.0;
         for(l=0;l<nspec;l++)
         {
            tmp = tmp+dh1dr[i][l]*drdc[l][k];
         }
         dh1dc[i][k] = tmp/(rhodelta*mudelta);
      }
   }
   for(i=0;i<nspec;i++) h1[i]=h1[i]/(rhodelta*mudelta)-Jiw[i]*fact/(rhodelta*mudelta)
                              -Jbw[i];

   free_vector(Jw); free_vector(rhoiw); free_vector(rhoiw1); free_vector(h1b);
   free_vector(Di); free_vector(x); free_matrix(drdc,nspec); free_vector(dri);
   free_matrix(dh1dr,nspec); free_vector(cw1);
}


void wall_stef(double *f1,double *g1,double *h1,double **dh1dc,double rhow,double *cw,
               double M,double Ms,double *Mw,double **Dij,double den)

                             
{
   void wall_flux(double *,double *,double *);
   double delta(int,int);
   extern double K_bl;

   int i,j,k,l,kder,kk;
   double M1,*Jw,fact,*rhoiw,*rhoiw1,*h1b,*dri,**dh1dr,**drdc,rhow1;
   const double eps=1.e-5,meps=1.e-5;

   Jw=vector(nspec); rhoiw=vector(nspec); rhoiw1=vector(nspec); h1b=vector(nspec);
   dri=vector(nspec); dh1dr=matrix(nspec,nspec); drdc=matrix(nspec,nspec);

   /* Computes derivatives of partial densities with respec to mass fractions */

   for(i=0;i<nspec;i++)
   {
      rhoiw[i] = rhow*cw[i];
      for(j=0;j<nspec;j++)
      {
//         drdc[i][j] = rhow*delta(i,j)+rhoiw[i]*(-M/Mw[j]+hspwall[j]/cpwall/twall);
       drdc[i][j]=rhow*delta(i,j)+rhoiw[i]*(-M/Mw[j]+hspwall[j]/cpwall/twalldelta);
      }
   }
   for(i=0;i<nspec;i++)
   {
      f1[i]=M/Mw[i]/deta*K_bl;
      g1[i]=Ms/Mw[i]*K_bl;
   }
   wall_flux(Jw,rhoiw,Mw);
   /*printf("\n 1 \n ");/*pietro*/
   /*printf(" %le %le %le %le %le \n",Jw[0],Jw[1],Jw[2],Jw[3],Jw[4]);/*pietro*/

   for(i=0;i<nspec;i++)
   {
      fact=0.0;
      for(j=0;j<nspec;j++)
      {
         fact=fact+cw[i]*M/Mw[i]*cw[j]*M/Mw[j]/Dij[i][j]*(Jw[j]/rhoiw[j]-Jw[i]/rhoiw[i]);
      }
      h1[i] = fact*den;
   }

   /* Computes derivatives of h1 with respect to partial densities */

   for(i=0;i<nspec;i++) rhoiw1[i]=rhoiw[i];

   for(kder=0;kder<nspec;kder++)
   {
      dri[kder] = meps*fabs(rhoiw[kder]);
      rhoiw1[kder] = rhoiw[kder]+dri[kder];
      if(fabs(rhoiw[kder])<=1.0)
      {
         dri[kder] = eps;
         rhoiw1[kder] = rhoiw[kder]+eps;
      }
      dri[kder] = rhoiw1[kder]-rhoiw[kder];

   /* Computes new density and new molecular weight */
      rhow1=0.0;
      for(i=0;i<nspec;i++) rhow1 = rhow1+rhoiw1[i];
      for(i=0;i<nspec;i++) cw[i] = rhoiw1[i]/rhow1;
      M1=0.0;
      for(i=0;i<nspec;i++) M1 = M1+cw[i]/Mw[i];
      M1=1.0/M1;
   /* Computes the new wall flux */   
      wall_flux(Jw,rhoiw1,Mw);
      /*printf("\n 2 \n ");/*pietro*/
      /*printf(" %le %le %le %le %le ",Jw[0],Jw[1],Jw[2],Jw[3],Jw[4]);/*pietro*/

   /* Computes the new h1 */
      for(i=0;i<nspec;i++)
      {
         fact = 0.0;
         for(j=0;j<nspec;j++)
         {
               fact=fact+cw[i]*cw[j]*M1*M1/Mw[i]/Mw[j]/Dij[i][j]*
                         (Jw[j]/rhoiw1[j]-Jw[i]/rhoiw1[i]);
         }
         h1b[i]=fact*den;
      }

      for(kk=0;kk<nspec;kk++) dh1dr[kk][kder] = (h1b[kk]-h1[kk])/dri[kder];

      rhoiw1[kder] = rhoiw[kder];

   }

   /* Computes the derivative of h1 with respect to mass fractions */

   for(i=0;i<nspec;i++)
   {
      for(k=0;k<nspec;k++)
      {
         fact = 0.0;
         for(l=0;l<nspec;l++)
         {
            fact = fact+dh1dr[i][l]*drdc[l][k];
         }
         dh1dc[i][k] = fact;
      }
   }


   free_vector(Jw); free_vector(rhoiw); free_vector(rhoiw1); free_vector(h1b);
   free_vector(dri); free_matrix(dh1dr,nspec); free_matrix(drdc,nspec);

}


void wall_ci(double *h1)
{
#include "peg_def.h"
   int i,flg_log,flg_guess,niter;
   double Tneq[4],*y,*y_guess,eps;
   void c_massfrac_(double *,double *,double *,double *,double *,int *,int *,int *,
                    int *,int *,int *,double *,int *);

   flg_log = 0;
   flg_guess = 0;
   niter=1000;
   eps=1.e-10;

   y=vector(nspec); y_guess=vector(nspec);
   for(i=0;i<nspec;i++) y[i]=0.0;
   for(i=0;i<nspec;i++) y_guess[i]=0.0;

   for(i=0;i<4;i++) Tneq[i]=twalldelta;

   c_massfrac_(y,&pdelta,&twalldelta,Tneq,&eps,&flg_log,&flg_anha,&flg_neq,&flg_stop,
               &flg_termo,&flg_guess,y_guess,&niter);

   for(i=0;i<nspec;i++) h1[i] = y[i];

   free_vector(y); free_vector(y_guess);

}

void wall_ram(double *f1,double *g1,double *h1,double **dh1dc,double rhow,double *cw,
              double *csw,double M,double Msw,double *Mw,double **Dij,double fact)
{
   void wall_flux(double *,double *,double *);
   double delta(int,int);
   extern double K_bl;

   int i,j,k,l,kk,kder;
   double *Jw,*rhoiw,*rhoiw1,*h1b,*Di,*x,tmp,tmp1,tmp2;
   double **drdc,*dri,**dh1dr,rhow1,M1,*cw1;
   const double eps=1.e-5,meps=1.e-5;

   Jw=vector(nspec); rhoiw=vector(nspec); rhoiw1=vector(nspec); h1b=vector(nspec);
   Di=vector(nspec); x=vector(nspec); drdc=matrix(nspec,nspec); dri=vector(nspec);
   dh1dr=matrix(nspec,nspec); cw1=vector(nspec);

   for(i=0;i<nspec;i++)
   {
      rhoiw[i] = rhow*cw[i];
      x[i] = cw[i]*M/Mw[i];
   }
   tmp1=0.0;
   tmp2=0.0;
   for(i=0;i<nspec;i++)
   {
      tmp1=tmp1+x[i];
      tmp2=tmp2+cw[i];
   }

   for(i=0;i<nspec;i++)
   {
      tmp=0.0;
      for(j=0;j<nspec;j++)
         if(j != i) tmp = tmp+x[j]/Dij[i][j];

      Di[i] = (1.0-x[i]/tmp1)/tmp;
   }

   wall_flux(Jw,rhoiw,Mw);
   /*printf("\n 3 \n ");/*pietro*/   
   /*printf(" %le %le %le %le %le ",Jw[0],Jw[1],Jw[2],Jw[3],Jw[4]);/*pietro*/

   for(i=0;i<nspec;i++) h1[i]=Jw[i]*fact;

   /*printf("\n h1 %le %le %le %le %le ",h1[0],h1[1],h1[2],h1[3],h1[4]);/*pietro*/
   for(i=0;i<nspec;i++)
   {
      f1[i] = -rhow*Di[i]/deta*K_bl;
      tmp=0.0;
      tmp1=0.0;
      for(j=0;j<nspec;j++)
      {
         tmp=tmp+rhow*Di[j]*(csw[j]+cw[j]*Msw/M)*K_bl;
         /*tmp=tmp+rhow*Di[j]*csw[j];
         tmp1=tmp1+cw[j];*/
      }
      /*g1[i] = tmp/tmp1;*/
      g1[i] = -rhow*Di[i]*Msw/M*K_bl+tmp/tmp2;
   }


   /* Computes derivatives of partial densities with respec to mass fractions */

   for(i=0;i<nspec;i++)
   {
      for(j=0;j<nspec;j++)
      {
//         drdc[i][j] = rhow*delta(i,j)+rhoiw[i]*(-M/Mw[j]+hspwall[j]/cpwall/twall);
         drdc[i][j] =rhow*delta(i,j)+rhoiw[i]*(-M/Mw[j]+hspwall[j]/cpwall/twalldelta);
      }
   }

   /* Computes derivatives of h1 with respect to partial densities */

   for(i=0;i<nspec;i++) rhoiw1[i]=rhoiw[i];

   for(kder=0;kder<nspec;kder++)
   {
      dri[kder] = meps*fabs(rhoiw[kder]);
      rhoiw1[kder] = rhoiw[kder]+dri[kder];
      if(fabs(rhoiw[kder])<=1.0)
      {
         dri[kder] = eps;
         rhoiw1[kder] = rhoiw[kder]+eps;
      }
      dri[kder] = rhoiw1[kder]-rhoiw[kder];

   /* Computes new density  */
      rhow1=0.0;
      for(i=0;i<nspec;i++) rhow1 = rhow1+rhoiw1[i];
      for(i=0;i<nspec;i++) cw1[i] = rhoiw1[i]/rhow1;
   /* Computes the new wall flux */   
      wall_flux(h1b,rhoiw1,Mw);
      /*printf("\n 4 \n ");/*pietro*/      
      /*printf(" %le %le %le %le %le ",Jw[0],Jw[1],Jw[2],Jw[3],Jw[4]);/*pietro*/

      for(i=0;i<nspec;i++) h1b[i]=h1b[i]*fact;

      for(kk=0;kk<nspec;kk++) dh1dr[kk][kder] = (h1b[kk]-h1[kk])/dri[kder];
      rhoiw1[kder] = rhoiw[kder];
   }

   for(i=0;i<nspec;i++)
   {
      for(k=0;k<nspec;k++)
      {
         tmp = 0.0;
         for(l=0;l<nspec;l++)
         {
            tmp = tmp+dh1dr[i][l]*drdc[l][k];
         }
         dh1dc[i][k] = tmp;
      }
   }

   free_vector(Jw); free_vector(rhoiw); free_vector(rhoiw1); free_vector(h1b);
   free_vector(Di); free_vector(x); free_matrix(drdc,nspec); free_vector(dri);
   free_matrix(dh1dr,nspec); free_vector(cw1);
}
