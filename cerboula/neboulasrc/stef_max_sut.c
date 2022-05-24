#include <math.h>
#include "externalvariables.h"

void ludcmp(double **,int,int *,double *);
void lubksb(double **,int,int *,double *);
double *vector(int);
double **matrix(int,int);
void free_vector(double *);
void free_matrix(double **,int);
void ableit(double *,double *,double,int);
void wall_flux(double *,double *,double *);


void stef_max(double **c,double **cs,double *rho,double ***Dij,double *Mw,double *M,
              double *Ms,double **J,double **Js,double duedx,double xit2,int jstart,
              void (*diffusion)(double ,double *,double *,double *,double **,
                                double ,double ,double *,double *,double) )
{
   extern int lew,stef,ram; /* Kind of b.c. to be applied */
   extern int axisymm,two_d,cone,flat_plate; /* Geometric informations */

   int i,j,k,kstart;
   double *Yi,*dYi,*Jw,*Jtmp,fact;
   double *rhoiw,*x,**Dbin;


   x=vector(nspec); Yi=vector(nspec); dYi=vector(nspec); Dbin=matrix(nspec,nspec);
   Jtmp=vector(nspec); rhoiw=vector(nspec); Jw=vector(nspec);

/*   printf("\n Begin stef_max \n"); */
   if(jstart != 0)
   {
      fact=udelta*rdelta/sqrt(xit2);
   }
   else
   {
      if(axisymm)
         if(cone)
            fact=1.0/(rhodelta*mudelta);
/* We modify to take into account the singularity at the cone tip. The factor
   above is multiplied by sqrt(2*csi)*rdelta/dcsidelta and the Ji are by consequence
   multiplied by this term. The same for the electric field */
         else
            fact=sqrt(2.0*duedx/rhodelta/mudelta);

      if(two_d)
         if(flat_plate)
            fact=1.0/(rhodelta*mudelta);
/* We modify to take into account the singularity at the flat plate tip. The factor
   above is multiplied by sqrt(2*csi)*rdelta/dcsidelta and the Ji are by consequence
   multiplied by this term. The same for the electric field */
         else
            fact=sqrt(duedx/rhodelta/mudelta);
   }


/* Apply the b.c. If stef kstart=0 and the Jwall is computed from Stefan-Maxwell.
   If ram or lew kstart=0 and the Jwall is computed from the Stefan-Maxwell equations */
   if(stef) kstart=0;
   else if(lew) kstart=0;
   else if(ram) kstart=0;
   
   for(k=kstart;k<neq;k++)
   {

      for(i=0;i<nspec;i++)
      {
         Yi[i]=c[i][k];
         dYi[i]=cs[i][k];
         x[i]=Yi[i]*M[k]/Mw[i];
         Jtmp[i]=0.0;
         for(j=0;j<nspec;j++) Dbin[i][j]=Dij[k][i][j];
      }

   /*printf("\n Before diffusion \n");*/
      diffusion(rho[k],Yi,x,dYi,Dbin,M[k],Ms[k],Mw,Jtmp,fact);
   /*printf("\n After diffusion \n");*/


      for(i=0;i<nspec;i++) J[i][k]=Jtmp[i];
   }

/* Wall flux, for lew and ram is already computed */
   /*Paolo*/
  /* if(stef)
   {
      for(i=0;i<nspec;i++) rhoiw[i]=rho[0]*c[i][0];

      wall_flux(Jw,rhoiw,Mw);

      for(i=0;i<nspec;i++) J[i][0] = Jw[i];

   }*/
/*   else if(ram)
   {
      double *Di,tmp,tmp1;
      Di=vector(nspec);

      for(i=0;i<nspec;i++)
         x[i] = c[i][0]*M[0]/Mw[i];

      for(i=0;i<nspec;i++)
      {
         tmp=0.0;
         for(j=0;j<nspec;j++)
            if(j != i) tmp = tmp+x[j]/Dij[0][i][j];

         Di[i] = (1.0-x[i])/tmp;
      }
      for(i=0;i<nspec;i++)
      {
         tmp=0.0;
         tmp1=0.0;
         for(j=0;j<nspec;j++)
         {
            tmp=tmp+rho[0]*Di[j]*(cs[j][0]+c[j][0]*Ms[0]/M[0]);
            tmp=tmp+rho[0]*Di[j]*cs[j][0];
            tmp1=tmp1+c[j][0];
         }
         J[i][0] = -rho[0]*Di[i]*cs[i][0] + tmp/tmp1;
        J[i][0] = -rho[0]*Di[i]*cs[i][0]+c[i][0]*(-rho[0]*Di[i]*Ms[0]/M[0]+tmp);
      }
      free_vector(Di);
   } */


   for(i=0;i<nspec;i++)
   {
      ableit(J[i],Js[i],deta,neq);
   }
   
   free_vector(x); free_vector(Yi); free_vector(dYi); free_matrix(Dbin,nspec);
   free_vector(Jtmp); free_vector(rhoiw); free_vector(Jw);

}




void stef_max_sutton(double rho,double *c,double *x,double *cs,double **Dbin,
                     double M,double Ms,double *Mw,double *Jtmp,double fact)
{
   extern double K_bl;

   int i,j,it;
   double *v1,*v2,*v3,*v2v4,*Jold,*Dim,s1,s2,s3,E,tmp,tmp1,tmp2,res;

   v1=vector(nspec); v3=vector(nspec); v2v4=vector(nspec); Jold=vector(nspec);
   Dim=vector(nspec); v2=vector(nspec);

   s1=0.0;
   E=0.0;
   s3=1.0;

   res=sutres*2.0;

   
   tmp2 = 0.0;
   for(i=0;i<nspec;i++)
      tmp2 = tmp2+c[i];

   for(i=0;i<nspec;i++)
   {
      tmp = 0.0;
      for(j=0;j<nspec;j++)
      {
         tmp = tmp+x[j]/Dbin[i][j];
      }
      Dim[i] = 1.0/tmp;
   }

   for(i=0;i<nspec;i++)   /* The K_bl factor takes into account finite thickness */
   {
      v1[i] = -rho*Mw[i]/M*Dim[i]*rho*fact*(cs[i]*M/Mw[i]+c[i]*Ms/Mw[i])*K_bl;
      v2[i] = c[i]*Dim[i]*M;
   }
   if(chflag)
   {
      for(i=0;i<nspec;i++)
         v3[i]=rho*rho*Mw[i]/M*Dim[i]/pdelta*charge[i]*c[i];
         /*v3[i]=charge[i]*(rho*Mw[i]/M)*(rho*Mw[i]/M)*Dim[i]*x[i]/pdelta;*/
      s1=0.0;
      s3=0.0;
      for(i=0;i<nspec;i++)
      {
         s1=s1+charge[i]*v1[i];
         s3=s3+charge[i]*v3[i];
      }
   }

   
   it = 0;
   while(1)
   {

      for(i=0;i<nspec;i++)
      {
         tmp1=0.0;
         for(j=0;j<nspec;j++) tmp1=tmp1+Jtmp[j]/(Mw[j]*Dbin[i][j]);
         v2v4[i] = tmp1*v2[i];
         Jold[i] = Jtmp[i];
      }


/*  Computes ambipolar E */
      s2=0.0;
      for(i=0;i<nspec;i++) s2=s2+charge[i]*v2v4[i];
      E = -(s1+s2)/s3;
  

/*  Uncorrected fluxes */
      for(i=0;i<nspec;i++) Jtmp[i] = v1[i]+v2v4[i]+v3[i]*E;


/*  Corrected fluxes */
      tmp1=0.0;
      for(i=0;i<nspec;i++) tmp1=tmp1+Jtmp[i];

      for(i=0;i<nspec;i++) Jtmp[i] = Jtmp[i]-c[i]/tmp2*tmp1;


      res=0.0;
      for(i=0;i<nspec;i++) res=res+(Jtmp[i]-Jold[i])*(Jtmp[i]-Jold[i]);


      res=sqrt(res);
      it=it+1;


      if((sutres==0.0)&&(it>=sutit)) break;
      if(res<sutres) break;
   }

   
   free_vector(v1); free_vector(v3); free_vector(v2v4); free_vector(Jold);
   free_vector(Dim); free_vector(v2);

}




void ramshaw(double rho,double *c,double *x,double *cs,double **Dbin,
                     double M,double Ms,double *Mw,double *Jtmp,double fact)
{
   extern double K_bl;

   int i,j,it;
   double *v1,*v2,*v3,*v2v4,*Jold,*Dim,s1,s2,s3,E,tmp,tmp1,tmp2,res;

   v1=vector(nspec); v3=vector(nspec); v2v4=vector(nspec); Jold=vector(nspec);
   Dim=vector(nspec); v2=vector(nspec);

   s1=0.0;
   E=0.0;
   s3=1.0;
   res=sutres*2.0;

   tmp2 = 0.0;
   for(i=0;i<nspec;i++)
      tmp2 = tmp2+x[i];

   for(i=0;i<nspec;i++)
   {
      tmp = 0.0;
      for(j=0;j<nspec;j++)
      {
         if(j!=i) tmp = tmp+x[j]/Dbin[i][j];
      }
      Dim[i] = (1.0-x[i]/tmp2)/tmp;
   }

   tmp2 = 0.0;
   for(i=0;i<nspec;i++)
      tmp2 = tmp2+c[i];

   for(i=0;i<nspec;i++)  /* K_bl takes into account finite thickness */
   {
      v1[i] = rho/M*Mw[i]*Dim[i]*rho*fact*(cs[i]*M/Mw[i]+c[i]*Ms/Mw[i])*K_bl;
   }

   tmp1=0.0;
   for(i=0;i<nspec;i++) tmp1=tmp1+v1[i];

/* Non corrected for the electrix field */
   for(i=0;i<nspec;i++) Jold[i]=-v1[i]+c[i]/tmp2*tmp1;

   if(chflag)
   {
      for(i=0;i<nspec;i++)
         v3[i]=rho*rho*Mw[i]/M*Dim[i]/pdelta*charge[i]*c[i];
      s1=0.0;
      s3=0.0;
      for(i=0;i<nspec;i++)
      {
         s1=s1+charge[i]*Jold[i];
         s3=s3+charge[i]*v3[i];
      }
      E=-s1/s3;
      for(i=0;i<nspec;i++) Jtmp[i]=Jold[i]+v3[i]*E;
   }
   else
      for(i=0;i<nspec;i++) Jtmp[i]=Jold[i];



   free_vector(v1); free_vector(v3); free_vector(v2v4); free_vector(Jold);
   free_vector(Dim); free_vector(v2);

}
