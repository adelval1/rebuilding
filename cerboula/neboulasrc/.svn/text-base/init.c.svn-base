#include "externalvariables.h"

void search(double,double *,int *,int *);

void inizialize_data(double *u,double *u1,double *u2,double *up1,double *h,double *h1,
                     double *h2,double *hp1,
                     double **c,double **c1,double **c2,double **cp1,
                     double *eta,double etastop,
                     double *xgrid,double *csi,double *xwr,double *csiwr)
{
   int i,i1,i2,j;
   if(only_stag)
   {
      csiwr[0]=0.0;
   }
   else
   {
      for(i=0;i<nwrite;i++)
      {
         search(xwr[i],xgrid,&i1,&i2);
         if(i1!=i2)
         {
            csiwr[i]=(xgrid[i2]-xwr[i])/(xgrid[i2]-xgrid[i1])*csi[i1]+(xwr[i]-xgrid[i1])/(xgrid[i2]-xgrid[i1])*csi[i2];
         }
         else
         {
            csiwr[i]=csi[i1];
         }
      }
   }

   for(i=0;i<neq;i++)
   {
      u[i]=eta[i]*(2.0-eta[i]/etastop)/etastop;
      u1[i]=u[i];		u2[i]=u[i];		up1[i]=0.0;
      for(j=0;j<nspec;j++)
      {
         c1[j][i]=c[j][i];	c2[j][i]=c[j][i];	cp1[j][i]=0.0;
      }
      h[i]=h[i]/h[neq-1];		
      h1[i]=h[i];		h2[i]=h[i];		hp1[i]=0.0;
   }
}
