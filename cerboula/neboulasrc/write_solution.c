#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "peg_def.h"
#define WL 80
#include "redefine.h"
#include "externalvariables.h"
/*---------------------------------------------------
c      subroutine write_solution
c---------------------------------------------------  */

void write_solution(double *eta,double *rho,double *Mw,
                    double *u,double *h,double *t,double **c,
                    double qw,double cf,double St,double xi,double x_wr,double dxi,
                    int it,int jstart,double duedx,
                    double qwc,double qwd,double *M,double *Jw,double **J,double
		    *V,int iwr)
{

   extern double R;
   extern double K_bl;
   double number_density_(double *,double *);
   double *vector(int);
   void free_vector(double *);
   void continuity(double *,double *);

   int i,j,k;
   char str_tmp[WL];
   FILE *fp;
   double p_sp,n_den,tmp,*Yi,*Y,*yc,*rhotmp,fact,el1,el2;
   char file_heat_skin[]="neboulaoutput/heatskin.dat";
   char file_prof[]="neboulaoutput/profile";
   char file_stag_sol[]="neboulaoutput/stagsol.dat";
   char file_mol_weight[]="neboulaoutput/molecular_weight_mix.dat";
   char file_flowfield_tecplot[]="neboulaoutput/flowfield.plt";
   
   Yi=vector(nspec);   Y=vector(nspec); yc=vector(neq); rhotmp=vector(neq);

   if(jstart==0)
   {
      if(axisymm)
      {
        if(cone) fact=1.0;
        else
          fact=sqrt(2.0*duedx/rhodelta/mudelta);
      }
      if(two_d)
      {
        if(flat_plate) fact=1.0;
        else fact=sqrt(duedx/rhodelta/mudelta);
      }
   }
   else
   {
      if(axisymm) for(k=0;k<n_elements;k++)
        fact=udelta*rdelta/sqrt(2.0*xi);
      if(two_d)
        fact=udelta*rdelta/sqrt(2.0*xi);
   }


   if(jstart==0)
   {
      if((fp = fopen(file_stag_sol,"w"))==NULL)
      {
         printf("Cannot open file %s \n",file_stag_sol);
         exit(1);
      }
      else
      {
         fprintf(fp," %s \n","T, h at the stagnation point");
         for(i=0;i<neq;i++)
         {
            fprintf(fp," %.8e %.8e \n",t[i],h[i]*hdelta);
         }
         fprintf(fp," %s \n","Species at the stagnation point");
         for(i=0;i<neq;i++)
         {
            for(j=0;j<nspec;j++)
            {
               fprintf(fp," %.8e",c[j][i]);
            }
            fprintf(fp," \n");
         }
      }
      fclose(fp);
   }
   


/* Gets the physical coordinate y from the transformed one */
   for(i=0;i<neq;i++) rhotmp[i]=1.0/rho[i];
   continuity(rhotmp,yc);
   for(i=0;i<neq;i++) yc[i]=yc[i]/fact/K_bl;


   if((fp = fopen(file_heat_skin,"a"))==NULL)
   {
      printf("Cannot open file %s\n",file_heat_skin);exit(1);
   }
   else
   {
      fprintf(fp,"%e %e %e %d %e %e %e %e %e %e %e %e %e\n",x_wr,xi,dxi,it,pdelta,tdelta,udelta,qw,cf,St,qwc,qwd,twalldelta);
   }
   fclose(fp);
  
  
    if(iwr==0)
    {
   	if((fp = fopen(file_flowfield_tecplot,"w"))==NULL)
   	{	printf("Cannot open file %s\n",file_flowfield_tecplot);exit(1);}
   	else
   	{      
   	fprintf(fp,"title = \"boundary layer flow field\"\n");
   	fprintf(fp,"variables =\n");
   	fprintf(fp,"y\n");
   	fprintf(fp,"x\n");
   	fprintf(fp,"u\n");
   	fprintf(fp,"rho\n");
   	fprintf(fp,"t\n");
   	fprintf(fp,"h\n");
	for(j=0;j<nspec;j++){fprintf(fp," c%d\n",j+1);}
	for(j=0;j<nspec;j++){fprintf(fp," x%d\n",j+1);}	
//	elemental fraction for co25 only, if other mixture is used comment the next two line	
   	fprintf(fp,"elem1\n");
   	fprintf(fp,"elem2\n");		
   	fprintf(fp,"zone t = \"boundary layer\" j= %d i = %d \n" ,nwrite, neq);
      	for(i=0;i<neq;i++)
      	{
         fprintf(fp," %.6e %.6e %.8e %.8e %.8e %.8e",yc[i],x_wr,u[i],rho[i],t[i],h[i]);
	tmp=0.0;
	for(j=0;j<nspec;j++){fprintf(fp," %.8e",c[j][i]);tmp=tmp+c[j][i]/Mw[j];Y[j]=c[j][i]/Mw[j];}
	for(j=0;j<nspec;j++){Y[j]=Y[j]/tmp;}	
	for(j=0;j<nspec;j++){fprintf(fp," %.8e",Y[j]);}
//	elemental fraction for co25 only, if other mixture is used comment the next three line
	el1=(2*Y[0]+2*Y[1]+Y[2]+Y[3])/(3*Y[0]+2*Y[1]+2*Y[2]+Y[3]+Y[4]);	
	el2=(Y[0]+Y[2]+Y[4])/(3*Y[0]+2*Y[1]+2*Y[2]+Y[3]+Y[4]);	
	fprintf(fp," %.8e %.8e",el1,el2);
	fprintf(fp,"\n");

      	}
   	}
   fclose(fp);
   }
   else
   {
   	if((fp = fopen(file_flowfield_tecplot,"a"))==NULL)
   	{	printf("Cannot open file %s\n",file_flowfield_tecplot);exit(1);}
   	else
   	{      

      	for(i=0;i<neq;i++)
      	{
         fprintf(fp," %.6e %.6e %.8e %.8e %.8e %.8e ",yc[i],x_wr,u[i],rho[i],t[i],h[i]);
	tmp=0.0;
	for(j=0;j<nspec;j++){fprintf(fp," %.8e",c[j][i]);tmp=tmp+c[j][i]/Mw[j];Y[j]=c[j][i]/Mw[j];}
	for(j=0;j<nspec;j++){Y[j]=Y[j]/tmp;}	
	for(j=0;j<nspec;j++){fprintf(fp," %.8e",Y[j]);}
//	elemental fraction for co25 only, if other mixture is used comment the next three line
	el1=(2*Y[0]+2*Y[1]+Y[2]+Y[3])/(3*Y[0]+2*Y[1]+2*Y[2]+Y[3]+Y[4]);	
	el2=(Y[0]+Y[2]+Y[4])/(3*Y[0]+2*Y[1]+2*Y[2]+Y[3]+Y[4]);	
	fprintf(fp," %.8e %.8e",el1,el2);
      	}
	}   
   fclose(fp);  
   }



   free_vector(Yi); free_vector(yc); free_vector(rhotmp);
   }


#undef WL
