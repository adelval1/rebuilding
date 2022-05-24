#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "peg_def.h"
#define WL 80
#include "redefine.h"
#include "externalvariables.h"


/*----------------------------------------------------------------------
     Reads the parameters needed for the memory allocation
     and for Pegase routines initialisation
-----------------------------------------------------------------------*/

void read_param()
{
   char dummystr[WL];
   int i,j;
   FILE *fp;
   char file_bl_mode[]="neboulainput/bl_mode.in";

   if((fp = fopen(file_bl_mode,"r"))==NULL)
   {
      printf("Cannot open file %s\n",file_bl_mode);exit(1);
   }
   else
   {
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d\n",&only_stag);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d\n",&only_stag_mod);
      fgets(dummystr,WL,fp);
/* If finite_thick is 1 the exact b.l. edge is taken into account. It can be used
   only in the stagnation point */
      fscanf(fp,"%d\n",&finite_thick);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d\n",&neq);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d\n",&ngrid);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d\n",&nwrite);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d\n",&nrwall);
      fgets(dummystr,WL,fp);
/* If diff_type is 1 Stefan Mawell equations are used to compute the diffusion
   fluxes, if diff_type is 2 Ramshaw formula is used  */
      fscanf(fp,"%d\n",&diff_type);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d\n",&axisymm);    /* Geometry informations */
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d\n",&cone);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d\n",&two_d);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d\n",&flat_plate);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d\n",&flg_anha);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d\n",&oper);
      fgets(dummystr,WL,fp);                      /*NANNI*/
      fscanf(fp,"%d\n",&total_therm_cond);        /*NANNI*/
      fgets(dummystr,WL,fp);                      /*NANNI*/
      fscanf(fp,"%d\n",&n_elements);              /*NANNI*/
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d\n",&oper_traco);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d\n",&mode);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d\n",&flg_neq);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d\n",&flg_termo);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d\n",&flg_traco);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d\n",&flg_stop);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d\n",&sonine);
      fgets(dummystr,WL,fp);             /*NANNI*/
      fscanf(fp,"%d\n",&eq_flg);         /*NANNI*/
      fgets(dummystr,WL,fp);             /*jan*/
      fscanf(fp,"%d\n",&choice);         /*jan*/
      fgets(dummystr,WL,fp);             /*jan*/
      fscanf(fp,"%d\n",&ga_choice);         /*jan*/
/* Read the kind of b.c. condition for the wall */
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d \n",&lew);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d \n",&stef);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d \n",&ram);
   }
   fclose(fp);
   if(axisymm)
   {
      if(two_d)
      {
         printf("\n The geometry cannot be axisymmetric and 2D at the same time \n");exit(-1);
      }
   }

   if(two_d)
   {
      if(axisymm)
      {
         printf("\n The geometry cannot be 2D and axisymmetric at the same time \n");exit(-1);
      }
   }

   if((!axisymm) && (!two_d))
   {
      printf("\n The geometry is not set \n");exit(-1);
   }

/* If only_stag is on ngrid and nwrite should be equal to one */
   if(only_stag)
   {
      if(ngrid!=1)
      {
         printf("\n When only_stag is on ngrid should be equal to one \n");exit(-1);
      }
      if(nwrite!=1)
      {
         printf("\n When only_stag is on nwrite should be equal to one \n");exit(-1);
      }
//      if(finite_thick)
//         printf("\n Finite thickness is taken into account \n");
   }
   else
   {
      if(only_stag_mod)
      {
         printf("\n When only_stag is off only_stag_mod should be equal to zero \n");exit(-1);
      }
      if(finite_thick)
      {
         printf("\n Finite thickness is allowed only if only_stag is on \n");exit(-1);
      }
   }

   if( (diff_type!=1)&&(diff_type!=2) )
   {
      printf("\n diff_type should be equal to 1 or 2 and not %d \n",diff_type);exit(-1);
   }


   if(only_stag)
   {
//      printf("\n The code will perform only the computation in the \n");
//      printf("stagnation point! Remember to give the right input data \n");
   }

/* Perform some tests to see if the b.c. is correct */
      if(lew)
      {
         if(stef)
         {
            printf("\n Error in the b.c. Specified lew and stef \n");exit(-1);
         }
         else if(ram)
         {
            printf("\n Error in the b.c. Specified lew and ram \n");exit(-1);
         }
      }
      else if(stef)
      {
         if(lew)
         {
            printf("\n Error in the b.c. Specified stef and lew \n");exit(-1);
         }
         else if(ram)
         {
            printf("\n Error in the b.c. Specified stef and ram \n");exit(-1);
         }
      }
      else if(ram)
      {
         if(lew)
         {
            printf("\n Error in the b.c. Specified ram and lew \n");exit(-1);
         }
         else if(stef)
         {
            printf("\n Error in the b.c. Specified ram and stef \n");exit(-1);
         }
      }
      else
      {
         printf("\n Error in the b.c. Nothing set for lew, stef, ram \n");exit(-1);
      }


}



/**********************************************
     reads the data   
     returns: all flow properties
************************************************/
void read_flow_conditions(double *xgrid,double *rhoext,double *twall,double *muext,double *uext,double *rbody,double *hext,double *pext,
	double **cext,double *h,double *t,double **c,double *eta,double *etastop,double *xwr,int *maxit,double *eps,double **nu,double ***mu,
	double **gam,double *rhoinf,double *uinf,double *under_eps,double *duedx,int **lambda_matrx){
#include "newtloop.h"

   int i,j,k;
   char file_num_par[]="neboulainput/numerical_parameter.in";
   char file_chem[]="neboulainput/wall_chemistry.in";
   char file_bc_ini[]="neboulainput/bc_ini.in";
   char file_bc_wall[]="neboulainput/bc_wall.in";
   char file_bc_out[]="neboulainput/bc_out.in";
   char file_lam[]="neboulainput/lambda_matrix_for_demixing.in";
   char dummystr[WL];
   FILE *fp;
/* In this file quantities like the number of points in eta direction, the
  integration tolerance, the maximum number of iterations for each integration 
  and conditions on the fluid (frozen,equilibrium) and on the wall catalicity are defined */

	
   if((fp = fopen(file_num_par,"r"))==NULL)
   {      printf("Cannot open file %s\n",file_num_par);exit(1);}
   else
   {
      fgets(dummystr,WL,fp);
      fscanf(fp,"%le\n",etastop);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%le\n",eps);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%le\n",under_eps);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d\n",maxit);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%le\n",&macheps);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d\n",&newton_it);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%le\n",&sutres);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%d\n",&sutit);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%le\n",&tdelta);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%le\n",rhoinf);
      fgets(dummystr,WL,fp);
      fscanf(fp,"%le\n",uinf);

      fgets(dummystr,WL,fp);
      for(i=0;i<nwrite;i++)
      {
         fscanf(fp,"%le\n",(xwr+i));
      }
   }
   fclose(fp);


   deta = *etastop/(double)(neq-1);

   for(i=0;i<neq;i++)
   {
      eta[i] = *etastop*(double)i/(double)(neq-1);
   }


/* Wall reaction rates are here */
   if((fp = fopen(file_chem,"r"))==NULL)
   {      printf("Cannot open file %s\n",file_chem);exit(1);}
   else
   {
/*    If it is a local equilibrium wall there is no need to read the reactiions */
      if(lew)
      {
      }
      else
      {
/* Matrices for reactions at the wall */
/* Read matrix nu */
         fgets(dummystr,WL,fp);
         fgets(dummystr,WL,fp);
         /*printf("%s",dummystr);*/

         for(i=0;i<nrwall;i++)
         {

            for(j=0;j<nspec;j++)
            {
               fscanf(fp,"%le",&nu[i][j]);
            }
            fgets(dummystr,WL,fp);
         }
         fgets(dummystr,WL,fp);
/*    read matrix mu  */

         for(i=0;i<nrwall;i++)
         {
            for(j=0;j<nspec;j++)
            {
               for(k=0;k<nspec;k++)
               {
                  fscanf(fp,"%le",&mu[i][j][k]);
               }
               fgets(dummystr,WL,fp);
            }
            fgets(dummystr,WL,fp);
         }
      }

   }
   fclose(fp);


   if((fp = fopen(file_lam,"r"))==NULL)
   {      printf("Cannot open file %s\n",file_lam);exit(1);}
   else
   {
      fgets(dummystr,WL,fp);
     for(i=0;i<nspec;i++)
        {  for(j=0;j<n_elements;j++){fscanf(fp,"%d",&lambda_matrx[i][j]);}        }
      }
	
	
/* initial conditions are written in this file */
   if((fp = fopen(file_bc_ini,"r"))==NULL)
   {
      printf("Cannot open file %s\n",file_bc_ini);
      exit(1);
   }
   else
   {
         fgets(dummystr,WL,fp);
         for(i=0;i<neq;i++)
         {fscanf(fp,"%le %le\n",&t[i],&h[i]);}
         fgets(dummystr,WL,fp);
         for(i=0;i<neq;i++)
         {
            for(j=0;j<nspec;j++)
            {fscanf(fp,"%le",&c[j][i]);}
            fgets(dummystr,WL,fp);
         }
   }
   fclose(fp);

/* outer conditions are written in this file */
   if((fp = fopen(file_bc_out,"r"))==NULL)
   {
      printf("Cannot open file %s\n",file_bc_out);
      exit(1);
   }
   else
   {
      if(only_stag)
      {
         fgets(dummystr,WL,fp);
         for(i=0;i<ngrid;i++){fscanf(fp,"%le %le %le %le\n",(xgrid+i),(rbody+i),(rhoext+i),(muext+i));}
         fgets(dummystr,WL,fp);
         for(i=0;i<ngrid;i++){fscanf(fp,"%le %le %le\n",(uext+i),(hext+i),(pext+i));}
         fgets(dummystr,WL,fp);
         for(i=0;i<ngrid;i++){fscanf(fp,"%le \n",(duedx+i));}
         fgets(dummystr,WL,fp);

         fscanf(fp,"%le %le %le \n",&v_e,&duedy,&delta_bl);
         fgets(dummystr,WL,fp);

         for(i=0;i<ngrid;i++)
 	 {
            for(j=0;j<nspec;j++)	{fscanf(fp,"%le",&cext[j][i]);}
         }
      
      }
      else
      {
         fgets(dummystr,WL,fp);
         for(i=0;i<ngrid;i++){fscanf(fp,"%le %le %le %le\n",(xgrid+i),(rbody+i),(rhoext+i),(muext+i));}
         fgets(dummystr,WL,fp);
         for(i=0;i<ngrid;i++){fscanf(fp,"%le %le %le\n",(uext+i),(hext+i),(pext+i));}
         fgets(dummystr,WL,fp);
         for(i=0;i<ngrid;i++)
         {
            for(j=0;j<nspec;j++)	{fscanf(fp,"%le",&cext[j][i]);}
         }
      
      }
   }
   fclose(fp);
/* wall conditions are written in this file */
 if((fp = fopen(file_bc_wall,"r"))==NULL)
   {
      printf("Cannot open file %s\n",file_bc_wall);
      exit(1);
   }
   else
   {
         fgets(dummystr,WL,fp);
         for(i=0;i<ngrid;i++){fscanf(fp,"%le  \n",(twall+i));   }
         fgets(dummystr,WL,fp);
         for(i=0;i<ngrid;i++){
		for(j=0;j<nrwall;j++)
			fscanf(fp,"%le ",&gam[j][i]);
  			
         	}
     		fscanf(fp,"\n ");
       
 	}
   fclose(fp);
   
   
}

#undef WL
