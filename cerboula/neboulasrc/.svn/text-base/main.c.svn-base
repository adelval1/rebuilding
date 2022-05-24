#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "prototip.h"
#include "peg_def.h"
#define WL 80
#include "redefine.h"
#include "externalvariables.h"
void MAIN__(void) 
{
//	nwrite is the number of points in x/csi direction where solution will be written
//	xwr[nwrite] is the array of actual points in x direction where the result is written
//	x_wr is the actual points in x direction where the result is written
//	iwr is the counter for xwr[]; it goes from 0 to nwrite

//	ngrid is the number of points where data is give at the outer edge and at the wall
//	xgrid[ngrid] is the actual position in x direction where outer bc are given
//	csi[ngrid] is the array of the transformed coordinate

//	neq is the number of points in y/eta direction

//	xi is the position in csi direction and dimensions where the code trys to compute and converge
//	xi1 is the last position in csi direction and in transformed units
//	xi2 is the second last position in csi direction and in transformed units
//	xit2 is xi times 2
   int maxit,i1,i2,jstart,it,i,j,iwr,k, stopnow, indx1,indx2,indx3,indx4;
   double *u,*u_old,*u1,*u2,*up1,*us;
   double *h,*h_old,*h1,*h2,*hp1;
   double **c,**c_old,**c1,**c2,**cp1,**cs,**css;
   double **csbis,**cssbis;
   double *V,*t,*ts,*y,*tmp;
   double *M,*Ms,*Mw,*Mw_scaled;
   double *l0,*l0s,***Dij,*l3,*l3s,**mapr,***dmapr,*rho;
   double **J,**Js,**hsp,**dhsp,*Jwall,*somma;
   double *xgrid,*csi,*dcsidx,*sdcsidx,*rbody,*srbody,*muext,*smuext,*uext,*suext;
   double *duedx,*sduedx,*pext,*spext,*hext,*shext,*dhdx,*sdhdx,*rhoext,*twall,**cext,**scext,**gam,**sgam;
   double *eta,etastop,eps,conv,x_wr,*xwr,*csiwr,*srhoext,*stwall;
   double xi,xit2,xi1,xi2,lam0,lam1,lam2;
   double rhoinf,uinf,qw,cf,St,hstag,under_eps,rhodelta2,qwc,qwd,cO,cN;;
   double *sum_ci_stef,*add_var;
   FILE *fp;
     
   void stef_max_sutton(double ,double *,double *,double *,double **,double ,double ,double *,double *,double);
   void ramshaw(double ,double *,double *,double *,double **,double ,double ,double *,double *,double);

  indx1=0;  indx2=0;  indx3=0; indx4=0;

/* Perform all the initializations needed to use the code and the Pegase libraries*/
   read_param();
//gives e.g. the number of species nspec
   set_1_(&flg_anha,&oper,&flg_termo,&flg_traco,&flg_stop,&sonine,&R,&pi,&kb,&Na,&nspec);
   Mw=vector(nspec);   Mw_scaled=vector(nspec);   charge=vector(nspec);
   set_2_(Mw,Mw_scaled,charge,&chflag);

/* Dynamic allocation of problem variables  */
   u=vector(neq); u_old=vector(neq); us=vector(neq); u1=vector(neq);   u2=vector(neq); up1=vector(neq);  h=vector(neq); h_old=vector(neq); h1=vector(neq); h2=vector(neq);   hp1=vector(neq); Jwall=vector(nspec); somma=vector(neq);
   c=matrix(nspec,neq); c_old=matrix(nspec,neq); cs=matrix(nspec,neq);   css=matrix(nspec,neq); c1=matrix(nspec,neq); c2=matrix(nspec,neq);
   cp1=matrix(nspec,neq); csbis=matrix(nspec,neq); cssbis=matrix(nspec,neq);   M=vector(neq); Ms=vector(neq);    J=matrix(nspec,neq); Js=matrix(nspec,neq); hsp=matrix(nspec,neq);   dhsp=matrix(nspec,neq);
   nu=matrix(nrwall,nspec);   mu=trimatrix(nrwall,nspec,nspec);
   gam=matrix(nrwall,ngrid); sgam=matrix(nrwall,ngrid);   V=vector(neq); t=vector(neq); y=vector(neq); ts=vector(neq);
   l0=vector(neq); l0s=vector(neq); Dij=trimatrix(neq,nspec,nspec);   l3=vector(neq); l3s=vector(neq);   mapr=matrix(neq,nspec); dmapr=trimatrix(neq,nspec,nspec); rho=vector(neq);   eta=vector(neq);
   xgrid=vector(ngrid); csi=vector(ngrid);   dcsidx=vector(ngrid); sdcsidx=vector(ngrid); rbody=vector(ngrid);   srbody=vector(ngrid); muext=vector(ngrid); smuext=vector(ngrid);
   uext=vector(ngrid); suext=vector(ngrid); duedx=vector(ngrid);   sduedx=vector(ngrid); pext=vector(ngrid); spext=vector(ngrid);    hext=vector(ngrid); shext=vector(ngrid); dhdx=vector(ngrid);
   sdhdx=vector(ngrid); rhoext=vector(ngrid); srhoext=vector(ngrid);   twall=vector(ngrid);   stwall=vector(ngrid);   cext=matrix(nspec,ngrid); scext=matrix(nspec,ngrid); rbody=vector(ngrid);
   cdelta=vector(nspec);   gamdelta=vector(nrwall); hspwall=vector(nspec);   xwr=vector(nwrite); csiwr=vector(nwrite);   sum_ci_stef=vector(neq); add_var=vector(neq);   lambda_matrx=matrix_int(nspec,n_elements);   tmp=vector(n_elements);
   

   read_flow_conditions(xgrid,rhoext,twall,muext,uext,rbody,hext,pext,cext,h,t,c,eta,&etastop,xwr,&maxit,&eps,nu,mu,gam,&rhoinf,&uinf,&under_eps,duedx,lambda_matrx);
   grid_gen(xgrid,rhoext,muext,uext,hext,rbody,csi,dcsidx,duedx,dhdx);   
   spline_gen(xgrid,rbody,uext,hext,duedx,dhdx,dcsidx,cext,pext,muext,srbody,suext,shext,sduedx,sdhdx,sdcsidx,scext,spext,smuext,rhoext,srhoext,twall,stwall,gam,sgam);


//  remove old files
    remove("neboulaoutput/convergence.dat");
    remove("neboulaoutput/heatskin.dat");
        
   for (i=0;i<nrwall;i++){
	   for (j=0;j<ngrid;j++){
		if(gam[i][j]>1){printf("warning this version of Neboula works with a gamma definition between 0 and 1");getchar();}
   		gam[i][j]=2*gam[i][j]/(2-gam[i][j]);
	   }
   }


   for(i=0;i<nspec;i++) {cdelta[i]=cext[i][0]; }
   for(i=0;i<nrwall;i++) {gamdelta[i]=gam[i][0]; }
   dhdelta=0.0; pdelta=pext[0];mudelta=muext[0]; udelta=uext[0]; hdelta=hext[0]; hstag=hext[0];twalldelta=twall[0];

   inizialize_data(u,u1,u2,up1,h,h1,h2,hp1,c,c1,c2,cp1,eta,etastop,xgrid,csi,xwr,csiwr);
      /* It is probably better to initialize data outside of the program with PEGASE */

   /* Initialise the coordinates */
   a_e=0.0;
   if(only_stag)
   {
      dcsi=1.0e-6;  xi=0.0; xit2=0.0; xi1=0.0; xi2=0.0; lam0=1.0/dcsi; lam1=-1.0/dcsi; lam2=0.0;
      if(only_stag_mod) a_e = v_e*duedy/(duedx[0]*duedx[0]);
   }
   else
   {
      dcsi=csi[1];  xi=0.0; xit2=0.0; xi1=0.0; xi2=0.0; lam0=1.0/dcsi; lam1=-1.0/dcsi; lam2=0.0;
   }
   
   if(axisymm) if(cone) beta=0.0; else beta=0.5;
   if(two_d) if(flat_plate) beta=0.0; else beta=1.0;

   rhodelta2=rhoext[0];

/*********************************************************************************/   
/*  Outer loop: integration in the csi direction   */
/*********************************************************************************/

   jstart=0;
   iwr=0;
   while(xi <= csi[ngrid-1])
   {
      it=0;
      conv=2.0*eps;
      if(jstart >= 1)
      {
 	var_interp(csi,rbody,uext,hext,duedx,dhdx,dcsidx,cext,pext,muext,srbody,suext,shext,sduedx,sdhdx,sdcsidx,scext,spext,smuext,xi,xgrid,xit2,rhoext,srhoext,twall,stwall,gam,sgam,&rhodelta2);
        make_csider(xit2,lam1,lam2,u1,u2,up1,h1,h2,hp1,c1,c2,cp1);
      }

      
   /**************************************************************************/
   /*    Inner loop: integration in the eta direction    */
   /**************************************************************************/
      printf("\n NEBOULA is iterating...");fflush(stdout);
      while((it<=maxit) && (conv>eps))
      {
//	printf("\niter %d %e",it,conv);fflush(stdout);
        for(i=0;i<neq;i++)
         {
            double dummy;
            u_old[i]=u[i];		h_old[i]=h[i];		dummy=0.0;
            for(j=0;j<nspec;j++)	{c_old[j][i]=c[j][i];	dummy=dummy+c[j][i];}
         }
         u[neq-1]=1.0;         		h[neq-1]=1.0;
         for(j=0;j<nspec;j++)		{c[j][neq-1]=cdelta[j];}
	 ent2temp(t,h_old,c_old,hdelta,twalldelta);
      	 der1(u_old,us,c_old,cs,css);

         for(i=0;i<nspec;i++)
            for(j=0;j<neq;j++)		{csbis[i][j]=cs[i][j];	cssbis[i][j]=css[i][j];}
	 msq_der(c_old,csbis,cssbis);
	
	coefficients(l0,Dij,l3,mapr,dmapr,hsp,rho,y,u_old,up1,c_old,cs,t,M,Ms,Mw,xit2,lam0,twalldelta);/*,total_therm_cond);*/

	 /* B.L. finite thickness */
	 if(finite_thick){finite_thickness_chi(rho,duedx[0],&chi,&K_bl);}
         else{K_bl=1.0;}

	 der2(l0,l0s,l3,l3s,hsp,dhsp);
      
	 if(diff_type==1)	{stef_max(c_old,cs,rho,Dij,Mw,M,Ms,J,Js,duedx[0],xit2,jstart,stef_max_sutton);}
         else if(diff_type==2)	{stef_max(c_old,cs,rho,Dij,Mw,M,Ms,J,Js,duedx[0],xit2,jstart,ramshaw);}

       	continuity(y,V);
	if (eq_flg==1)  	{eq_conc_cfe(c,t);}
        /*Paolo*/
	if(eq_flg==2)
        {
           species(c_old,csbis,css,cp1,u_old,V,J,Js,Dij,l0,l0s,mapr,dmapr,M,Ms,rho,Mw,
                   c,xit2,duedx[0],jstart,lam0,Mw_scaled,csbis,Jwall,iwr);
           if(chflag) neutrality(c);
        }
	if (eq_flg==3) 		{elements_eq(c_old,csbis,css,cp1,u_old,V,J,Js,Dij,l0,l0s,mapr,dmapr,M,Ms,rho,Mw,c,xit2,duedx[0],jstart,lam0,Mw_scaled,csbis,t,add_var,it);/*,Jel);*/}	 
	momentum(u_old,up1,rho,l0,l0s,V,u,xit2,lam0);
        energy(u_old,us,h_old,hp1,cs,css,V,rho,l0,l3,l3s,hsp,dhsp,J,Js,h,xit2,duedx[0],jstart,lam0,total_therm_cond);
	
	 conv = convergence(neq,nspec,u,u_old,h,h_old,c,c_old);
         it++;
	 
	 if((fp = fopen("neboulaoutput/convergence.dat","a"))==NULL)	{printf("Cannot open file %s\n","neboulaoutput/convergence.dat");exit(1);}
	 else           					{fprintf(fp,"%d %le \n",it,log(conv));}
	 fclose(fp);
	 
         under_relax(neq,nspec,u,u_old,h,h_old,c,c_old,under_eps);

         if(chflag) neutrality(c);
       }
   /***********************************************************************/
   /*     End of inner loop: eta integration     */
   /***********************************************************************/
//----------------------
      if(it <= maxit+1) //convergence reached within limit(maxit)?
      {
//////////////
         if((xi>=csiwr[iwr]) && (iwr < nwrite))
         {
            if(xi==csi[iwr]) x_wr=xwr[iwr];//here should be a bug - it is pointless to use the counter iwr with csi
            else
            {
               search(xi,csi,&i1,&i2);
               if(i1!=i2){x_wr = (csi[i2]-xi)/(csi[i2]-csi[i1])*xgrid[i1]+(xi-csi[i1])/(csi[i2]-csi[i1])*xgrid[i2];}
               else{x_wr = xgrid[i1];}
            }
	    ent2temp(t,h,c,hdelta,twalldelta);
            ableit(t,ts,deta,neq);            
	    der1(u,us,c,cs,css);
            tdelta=t[neq-1];
            if(diff_type==1)stef_max(c,cs,rho,Dij,Mw,M,Ms,J,Js,duedx[0],xit2,jstart,stef_max_sutton);
            else if(diff_type==2)stef_max(c,cs,rho,Dij,Mw,M,Ms,J,Js,duedx[0],xit2,jstart,ramshaw);
            coefficients(l0,Dij,l3,mapr,dmapr,hsp,rho,y,u_old,up1,c,cs,t,M,Ms,Mw,xit2,lam0,twalldelta);/*,total_therm_cond/*NANNI);*/
	    if(finite_thick) finite_thickness_chi(rho,duedx[0],&chi,&K_bl);
            else K_bl=1.0;   
	    skin_and_heat(xit2,duedx[0],jstart,ts,us,rho,hsp,J,&qw,&cf,&St,rhoinf,uinf,hstag,nspec,&qwc,&qwd,Jwall);
            write_solution(eta,rho,Mw,u,h,t,c,qw,cf,St,xi,x_wr,dcsi,it,jstart,duedx[0],qwc,qwd,M,Jwall,J,V,iwr);
            if((eq_flg==3)&&(stef)){ for(i=0;i<neq;i++){for(j=0;j<nspec;j++)sum_ci_stef[i]=sum_ci_stef[i]+c[j][i];} }
	    if (eq_flg==3)
	      {
	      for(i=0;i<neq;i++)
	        {
		for(j=0;j<n_elements;j++)
		  {
		   tmp[j]=0.0; 
		   for(k=0;k<nspec;k++){tmp[j]=tmp[j]+lambda_matrx[k][j]*J[k][i]/Mw[k];}
		  }
	        }
	      }
	  
	  if (eq_flg==3)
	     {
	     for(i=0;i<neq;i++)	{for(j=0;j<nspec;j++)sum_ci_stef[i]=sum_ci_stef[i]+c[j][i];}  
	     }
	  iwr++;
         }
//////////////

         if(jstart==0)
         {
            if(only_stag) {
			      printf("\n NEBOULA has finished after %d iterations to a residual of:%2.2e \n",it,conv);
	    		exit(1);
			}
            xi = csi[1];		dcsi = xi;
            xit2 = 2.0*xi;
            for(i=0;i<neq;i++)
            {
               u2[i]=u[i];		u1[i]=u[i];
               h2[i]=h[i];		h1[i]=h[i];
               for(j=0;j<nspec;j++)
               {
                  c2[j][i]=c[j][i];	c1[j][i]=c[j][i];
                }
             }
            jstart=1;
          }
         else if(jstart==1)
         {
            stepsize_control(it,8,12,&dcsi);
            new_coordinate(dcsi,&xi,&xi1,&xi2,&lam0,&lam1,&lam2,&xit2);
            for(i=0;i<neq;i++)
            {
               u2[i]=u1[i];		u1[i]=u[i];
               h2[i]=h1[i];		h1[i]=h[i];
               for(j=0;j<nspec;j++)
               {
                  c2[j][i]=c1[j][i];	c1[j][i]=c[j][i];
               }
            }
            jstart=2;
         }
         else
         {
            stepsize_control(it,8,12,&dcsi);
            new_coordinate(dcsi,&xi,&xi1,&xi2,&lam0,&lam1,&lam2,&xit2);
            for(i=0;i<neq;i++)
            {
               u2[i]=u1[i];		u1[i]=u[i];
               h2[i]=h1[i];		h1[i]=h[i];
               for(j=0;j<nspec;j++)
               {
                  c2[j][i]=c1[j][i];	c1[j][i]=c[j][i];
               }
            }
         }

      }
   //----------------------
      else //to many iterations
      {
         if(jstart==0);
         else if(jstart==1)
         {
           /* printf("\n Reduce stepsize jstart=1 \n"); */
            reduce_stepsize(&dcsi,&xi);
            for(i=0;i<neq;i++)
            {
               u[i]=u1[i];
               h[i]=h1[i];
               for(j=0;j<nspec;j++)
               {
                  c[j][i]=c1[j][i];
               }
            }
            lam0 = 1.0/dcsi;lam1 = -1.0/dcsi;lam2 = 0.0;
            xit2=2.0*xi;
	 }
         else
         {
            reduce_stepsize(&dcsi,&xi);
            for(i=0;i<neq;i++)
            {
               u[i]=u1[i];
               h[i]=h1[i];
               for(j=0;j<nspec;j++)
               {
                  c[j][i]=c1[j][i];
               }
            }
            lam0 = 1.0/(xi-xi1)+1.0/(xi-xi2);lam1 = (xi-xi2)/((xi1-xi)*(xi1-xi2));lam2 = (xi-xi1)/((xi2-xi)*(xi2-xi1));
            xit2=2.0*xi;
         }
      }
         }  
   /***********************************************************************/	 
	 /*****  End of outer loop: integration in csi direction  *****/
   /***********************************************************************/
    printf("\n-------------------------------\n");
    printf("VKI Boundary Layer Code has finished\n") ;   
    printf("-------------------------------\n");
   }


#undef WL
