#include "externalvariables.h"

void spline(double *,double *,int,double,double,double *);
void splint(double *,double *,double *,double,double *,int,int);


/*--------------------------------------------------------------
c subroutine grid_gen
c-----------------------------------------------------------*/

void grid_gen(double *xgrid,double *rhoext,double *muext,double *uext,double *text,
              double *rbody,double *csi,double *dcsidx,double *duedx,double *dtedx)
{
   int i;
   double dx;
   void ableit(double *,double *,double,int);
               
   if(only_stag)
   {
      csi[0] = 0.0;
      dcsidx[0] = rhoext[0]*muext[0]*uext[0]*rbody[0]*rbody[0];
   }
   else
   {
      dx = xgrid[1]-xgrid[0];
        
      csi[0] = 0.0;
      csi[1] = (17.0*rhoext[0]*muext[0]*uext[0]*rbody[0]*rbody[0]+
                42.0*rhoext[1]*muext[1]*uext[1]*rbody[1]*rbody[1]-
                16.0*rhoext[2]*muext[2]*uext[2]*rbody[2]*rbody[2]+
                6.0*rhoext[3]*muext[3]*uext[3]*rbody[3]*rbody[3]-
                rhoext[4]*muext[4]*uext[4]*rbody[4]*rbody[4])*dx/48.0;
        
      for(i=2;i < ngrid;i++)
      {
         csi[i] = csi[i-2]+
                  (rhoext[i-2]*muext[i-2]*uext[i-2]*rbody[i-2]*rbody[i-2]+
                   4.0*rhoext[i-1]*muext[i-1]*uext[i-1]*rbody[i-1]*rbody[i-1]+
                   rhoext[i]*muext[i]*uext[i]*rbody[i]*rbody[i])*dx/3.0;
      }

      for(i=0;i < ngrid;i++)
      {
         dcsidx[i] = rhoext[i]*muext[i]*uext[i]*rbody[i]*rbody[i];
      }

      ableit(uext,duedx,dx,ngrid);
      ableit(text,dtedx,dx,ngrid);
   }

}

/*----------------------------------------------------------------------------------
 subroutine spline_gen
------------------------------------------------------------------------------*/
void spline_gen(double *xgrid,double *rbody,double *uext,double *text,
                double *duedx,double *dtedx,double *dcsidx,double **cext,
                double *pext,double *muext,double *srbody,double *suext,
                double *stext,double *sduedx,double *sdtedx,double *sdcsidx,
                double **scext,double *spext,double *smuext,double *rhoext
		,double *srhoext,double *twall, double *stwall,double **gam,double **sgam)
{
   int i;
   if(only_stag)
   {
   }
   else
   {
      spline(xgrid,rbody,ngrid,2.e30,2.e30,srbody);
      spline(xgrid,uext,ngrid,2.e30,2.e30,suext);
      spline(xgrid,text,ngrid,2.e30,2.e30,stext);
      spline(xgrid,duedx,ngrid,2.e30,2.e30,sduedx);
      spline(xgrid,dtedx,ngrid,2.e30,2.e30,sdtedx);
      spline(xgrid,dcsidx,ngrid,2.e30,2.e30,sdcsidx);
      for(i=0;i<nspec;i++){spline(xgrid,cext[i],ngrid,2.e30,2.e30,scext[i]);}
      spline(xgrid,pext,ngrid,2.e30,2.e30,spext);
      spline(xgrid,muext,ngrid,2.e30,2.e30,smuext);
      spline(xgrid,rhoext,ngrid,2.e30,2.e30,srhoext);
      spline(xgrid,twall,ngrid,2.e30,2.e30,stwall);
      for(i=0;i<nrwall;i++){spline(xgrid,gam[i],ngrid,2.e30,2.e30,sgam[i]);}

   }
}
/*---------------------------------------------------------------------
  subroutine search
---------------------------------------------------------------------*/
void search(double xi,double *csi,int *index1,int *index2)
{
   int i;
   for(i=0;i < ngrid; i++)
   {
      if(xi > csi[i]){*index1=i;}
      else if(xi < csi[i]){*index2=i;return;}
      else if(xi == csi[i]){*index1=i;*index2=i;return;}
   }
}
/*------------------------------------------------------------------------------
  subroutine var_interp
-----------------------------------------------------------------------------*/
void var_interp(double *csi,double *rbody,double *uext,double *hext,double *duedx,
                double *dhdx,double *dcsidx,double **cext,double *pext,double *muext,
                double *srbody,double *suext,double *shext,double *sduedx,double *sdhdx,
                double *sdcsidx,double **scext,double *spext,double *smuext,
                double xi,double *xgrid,double xit2,double *rhoext, double *srhoext,double
		*twall,double *stwall,double **gam,double **sgam,double *rhodelta2)
{
   int i,index1,index2;
   double xix;
   void search(double,double *,int *,int *);

   search(xi,csi,&index1,&index2);
   if(index1 != index2)
   {
      xix=(csi[index2]-xi)/(csi[index2]-csi[index1])*xgrid[index1]+(xi-csi[index1])/(csi[index2]-csi[index1])*xgrid[index2];
      splint(xgrid,uext,suext,xix,&udelta,index1,index2);
      splint(xgrid,hext,shext,xix,&hdelta,index1,index2);
      for(i=0;i<nspec;i++){splint(xgrid,cext[i],scext[i],xix,&cdelta[i],index1,index2);}
      splint(xgrid,pext,spext,xix,&pdelta,index1,index2);
      splint(xgrid,twall,stwall,xix,&twalldelta,index1,index2);
      splint(xgrid,duedx,sduedx,xix,&duedelta,index1,index2);
      splint(xgrid,dhdx,sdhdx,xix,&dhdelta,index1,index2);
      splint(xgrid,dcsidx,sdcsidx,xix,&dcsidelta,index1,index2);
      duedelta = duedelta/dcsidelta;
      dhdelta = dhdelta/dcsidelta;
      splint(xgrid,rbody,srbody,xix,&rdelta,index1,index2);
      beta = xit2*duedelta/udelta;
      splint(xgrid,muext,smuext,xix,&mudelta,index1,index2);
      splint(xgrid,rhoext,srhoext,xix,rhodelta2,index1,index2);
      for(i=0;i<nrwall;i++){splint(xgrid,gam[i],sgam[i],xix,&gamdelta[i],index1,index2);}
   }
   else
   {
      rdelta=rbody[index1];
      udelta=uext[index1];
      hdelta=hext[index1];
      for(i=0;i<nspec;i++) cdelta[i]=cext[i][index1];
      pdelta=pext[index1];
      twalldelta=twall[index1];
      duedelta=duedx[index1]/dcsidx[index1];
      dhdelta=dhdx[index1]/dcsidx[index1];
      dcsidelta=dcsidx[index1];
      mudelta=muext[index1];
      beta=xit2*duedelta/udelta;
      *rhodelta2 = rhoext[index1];
      for(i=0;i<nrwall;i++) gamdelta[i]=gam[i][index1];

   }

}
/*------------------------------------------------------      
c     subroutine new_coordinates
c     Calculates the new coordinates and other geometrical parameters
c------------------------------------------------------  */
void new_coordinate(double dxi,double *xi,double *xi1,double *xi2,double *lam0,double *lam1,double *lam2,double *xit2)
{
   *xi2 = *xi1;
   *xi1 = *xi;
   *xi = *xi + dxi;

   *lam0 = 1.0/(*xi - *xi1)+1.0/(*xi - *xi2);
   *lam1 = (*xi - *xi2)/((*xi1 - *xi)*(*xi1 - *xi2));
   *lam2 = (*xi - *xi1)/((*xi2 - *xi)*(*xi2 - *xi1));

   *xit2 = 2.0*(*xi);
}

/*------------------------------------------------------      
c     subroutine reduce_stepsize
c     Reduces the stepsize delta xi if there is no
c     convergence. The reduction is 50%
c------------------------------------------------------ */

void reduce_stepsize(double *dxi,double *xi)
{
   *xi = *xi - *dxi;
   *dxi = 0.5*(*dxi);
   *xi = *xi + *dxi;
}


/*-------------------------------------------------------
c     subroutine stepsize_control
c     Adjusts the stepsize by means of number of 
c     iterations
c-------------------------------------------------------*/

void stepsize_control(int nn,int lower,int upper,double *dxi)
{
   if(nn<=lower)
   {
      *dxi = *dxi*(1.0+1.0/(double)nn);
   }
   else
   {
      if(nn>=upper)
      {
         *dxi = *dxi*(double)upper/(double)nn;
      }
   }
}


/**********************************************************
  subroutine make_csider
    computes the csi derivative of the unknown quantities
***********************************************************/

void make_csider(double xit2,double lam1,double lam2,double *u1,double *u2,
                 double *up1,double *h1,double *h2,double *hp1,double **c1,
                 double **c2,double **cp1)
{
   int i,j;
   double dummy;

   /*printf("\n xit2*lam1=%e xit2*lam2=%e\n",xit2*lam1,xit2*lam2);*/
   for(i=0;i<neq;i++)
   {
      up1[i] = (xit2*lam1)*u1[i]+(xit2*lam2)*u2[i];
      hp1[i] = (xit2*lam1)*h1[i]+(xit2*lam2)*h2[i];
      dummy = 0.0;
      for(j=0;j<nspec;j++)
      {
         cp1[j][i] = (xit2*lam1)*c1[j][i]+(xit2*lam2)*c2[j][i];
         /*dummy=dummy+cp1[j][i];
         printf(" %e  %e  %e  %e",c1[j][i],(xit2*lam1)*c1[j][i],(xit2*lam2)*c2[j][i],
                                 cp1[j][i]);*/
      }
      /*printf(" %e \n",dummy);*/
   }
}
