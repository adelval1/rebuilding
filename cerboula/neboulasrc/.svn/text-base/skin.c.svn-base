#include "externalvariables.h"
#include <math.h>

void skin_and_heat(double xit2,double duedx,int jstart,double *ts,double *us,
                   double *rho,
                   double **hsp,double **J,double *qw,double *cf,double *St,
                   double rhoinf,double uinf,double hstag,int nspec,
                   double *qwc,double *qwd,double *Jwall/*NANNI*/)
{
   extern int axisymm,two_d,cone,flat_plate;
   extern double K_bl;

   int i;
   double aux,dummy1,dummy2;

   if(jstart==0)
   {
      if(axisymm)
         if(cone)
         {
            *cf = 0.0;  /* Are fictious values because in reality it goes to infinity */
            *St = 0.0; *qw=0.0; *qwc=0.0; *qwd=0.0;
         }
         else
         {
            aux=sqrt(2.0*duedx/(rhodelta*mudelta))*rho[0]*K_bl;
            dummy1=0.0;
            for(i=0;i<nspec;i++) dummy1=dummy1+J[i][0]*hsp[i][0];
/*            printf("\n ***** kwall=%e, dummy=%e, dtdn=%e, aux*dtdy=%e, rho[0]=%e ***\n"
                   ,kwall,dummy1,ts[0],aux*ts[0],rho[0]);
            printf("\n ***** rhodelta=%e, mudelta=%e, duedx=%e *** \n",
                    rhodelta,mudelta,duedx); */
                    
            *qw=-aux*kwall*ts[0]+dummy1;
            *qw=-*qw;
            *cf = 0.0;
            *qwc = aux*kwall*ts[0];
            *qwd = -dummy1;
	    dummy2=0;
	    for(i=0;i<nspec;i++) dummy2=dummy2+Jwall[i]*hsp[i][0];
//	    printf("\nqwd %le, qwc=%le, qw=%le \n",*qwd,*qwc,*qw);/*NANNI*/
            /*for (i=0;i<nspec;i++) printf("hsp[%d]=%le",i,hsp[i]);*/
	    /* *St = *qw/(rhoinf*uinf*(hstag-hwall));*/
            *St = *qw/(0.5*rhoinf*uinf*uinf*uinf);
         }

      if(two_d)
         if(flat_plate)
         {
            *cf = 0.0;  /* Are fictious values because in reality it goes to infinity */
            *St = 0.0; *qw=0.0; *qwc=0.0; *qwd=0.0;
         }
         else
         {
            aux=sqrt(duedx/(rhodelta*mudelta))*rho[0];
            dummy1=0.0;
            for(i=0;i<nspec;i++) dummy1=dummy1+J[i][0]*hsp[i][0];
            *qw=-aux*kwall*ts[0]+dummy1;
            *qw=-*qw;
            *cf = 0.0;
            *qwc = aux*kwall*ts[0];
            *qwd = -dummy1;
 	    /* *St = *qw/(rhoinf*uinf*(hstag-hwall)); */
            *St = *qw/(0.5*rhoinf*uinf*uinf*uinf);
         }
   }
   else
   {
      aux=rho[0]*rdelta*udelta/sqrt(xit2);
      dummy1=0.0;
      for(i=0;i<nspec;i++) dummy1=dummy1+J[i][0]*hsp[i][0];
      *qw=-aux*kwall*ts[0]+dummy1;
      *qw=-*qw;
      *qwc = aux*kwall*ts[0];
      *qwd = -dummy1;
      *cf=rho[0]*muwall*rdelta*udelta*udelta/sqrt(xit2)*us[0]*2./
          rhoinf/uinf/uinf;
      /* *St = *qw/(rhoinf*uinf*(hstag-hwall)); */
      *St = *qw/(0.5*rhoinf*uinf*uinf*uinf);
   }
}
