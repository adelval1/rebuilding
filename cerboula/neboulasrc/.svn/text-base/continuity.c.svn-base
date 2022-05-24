#include <math.h>


/********************************************************************
     subroutine continuity
     Integrating the continuity equation, by mean of
     Simpson rule
     Returns: V, dimensionless normal velocity
********************************************************************/

void continuity(double *y,double *result)
{
   extern int neq;
   extern double deta;

   int i;


   result[0] = 0.0;

   result[1] = (17.0*y[0]+42.0*y[1]-16.0*y[2]+6.0*y[3]-y[4])*deta/48.0;

   for(i=2;i<neq;i++)
   {
      result[i] = result[i-2]+(y[i-2]+4.0*y[i-1]+y[i])*deta/3.0;
   }
}


/********************************************************************
   subroutine finite_thickness_chi
   computes the normalisation factor used when the real thickness
   of the boundary layer at the stagnation point is taken into account
********************************************************************/

void finite_thickness_chi(double *rho,double duedx,double *chi,double *K_bl)
{
   extern int neq;
   extern double rhodelta,mudelta;
   extern double delta_bl;  /* boundary layer real thickness */
   double *vector(int);
   void free_vector(double *);
   void continuity(double *,double *);

   int i;
   double fact,*yc,*rhotmp;

   yc=vector(neq);
   rhotmp=vector(neq);

   fact = sqrt( (rhodelta*mudelta)/(2.0*duedx) );
   for(i=0;i<neq;i++) rhotmp[i]=1.0/rho[i];
   continuity(rhotmp,yc);

   *chi = yc[neq-1];
   *K_bl = fact*(*chi)/delta_bl;

   free_vector(yc); free_vector(rhotmp);
}
