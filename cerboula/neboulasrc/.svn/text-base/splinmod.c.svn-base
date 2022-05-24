#include <stdio.h>
#include <stdlib.h>
extern double *vector(int);
extern void free_vector(double *);


/*-----------------------------------------------------------------------------
     subroutine spline
-----------------------------------------------------------------------------*/

void spline(double *x,double *y,int n,double yp1,double ypn,double *y2)
{
   int i;
   double p,qn,sig,un,*u;
   u=vector(n);
   if(!u)
   {
      printf("Dynamic allocation failed in spline");
      exit(1);
   }
   if(yp1 > .99e30) y2[0]=u[0]=0.0;
   else
   {
      y2[0]=-0.5;
      u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
   }

   for(i=1;i < n-1; i++)
   {
      sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
      p=sig*y2[i-1]+2.0;
      y2[i]=(sig-1.0)/p;
      u[i]=(6.0*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]))/
	   (x[i+1]-x[i-1])-sig*u[i-1])/p;
   }

   if(ypn > .99e30) qn=un=0.0;
   else
   {
      qn=0.5;
      un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
   }
   y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
   for(i=n-2;i>=0;i--)
   {
      y2[i]=y2[i]*y2[i+1]+u[i];
   }
   free_vector(u);
}


/*-----------------------------------------------------------------------------
      subroutine splint
-----------------------------------------------------------------------------*/

void splint(double *xa,double *ya,double *y2a,double x,double *y,int klo,int khi)
{
   double a,b,h;

   h=xa[khi]-xa[klo];
   if(h==0.0){ printf("bad xa input in splint");exit(1);}
   a=(xa[khi]-x)/h;
   b=(x-xa[klo])/h;
   *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;     
}
