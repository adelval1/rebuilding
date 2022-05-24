#include <stdio.h>
#include "externalvariables.h"

void msq_der(double **c,double **cs,double **css)
{
   int i,j;
   double tmp1,tmp2,tmp3;

   for(i=0;i<neq;i++)
   {
   /* Computes the sum of the ci and of the csi and of the cssi */

      tmp1 = 0.0;
      tmp2 = 0.0;
      tmp3 = 0.0;
      for(j=0;j<nspec;j++)
      {
         tmp1 = tmp1+c[j][i];
         tmp2 = tmp2+cs[j][i];
         tmp3 = tmp3+css[j][i];
      }

    /* Computes the corrected values of csi and cssi */

      for(j=0;j<nspec;j++)
      {
         cs[j][i] = cs[j][i]-c[j][i]/tmp1*tmp2;
         css[j][i] = css[j][i]-c[j][i]/tmp1*tmp3;
      }

   }


}
