#include "stdio.h"
#include "stdlib.h"
/*#include "strdef.h"*/

void nerror(char *text)
{
   printf("%s\n",text);
   exit(1);
}


double *vector(int dim)
{
   double *v;
   v = (double*) calloc((size_t)dim,(size_t)sizeof(double));
   if(!v) nerror("Dynamic allocation failed in vector()");
   return v;
}



double **matrix(int isize,int jsize) /*NANNI*/
{
   int i;
   double **m;
   m = (double **) calloc((size_t)isize,(size_t)sizeof(double*));
   if(!m) nerror("Dynamic allocation failure in matrix()");
   for(i=0;i<isize;i++)
   {
      m[i] = (double*) calloc((size_t)jsize,(size_t)sizeof(double));
      if(!m[i]) nerror("Dynamic allocation failure in matrix()");
   }
   return m;
}


int **matrix_int(int isize,int jsize)
{
   int i;
   int **m;
   m = (int **) calloc((size_t)isize,(size_t)sizeof(int*));
   if(!m) nerror("Dynamic allocation failure in matrix()");
   for(i=0;i<isize;i++)
   {
      m[i] = (int*) calloc((size_t)jsize,(size_t)sizeof(int));
      if(!m[i]) nerror("Dynamic allocation failure in matrix()");
   }
   return m;
}



double ***trimatrix(int isize,int jsize,int ksize)
{
   int i,j;
   double ***m;
   m = (double ***) calloc((size_t)isize,(size_t)sizeof(double**));
   if(!m) nerror("Dynamic allocation failure in trimatrix()");
   for(i=0;i<isize;i++)
   {
      m[i] = (double **) calloc((size_t)jsize,(size_t)sizeof(double*));
      if(!m[i]) nerror("Dynamic allocation failure in trimatrix()");
   }
   for(i=0;i<isize;i++)
   {
      for(j=0;j<jsize;j++)
      {
         m[i][j] = (double *) calloc((size_t)ksize,(size_t)sizeof(double));
         if(!m[i][j]) nerror("Dynamic allocation failure in trimatrix()");
      }
   }
   return m;
}



int *intvector(int dim)
{
   int *v;
   v = (int*) calloc((size_t)dim,(size_t)sizeof(int));
   if(!v) nerror("Dynamic allocation failed in intrvector()");
   return v;
}


/*struct chem *struc_alloc(int dim)
{
   struct chem *v;
   v = (struct chem *) calloc((size_t)dim,(size_t)sizeof(struct chem));
   if(!v) nerror("Dinamic allocation failed in struc_alloc()");
   return v;
}*/




void free_vector(double *v)
{
   free(v);
}


 
void free_matrix(double **m,int isize)
{
   int i;
   for(i=isize-1;i>=0;i--) free(m[i]);
   free(m);
}



void free_trimatrix(double ***m,int isize,int jsize)
{
   int i,j;
   for(i=isize-1;i>=0;i--)
   {
      for(j=jsize-1;j>=0;j--) free(m[i][j]);
   }
   for(i=isize-1;i>=0;i--) free(m[i]);
   free(m);
}


void free_intvector(int *v)
{
   free(v);
}



char *stringa(int dim)
{
   char *v;
   v = (char*) calloc(dim,(size_t)sizeof(char));
   if(!v) nerror("Dynamic allocation failed in stringa()");
   return v;
}
