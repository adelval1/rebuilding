#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define TINY 1.0e-20
double *vector(int);
double **matrix(int,int);
double ***trimatrix(int,int,int);
int *intvector(int);
void nerror(char *);
void free_vector(double *);
void free_intvector(int *);
void free_matrix(double **,int);
void free_trimatrix(double ***,int,int);
void inv(double **,double **,int);
void inv_prod(double **,double **,double **,double *,double *,int);
void inv_vec(double **,double *,double *,int);
void ludcmp(double **,int,int *,double *);
void lubksb(double **,int,int *,double *);


void trisol(double ***ca,double ***cb,double ***cc,double **cr,int nstart,int neq,
            int nspec)
{
   int i,j,k,l;
   double ***beta,**cy,d;
   
   beta=trimatrix(neq,nspec,nspec);
   cy=matrix(neq,nspec);
   
   inv_prod(cb[nstart],cc[nstart],beta[nstart],cr[nstart],cy[nstart],nspec);
   
   
   for(i=nstart+1;i<neq;i++)
   {
      
      for(j=0;j<nspec;j++)
      {
         for(k=0;k<nspec;k++)
         {
            d=0.0;
            for(l=0;l<nspec;l++)
            {
               d=d+ca[i][j][l]*beta[i-1][l][k];
            }
            cb[i][j][k]=cb[i][j][k]-d;
         }
      }


      for(k=0;k<nspec;k++)
      {
         d=0.0;
         for(l=0;l<nspec;l++)
         {
            d=d+ca[i][k][l]*cy[i-1][l];
         }
         cr[i][k]=cr[i][k]-d;
      }

      
      if(i != (neq-1))
      {
         inv_prod(cb[i],cc[i],beta[i],cr[i],cy[i],nspec);
      }
      else
      {
         inv_vec(cb[i],cr[i],cy[i],nspec);
      }
   }
   
   for(k=0;k<nspec;k++)
   {
      cr[neq-1][k]=cy[neq-1][k];
      /*printf("cr[%d][%d]=%le              cy[%d][%d]=%le \n",
              neq-1,k,cr[neq-1][k],neq-1,k,cy[neq-1][k]);/*NANNI*/
   }
   for(i=neq-2;i>=nstart;i--)
   {
      for(j=0;j<nspec;j++)
      {
         d=0.0;
         for(k=0;k<nspec;k++)
         {
            d=d+beta[i][j][k]*cr[i+1][k];
         }
         cr[i][j]=cy[i][j]-d;
      }
   }
   
   free_trimatrix(beta,neq,nspec);
   free_matrix(cy,neq);
}


void ludcmp(double **a,int n,int *indx,double *d)
{
   int i,imax,j,k;
   double big,dum,sum,temp,*vv;
   vv=vector(n);
   *d=1.0;
   for(i=0;i<n;i++)
   {
      big=0.0;
      for(j=0;j<n;j++)
         if((temp=fabs(a[i][j])) > big) big=temp;
      if(big == 0.0) nerror("Singular matrix in routine LUDCMP");
      vv[i]=1.0/big;
   }
   for(j=0;j<n;j++)
   {
      for(i=0;i<j;i++)
      {
         sum=a[i][j];
         for(k=0;k<i;k++) sum -= a[i][k]*a[k][j];
         a[i][j]=sum;
      }
      big=0.0;
      for(i=j;i<n;i++)
      {
         sum=a[i][j];
         for (k=0;k<j;k++) sum -= a[i][k]*a[k][j];
         a[i][j]=sum;
         if ( (dum=vv[i]*fabs(sum)) >= big)
         {
            big=dum;
            imax=i;
         }
      }
      if(j != imax)
      {
         for(k=0;k<n;k++)
         {
            dum=a[imax][k];
            a[imax][k]=a[j][k];
            a[j][k]=dum;
         }
         *d = -(*d);
         vv[imax]=vv[j];
      }
      indx[j]=imax;
      if (a[j][j] == 0.0) a[j][j]=TINY;
      if (j != (n-1))
      {
         dum=1.0/(a[j][j]);
         for (i=j+1;i<n;i++) a[i][j] *= dum;
      }
   }
   free_vector(vv);
}


void lubksb(double **a,int n,int *indx,double *b)
{
   int i,ii=-1,ip,j;
   double sum;

   for (i=0;i<n;i++)
   {
      ip=indx[i];
      sum=b[ip];
      b[ip]=b[i];
      if(ii!=-1) for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
      else if (sum) ii=i;
      b[i]=sum;
   }
   for(i=n-1;i>=0;i--)
   {
      sum=b[i];
      for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
      b[i]=sum/a[i][i];
   }
}


void inv(double **a,double **y,int n)
{
   int i,j,*indx;
   double d,*col;
   col=vector(n);
   indx=intvector(n);
   
   ludcmp(a,n,indx,&d);
   for(j=0;j<n;j++)
   {
      for(i=0;i<n;i++) col[i]=0.0;
      col[j]=1.0;
      lubksb(a,n,indx,col);
      for(i=0;i<n;i++) y[i][j]=col[i];
   }
   free_vector(col);
   free_intvector(indx);
}



void inv_prod(double **a,double **b,double **y,double *r,double *x,int n)
{                     /* computes the product a^-1*b */
   int i,j,*indx;
   double d,*col;
   col=vector(n);
   indx=intvector(n);
   
   ludcmp(a,n,indx,&d);
   for(j=0;j<n;j++)
   {
      for(i=0;i<n;i++) col[i]=b[i][j];
      lubksb(a,n,indx,col);
      for(i=0;i<n;i++) y[i][j]=col[i];
   }
   for(i=0;i<n;i++) col[i]=r[i];
   lubksb(a,n,indx,col);
   for(i=0;i<n;i++) x[i]=col[i];
   free_vector(col);
   free_intvector(indx);
}



void inv_vec(double **a,double *b,double *y,int n)
{
   int i,*indx;
   double d,*col;
   col=vector(n);
   indx=intvector(n);
   
   ludcmp(a,n,indx,&d);
   for(i=0;i<n;i++) col[i]=b[i];
   lubksb(a,n,indx,col);
   for(i=0;i<n;i++) y[i]=col[i];

   free_vector(col);
   free_intvector(indx);
}


void invbis(double **a,double **y,int n,int *ind)
{
   int i,j,*indx;
   double d,*col;
   col=vector(n);
   indx=intvector(n);
   
   ludcmp(a,n,indx,&d);
   for(i=0;i<n;i++)
   {
      ind[i]=indx[i];
   }
   for(j=0;j<n;j++)
   {
      for(i=0;i<n;i++) col[i]=0.0;
      col[j]=1.0;
      lubksb(a,n,indx,col);
      for(i=0;i<n;i++) y[i][j]=col[i];
   }
   free_vector(col);
   free_intvector(indx);
}



void thomas(int ii,double *a,double *b,double *c,double *d,double *result)
{
   int i;
   
/* Hirsch version   */

   b[0] = 1.0/b[0];
   a[0] = d[0]*b[0];
   
   for(i=1;i < ii-1;i++)
   {
      c[i-1] = c[i-1]*b[i-1];
      b[i] = b[i]-a[i]*c[i-1];
      b[i]=1.0/b[i];
      a[i]=(d[i]-a[i]*a[i-1])*b[i];
   }

   result[ii-2] = a[ii-2];
   for(i=ii-3;i >= 0;i--)
   {
      result[i] = a[i]-c[i]*result[i+1];
   }
}

/*c  Fletcher version
c      c(1)=c(1)/b(1)
c      d(1)=d(1)/b(1)
c      do i=2,ngm1
c         c(i)=c(i)/(b(i)-a(i)*c(i-1))
c         d(i)=(d(i)-a(i)*d(i-1))/(b(i)-a(i)*c(i-1))
c      end do
c      result(ii+ngm1)=d(ngm1)
c      do i=(ngm1-1),1,-1
c         result(ii+i)=d(i)-result(ii+i+1)*c(i)
c      end do

c Anderson,Tannehill version

c      do i=2,ngm1
c         b(i)=b(i)-a(i)/b(i-1)*c(i-1)
c         d(i)=d(i)-a(i)/b(i-1)*d(i-1)
c      end do
c      result(ii+ngm1)=d(ngm1)/b(ngm1)
c      do i=(ngm1-1),1,-1
c         result(ii+i)=(d(i)-c(i)*result(ii+i+1))/b(i)
c      end do
c      return
      end   */


void matmat(double **a,int ia,int ja,double **b,int ib,int jb,double **c)
{
   int i,j,k;
   double tmp;

   if(ja!=ib)
   {
      printf("\n Error in matmat ja=%d and ib=%d are not equal \n",ja,ib);
      exit(1);
   }

   for(i=0;i<ia;i++)
   {
      for(j=0;j<jb;j++)
      {
         tmp=0.0;
         for(k=0;k<ja;k++) tmp=tmp+a[i][k]*b[k][j];
         c[i][j]=tmp;
      }
   }

}



void matvec(double **a,int ia,int ja,double *b,int ib,double *c)
{
   int i,j;
   double tmp;

   if(ja!=ib)
   {
      printf("\n Error in matvec ja=%d and ib=%d are not equal \n",ja,ib);
      exit(1);
   }

   for(i=0;i<ia;i++)
   {
      tmp=0.0;
      for(j=0;j<ja;j++)
         tmp=tmp+a[i][j]*b[j];

      c[i]=tmp;
   }

}


void invp(double **a,double **b,double **y,int *indx,int n)
{                     /* computes the product a^-1*b */
   int i,j;
   double *col;
   col=vector(n);
   
   for(j=0;j<n;j++)
   {
      for(i=0;i<n;i++) col[i]=b[i][j];
      lubksb(a,n,indx,col);
      for(i=0;i<n;i++) y[i][j]=col[i];
   }

   free_vector(col);
}



void invv(double **a,double *b,double *y,int *indx,int n)
{
   int i;
   double *col;
   col=vector(n);
   
   for(i=0;i<n;i++) col[i]=b[i];
   lubksb(a,n,indx,col);
   for(i=0;i<n;i++) y[i]=col[i];

   free_vector(col);
}




#undef TINY
