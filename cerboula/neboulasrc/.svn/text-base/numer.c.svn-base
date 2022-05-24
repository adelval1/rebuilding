#include <math.h>
#define max(x,y)  (((x) >= (y)) ? (x) : (y))
#include <stdio.h> 

void thomas(int,double *,double *,double *,double *,double *);
double *vector(int);
double **matrix(int,int);
double ***trimatrix(int,int,int);
int *intvector(int);
void free_vector(double *);
void free_intvector(int *);
void free_matrix(double **,int);
void free_trimatrix(double ***,int,int);
void trisol(double ***,double ***,double ***,double **,int,int,int);
void lubksb(double **,int,int *,double *);
void invbis(double **,double **,int,int *);
void norma(double ***,double ***,double ***,int,int,int);


/*--------------------------------------------------------
      Subroutine ableit
      Calculates the derivative of bb
 -------------------------------------------------------*/

void ableit(double *aa,double *bb,double deta,int nn)
{
   int i;
   double deta12;
   deta12 = deta*12.0;
   
   
   for(i=2;i<nn-2;i++)
   {
      bb[i]=(aa[i-2]-8.0*aa[i-1]+8.0*aa[i+1]-aa[i+2])/deta12;
   }

   bb[0]=(-25.0*aa[0]+48.0*aa[1]-36.0*aa[2]+16.0*aa[3]-
            3.0*aa[4])/deta12;

   bb[1]=(-3.0*aa[0]-10.0*aa[1]+18.0*aa[2]-6.0*aa[3]+
             aa[4])/deta12;

   bb[nn-2]=(3.0*aa[nn-6]-16.0*aa[nn-5]+36.0*aa[nn-4]
             -48.0*aa[nn-3]+25.0*aa[nn-2])/deta12;

   bb[nn-1]=(3.0*aa[nn-5]-16.0*aa[nn-4]+36.0*aa[nn-3]
             -48.0*aa[nn-2]+25.0*aa[nn-1])/deta12;

}



void tridiag_inv(int ii,double *wrhs,double f1,double g1,double h1,
                 double *a,double *b,double *c,double *d,double *F,double *dF,
                 double *result)
{
   int i;
   double delnp1,deln,delnm1,r,r0,r1,r2;
   double *A,*B,*C,*D;
   double mnp1m1,mnp1al,mnp1be,mnp1,mn0,mnm1,mrp1,mr0,mrm1,mral,mrbe;
   double mnal,mnbe,mnm1p1,mnm10,mnm1m1,mnm1al,mnm1be,mnp1p1,mnp10;
   double m1p1,m10,m1m1,m1al,m1be,m2p1,m20,m2m1,m2al,m2be;

   A=vector(ii); B=vector(ii); C=vector(ii); D=vector(ii);

/*    Coefficients of the system         */
/* At the interior points */

   for(i=1;i < ii-1;i++)
   {
      mnp1p1 = a[i+1]+1.5*b[i+1]+c[i+1];
      mnp10 = -(2.0*a[i+1]+2.0*b[i+1]);
      mnp1m1 = a[i+1]+.5*b[i+1];
      mnp1al = -(6.0*a[i+1]+2.0*b[i+1]);
      mnp1be = -(10.0*a[i+1]+2.0*b[i+1]);
      mnp1 = a[i]+.5*b[i];
      mn0 = -2.0*a[i]+c[i];
      mnm1 = a[i]-.5*b[i];
      mnal = b[i];
      mnbe = 2.0*a[i];
      mnm1p1 = a[i-1]-.5*b[i-1];
      mnm10 = 2.0*(-a[i-1]+b[i-1]);
      mnm1m1 = a[i-1]-1.5*b[i-1]+c[i-1];
      mnm1al = 6.0*a[i-1]-2.0*b[i-1];
      mnm1be = -10.0*a[i-1]+2.0*b[i-1];

/*    Determinants to eliminate alfa and beta    */
      delnp1 = mnal*mnp1be-mnp1al*mnbe;
      deln = mnp1al*mnm1be-mnm1al*mnp1be;
      delnm1 = mnm1al*mnbe-mnal*mnm1be;

/*    Coefficients of the system    */
      A[i] = (mnm1m1-dF[i-1])*delnp1+mnm1*deln+mnp1m1*delnm1;
      B[i] = mnm10*delnp1+(mn0-dF[i])*deln+mnp10*delnm1;
      C[i] = mnm1p1*delnp1+mnp1*deln+(mnp1p1-dF[i+1])*delnm1;
      D[i] = (d[i-1]+F[i-1]-dF[i-1]*wrhs[i-1])*delnp1
            +(d[i]+F[i]-dF[i]*wrhs[i])*deln
            +(d[i+1]+F[i+1]-dF[i+1]*wrhs[i+1])*delnm1;
   }
 
/*   At the wall */
   mrp1=-.5*f1;
   mr0=2.0*f1;
   mrm1=-1.5*f1+g1;
   mral=-2.0*f1;
   mrbe=2.0*f1;
   m1p1=a[0]-.5*b[0];
   m10=-2.0*a[0]+2.0*b[0];
   m1m1=a[0]-1.5*b[0]+c[0];
   m1al=6.0*a[0]-2.0*b[0];
   m1be=-10.0*a[0]+2.0*b[0];
   m2p1=a[1]+.5*b[1];
   m20=-2.0*a[1]+c[1];
   m2m1=a[1]-.5*b[1];
   m2al=b[1];
   m2be=2.0*a[1];

/* Determinants to eliminate alpha and beta   */
   delnp1 = m1al*m2be-m2al*m1be;
   deln = m2al*mrbe-mral*m2be;
   delnm1 = mral*m1be-m1al*mrbe;

/* Coefficients of the matrix at the wall  */
   r2 = mrp1*delnp1+m1p1*deln+m2p1*delnm1;
   r1 = mr0*delnp1+m10*deln+(m20-dF[1])*delnm1;
   r0 = mrm1*delnp1+(m1m1-dF[0])*deln+m2m1*delnm1;
   r = h1*delnp1+(d[0]+F[0]-dF[0]*wrhs[0])*deln
      +(d[1]+F[1]-dF[1]*wrhs[1])*delnm1;

 
/*   Making tridiagonal  */
/*   D[ii-2] = D[ii-2]-C[ii-2]*wrhs;
   B[0] = r0*C[1]-r2*A[1];
   C[0] = r1*C[1]-r2*B[1];
   D[0] = r*C[1]-r2*D[1];
   A[0] = 0.0;
   C[ii-2] = 0.0; */

   D[ii-2] = D[ii-2]-C[ii-2]*wrhs[ii-1];
   B[0] = r0*C[1]-r2*A[1];
   C[0] = r1*C[1]-r2*B[1];
   D[0] = r*C[1]-r2*D[1];
   A[0] = 0.0;
   C[ii-2] = 0.0;

   thomas(ii,A,B,C,D,result);

   free_vector(A); free_vector(B); free_vector(C); free_vector(D);
}



void block_tridiag(double **a,double **b,double **c,double **d,double **F,
                   double ***dF,double *f1,double *g1,double *h1,double **dh1dc,
                   double **w,double **res,int neq,int nspec)
{
   int i,j,k,*ind,ii,kk,jj;
   double db;
   double ***A,***B,***C,**mattmp,**A1,**B1,**C1,*D1,**mt0,**mt1,*vt;
   double *tmp1,*r,**r0,**r1,**r2,**C_1,tmp2;

   double *mnp1m1,*mnp1al,*mnp1be,*mnp1,*mn0,*mnm1,*mrp1,*mr0,*mrm1,*mral,*mrbe;
   double *mnal,*mnbe,*mnm1p1,*mnm10,*mnm1m1,*mnm1al,*mnm1be,*mnp1p1,*mnp10;
   double *m1p1,*m10,*m1m1,*m1al,*m1be,*m2p1,*m20,*m2m1,*m2al,*m2be;
   double *delnp1,*deln,*delnm1;

   void matmat(double **,int,int,double **,int,int,double **);
   void matvec(double **,int,int,double *,int,double *);
   void ludcmp(double **,int,int *,double *);
   void invp(double **,double **,double **,int *,int);
   void invv(double **,double *,double *,int *,int);


   ind=intvector(nspec);
   A=trimatrix(neq,nspec,nspec); B=trimatrix(neq,nspec,nspec);
   C=trimatrix(neq,nspec,nspec);
   A1=matrix(nspec,nspec); B1=matrix(nspec,nspec); C1=matrix(nspec,nspec);
   D1=vector(nspec);
   mt0=matrix(nspec,nspec); mt1=matrix(nspec,nspec); vt=vector(nspec);
   tmp1=vector(nspec); r=vector(nspec); r0=matrix(nspec,nspec);
   r1=matrix(nspec,nspec); r2=matrix(nspec,nspec); 
   mattmp=matrix(nspec,nspec);

   mnp1m1=vector(nspec); mnp1al=vector(nspec); mnp1be=vector(nspec); mnp1=vector(nspec);
   mn0=vector(nspec); mnm1=vector(nspec); mrp1=vector(nspec); mr0=vector(nspec);
   mrm1=vector(nspec); mral=vector(nspec); mrbe=vector(nspec); mnal=vector(nspec);
   mnbe=vector(nspec); mnm1p1=vector(nspec); mnm10=vector(nspec); mnm1m1=vector(nspec);
   mnm1al=vector(nspec); mnm1be=vector(nspec); mnp1p1=vector(nspec); mnp10=vector(nspec);
   m1p1=vector(nspec); m10=vector(nspec); m1m1=vector(nspec); m1al=vector(nspec);
   m1be=vector(nspec); m2p1=vector(nspec); m20=vector(nspec); m2m1=vector(nspec);
   m2al=vector(nspec); m2be=vector(nspec);
   delnp1=vector(nspec); deln=vector(nspec); delnm1=vector(nspec);

/* Initialises the matrix of coefficients of the unknowns */
   
   for(i=0;i<nspec;i++)
   {
      r[i] = 0.0;
      for(j=0;j<nspec;j++)
      {
         r0[i][j] = 0.0;
         r1[i][j] = 0.0;
         r2[i][j] = 0.0;
      }
   }

   for(i=0;i<neq;i++)
   {
      for(j=0;j<nspec;j++)
      {
         for(k=0;k<nspec;k++)
         {
            A[i][j][k] = 0.0;
            B[i][j][k] = 0.0;
            C[i][j][k] = 0.0;
         }
      }
   }



/* At the interior points */

   for(i=1;i<neq-1;i++)
   {
      for(j=0;j<nspec;j++)
      {
         mnp1p1[j] = a[i+1][j]+1.5*b[i+1][j]+c[i+1][j];
         mnp10[j] = -(2.0*a[i+1][j]+2.0*b[i+1][j]);
         mnp1m1[j] = a[i+1][j]+.5*b[i+1][j];
         mnp1al[j] = -(6.0*a[i+1][j]+2.0*b[i+1][j]);
         mnp1be[j] = -(10.0*a[i+1][j]+2.0*b[i+1][j]);
         mnp1[j] = a[i][j]+.5*b[i][j];
         mn0[j] = -2.0*a[i][j]+c[i][j];
         mnm1[j] = a[i][j]-.5*b[i][j];
         mnal[j] = b[i][j];
         mnbe[j] = 2.0*a[i][j];
         mnm1p1[j] = a[i-1][j]-.5*b[i-1][j];
         mnm10[j] = 2.0*(-a[i-1][j]+b[i-1][j]);
         mnm1m1[j] = a[i-1][j]-1.5*b[i-1][j]+c[i-1][j];
         mnm1al[j] = 6.0*a[i-1][j]-2.0*b[i-1][j];
         mnm1be[j] = -10.0*a[i-1][j]+2.0*b[i-1][j];

/*    Determinants to eliminate alfa and beta    */

         delnp1[j] = mnal[j]*mnp1be[j]-mnp1al[j]*mnbe[j];
         deln[j] = mnp1al[j]*mnm1be[j]-mnm1al[j]*mnp1be[j];
         delnm1[j] = mnm1al[j]*mnbe[j]-mnal[j]*mnm1be[j];

/*    Coefficients of the system    */
          
         A[i][j][j] = mnm1m1[j]*delnp1[j]+mnm1[j]*deln[j]+mnp1m1[j]*delnm1[j];
         B[i][j][j] = mnm10[j]*delnp1[j]+mn0[j]*deln[j]+mnp10[j]*delnm1[j];
         C[i][j][j] = mnm1p1[j]*delnp1[j]+mnp1[j]*deln[j]+mnp1p1[j]*delnm1[j];
         res[i][j] = (d[i-1][j]+F[i-1][j])*delnp1[j]+(d[i][j]+F[i][j])*deln[j]
                    +(d[i+1][j]+F[i+1][j])*delnm1[j];
      }
      for(j=0;j<nspec;j++)
      {
         for(k=0;k<nspec;k++)
         {
            A[i][j][k] = A[i][j][k]-delnp1[j]*dF[i-1][j][k];
            B[i][j][k] = B[i][j][k]-deln[j]*dF[i][j][k];
            C[i][j][k] = C[i][j][k]-delnm1[j]*dF[i+1][j][k];
         }
      }
      for(j=0;j<nspec;j++)
      {
         tmp1[j]=0.0;
         for(k=0;k<nspec;k++)
         {
            tmp1[j]=tmp1[j]+delnm1[j]*dF[i+1][j][k]*w[k][i+1]+deln[j]*dF[i][j][k]*w[k][i]
                           +delnp1[j]*dF[i-1][j][k]*w[k][i-1];
         }
      }
      for(j=0;j<nspec;j++)
      {
         res[i][j] = res[i][j]-tmp1[j]; 
      }
   }


/*   At the wall */
 
   for(j=0;j<nspec;j++)
   {
      mrp1[j] = -.5*f1[j];
      mr0[j] = 2.0*f1[j];
      mrm1[j] = -1.5*f1[j]+g1[j];
      mral[j] = -2.0*f1[j];
      mrbe[j] = 2.0*f1[j];/* occhio */
      m1p1[j] = a[0][j]-.5*b[0][j];
      m10[j] = -2.0*a[0][j]+2.0*b[0][j];
      m1m1[j] = a[0][j]-1.5*b[0][j]+c[0][j];
      m1al[j] = 6.0*a[0][j]-2.0*b[0][j];
      m1be[j] = -10.0*a[0][j]+2.0*b[0][j];
      m2p1[j] = a[1][j]+.5*b[1][j];
      m20[j] = -2.0*a[1][j]+c[1][j];
      m2m1[j] = a[1][j]-.5*b[1][j];
      m2al[j] = b[1][j];
      m2be[j] = 2.0*a[1][j];

/* Determinants to eliminate alpha and beta   */

      delnp1[j] = m1al[j]*m2be[j]-m2al[j]*m1be[j];
      deln[j] = m2al[j]*mrbe[j]-mral[j]*m2be[j];
      delnm1[j] = mral[j]*m1be[j]-m1al[j]*mrbe[j];

/* Coefficients of the matrix at the wall  */

      r2[j][j] = mrp1[j]*delnp1[j]+m1p1[j]*deln[j]+m2p1[j]*delnm1[j];
      r1[j][j] = mr0[j]*delnp1[j]+m10[j]*deln[j]+m20[j]*delnm1[j];
      r0[j][j] = mrm1[j]*delnp1[j]+m1m1[j]*deln[j]+m2m1[j]*delnm1[j];
      r[j] = h1[j]*delnp1[j]+(d[0][j]+F[0][j])*deln[j]+(d[1][j]+F[1][j])*delnm1[j];

   }

   for(i=0;i<nspec;i++)
   {
      for(j=0;j<nspec;j++)
      {
         r0[i][j] = r0[i][j]-delnp1[i]*dh1dc[i][j]-deln[i]*dF[0][i][j];
         r1[i][j] = r1[i][j]-delnm1[i]*dF[1][i][j];
      }
   }

   for(i=0;i<nspec;i++)
   {
      tmp1[i] = 0.0;
      for(j=0;j<nspec;j++)
      {
         tmp1[i] = tmp1[i]+delnp1[i]*dh1dc[i][j]*w[j][0]+deln[i]*dF[0][i][j]*w[j][0]
                          +delnm1[i]*dF[1][i][j]*w[j][1];
      }
   }
   for(i=0;i<nspec;i++) r[i] = r[i]-tmp1[i];


   for(i=0;i<nspec;i++)
   {
      D1[i] = res[1][i];
      for(j=0;j<nspec;j++)
      {
         A1[i][j] = A[1][i][j];
         B1[i][j] = B[1][i][j];
         C1[i][j] = C[1][i][j];
      }
   }

   ludcmp(C1,nspec,ind,&db);
   invp(C1,A1,mt0,ind,nspec);
   invp(C1,B1,mt1,ind,nspec);
   invv(C1,D1,vt,ind,nspec);

   matmat(r2,nspec,nspec,mt0,nspec,nspec,B[0]);
   matmat(r2,nspec,nspec,mt1,nspec,nspec,C[0]);
   matvec(r2,nspec,nspec,vt,nspec,res[0]);

   for(i=0;i<nspec;i++)
   {
      res[0][i] = r[i]-res[0][i];
      for(j=0;j<nspec;j++)
      {
         B[0][i][j] = r0[i][j]-B[0][i][j];
         C[0][i][j] = r1[i][j]-C[0][i][j];
      }
   }

  /* for(kk=0;kk<100;kk++)
    {
    for(ii=0;ii<nspec;ii++)
     for(jj=0;jj<nspec;jj++)
       printf("\n A[%d][%d][%d]=%le \n B[%d][%d][%d]=%le  \n C[%d][%d][%d]=%le",
                       kk,ii,jj,A[kk][ii][jj],kk,ii,jj,B[kk][ii][jj],kk,ii,jj,C[kk][ii][jj]);
		 
   
    getchar();
    } */
   
   /*printf(" \n");
   for(i=0;i<neq;i++)
   {
      for(j=0;j<nspec;j++)
      {
         printf(" %le",res[i][j]);
      }
      printf(" \n");
   }*/

   for(i=0;i<nspec;i++)
   {
      tmp2 = 0.0;
      for(j=0;j<nspec;j++) tmp2=tmp2+C[neq-2][i][j]*w[j][neq-1];
      /*res[neq-2][i]=res[neq-2][i]-C[neq-2][i][i]*w[i][neq-1];*/
      res[neq-2][i]=res[neq-2][i]-tmp2;
   }
   for(i=0;i<nspec;i++)
   {
      for(j=0;j<nspec;j++)
      {
         C[neq-2][i][j]=0.0;
      }
   }

   /*printf(" \n");
   for(i=0;i<neq;i++)
   {
      for(j=0;j<nspec;j++)
      {
         printf(" %le",res[i][j]);
      }
      printf(" \n");
   }*/

   trisol(A,B,C,res,0,neq-1,nspec);


   free_intvector(ind); free_vector(tmp1); free_vector(r); free_matrix(r0,nspec);
   free_matrix(r1,nspec); free_matrix(r2,nspec);
   free_matrix(mattmp,nspec);
   free_trimatrix(A,neq,nspec); free_trimatrix(B,neq,nspec);
   free_trimatrix(C,neq,nspec);

   free_matrix(A1,nspec); free_matrix(B1,nspec); free_matrix(C1,nspec);
   free_vector(D1);
   free_matrix(mt0,nspec); free_matrix(mt1,nspec); free_vector(vt);

   free_vector(mnp1m1); free_vector(mnp1al); free_vector(mnp1be); free_vector(mnp1);
   free_vector(mn0); free_vector(mnm1); free_vector(mrp1); free_vector(mr0);
   free_vector(mrm1); free_vector(mral); free_vector(mrbe); free_vector(mnal);
   free_vector(mnbe); free_vector(mnm1p1); free_vector(mnm10); free_vector(mnm1m1);
   free_vector(mnm1al); free_vector(mnm1be); free_vector(mnp1p1); free_vector(mnp10);
   free_vector(m1p1); free_vector(m10); free_vector(m1m1); free_vector(m1al);
   free_vector(m1be); free_vector(m2p1); free_vector(m20); free_vector(m2m1);
   free_vector(m2al); free_vector(m2be);
   free_vector(delnp1); free_vector(deln); free_vector(delnm1);
}



double convergence(int neq,int nspec,double *u,double *u_old,double *h,
                   double *h_old,double **c,double **c_old)
{
   int i,j,ind1,ind2,ind3,ind4;
   double diff,max;
   FILE *fp;

   ind1=0;
   ind2=0;
   ind3=0;
   ind4=0;
   
   max = 0.0;      
   for(i=0;i<neq;i++)
   {
      diff = fabs(u_old[i]-u[i]);
      if(diff >= max)
      {
         max = diff;
         ind1=1; ind2=i;
      }
      diff = fabs(h_old[i]-h[i]);
      if(diff >= max)
      {
         max = diff;
         ind1=2; ind2=i;
       }

      for(j=0;j<nspec;j++)
      {
         diff = fabs(c_old[j][i]-c[j][i]);
         if(diff >= max)
         {
            max = diff;
            ind1=3; ind2=i; ind3=j;
         }
      }
   }
//   	printf("\nind1 %d, ind2 %d,ind3 %d ",ind1, ind2,ind3);
     return max;
}


void under_relax(int neq,int nspec,double *u,double *u_old,double *h,double *h_old,
                 double **c,double **c_old,double under_eps)
{
   int i,j;

   for(i=0;i<neq;i++)
   {
      u[i] = under_eps*u[i]+(1.0-under_eps)*u_old[i];
      h[i] = under_eps*h[i]+(1.0-under_eps)*h_old[i];
      for(j=0;j<nspec;j++)
      {
         c[j][i] = under_eps*c[j][i]+(1.0-under_eps)*c_old[j][i];
        /* c[j][i] = max(0.0,c[j][i]);*/
        /* if(c[j][i]<0.0) c[j][i]=c_old[j][i]; */
      }
   }

/*   for(i=0;i<neq;i++)
      for(j=0;j<nspec;j++)
         c[j][i] = fabs(c[j][i]); */

}
