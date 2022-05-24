#include <stdlib.h> 
#include <math.h>
#include<stdio.h>
#define TINY 1.0e-20
#define NR_END 1
#define  FREE_ARG char*

double *CATC_vector(long nl, long nh){
	double *v;
	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) CATC_nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

double **CATC_matrix(long nrl, long nrh, long ncl, long nch){
 	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
 	double **m;
 	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
 	if (!m) CATC_nrerror("allocation failure 1 in matrix()");
 	m += NR_END;
 	m -= nrl;
 	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
 	if (!m[nrl]) CATC_nrerror("allocation failure 2 in matrix()");
 	m[nrl] += NR_END;
 	m[nrl] -= ncl;
 	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
 	/* return pointer to array of pointers to rows */
 	return m;
}
void CATC_free_vector(double *v, long nl, long nh){
 	free((FREE_ARG) (v+nl-NR_END));
}
 void CATC_free_matrix(double **m, long nrl, long nrh, long ncl, long nch) {
 	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}



CATC_ludecomposition(int n, double *d,double b[n][n], int indx[n]){
int i,imax,j,k;
double big,dum,sum,temp;
double *vv,**a;

vv=CATC_vector(1,n);
a=CATC_matrix(1,n,1,n);
    	for (i=1;i<=n;i++) {
	    	for (j=1;j<=n;j++) {
			a[i][j]=b[i-1][j-1];
		}
	}


*d=1.0;

for (i=1;i<=n;i++) {
	big=0.0;
    	for (j=1;j<=n;j++) {

		if ((temp=fabs(a[i][j])) > big) big=temp;
	};
	if (big == 0.0) {
	    	for (i=1;i<=n;i++) {
		    	for (j=1;j<=n;j++) {
				printf(" %e ",a[i][j]);
			}
			printf("\n");
		}

	CATC_nrerror("Singular matrix in routine ludcmp");
	}
 	vv[i]=1.0/big;
}

for (j=1;j<=n;j++) {
	for (i=1;i<j;i++) {
    		sum=a[i][j];
		for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
  		a[i][j]=sum;
  	}
 	big=0.0;
	for (i=j;i<=n;i++) {
   		sum=a[i][j];
 		for (k=1;k<j;k++)sum -= a[i][k]*a[k][j];
		a[i][j]=sum;
    		if ( (dum=vv[i]*fabs(sum)) >= big) {
    			big=dum;
 	   		imax=i;
    		}
	}
	if (j != imax) {
    		for (k=1;k<=n;k++) {
    			dum=a[imax][k];
	    		a[imax][k]=a[j][k];
    			a[j][k]=dum;
		}
		*d = -(*d);
		vv[imax]=vv[j];
	}
	indx[j]=imax;
	if (a[j][j] == 0.0) a[j][j]=TINY;
	if (j != n) {
		dum=1.0/(a[j][j]);
	    	for (i=j+1;i<=n;i++) a[i][j] *= dum;
	}
}
    	for (i=1;i<=n;i++) {
	    	for (j=1;j<=n;j++) {
			b[i-1][j-1]=a[i][j];
		}
	}

CATC_free_vector(vv,1,n);
CATC_free_matrix(a,1,n,1,n);
}

 CATC_lubksb( int n, int indx[n], double bb[n],double ab[n][n]){
 	int i,ii=0,ip,j;
 	double sum;
	double **a,*b;

	a=CATC_matrix(1,n,1,n);
	b=CATC_vector(1,n);
	
    	for (i=1;i<=n;i++) {
	    	for (j=1;j<=n;j++) {
			a[i][j]=ab[i-1][j-1];
		}
	b[i]=bb[i-1];
	}	

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii) for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
 		else if (sum) ii=i;
 		b[i]=sum;
	}
	 for (i=n;i>=1;i--) {
 		sum=b[i];
		 for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
 		b[i]=sum/a[i][i];
	}
    	for (i=1;i<=n;i++) {
			bb[i-1]=b[i];
	}	
	CATC_free_matrix(a,1,n,1,n);
	CATC_free_vector(b,1,n);
}



