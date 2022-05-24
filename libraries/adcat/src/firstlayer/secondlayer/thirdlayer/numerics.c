#include<stdio.h>
#include <math.h>
#include <stdlib.h> 
#define underrelaxation 0.9

void catc_mnewt_(int nb_reactions,int nb_adspecies,int nb_species,int ntrial,double tolx,double tolf,double rate[nb_reactions],double concentration[nb_species+nb_adspecies],int actual_species[nb_species+nb_adspecies],int stochio[nb_reactions][nb_species+nb_adspecies],int stochiodiff[nb_reactions][nb_species+nb_adspecies],double max_concentration){
 int n,NP;
 double x[nb_adspecies-1],sum;
//!achtung---------------------------------------------------------------------------------------------------
	NP=nb_adspecies-1;
	n=nb_adspecies-1;
 
//! 	Given an initial guess x for a root in n dimensions, take ntrial
//!	Newton-	Raphson steps to improve the root. Stop if the root
//!	converges in 	either summed absolute variable increments
//!	tolx or summed absolute function values tolf.
  	int i,k,indx[NP],j,l;
  	double d,errf,errx,fjac[NP][NP],fvec[NP],p[NP];

	for(l=0;l<(nb_adspecies-1);l++){
		x[l]=concentration[l+nb_species];
	}

	for(k=0;k<ntrial;k++){
		sum=0.;
		for(l=0;l<nb_adspecies-1;l++){
			concentration[l+nb_species]=x[l];
			sum=sum+x[l];
		}
		concentration[nb_species+nb_adspecies-1]=max_concentration-sum;

	       	CATC_usrfun(nb_reactions,nb_species,nb_adspecies,fvec,fjac,rate,concentration,stochio,stochiodiff,actual_species);

/*		printf("\nfjac");
		for(j=0;j<5;j++){		
		for(i=0;i<5;i++){

 			printf(" %e ",fjac[i][j]);
	  	}
		printf("\n");
		}
*/		
		
//!		User subroutine supplies function values at x in fvec
//!   		and Jacobian matrix in fjac.
		errf=0.;
//!   		Check function convergence.
		for(i=0;i<n;i++){
 			errf=errf+fabs(fvec[i]);
	  	}
 		if(errf<=tolf){ 
//			printf("returned at position1 k= %d errf = %f \n",k,errf);
			sum=0.;
			for(l=0;l<nb_adspecies-1;l++){
				concentration[l+nb_species]=x[l];
				sum=sum+x[l];
			}	
			concentration[nb_species+nb_adspecies-1]=max_concentration-sum;
			return;
		}

//!		Right-hand side of linear equations.
	  	for(i=0;i<n;i++){
			p[i]=-fvec[i];
	  	}
/*
		printf("\nfjac");
		for(j=0;j<5;j++){		
		for(i=0;i<5;i++){

 			printf(" %e ",fjac[i][j]);
	  	}
		printf("\n");
		}
*/
		CATC_ludecomposition(n, &d,fjac,indx);
 		
//		CATC_ludecomposition(n,NP,indx,d,fjac);

/*
		printf("\nfjac");
		for(j=0;j<5;j++){		
		for(i=0;i<5;i++){

 			printf(" %e ",fjac[i][j]);
	  	}
		}
		getchar();
*/		
		
//! 		Solve linear equations using LU decomposition.
//	  	CATC_lubacksubstitution(n,NP,indx,p,fjac);
		CATC_lubksb(n,indx, p, fjac);

//!		Check root convergence.
 		errx=0.;
//!  		Update solution.

		for(i=0;i<n;i++){
			errx=errx+fabs(p[i]);
//!		ACHTUNG UNTERRELAXIERT
			x[i]=x[i]+underrelaxation*p[i];
		}

		if(errx<=tolx) {
//			printf("returned at position2 k= %d",k);		
			sum=0.;
			for(l=0;l<(nb_adspecies-1);l++){
				concentration[l+nb_species]=x[l];
				sum=sum+x[l];
			}
			concentration[nb_species+nb_adspecies-1]=max_concentration-sum;
			return;
		}
	} 

	printf("returned at position3 k= %d",k);
	sum=0.;
	for(l=0;l<(nb_adspecies-1);l++){
		concentration[l+nb_species]=x[l];
		sum=sum+x[l];
		printf("\n x( %d ) %e ",l,x[l]);
	}
		printf("\n");fflush(stdout);
	concentration[nb_species+nb_adspecies-1]=max_concentration-sum;


	return;

}



//!********************************************************************************************************************
//!********************************************************************************************************************

CATC_usrfun(int nb_reactions,int nb_species,int nb_adspecies,double fvec[nb_adspecies-1],double fjac[nb_adspecies-1][nb_adspecies-1],double rate[nb_reactions],double concentration[nb_species+nb_adspecies],int stochio[nb_reactions][(nb_species+nb_adspecies)],int stochiodiff[nb_reactions][(nb_species+nb_adspecies)],int actual_species[nb_species+nb_adspecies]){
	int h,i,j,k,l,r;
	double temp,temp2,schalter;


	for(i=0;i<(nb_adspecies-1);i++){
		fvec[i]=0.;
	}
	for(i=0;i<nb_adspecies-1;i++){
		for(j=0;j<nb_reactions;j++){
			temp=1.;
			for(k=0;k<(nb_species+nb_adspecies);k++){
				temp=temp*pow(concentration[k],(double)stochio[j][k]);
//					printf("\n adspecies %d reaction %d  species %d temp %e",i,j,k,pow(concentration[k],(double)stochio[j][k]));
			}
		fvec[i]=fvec[i]+rate[j]*stochiodiff[j][(i+nb_species)]*temp;
//		printf("\n fvec %e",rate[j]*stochiodiff[j][(i+nb_species)]*temp);
		}
//		getchar();
	}

	
	for(i=0;i<(nb_adspecies-1);i++){
		for(l=0;l<(nb_adspecies-1);l++){
			fjac[i][l]=0.;
		}
	}			

	for(i=0;i<(nb_adspecies-1);i++){
		for(l=0;l<(nb_adspecies-1);l++){
			for(r=0;r<nb_reactions;r++){
				temp=1.;
				for(j=0;j<(nb_adspecies+nb_species);j++){
					if(j!=(l+nb_species)){
						temp=temp*pow(concentration[j],stochio[r][j]);
					}
					else{
						if(stochio[r][(l+nb_species)]!=0){
							temp=temp*stochio[r][j]*pow(concentration[j],(stochio[r][j]-1));
						}
						else{
							temp=0.;
						}
					}	


					schalter=0.;
					if(stochio[r][(nb_adspecies+nb_species-1)]==1){
						schalter=1.;
					}
					else{
						if(stochio[r][(nb_adspecies+nb_species-1)]==2){
							schalter=2.;
						}
						else{
							if(stochio[r][(nb_adspecies+nb_species-1)]==0){
//								!do nothing
							}
							else{
								printf("error in differentiation");
							}
						}
					}
				}
				if(schalter==1.){
					temp2=1.;
					for(j=0;j<(nb_adspecies+nb_species-1);j++){
						temp2=temp2*pow(concentration[j],stochio[r][j]);
					}
					temp=temp-temp2;
				}
				else{
					if(schalter==2.){
						temp2=1.;
						for(j=0;j<(nb_adspecies+nb_species-1);j++){
							temp2=temp2*pow(concentration[j],stochio[r][j]);
						}
						temp=temp-2*temp2*concentration[nb_adspecies+nb_species-1];
					}


					
				}

			fjac[i][l]=fjac[i][l]+rate[r]*stochiodiff[r][(i+nb_species)]*temp;

			}
		}
	}

//????????????????
	for(k=(nb_species);k<(nb_species+nb_adspecies-1);k++){
		if(actual_species[k]==0){
			for(i=0;i<nb_adspecies-1;i++){
				fjac[i][k-nb_species]=0.;
				fjac[k-nb_species][i]=0.;
			}
			fjac[(k-nb_species)][(k-nb_species)]=1.;
			fvec[(k-nb_species)]=0.;
		}
	}

/*
printf("\n");	
for(i=0;i<(nb_adspecies-1);i++){
	for(l=0;l<(nb_adspecies-1);l++){
		printf("%e ",fjac[i][l]);
	}
	printf("\n");
}	
//getchar();

for(l=0;l<(nb_adspecies-1);l++){
	printf("%e ",fvec[l]);
}
printf("\n");
//getchar();

*/
}


/*
CATC_ludecomposition(int n,int np,int indx[n],double d,double a[np][np]){
int NMAX,aaa,sss;
double TINY;
NMAX=500;
TINY=1.0e-20;
/*!	Largest expected n, and a small number.
!	Given a matrix a(1:n,1:n), with physical dimension np by np, this routine replaces it by
! 	the LU decomposition of a rowwise permutation of itself. a and n are input. a is output,
!  	arranged as in equation (2.3.14) above; indx(1:n) is an output vector that records the
!   	row permutation e\ufb00ected by the partial pivoting; d is output as ±1 depending on whether
!  	the number of row interchanges was even or odd, respectively. This routine is used in
! 	combination with lubksb to solve linear equations or invert a matrix.
	int i,imax,j,k;
//! 	vv stores the implicit scaling of each row.
   	double aamax,dum,sum,vv[NMAX];
//!	No row interchanges yet.

	
   	d=1.;
//! 	Loop over rows to get the implicit scaling information.
 	for( i=0;i<n;i++){
		aamax=0.;
	  	for(j=0;j<n;j++){
 			if (fabs(a[i][j])>=aamax){ aamax=fabs(a[i][j]);}
		} 
	  	if (aamax==0.){ printf("numerics.c:line233");getchar();}
		vv[i]=1./aamax;
	} 
	
//!		This is the loop over columns of Crout\u2019s method.
 	for(j=0;j<n;j++){
  		for(i=0;i<j-1;i++){
  			sum=a[i][j];
	  		for(k=0;k<i-1;k++){
				sum=sum-a[i][k]*a[k][j];
			} 
			a[i][j]=sum;
		} 
//!	 	Initialize for the search for largest pivot element.
 		aamax=0.;
//!	 	This is i = j of equation (2.3.12) and i = j + 1 . . . Nof equation (2.3.13).
		for(i=j;i<n;i++){		
	  		sum=a[i][j];
			for(k=0;k<j-1;k++){
	 			sum=sum-a[i][k]*a[k][j];
		 	} 
	   		a[i][j]=sum;
//!			Figure of merit for the pivot.
			dum=vv[i]*fabs(sum);
//!			Is it better than the best so far?
 	   		if (dum>=aamax){
  				imax=i;
				aamax=dum;
			}
  		} 
//!	  	do we need to interchange rows?
 		if (j!=imax){
//!  			Yes,do so...
			for(k=0;k<n;k++){
 				dum=a[imax][k];
			 	a[imax][k]=a[j][k];
				a[j][k]=dum;
	   		} 
//!	  		...and change the parity of d.
			d=-d;
//!  			Also interchange the scale factor.
	 		vv[imax]=vv[j];
		}
		indx[j]=imax;
		if(a[j][j]==0.){a[j][j]=TINY;}
//!		If the pivot element is zero the matrix is singular (at least to the precision of the al-
//!	  	gorithm). For some applications on singular matrices, it is desirable to substitute TINY
//!		for zero.
//!		finally, divide by the pivot element.
		if(j!=n){
  			dum=1./a[j][j];
	  		for(i=j+1;i<n;i++){
					a[i][j]=a[i][j]*dum;
	  			} 
			}
		} 



return;
}


*/



/*
CATC_lubacksubstitution(int n,int np,int indx[n],double b[n],double a[np][np]){

	Solves the set of n linear equations A · X = B . Here a is input, not as the matrix A but
! 	rather as its LU decomposition, determined by the routine ludcmp. indx is input as the
!	permutation vector returned by ludcmp. b(1:n) is input as the right-hand side vector B ,
!	and returns with the solution vector X . a, n, np, and indx are not modi\ufb01ed by this routine
!	and can be left in place for successive calls with di\ufb00erent right-hand sides b. This routine
!	takes into account the possibility that b will begin with many zero elements, so it is e\ufb03cient
!	for use in matrix inversion.
  	int i,ii,j,ll,aaa,sss;
  	double sum;
! 	When ii is set to a positive value, it will become the in-
! 	dex of the \ufb01rst nonvanishing element of b. We now do
!	the forward substitution, equation (2.3.6). The only new
!	wrinkle is to unscramble the permutation as we go.
	ii=0;
	
	for(i=0;i<n;i++){
 		ll=indx[i];
		sum=b[ll];
 		b[ll]=b[i];
 		if (ii!=0){
			for(j=ii;j<i-1;j++){
				sum=sum-a[i][j]*b[j];
			}
		}
 		else{
			 if (sum!=0.){
//!	  			A nonzero element was encountered, so from now on we will
//!				have to do the sums in the loop above.
				ii=i;
			}
  		}
		b[i]=sum;
 	}
	
	

//!  	Now we do the backsubstitution, equation (2.3.7).
 	for(i=n-1;i>=0;i--){
 		sum=b[i];
  		for(j=i+1;j<n;j++){
 			sum=sum-a[i][j]*b[j];
 		}
//! 		Store a component of the solution vector X .
		b[i]=sum/a[i][i];
 	}

  return;
}
*/

void CATC_nrerror(char error_text[]){
  	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
  	fprintf(stderr,"...now exiting to system...\n");
  	exit(1);
}


