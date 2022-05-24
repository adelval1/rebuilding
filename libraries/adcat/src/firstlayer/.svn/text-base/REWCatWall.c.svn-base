#include<stdio.h>
void catc_invstefanmaxwelldiffusion_(int *nb_species, double temperaturegradient[*nb_species], double *walltemperature, double thermaldiffusioncoefficient[*nb_species],double moleproduction[*nb_species],double species_concentration[*nb_species],double D[*nb_species][*nb_species],double speciesconcentrationgradient[*nb_species]){
/*

status:
PRELIMINARY IMPLEMENTED!!! DO NOT USE!

input
-number of species
-species concentration
-mass priduction of catalysis

output
-species gradients
*/

	int i,j;
	double sum,mixture_concentration;
	
	mixture_concentration=0.;
	for(i=0;i<*nb_species;i++){
		mixture_concentration=mixture_concentration+species_concentration[i];
	}

// calculate species gradient

	for(i=0;i<*nb_species;i++){
		sum=0;
		for(j=0;j<*nb_species;j++){
			sum=sum+species_concentration[i]*moleproduction[j]/D[i][j]-species_concentration[j]*moleproduction[i]/D[i][j];
		}
//	speciesconcentrationgradient[i]=sum/mixture_concentration-thermaldiffusioncoefficient[i]/ *walltemperature*temperaturegradient[i];	
	speciesconcentrationgradient[i]=-sum/mixture_concentration;	
	}
}
