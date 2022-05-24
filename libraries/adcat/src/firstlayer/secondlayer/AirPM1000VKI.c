#include<stdio.h>
#include<math.h>
void catc_airpm1000vki_(int nb_species, double concentration[nb_species],double  molarmass[nb_species],double walltemperature,double wallmoleproduction[nb_species]){
//status:
//	not tested

/*
Source:
\bibitem{b:vki1} Chazot, O., Collin, P. , Asma, C.,
"TPS Plasma Tests for Analysis of the Gas ? Surface Interaction: Steady State Testing of Pre-Oxidized PM1000 Material ? 3rd Campaign,"
VKI Contract Report 2005-26


to do:
use only pointers
*/

	double catalycity,reactantconsumption;

	catalycity=2.3191e-6*exp(7.1554e-3*walltemperature);

	CATC_CatalycityMoleProduction(concentration[3],molarmass[3],walltemperature,catalycity,&reactantconsumption);
	wallmoleproduction[0]=-reactantconsumption*1/2;
	wallmoleproduction[3]=reactantconsumption;

	CATC_CatalycityMoleProduction(concentration[4],molarmass[4],walltemperature,catalycity,&reactantconsumption);
	wallmoleproduction[1]=-reactantconsumption*1/2;
	wallmoleproduction[4]=reactantconsumption;

	if(walltemperature<1000||walltemperature>1500){
		printf("\n error:");
		printf("\n the catalysis model is valid in the range between 1000 and 1500K");
		printf("\n currently the wall temperature is about %e", walltemperature);
		getchar();
	}
}
