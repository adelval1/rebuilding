#include<stdio.h>
#include<math.h>
void catc_aircsicvki_(int nb_species, double concentration[nb_species],double  molarmass[nb_species],double walltemperature,double wallmoleproduction[nb_species]){
//status:
//	not tested

/*
Source:
\bibitem{b:vki2} Chazot, O., Collin, P. , Asma, C.,
"TPS Plasma Tests for Analysis of the Gas ? Surface Interaction: Steady State Testing of C/C-SiC Material,"
VKI Contract Report 2005-13


to do:
use only pointers
*/

	double catalycity,reactantconsumption;

	catalycity=6.4568E-06*exp(3.5710E-03*walltemperature); 

	CATC_CatalycityMoleProduction(concentration[3],molarmass[3],walltemperature,catalycity,&reactantconsumption);
	wallmoleproduction[0]=-reactantconsumption*1/2;
	wallmoleproduction[3]=reactantconsumption;

	CATC_CatalycityMoleProduction(concentration[4],molarmass[4],walltemperature,catalycity,&reactantconsumption);
	wallmoleproduction[1]=-reactantconsumption*1/2;
	wallmoleproduction[4]=reactantconsumption;

	if(walltemperature<1400||walltemperature>2000){
		printf("\n error:");
		printf("\n the catalysis model is valid in the range between 1400 and 2000K");
		printf("\n currently the wall temperature is about %e", walltemperature);
		getchar();
	}
}
