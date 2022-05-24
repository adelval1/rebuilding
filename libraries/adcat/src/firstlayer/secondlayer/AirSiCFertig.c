#include<stdio.h>
#include<math.h>
void catc_airsicfertig_(int nb_species, double concentration[nb_species],double  molarmass[nb_species],double walltemperature,double wallmoleproduction[nb_species]){
//status:
//	not double checked
//	not tested

/*
Source:
\bibitem{b:mf}  Fertig, M.,
"Modellierung reaktiver Prozesse auf Siliziumkarbid-Oberfl\"achen on verd\"unten Nichtgleichgewichts-Luftstr\"omungen,"
PhD thesis, University of Stuttgart, Germany, November 2005


to do:
use only pointers
*/



	double catalycity,reactantconsumption,coefficient[3][7];
//Nitrogen Data
	coefficient[0][0]=-8.991347005726e0;
	coefficient[0][1]= 1.343438622951e-3;
	coefficient[0][2]= 2.942582244104e-5;
	coefficient[0][3]=-5.498614068676e-8;
	coefficient[0][4]= 3.635962640421e-11;
	coefficient[0][5]=-8.018402930175e-15;
	coefficient[0][6]=1.;
//Oxygen Data
	coefficient[1][0]=-1.074716503382e1;
	coefficient[1][1]= 7.032362079743e-3;
	coefficient[1][2]= 1.078667217413e-5;
	coefficient[1][3]=-9.984140196274e-9;
	coefficient[1][4]= 1.492827945385e-12;
	coefficient[1][5]=0;
	coefficient[1][6]=1.;
//Nitric Oxide Data
	coefficient[2][0]=-8.991347005726e0;
	coefficient[2][1]= 1.343438622951e-3;
	coefficient[2][2]= 2.942582244104e-5;
	coefficient[2][3]=-5.498614068676e-8;
	coefficient[2][4]= 3.635962640421e-11;
	coefficient[2][5]=-8.018402930175e-15;
	coefficient[2][6]=1.;
	

	if(walltemperature<300||walltemperature>2100){
		printf("\n error:");
		printf("\n the catalysis model is valid in the range between 300 and 2100K");
		printf("\n currently the wall temperature is about %e", walltemperature);
		getchar();
	}

//Nitrogen recombination
	catalycity=coefficient[0][6]*exp(coefficient[0][0]+coefficient[0][1]*walltemperature+coefficient[0][2]*walltemperature*walltemperature+coefficient[0][3]*walltemperature*walltemperature*walltemperature+coefficient[0][4]*walltemperature*walltemperature*walltemperature*walltemperature+coefficient[0][5]*walltemperature*walltemperature*walltemperature*walltemperature*walltemperature);
	CATC_CatalycityMoleProduction(concentration[4],molarmass[4],walltemperature,catalycity,&reactantconsumption);
	wallmoleproduction[1]=-reactantconsumption*1/2;
	wallmoleproduction[4]=reactantconsumption;

//Oxygen recombination
	catalycity=coefficient[1][6]*exp(coefficient[1][0]+coefficient[1][1]*walltemperature+coefficient[1][2]*walltemperature*walltemperature+coefficient[1][3]*walltemperature*walltemperature*walltemperature+coefficient[1][4]*walltemperature*walltemperature*walltemperature*walltemperature+coefficient[1][5]*walltemperature*walltemperature*walltemperature*walltemperature*walltemperature);
	CATC_CatalycityMoleProduction(concentration[3],molarmass[3],walltemperature,catalycity,&reactantconsumption);
	wallmoleproduction[0]=-reactantconsumption*1/2;
	wallmoleproduction[3]=reactantconsumption;

//Nitric Oxide Dissociation
	catalycity=coefficient[2][6]*exp(coefficient[2][0]+coefficient[2][1]*walltemperature+coefficient[2][2]*walltemperature*walltemperature+coefficient[2][3]*walltemperature*walltemperature*walltemperature+coefficient[2][4]*walltemperature*walltemperature*walltemperature*walltemperature+coefficient[2][5]*walltemperature*walltemperature*walltemperature*walltemperature*walltemperature);
	CATC_CatalycityMoleProduction(concentration[2],molarmass[2],walltemperature,catalycity,&reactantconsumption);
	wallmoleproduction[0]=reactantconsumption;
	wallmoleproduction[3]=-reactantconsumption;
	wallmoleproduction[4]=-reactantconsumption;
}
