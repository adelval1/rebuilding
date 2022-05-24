#include<stdio.h>
#include<math.h>
void catc_airpm1000preci_(int nb_species, double concentration[nb_species],double  molarmass[nb_species],double walltemperature,double wallmoleproduction[nb_species]){
//status:
//	not tested

/*
Source:
\bibitem{b:preci}   Preci,A., Lein, S., Schuessler, M., Auweter-Kurtz, M., Fertig, M., Herdrich, G., and M. Winter,
"Numerical Simulation and IRS Instrumentation Design for EXPERT,"
\textit{25th AIAA Aerodynamic Measurement Technology and Ground Testing Conference}, 5 - 8 June 2006, San Francisco, California

to do:
use only pointers
*/



	double catalycity,reactantconsumption,coefficient[2][7];
//Nitrogen Data
	coefficient[0][0]=-107.53;
	coefficient[0][1]=0.58;
	coefficient[0][2]=-9.71e-4;
	coefficient[0][3]=1.13e-7;
	coefficient[0][4]=1.04e-9;
	coefficient[0][5]=-6.68e-13;
	coefficient[0][6]=1.;
//Oxygen Data
	coefficient[1][0]=19.34;
	coefficient[1][1]=-0.17;
	coefficient[1][2]=4.08e-4;
	coefficient[1][3]=-4.21e-7;
	coefficient[1][4]=2.01e-10;
	coefficient[1][5]=-3.64e-14;
	coefficient[1][6]=1.;
	

	if(walltemperature<400||walltemperature>1900){
		printf("\n error:");
		printf("\n the catalysis model is valid in the range between 400 and 1900K");
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
}
