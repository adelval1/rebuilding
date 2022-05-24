#include<stdio.h>
#include<math.h>
void catc_airbetacristobalitebalat_(int nb_species, double concentration[nb_species],double  molarmass[nb_species],double walltemperature,double wallmoleproduction[nb_species]){
/*
Source:
\bibitem{b:b}  Balat, M.,
"Interaction of Reactive Gas Flows and Ceramics at High Temperature-Experimental Methods for the Measurement of Species Recombination during Planetary Entry,"
\textit{RTO-AVT-VKI Lectures Series 2006: Experiment, Modelling and Simulation of Gas-Surface Interactions for Reactive Flows in Hypersonic Flights}
\bibitem{b:b}  Balat, M.,
"Interaction of Reactive Gas Flows and Ceramics at High Temperature-Experimental Methods for the Measurement of Species Recombination during Planetary Entry,"
\textit{RTO-AVT-VKI Lectures Series 2006: Experiment, Modelling and Simulation of Gas-Surface Interactions for Reactive Flows in Hypersonic Flights}

to do:
use only pointers
*/

	double catalycity,reactantconsumption;
	
	catalycity=0.6382*exp(-3374./walltemperature);
	catc_catalycitymoleproduction_(concentration[3],molarmass[3],walltemperature,catalycity,&reactantconsumption);
	wallmoleproduction[0]=-reactantconsumption*1/2;
	wallmoleproduction[3]=reactantconsumption;
	if(walltemperature<800||walltemperature>1830){
		printf("\n error:");
		printf("\n the catalysis model is valid in the range between 800 and 1830K");
		printf("\n currently the wall temperature is about %e", walltemperature);
		getchar();
	}
}
