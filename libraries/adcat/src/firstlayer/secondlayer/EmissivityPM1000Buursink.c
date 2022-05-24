#include<stdio.h>
void catc_emissivitypm1000buursink_(double *walltempperature,double  *emissivity ){
/*
Source:
\bibitem{b:buursinkphd} Buursink, J.,
"On the Development of a Cooled Metallic Thermal Protection System for Spacecraft,"
\textit{PhD thesis, TU-Delft}, 2005


not double checked!
*/

*emissivity=-1.49e-7* *walltempperature* *walltempperature+4.89e-4* *walltempperature+0.456;

	if(*walltempperature<950||*walltempperature>1600){
		printf("\n error:");
		printf("\n the Emissivity model is valid in the range between 950 and 1600K");
		printf("\n currently the wall temperature is about %e", *walltempperature);
		getchar();
	}
}
