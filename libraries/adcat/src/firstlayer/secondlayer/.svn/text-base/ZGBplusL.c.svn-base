#include <math.h>
#include <stdio.h>
#include <physicalconstants.h>
//basic catalysis model specific data
#define nb_adspecies 6					//number of adsorbed species(including empty sites)
#define nb_reactions 10					//number of gas surface interactio reactions
#define max_concentration 2.16E-05			//maximum concentration on the surface in [mol/m^2]
#define tolf 1.e-15					//convergence criterium for newton method
#define tolx 1.e-15					//convergence criterium for newton method
#define nb_iterations 1000				//convergence criterium for newton method

void catc_zgbplusl_(int *nb_species,double catalycities[*nb_species],double adsorbed_concentration[nb_adspecies],double *walltemperature,double gas_concentration[*nb_species],double molarmass[*nb_species],double wallmoleproduction[*nb_species]){

/*-------------------------------------------------------------------------------
This routines provides the boundary condition for the species equation in a 
chemical non equilibrium computation. It computes the mass production rate at
the wall applying a reaction rate based model for CO2.

Name of the Model:
ZGB plus L (Ziff Gulari Barshad extented low temperature)

Publication:
(not pubished yet)

status:
	not tested

output variables:
	double wallmoleproduction[*nb_species]				vector of double of size 5(CO2,O2,CO,O,C),wallproduction rate in [mol/m^3]
in/output variables
	double adsorbed_concentration[nb_adspecies]	vector of double of size 6(CO2s,O2s,COs,Os,Cs,V), concentration on the surface in [mole/m^2], serves only as initial condition for efficiency, dummy values are allowed
input variables:
	double T					double, walltemperature in [K]
	double gas_concentration[*nb_species]		vector of double of size 5(CO2,O2,CO,O,C), concentration in the gas phase in [mol/m^3]
	double molarmass[*nb_species]			vector of double of size 5(CO2,O2,CO,O,C), molar mass of species in [kg/mol]

implementation done by: jan thoemel, November 07

to do:
make working for arcitrary co
-------------------------------------------------------------------------------*/


//1.Declaration of local variables
      	int i,j,sw,sb,r;
      	double alpha1,alpha2;
	double concentration[*nb_species+nb_adspecies];
	int actual_species[*nb_species+nb_adspecies];
	int stochio[nb_reactions][*nb_species+nb_adspecies];
	int stochiodiff[nb_reactions][*nb_species+nb_adspecies];
	double param_arrh[nb_reactions][3];
	double rate[nb_reactions];
	double production[nb_reactions][*nb_species+nb_adspecies];
	double M_inc[*nb_species];


//CHECK
       if(*nb_species!=5){printf("\n\nstop this version of the catalysis library works only with 5 gas species\n\n");fflush(stdout);getchar();}


//2.Initilization
	for(i=0;i<nb_reactions;i++){	
		for(j=0;j<*nb_species+nb_adspecies;j++){	
        	 	stochio[i][j]=0;
         		stochiodiff[i][j]=0;
		}
	}
	for(i=0;i<(*nb_species+nb_adspecies);i++){actual_species[i]=1;}
        for(i=0;i<*nb_species;i++){concentration[i]=gas_concentration[i];}
        for(i=0;i<nb_adspecies;i++){concentration[i+*nb_species]=adsorbed_concentration[i]; }

   	for(i=0;i<*nb_species;i++){
      		M_inc[i]=concentration[i]*sqrt(8.314* *walltemperature/(2.0*catc_pi*molarmass[i]));
   	}	
	
//3.Further catalysis model specific data

	actual_species[4]=0;
        actual_species[5]=0;
        actual_species[6]=0;
        actual_species[9]=0;

//reaction 1
//Dissoziative Adsorption of O2
//O2+2s->2Os
         stochio[0][1]=1;
         stochio[0][10]=2;
         stochiodiff[0][1]=-1;
         stochiodiff[0][10]=-2  ;       
         stochiodiff[0][8]=2 ;        
//reaction 2
//Adsorption of CO
//CO+s->COs
         stochio[1][2]=1;
         stochio[1][10]=1;
         stochiodiff[1][2]=-1;
         stochiodiff[1][10]=-1;
         stochiodiff[1][7]=1;
//reaction 3
//Adsorption of O
//O+s->Os
         stochio[2][3]=1;
         stochio[2][10]=1;
         stochiodiff[2][3]=-1;         
         stochiodiff[2][10]=-1;         
         stochiodiff[2][8]=1;         
//reaction 4
//Desorption of CO
//COs->CO+s
         stochio[3][7]=1;
         stochiodiff[3][7]=-1;         
         stochiodiff[3][2]=1;         
         stochiodiff[3][10]=1;         
//!reaction 5
//Desorption of O
//Os->O+s
         stochio[4][8]=1;
         stochiodiff[4][8]=-1;         
         stochiodiff[4][3]=1;         
         stochiodiff[4][10]=1;         
//reaction 6
//Eley Rideal of CO with O
//COs+O->CO2 +s
         stochio[5][7]=1;
         stochio[5][3]=1;
         stochiodiff[5][7]=-1;         
         stochiodiff[5][3]=-1;         
         stochiodiff[5][0]=1;         
         stochiodiff[5][10]=1;         
//reaction 7
//Eley Rideal of O with CO
//Os+CO->CO2+s
         stochio[6][8]=1;
         stochio[6][2]=1;
         stochiodiff[6][8]=-1;         
         stochiodiff[6][2]=-1;         
         stochiodiff[6][0]=1;         
         stochiodiff[6][10]=1;         
//reaction 8
//Eley Rideal of O with O
//O(s)+O-> O2 +s
         stochio[7][8]=1;
         stochio[7][3]=1;
         stochiodiff[7][8]=-1;         
         stochiodiff[7][3]=-1;         
         stochiodiff[7][1]=1;         
         stochiodiff[7][10]=1;         
//reaction 9
//Lamgmuir Hinshelwood of CO and O
//CO(s)+O(s)->CO2 +2s
         stochio[8][7]=1;
         stochio[8][8]=1;
         stochiodiff[8][7]=-1;         
         stochiodiff[8][8]=-1;         
         stochiodiff[8][0]=1;         
         stochiodiff[8][10]=2;       
//!reaction 10
//Langmuir Hinshelwood of O
//2O(s)->O2 + 2s
         stochio[9][8]=2;
         stochiodiff[9][8]=-2;         
         stochiodiff[9][1]=1;         
         stochiodiff[9][10]=2;         
/*

//Dissoziative Adsorption of O2
	param_arrh[0][0]=5.64948E+11;
	param_arrh[0][1]=0.00;
	param_arrh[0][2]=3387.83928;
//Adsorption of CO	
	param_arrh[1][0]=13592109.52;
	param_arrh[1][1]=0.0;
	param_arrh[1][2]=3798.23472;
//Adsorption of O	
	param_arrh[2][0]=17894429.12;
	param_arrh[2][1]=0.0;
	param_arrh[2][2]=3775.8042;
//Desorption of CO	
	param_arrh[3][0]=7.56709E+14;
	param_arrh[3][1]=0.0;
	param_arrh[3][2]=171551.94;
//Desorption of O	
	param_arrh[4][0]=5.1191E+12;
	param_arrh[4][1]=0.0;
	param_arrh[4][2]=250415.9868;
//Eley Rideal of CO with O	
	param_arrh[5][0]=21252816.8;
	param_arrh[5][1]=0.0;
	param_arrh[5][2]=5640.02964;
//Eley Rideal of O with CO	
	param_arrh[6][0]=11524664.09;
	param_arrh[6][1]=0.0;
	param_arrh[6][2]=3721.8048;
//Eley Rideal of O with O
	param_arrh[7][0]=15248626.63;
	param_arrh[7][1]=0.0;
	param_arrh[7][2]=3723.05094;
//Lamgmuir Hinshelwood of CO and O	
	param_arrh[8][0]=2.8899E+15;
	param_arrh[8][1]=0.0;
	param_arrh[8][2]=84239.064;
//Langmuir Hinshelwood of O	
	param_arrh[9][0]=1.16326E+15;
	param_arrh[9][1]=0.0;
	param_arrh[9][2]=145167.0024;

*/

rate[0]=0.0000e+00;
rate[1]=1.7275e+06;
rate[2]=5.9089e+06;
rate[3]=0.0000e+00;
rate[4]=0.0000e+00;
rate[5]=4.9863e+06;
rate[6]=3.7740e+06;
rate[7]=5.0245e+06;
rate[8]=0.0000e+00;
rate[9]=0.0000e+00;


//4.Computation of rates
/*
	for(i=0;i<nb_reactions;i++){
         	rate[i]=param_arrh[i][0]*pow(*walltemperature,param_arrh[i][1])*exp(-param_arrh[i][2]/catc_Runiv/ *walltemperature);
   	}
*/

//initial condition for newton raphson

              concentration[5]=0.;
              concentration[6]=0.;
              concentration[7]=0.028*max_concentration;
              concentration[8]=0.018*max_concentration;
              concentration[9]=0.;
              concentration[10]=0.8*max_concentration;

//5.Compute surface coverage
/*	
	printf("\n nb_reactions %d, nb_adspecies %d , *nb_species %d, nb _iterations%d, tolx %e ,tolf %e\n",nb_reactions,nb_adspecies,*nb_species,nb_iterations, tolx, tolf);
for (i=0;i<11;i++){
	printf(" %e",concentration[i]);
}
printf("\n ");

for (i=0;i<10;i++){
	for (j=0;j<3;j++){
		printf("sd %d",stochiodiff[i][j]);
	}
		printf("\n ");
	
}
for (i=0;i<10;i++){
	for (j=0;j<3;j++){
		printf("s %d",stochio[i][j]);
	}
		printf("\n ");
	
}
for (i=0;i<10;i++){
		printf("s %e",rate[i]);	
}

//printf("\n h1");getchar();fflush(stdout);
*/
      	catc_mnewt_(nb_reactions,nb_adspecies,*nb_species,nb_iterations,tolx,tolf,rate,concentration,actual_species,stochio,stochiodiff,max_concentration);
//printf("\n h2");getchar();fflush(stdout);
//6.computation of surface production rate
	for(i=0;i<nb_reactions;i++){
//corrected		for(j=0;j<(*nb_species+*nb_species);j++){
		for(j=0;j<(*nb_species+nb_adspecies);j++){
         		production[i][j]=1.;			
	 	}
	}



	for(sb=0;sb<*nb_species;sb++){
		for(r=0;r<nb_reactions;r++){
			for(sw=0;sw<(*nb_species+nb_adspecies);sw++){
                                             production[r][sb]=production[r][sb]*pow(concentration[sw],stochio[r][sw]);

                 	}
			production[r][sb]=stochiodiff[r][sb]*rate[r]*production[r][sb];
        	}

    	}
//CO catalycity based on impinging CO
// LH
catalycities[0]=production[8][0]/M_inc[2];
// ER
catalycities[1]=(production[5][0]+production[6][0])/M_inc[2];

//O catalycity based on sum of impinging O + 2 x O2
//LH
catalycities[2]=production[9][1]/(2*M_inc[1]+M_inc[3]);
//ER
catalycities[3]=production[7][1]/(2*M_inc[1]+M_inc[3]);



	for(i=0;i<*nb_species;i++){
         	wallmoleproduction[i]=0.0;
	}

	for(sb=0;sb<*nb_species;sb++){
		for(r=0;r<nb_reactions;r++){
                	wallmoleproduction[sb]=wallmoleproduction[sb]+production[r][sb];
		}
//               printf("\n wall %e",wallmoleproduction[sb]);
	}

//7.transfer surface coverage to in/output vector
	for(i=0;i<nb_adspecies;i++){
                    adsorbed_concentration[i]=concentration[i+*nb_species];
	}

return;
}
