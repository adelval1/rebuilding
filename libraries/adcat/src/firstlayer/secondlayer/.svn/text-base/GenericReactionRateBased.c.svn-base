#include <math.h>
#include <stdio.h>

//basic definitions for the catalysis model
#define tolf 1.e-10					//convergence criterium for newton method
#define tolx 1.e-10					//convergence criterium for newton method
#define nb_iterations 10000				//convergence criterium for newton method

void catc_genericreactionratebased_(double T,int nb_species,int nb_adspecies,int nb_reactions,double max_concentration,int stochio[nb_reactions][nb_species+nb_adspecies],int stochiodiff[nb_reactions][nb_species+nb_adspecies],double WALL[nb_species],double adsorbed_concentration[nb_adspecies],double gas_concentration[nb_species],double param_arrh[nb_reactions][3],int actual_species[nb_species+nb_adspecies],double molarmass[nb_species]){
/*-------------------------------------------------------------------------------
This routines provides the boundary condition for the species equation in a 
chemical non equilibrium computation. It computes the mass production rate at
the wall applying a reaction rate based model for an arbitrary mixture provided suitable input.

Name of the Model:
GenericReactionRatebased

Publication:
(Generic)
status:
	not checked yet

output variables:
	double WALL[nb_species]				vector of double of size 5(CO2,O2,CO,O,C),wallproduction rate in [mol/m^3]
in/output variables
	double adsorbed_concentration[nb_adspecies]	vector of double of size 6(CO2s,O2s,COs,Os,Cs,V), concentration on the surface in [mole/m^2], serves only as initial condition for efficiency, dummy values are allowed
input variables:
	double T					double, walltemperature in [K]
	double gas_concentration[nb_species]		vector of double of size 5(CO2,O2,CO,O,C), concentration in the gas phase in [mol/m^3]
	double molarmass[nb_species]			vector of double of size 5(CO2,O2,CO,O,C), molar mass of species in [kg/mol]
	int nb_species
	int nb_adspecies
	int nb_reactions,
	double max_concentration
	int stochio[nb_reactions][nb_species+nb_adspecies]
	int stochiodiff[nb_reactions][nb_species+nb_adspecies]
	double param_arrh[nb_reactions][3]
	int actual_species[nb_species+nb_adspecies]

implementation done by: jan thoemel, March 07

to do:
convert that variables are all pointers

-------------------------------------------------------------------------------*/


//1. Declare Local variables
      	int i,j,sw,sb,r;
      	double alpha1,alpha2;
	double concentration[nb_species+nb_adspecies];
	
	double rate[nb_reactions];
	double production[nb_reactions][nb_species+nb_adspecies];

//2. Initilization of local variables
         for(i=0;i<nb_species;i++){concentration[i]=gas_concentration[i];}
         for(i=0;i<nb_adspecies;i++){concentration[i+nb_species]=adsorbed_concentration[i]; }

	for(i=0;i<nb_reactions;i++){	
		for(j=0;j<nb_species+nb_adspecies;j++){	
        	 	stochio[i][j]=0;
         		stochiodiff[i][j]=0;
		}
	}
	for(i=0;i<nb_reactions;i++){
		for(j=0;j<(nb_species+nb_species);j++){
         		production[i][j]=1.;			
	 	}
	}
	for(i=0;i<nb_species;i++){
         	WALL[i]=0.0;
	}
	


//3. Computation of reaction rates fromthe Arrhenius parameters
	for(i=0;i<nb_reactions;i++){
         	rate[i]=param_arrh[i][0]*pow(T,param_arrh[i][1])*exp(-param_arrh[i][2]/8.314/T);
   	}

//4. Compute surface coverage
      	CATC_mnewt(nb_reactions,nb_adspecies,nb_species,nb_iterations,tolx,tolf,rate,concentration,actual_species,stochio,stochiodiff,max_concentration);

//5. Compute surface production rate per reaction and per species
	for(sb=0;sb<nb_species;sb++){
		for(r=0;r<nb_reactions;r++){
			for(sw=0;sw<(nb_species+nb_adspecies);sw++){
                                             production[r][sb]=production[r][sb]*pow(concentration[sw],stochio[r][sw]);
                 	}
                                production[r][sb]=stochiodiff[r][sb]*rate[r]*production[r][sb];
        	}
    	}

//6. sum up to obtain surface production rate per species
	for(sb=0;sb<nb_species;sb++){
		for(r=0;r<nb_reactions;r++){
                	WALL[sb]=WALL[sb]+production[r][sb];
		}
               
	}     

//7. transfer surface concentration to in/output vector
	for(i=0;i<nb_adspecies;i++){
                    adsorbed_concentration[i]=concentration[i+nb_species];
	}
return;
}
