#include <math.h>
#include <stdio.h>

#define nb_adspecies 6
#define nb_reactions 10
#define max_concentration 2.24E-05

void catc_drw_(int *nb_species,double adsorbed_concentration[nb_adspecies],double *walltemperature,double gas_concentration[*nb_species],double molarmass[*nb_species],double wallmoleproduction[*nb_species]){

/*-------------------------------------------------------------------------------
This routines provides the boundary condition for the species equation in a 
chemical non equilibrium computation. It computes the mass production rate at
the wall applying a reaction rate based model for Air on RCG.

Name of the Model:
DRW (Deutschmann Riedel Warnatz)

Publication:
"Modeling of Nitrogen and Oxygen Recombination on Partial Catalytic Surfaces," Journal of Heat Transfer, 1995

status:
	tested

output variables:
	double wallmoleproduction[*nb_species]				vector of double of size 5(CO2,O2,CO,O,C),wallproduction rate in [mol/m^3]
in/output variables
	double adsorbed_concentration[nb_adspecies]	vector of double of size 6(CO2s,O2s,COs,Os,Cs,V), concentration on the surface in [mole/m^2], serves only as initial condition for efficiency, dummy values are allowed
input variables:
	double *walltemperature					double, walltemperature in [K]
	double gas_concentration[*nb_species]		vector of double of size 5(CO2,O2,CO,O,C), concentration in the gas phase in [mol/m^3]
	double molarmass[*nb_species]			vector of double of size 5(CO2,O2,CO,O,C), molar mass of species in [kg/mol]

implementation done by: jan thoemel, March 07
to do:
make working for arcitrary co
-------------------------------------------------------------------------------*/




//   Local variables
//   ^^^^^^^^^^^^^^^
      	int i,j,sw,sb,r;
      	double alpha1,alpha2;
	double concentration[*nb_species+nb_adspecies];
	int actual_species[*nb_species+nb_adspecies];
	double param_arrh[nb_reactions][3];
	double rate[nb_reactions];
	double production[nb_reactions][*nb_species+nb_adspecies];//,sum,adsorbed_concentration(nb_adspecies),gas_concentration(*nb_species)
		
	double tolf,tolx;

	double stickingcoefficientO;
	double stickingcoefficientN ;     
	int stochio[nb_reactions][*nb_species+nb_adspecies];
	int stochiodiff[nb_reactions][*nb_species+nb_adspecies];

//CHECK
       if(*nb_species!=5){printf("\n\nstop this version of the catalysis library works only with 5 gas species\n\n");fflush(stdout);getchar();}



//   Initilization
//   ^^^^^^^^^^^^^
         for(i=0;i<*nb_species;i++){concentration[i]=gas_concentration[i];}
         for(i=0;i<nb_adspecies;i++){concentration[i+*nb_species]=adsorbed_concentration[i]; }


         stickingcoefficientO=0.1;
         stickingcoefficientN=0.1;
 
//order for species O2,N2,NO,O,N,O2s,N2s,NOs,Os,Ns,V
//                          1  2  3  4 5 6   7   8   9  10 11


	for(i=0;i<nb_reactions;i++){	
		for(j=0;j<*nb_species+nb_adspecies;j++){	
        	 	stochio[i][j]=0;
         		stochiodiff[i][j]=0;
		}
	}


//reaction 1
         stochio[0][3]=1;
         stochio[0][10]=1;
         stochiodiff[0][3]=-1;
         stochiodiff[0][10]=-1;
         stochiodiff[0][8]=1;
//reaction 2
         stochio[1][4]=1;
         stochio[1][10]=1;
         stochiodiff[1][4]=-1;
         stochiodiff[1][10]=-1  ;       
         stochiodiff[1][9]=1 ;        
//!reaction 3
         stochio[2][8]=1;
         stochiodiff[2][8]=-1;         
         stochiodiff[2][3]=1;         
         stochiodiff[2][10]=1;         
//!reaction 4
         stochio[3][9]=1;
         stochiodiff[3][9]=-1;
         stochiodiff[3][4]=1;
         stochiodiff[3][10]=1;         
//!reaction 5
         stochio[4][5]=1;
         stochiodiff[4][5]=-1;
         stochiodiff[4][0]=1;         
         stochiodiff[4][10]=1 ;        
//!reaction 6
         stochio[5][6]=1;
         stochiodiff[5][6]=-1;
         stochiodiff[5][1]=1;
         stochiodiff[5][10]=1;                               
//!reaction 7
         stochio[6][3]=1;
         stochio[6][8]=1;
         stochiodiff[6][3]=-1 ;        
         stochiodiff[6][8]=-1;
         stochiodiff[6][0]=1 ;        
         stochiodiff[6][10]=1;         
//!reaction 8
         stochio[7][4]=1;
         stochio[7][9]=1;
         stochiodiff[7][4]=-1;
         stochiodiff[7][9]=-1;
         stochiodiff[7][1]=1;
         stochiodiff[7][10]=1 ;                                         
//!reaction 9
         stochio[8][8]=2;
         stochiodiff[8][8]=-2;         
         stochiodiff[8][5]=1;
         stochiodiff[8][10]=1;                    
//!reaction 10
         stochio[9][9]=2;
         stochiodiff[9][9]=-2;         
         stochiodiff[9][6]=1;
         stochiodiff[9][10]=1 ;                   


	param_arrh[2][0]=5.100e+11;
	param_arrh[2][1]=0.00;
	param_arrh[2][2]=200.E3;
	
	param_arrh[3][0]=7.300e+11;
	param_arrh[3][1]=0.0;
	param_arrh[3][2]=215.E3;
	
	param_arrh[4][0]=1.0e+12;
	param_arrh[4][1]=0.0;
	param_arrh[4][2]=10.E3;

	param_arrh[5][0]=1.0e+12;
	param_arrh[5][1]=0.00;
	param_arrh[5][2]=10.E3;
	
	param_arrh[6][0]=6.0e+13;
	param_arrh[6][1]=0.0;
	param_arrh[6][2]=60.E3;
	
	param_arrh[7][0]=2.5e+13;
	param_arrh[7][1]=0.0;
	param_arrh[7][2]=40.E3;
	
	param_arrh[8][0]=2.0e+19;
	param_arrh[8][1]=0.0E3;
	param_arrh[8][2]=160.E3;
	
	param_arrh[9][0]=7.0e+17;
	param_arrh[9][1]=0.0;
	param_arrh[9][2]=130.E3;



	
//   Transformation of preexponetial factors in SI units
        for(i=2;i<6;i++){
              param_arrh[i][0] = param_arrh[i][0]/100.;
        }
        for(i=6;i<8;i++){
              param_arrh[i][0] = param_arrh[i][0]/1000000.;
	}
	for(i=8;i<10;i++){
              param_arrh[i][0] = param_arrh[i][0]/10000.;
	}

//   Computation of rates
//   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   	rate[0]=stickingcoefficientO/(1-stickingcoefficientO/2)/max_concentration*sqrt(8.314* *walltemperature/2./3.14/molarmass[0]);
	rate[1]=stickingcoefficientN/(1-stickingcoefficientN/2)/max_concentration*sqrt(8.314* *walltemperature/2./3.14/molarmass[1]);

	for(i=2;i<10;i++){
         	rate[i]=param_arrh[i][0]*pow(*walltemperature,param_arrh[i][1])*exp(-param_arrh[i][2]/8.314/ *walltemperature);
   	}


//   Compute surface coverage
//   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

//order for species O2,N2,NO,O,N,O2s,N2s,NOs,Os,Ns,V
//                          1  2  3  4 5 6   7   8   9  10 11

         
         tolf=1.e-10;
         tolx=1.e-10;

	for(i=0;i<(*nb_species+nb_adspecies);i++){actual_species[i]=1;}
         actual_species[2]=0;
         actual_species[7]=0;

      	catc_mnewt_(nb_reactions,nb_adspecies,*nb_species,100,tolx,tolf,rate,concentration,actual_species,stochio,stochiodiff,max_concentration);

//!   surface production rate
//!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	for(i=0;i<nb_reactions;i++){
		for(j=0;j<(*nb_species+*nb_species);j++){
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

	for(i=0;i<*nb_species;i++){
         	wallmoleproduction[i]=0.0;
	}

	for(sb=0;sb<*nb_species;sb++){
		for(r=0;r<nb_reactions;r++){
                	wallmoleproduction[sb]=wallmoleproduction[sb]+production[r][sb];
		}
               
	}     

	for(i=0;i<nb_adspecies;i++){
                    adsorbed_concentration[i]=concentration[i+*nb_species];
	}

return;
}
