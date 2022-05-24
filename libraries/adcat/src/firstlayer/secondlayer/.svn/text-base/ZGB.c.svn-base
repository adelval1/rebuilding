#include <math.h>
#include <stdio.h>


//basic catalysis model specific data
#define nb_adspecies 6					//number of adsorbed species(including empty sites)
#define nb_reactions 3					//number of gas surface interactio reactions
#define max_concentration 1.8457E-05			//maximum concentration on the surface in [mol/m^2]
#define tolf 1.e-20					//convergence criterium for newton method
#define tolx 1.e-20					//convergence criterium for newton method
#define nb_iterations 200				//convergence criterium for newton method

void catc_zgb_(int *nb_species,double adsorbed_concentration[nb_adspecies],double *walltemperature,double gas_concentration[*nb_species],double molarmass[*nb_species],double wallmoleproduction[*nb_species]){

/*-------------------------------------------------------------------------------
This routines provides the boundary condition for the species equation in a 
chemical non equilibrium computation. It computes the mass production rate at
the wall applying a reaction rate based model for CO2.

Name of the Model:
ZGB(Ziff Gulari Barshad)

Publication:
(not pubished yet)

status:
	tested, published

output variables:
	double wallmoleproduction[*nb_species*nb]				vector of double of size 5(CO2,O2,CO,O,C),wallproduction rate in [mol/m^3]
in/output variables
	double adsorbed_concentration[nb_adspecies]	vector of double of size 6(CO2s,O2s,COs,Os,Cs,V), concentration on the surface in [mole/m^2], serves only as initial condition for efficiency, dummy values are allowed
input variables:
	double *walltemperature				double, walltemperature in [K]
	double gas_concentration[*nb_species]		vector of double of size 5(CO2,O2,CO,O,C), concentration in the gas phase in [mol/m^3]
	double molarmass[*nb_species]			vector of double of size 5(CO2,O2,CO,O,C), molar mass of species in [kg/mol]

implementation done by: jan thoemel, March 07


watch out: seems to be buggy

to do:
make working for arcitrary co

-------------------------------------------------------------------------------*/



//1.Declaration of local variables
      	int i,j,sw,sb,r;
      	double alpha1,alpha2,Y;
	double concentration[*nb_species+nb_adspecies];
	int actual_species[*nb_species+nb_adspecies];
	int stochio[nb_reactions][*nb_species+nb_adspecies];
	int stochiodiff[nb_reactions][*nb_species+nb_adspecies];
	double param_arrh[nb_reactions][3];
	double rate[nb_reactions];
	double production[nb_reactions][*nb_species+nb_adspecies];
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

	
//3.Further catalysis model specific data
//        actual_species[0]=0;
        actual_species[3]=0;
	actual_species[4]=0;
	
        actual_species[5]=0;
        actual_species[6]=0;
        actual_species[9]=0;


//!reaction 1
//!      CO+s->COs


         stochio[0][2]=1;
         stochio[0][10]=1;
         stochiodiff[0][2]=-1;
         stochiodiff[0][10]=-1;
         stochiodiff[0][7]=1;

//!reaction 2
//!       O2+2s->2Os

         stochio[1][1]=1;
         stochio[1][10]=2;
         stochiodiff[1][1]=-1;
         stochiodiff[1][10]=-2  ;       
         stochiodiff[1][8]=2 ;        



//!reaction 3
//!       COs+Os->CO2+2s
         stochio[2][7]=1;
         stochio[2][8]=1;
         stochiodiff[2][7]=-1;         
         stochiodiff[2][8]=-1;         
         stochiodiff[2][0]=1;         
         stochiodiff[2][10]=2;         
	
	
	param_arrh[0][0]=23063898.1;
	param_arrh[0][1]=0.00;
	param_arrh[0][2]=5.71E+003;

	param_arrh[1][0]=2.64748E+12;
	param_arrh[1][1]=0.0;
	param_arrh[1][2]=5.71E+003;

	param_arrh[2][0]=6.26651E+11;
	param_arrh[2][1]=0.0;
	param_arrh[2][2]=-5.71E+003;

/*
	param_arrh[0][0]=1.1751E+07;
	param_arrh[0][1]=0.00;
	param_arrh[0][2]=0.00;
	
	param_arrh[1][0]=1.3492E+12;
	param_arrh[1][1]=0.0;
	param_arrh[1][2]=0.0;
	
	param_arrh[2][0]=1.2296E+12;
	param_arrh[2][1]=0.0;
	param_arrh[2][2]=0.0;

*/

//4.Computation of rates
	for(i=0;i<nb_reactions;i++){
         	rate[i]=param_arrh[i][0]*pow(*walltemperature,param_arrh[i][1])*exp(-param_arrh[i][2]/8.314/ *walltemperature);
   	}

//initial condition for newton raphson
       Y=concentration[2]/(concentration[2]+concentration[1]);

       if(Y<=0.4){
              concentration[5]=0.;
              concentration[6]=0.;
              concentration[7]=0.028874593*max_concentration;	
              concentration[8]=0.771125407*max_concentration;
              concentration[9]=0.;
              concentration[10]=max_concentration-concentration[7]-concentration[8];
	}
       else if(Y<=0.68){
              concentration[5]=0.;
              concentration[6]=0.;
              concentration[7]=0.0002*exp(12.431*Y)*max_concentration;
              concentration[8]=(-0.0002*exp(12.431*Y)+0.8)*max_concentration;
              concentration[9]=0.;
              concentration[10]=max_concentration-concentration[7]-concentration[8];
	}

       else if (Y<=1.0){
              concentration[5]=0.;
              concentration[6]=0.;
              concentration[7]=0.9*max_concentration;	
              concentration[8]=0.1*max_concentration;
              concentration[9]=0.;
              concentration[10]=max_concentration-concentration[7]-concentration[8];
       }
       else{printf("error in CATC ZGB.c line 74");getchar();}
       
       
//5.Compute surface coverage

      	catc_mnewt_(nb_reactions,nb_adspecies,*nb_species,nb_iterations,tolx,tolf,rate,concentration,actual_species,stochio,stochiodiff,max_concentration);


	for(i=0;i<nb_adspecies;i++){
		if(concentration[i+*nb_species]<1.E-6*max_concentration){concentration[i+*nb_species]=0.;}
	}
	
//6.computation of surface production rate
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

//7.transfer surface coverage to in/output vector
	for(i=0;i<nb_adspecies;i++){
                    adsorbed_concentration[i]=concentration[i+*nb_species];
	}
return;
}
