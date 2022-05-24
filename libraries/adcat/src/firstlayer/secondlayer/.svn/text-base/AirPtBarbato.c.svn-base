#include<stdio.h>
#include <math.h>
#include <physicalconstants.h>

//basic catalysis model specific data
#define nb_adspecies 6					//number of adsorbed species(including empty sites)
#define nb_reactions 10                               //number of gas surface interactio reactions
				
#define max_concentration 3.32107e-6			//maximum concentration on the surface in [mol/m^2]
#define tolf 1.e-20					//convergence criterium for newton method
#define tolx 1.e-20					//convergence criterium for newton method
#define nb_iterations 2000				//convergence criterium for newton method
#define S0_o 1.47
#define s0abar_o 0.8
#define ba_o 0.00045
#define T0a_o 300
#define s0mbar_o 0.4
#define bm_o 0.0006
#define T0m_o 300
#define ar_o 0.085
#define	br_o 0.0002
#define T0r_o 300
#define E4_o 120
#define E5_o 309
#define S0_n 147
#define s0abar_n 0.55
#define ba_n 0.25
#define T0a_n 300
#define s0mbar_n 000000000000000
#define bm_n 0000000000000000
#define T0m_n 000000000000000
#define ar_n 0.2
#define	br_n 0.0001
#define T0r_n 300
#define E4_n 175
#define E5_n 560
#define d_o 0000000000000000
#define d_n 0000000000000000


void catc_airptbarbato_(int *nb_species, double moleproduction[*nb_species],double catalycities[4],double adsorbed_concentration[nb_adspecies],double *walltemperature,double gas_concentration[*nb_species],double molarmass[*nb_species]){
/*-------------------------------------------------------------------------------
This routines provides the boundary condition for the species equation in a 
chemical non equilibrium computation. It computes the mass production rate at
the moleproduction applying a reaction rate based model for Air based on Brabato

source:

\bibitem{b:brbm}  Barbato, M., Reggiani, S., Bruno, C., Muylaert, J.,
"Model for Heterogenous Catalysis on Metal Surfaces with Application to Hypersonic Flows,"
\textit{Journal of Thermophysics and Heat Transfer}, Vol. 14, No. 3 July-September 2000, pp. 412-420

status:
NOT IMPLEMENTED!



output variables:
	double moleproduction[*nb_species]				vector of double of size 5(CO2,O2,CO,O,C),moleproductionproduction rate in [mol/m^3]
in/output variables
	double adsorbed_concentration[nb_adspecies]	vector of double of size 6(CO2s,O2s,COs,Os,Cs,V), concentration on the surface in [mole/m^2], serves only as initial condition for efficiency, dummy values are allowed
input variables:
	double T					double, *walltemperature in [K]
	double gas_concentration[*nb_species]		vector of double of size 5(CO2,O2,CO,O,C), concentration in the gas phase in [mol/m^3]
	double molarmass[*nb_species]			vector of double of size 5(CO2,O2,CO,O,C), molar mass of species in [kg/mol]

implementation done by: jan thoemel, August 07


-------------------------------------------------------------------------------*/


//1.Declaration of local variables
      	int i,j,sw,sb,r;
      	double alpha1,alpha2;
	double concentration[*nb_species+nb_adspecies];
	int actual_species[*nb_species+nb_adspecies];
	int stochio[nb_reactions][*nb_species+nb_adspecies];
	int stochiodiff[nb_reactions][*nb_species+nb_adspecies];
	double rate[nb_reactions];
	double production[nb_reactions][*nb_species+nb_adspecies];
	double M_inc[*nb_species];
	double s0m_o,s0m_n,s0a_o,s0a_n,P0r_o,P0r_n;


printf("STOP, this routine is not implemented");fflush(stdout);getchar();

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
      		M_inc[i]=concentration[i]*sqrt(catc_Runiv* *walltemperature/(2.0*catc_pi*molarmass[i]));
   	}	
	
//3.Further catalysis model specific data

//order for species O2,N2,NO,O,N,O2s,N2s,NOs,Os,Ns,V
//                  0  1  2  3 4 5    6  7   8   9  10 

	actual_species[2]=0;
	actual_species[5]=0;
	actual_species[6]=0;
	actual_species[7]=0;


//not double checked!!!!!!!!!!!!!!!!!!!
//reaction 1,Adsorption of O,O+s->Os
         stochio[0][3]=1;
         stochio[0][10]=1;
         stochiodiff[0][3]=-1;
         stochiodiff[0][10]=-1;       
         stochiodiff[0][8]=1;        
//reaction 2,Desorption of O,Os->O+s
         stochio[1][8]=1;
         stochiodiff[1][8]=-1;
         stochiodiff[1][3]=1;
         stochiodiff[1][10]=1;
//reaction 3,Adsorption of N,N+s->Ns
         stochio[2][4]=1;
         stochio[2][10]=1;
         stochiodiff[2][4]=-1;
         stochiodiff[2][10]=-1;         
         stochiodiff[2][9]=1; 
//reaction 4,Desorption of N,Ns->N+s
         stochio[3][8]=1;
         stochiodiff[3][8]=-1;         
         stochiodiff[3][4]=1;         
         stochiodiff[3][10]=1;         	         
//reaction 5,ER of O and Os,O+Os->O2+s
         stochio[4][3]=1;
         stochio[4][8]=1;
         stochiodiff[4][3]=-1;         
         stochiodiff[4][8]=-1;         
         stochiodiff[4][0]=1;         
         stochiodiff[4][10]=1;         
//reaction 6,ER of N with N,N+Ns->N2+s
         stochio[5][4]=1;
         stochio[5][9]=1;
         stochiodiff[5][4]=-1;
         stochiodiff[5][9]=-1  ;       
         stochiodiff[5][1]=2 ;        
         stochiodiff[5][10]=2 ;        
//reaction 7,LH of O,2Os->O2 +2s
         stochio[6][8]=2;
         stochiodiff[6][8]=-2;         
         stochiodiff[6][0]=1;
         stochiodiff[6][10]=2;         
//reaction 8,Dissociative Adsorption of O2,O2+2s->2Os
         stochio[7][0]=1;
         stochio[7][10]=2;
         stochiodiff[7][0]=-1;         
         stochiodiff[7][10]=-2;         
         stochiodiff[7][8]=2;         
//reaction 9,LH of N,2Ns->N2+2s
         stochio[8][9]=2;
         stochiodiff[8][9]=-2;
         stochiodiff[8][1]=1;
         stochiodiff[8][10]=2;
//reaction 10,Dissociative Adsorption of N2,N2+2s->2Ns
         stochio[9][8]=1;
         stochio[9][10]=2;	 
         stochiodiff[9][8]=-1; 
         stochiodiff[9][10]=-2;         
         stochiodiff[9][9]=2;        
	 
//3.Computation of rates

//reaction 1,Adsorption of O,O+s->Os
	if(*walltemperature<T0a_o){
		s0a_o=s0abar_o;
	}
	else{
		s0a_o=s0abar_o*exp(-ba(*walltemperature-T0a_o));
	}		
	rate[0]=s0a_o*sqrt(catc_Runiv**walltemperature/(2*catc_pi*molarmass[3]))*(1/S0_o);
//reaction 2,Desorption of O,Os->O+s
	rate[1]=catc_Runiv**walltemperature/catc_planckconst*exp(-E5_o/(catc_Runiv**walltemperature));
	
//reaction 3,Adsorption of N,N+s->Ns
	if(*walltemperature<T0a_n){
		s0a_n=s0abar_n;
	}
	else{
		s0a_n=s0abar_n*exp(-ba(*walltemperature-T0a_n));
	}		
	rate[2]=s0a_n*sqrt(catc_Runiv**walltemperature/(2*catc_pi*molarmass[3]))*(1/S0_n);

//reaction 4,Desorption of O,Os->O+s
	rate[3]=catc_Runiv**walltemperature/catc_planckconst*exp(-E5_n/(catc_Runiv**walltemperature));
	
//reaction 5,ER of O and Os,O+Os->O2+s
	if(*walltemperature<T0r_o){
		P0r_o=ar_o;
	}
	else{
		P0r_o=ar_o*exp(br_o*(*walltemperature-T0r_o));
	}
	rate[4]=P0r_o*sqrt(catc_Runiv**walltemperature/(2*catc_pi*molarmass[3]))*(1/S0_o);

//reaction 6,ER of N with N,N+Ns->N2+s
	if(*walltemperature<T0r_n){
		P0r_n=ar_n;
	}
	else{
		P0r_n=ar_n*exp(br_n*(*walltemperature-T0r_n));
	}
	rate[5]=P0r_n*sqrt(catc_Runiv**walltemperature/(2*catc_pi*molarmass[3]))*(1/S0_n);

//reaction 7,LH of O,2Os->O2 +2s
	rate[6]=sqrt(catc_pi*catc_Runiv**walltemperature/(2*molarmass[3]))*d_o*exp(-E4_o/(catc_Runiv**walltemperature));

//reaction 8,Dissociative Adsorption of O2,O2+2s->2Os
	if(*walltemperature<T0m_o){
		s0m_o=s0mbar_o;
	}
	else{
		s0m_o=s0mbar_o*exp(-ba_o*(*walltemperature-T0m_o));
	}
	rate[7]=s0m_o*sqrt(catc_Runiv**walltemperature/(2*catc_pi*molarmass[0]))*(1/S0_o*S0_o);
	
//reaction 9,LH of N,2Ns->N2+2s
	rate[8]=sqrt(catc_pi*catc_Runiv**walltemperature/(2*molarmass[3]))*d_n*exp(-E4_n/(catc_Runiv**walltemperature));


//reaction 10,Dissociative Adsorption of N2,N2+2s->2Ns
	if(*walltemperature<T0m_n){
		s0m_n=s0mbar_n;
	}
	else{
		s0m_n=s0mbar_n*exp(-ba_n*(*walltemperature-T0m_n));
	}
	rate[9]=s0m_n*sqrt(catc_Runiv**walltemperature/(2*catc_pi*molarmass[0]))*(1/S0_n*S0_n);




//initial condition for newton raphson

              concentration[5]=0.;
              concentration[6]=0.;
              concentration[7]=0.1*max_concentration;
              concentration[8]=0.1*max_concentration;
              concentration[9]=0.;
              concentration[10]=0.8*max_concentration;

//4.Compute surface coverage
      	catc_mnewt_(nb_reactions,nb_adspecies,*nb_species,nb_iterations,tolx,tolf,rate,concentration,actual_species,stochio,stochiodiff,max_concentration);

//5.computation of surface production rate
	for(i=0;i<nb_reactions;i++){
		for(j=0;j<(*nb_species+*nb_species);j++){
         		production[i][j]=1.;			
	 	}
	}

	for(sw=0;sw<(*nb_species+nb_adspecies);sw++){
//	printf("\ngasc %e",concentration[sw]);
	}
	


	for(sb=0;sb<*nb_species;sb++){
		for(r=0;r<nb_reactions;r++){
			for(sw=0;sw<(*nb_species+nb_adspecies);sw++){
                                             production[r][sb]=production[r][sb]*pow(concentration[sw],stochio[r][sw]);

                 	}
			production[r][sb]=stochiodiff[r][sb]*rate[r]*production[r][sb];
        	}

    	}
//CO catalycity bosed on impinging CO
// LH
catalycities[0]=production[8][0]/M_inc[2];
// ER
catalycities[1]=(production[5][0]+production[6][0])/M_inc[2];

//O catalycity bosed on sum of impinging O + 2 x O2
//LH
catalycities[2]=production[9][1]/(2*M_inc[1]+M_inc[3]);
//ER
catalycities[3]=production[7][1]/(2*M_inc[1]+M_inc[3]);



	for(i=0;i<*nb_species;i++){
         	moleproduction[i]=0.0;
	}

	for(sb=0;sb<*nb_species;sb++){
		for(r=0;r<nb_reactions;r++){
                	moleproduction[sb]=moleproduction[sb]+production[r][sb];
		}
//               printf("\n moleproduction %e",moleproduction[sb]);
	}

//6.transfer surface coverage to in/output vector
	for(i=0;i<nb_adspecies;i++){
                    adsorbed_concentration[i]=concentration[i+*nb_species];
	}
}


