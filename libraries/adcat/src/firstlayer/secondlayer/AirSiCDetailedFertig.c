#include<stdio.h>
#include <math.h>
#include <physicalconstants.h>

//basic catalysis model specific data
#define nb_adspecies 6					//number of adsorbed species(including empty sites)
#define nb_reactions 20                                //number of gas surface interactio reactions
				
#define max_concentration 3.32107e-6			//maximum concentration on the surface in [mol/m^2]
#define tolf 1.e-20					//convergence criterium for newton method
#define tolx 1.e-20					//convergence criterium for newton method
#define nb_iterations 2000				//convergence criterium for newton method
double catc_vibrationalpartitionfunction(double ,double,double);
double catc_translationalpartitionfunction(double , double);
double catc_equilibriumconstant(double,double,double,double,double,double );
double catc_equilibriumconstantadsorption(double,double,double,double );

void catc_airsicdetailedfertig_(int *nb_species, double moleproduction[*nb_species],double catalycities[4],double adsorbed_concentration[nb_adspecies],double *walltemperature,double gas_concentration[*nb_species],double molarmass[*nb_species]){
/*-------------------------------------------------------------------------------
This routines provides the boundary condition for the species equation in a 
chemical non equilibrium computation. It computes the mass production rate at
the moleproduction applying a reaction rate based model for Air based on Fertig.

source:
\bibitem{b:mf}  Fertig, M.,
"Modellierung reaktiver Prozesse auf Siliziumkarbid-Oberfl\"achen on verd\"unten Nichtgleichgewichts-Luftstr\"omungen,"
PhD thesis, University of Stuttgart, Germany, November 2005
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

to do:
make working for arcitrary air
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
//reaction 6,Dissociative Adsorption of O,O2+s->O + Os
         stochio[5][0]=1;
         stochio[5][10]=1;
         stochiodiff[5][0]=-1;         
         stochiodiff[5][10]=-1;         
         stochiodiff[5][3]=1;         
         stochiodiff[5][8]=1;         
//reaction 7,ER of O and Ns,O+Ns->NO+s
         stochio[6][4]=1;
         stochio[6][9]=1;
         stochiodiff[6][4]=-1;         
         stochiodiff[6][9]=-1;         
         stochiodiff[6][2]=1;         
         stochiodiff[6][10]=1;         
//reaction 8,Dissociative Adsorption of NO,NO+s-> O+Ns
         stochio[7][2]=1;
         stochio[7][10]=1;
         stochiodiff[7][2]=-1;         
         stochiodiff[7][10]=-1;         
         stochiodiff[7][2]=1;         
         stochiodiff[7][9]=1;         
//reaction 9,ER of N and Os,N+Os->NO+s
         stochio[8][1]=1;
         stochio[8][8]=1;
         stochiodiff[8][1]=-1;         
         stochiodiff[8][8]=-1;         
         stochiodiff[8][2]=1;         
         stochiodiff[8][10]=1;    
//reaction 10,Dissociative Adsorption of NO,NO+s->N+Os
         stochio[9][2]=1;
         stochio[9][10]=1;
         stochiodiff[9][2]=-1;         
         stochiodiff[9][10]=-1;         
         stochiodiff[9][2]=1;         
         stochiodiff[9][8]=1;         
//reaction 11,ER of N with N,N+Ns->N2+s
         stochio[10][4]=1;
         stochio[10][9]=1;
         stochiodiff[10][4]=-1;
         stochiodiff[10][9]=-1  ;       
         stochiodiff[10][1]=2 ;        
         stochiodiff[10][10]=2 ;        
//reaction 12,Dissociative Adsorption of N2,N2+s->N+Ns
         stochio[11][1]=1;
         stochio[11][10]=1;
         stochiodiff[11][1]=-1;
         stochiodiff[11][10]=-1;
         stochiodiff[11][4]=1;
         stochiodiff[11][9]=1;
//reaction 13,LH of O,2Os->O2 +2s
         stochio[12][8]=2;
         stochiodiff[12][8]=-2;         
         stochiodiff[12][0]=1;
         stochiodiff[12][10]=2;         
//reaction 14,Dissociative Adsorption of O2,O2+2s->2Os
         stochio[13][0]=1;
         stochio[13][10]=2;
         stochiodiff[13][0]=-1;         
         stochiodiff[13][10]=-2;         
         stochiodiff[13][8]=2;         
//reaction 15,LH of O and N,Os+Ns->NO+2s
         stochio[14][8]=1;
         stochio[14][9]=1;
         stochiodiff[14][8]=-1;         
         stochiodiff[14][9]=1;         
         stochiodiff[14][2]=1;         
         stochiodiff[14][10]=2;         
//reaction 16,Dissociative Adsorption of NO,NO+2s->Os + Ns
         stochio[15][2]=1;
         stochio[15][10]=2;
         stochiodiff[15][2]=-1;         
         stochiodiff[15][10]=-2;         
         stochiodiff[15][8]=1;         
         stochiodiff[15][9]=1;         
//reaction 17,LH of N and O,N+Os->NO+2s
         stochio[16][4]=1;
         stochio[16][8]=1;
         stochiodiff[16][4]=-1;         
         stochiodiff[16][8]=-1;         
         stochiodiff[16][2]=1;         
         stochiodiff[16][10]=2;         
//reaction 18,Dissociative Adsorption of NO,NO+2s-> Os+Ns
         stochio[17][2]=1;
         stochio[17][10]=2;
         stochiodiff[17][2]=-1;         
         stochiodiff[17][10]=-2;         
         stochiodiff[17][8]=1;         
         stochiodiff[17][9]=1;         
//reaction 19,LH of N,2Ns->N2+2s
         stochio[18][9]=2;
         stochiodiff[18][9]=-2;
         stochiodiff[18][1]=1;
         stochiodiff[18][10]=2;
//reaction 20,Dissociative Adsorption of N2,N2+2s->2Ns
         stochio[19][8]=1;
         stochio[19][10]=2;	 
         stochiodiff[19][8]=-1; 
         stochiodiff[19][10]=-2;         
         stochiodiff[19][9]=2;        
/*
//reaction 21,mediated Dissociation of O2,O2->2O
         stochio[20][0]=1;
         stochiodiff[20][0]=-1;         
         stochiodiff[20][3]=2;         
//reaction 22,mediated Rekombination of O,2O->O2
         stochio[21][3]=2;
         stochiodiff[21][3]=-2;         
         stochiodiff[21][0]=1;         
//reaction 23,mediated Dissociation of N2,N2->2N
         stochio[22][1]=1;
         stochiodiff[22][1]=-1;
         stochiodiff[22][4]=2;
//reaction 24,mediated Recombination of N,2N-> 2N
         stochio[23][4]=2;
         stochiodiff[23][4]=-2;         
         stochiodiff[23][1]=1;         
//reaction 25,mediated Dissociation of NO,NO->N+O
         stochio[24][2]=1;
         stochiodiff[24][2]=-1;         
         stochiodiff[24][0]=1;         
         stochiodiff[24][1]=1;       
//reaction 26,mediated Combination of NO,N+O->NO
         stochio[25][0]=1;
         stochio[25][1]=1;
         stochiodiff[25][0]=-1;         
         stochiodiff[25][1]=-1;         
         stochiodiff[25][2]=1;         
*/	 
	 
	 
//3.Computation of rates
//reaction 1,Adsorption of O,O+s->Os
rate[0]=0.2*0.8*(1-exp(-2887/ *walltemperature))/max_concentration*sqrt(catc_Runiv* *walltemperature/2./catc_pi/molarmass[3]);
//reaction 2,Desorption of O,Os->O+s
rate[1]=1.*0.2*0.8*(1-exp(-2887/ *walltemperature))*catc_boltzconst/catc_planckconst* *walltemperature/vibrationalpartitionfunction( *walltemperature,1714,28866)*max_concentration/catc_avogadro*exp(-28866/ *walltemperature)/max_concentration; 
//reaction 3,Adsorption of N,N+s->Ns
rate[2]=0.03*0.8*(1-exp(-3007/ *walltemperature))/max_concentration*sqrt(catc_Runiv* *walltemperature/2./catc_pi/molarmass[4]);;
//reaction 4,Desorption of N,Ns->N+s
rate[3]=1.*0.2*0.8*(1-exp(-2887/ *walltemperature))*catc_boltzconst/catc_planckconst* *walltemperature/vibrationalpartitionfunction( *walltemperature,1589,30068)*max_concentration/catc_avogadro*exp(-30068/ *walltemperature)/max_concentration; 

//reaction 5,ER of O and Os,O+Os->O2+s
rate[4]=0.002*(exp(-2045/ *walltemperature)-exp(28886/ *walltemperature))/max_concentration*sqrt(catc_Runiv* *walltemperature/2./catc_pi/molarmass[3]);
//reaction 6,Dissociative Adsorption of O,O2+s->O + Os
rate[5]=0.002*2/2*sqrt(reducedmass(molarmass[3],molarmass[4])/molarmass[3])*catc_equilibriumconstantO(*walltemperature)*catc_equilibriumconstantadsorption(*walltemperature,28866.,1714.,molarmass[3])*(exp(-2045/ *walltemperature)-exp(28886/ *walltemperature))/max_concentration*sqrt(catc_Runiv* *walltemperature/2./catc_pi/molarmass[0]);

//reaction 7,ER of O and Ns,O+Ns->NO+s
rate[6]=0.001125*(exp(-1864/ *walltemperature)-exp(30068/ *walltemperature))/max_concentration*sqrt(catc_Runiv* *walltemperature/2./catc_pi/molarmass[3]);
//reaction 8,Dissociative Adsorption of NO,NO+s-> O+Ns
rate[7]=0.001125*1/2*sqrt(reducedmass(molarmass[3],molarmass[4])/molarmass[3])*catc_equilibriumconstantNO(*walltemperature)*catc_equilibriumconstantadsorption(*walltemperature,28866.,1714.,molarmass[3])*(exp(-1864/ *walltemperature)-exp(30068/ *walltemperature))/max_concentration*sqrt(catc_Runiv* *walltemperature/2./catc_pi/molarmass[2]);

//reaction 9,ER of N and Os,N+Os->NO+s
rate[8]=0.001125*(exp(-1864/ *walltemperature)-exp(28866/ *walltemperature))/max_concentration*sqrt(catc_Runiv* *walltemperature/2./catc_pi/molarmass[3]);
//reaction 10,Dissociative Adsorption of NO,NO+s->N+Os
rate[9]=0.001125*1/2*sqrt(reducedmass(molarmass[3],molarmass[4])/molarmass[4])*catc_equilibriumconstantNO(*walltemperature)*catc_equilibriumconstantadsorption(*walltemperature,30068,1589.,molarmass[4])*(exp(-1864/ *walltemperature)-exp(28866/ *walltemperature))/max_concentration*sqrt(catc_Runiv* *walltemperature/2./catc_pi/molarmass[2]);

//reaction 11,ER of N with N,N+Ns->N2+s
rate[10]=0.001125*(exp(-1864/ *walltemperature)-exp(28866/ *walltemperature))/max_concentration*sqrt(catc_Runiv* *walltemperature/2./catc_pi/molarmass[3]);
//reaction 12,Dissociative Adsorption of N2,N+s->N+Ns
rate[11]=0.001125*1/2*sqrt(reducedmass(molarmass[3],molarmass[4])/molarmass[4])*catc_equilibriumconstantNO(*walltemperature)*catc_equilibriumconstantadsorption(*walltemperature,30068,1589.,molarmass[4])*(exp(-1864/ *walltemperature)-exp(28866/ *walltemperature))/max_concentration*sqrt(catc_Runiv* *walltemperature/2./catc_pi/molarmass[2]);

//reaction 13,LH of O,2Os->O2 +2s
rate[12]=1.*catc_boltzconst/catc_planckconst/catc_vibrationalpartitionfunction(*walltemperature,0.048,14433)*max_concentration/catc_avogadro*exp(-(14433+2045)/ *walltemperature); 
//reaction 14,Dissociative Adsorption of O2,O2+2s->2Os
rate[13]=1.*2/2*catc_boltzconst/catc_planckconst/catc_vibrationalpartitionfunction(*walltemperature,0.048,14433)*catc_equilibriumconstantO(*walltemperature)*pow(catc_equilibriumconstantadsorption(*walltemperature,28866.,1714.,molarmass[3]),2)*max_concentration/catc_avogadro*exp(-(14433+2045)/ *walltemperature)/max_concentration/max_concentration;

//reaction 15,LH of O and N,Os+Ns->NO+2s
rate[14]=1.*catc_boltzconst/catc_planckconst/catc_vibrationalpartitionfunction(*walltemperature,0.048,14433)*max_concentration/catc_avogadro*exp(-(14433+2045)/ *walltemperature); ;
//reaction 16,Dissociative Adsorption of NO,NO+2s->Os + Ns
rate[15]=1.*1/2*catc_boltzconst/catc_planckconst/catc_vibrationalpartitionfunction(*walltemperature,0.048,14433)*catc_equilibriumconstantNO(*walltemperature)*catc_equilibriumconstantadsorption(*walltemperature,28866.,1714.,molarmass[3])*catc_equilibriumconstantadsorption(*walltemperature,30068.,1714.,molarmass[4])*max_concentration/catc_avogadro*exp(-(14433+1864)/ *walltemperature)/max_concentration/max_concentration;

//reaction 17,LH of N and O,N+Os->NO+2s
rate[16]=1.*catc_boltzconst/catc_planckconst/catc_vibrationalpartitionfunction(*walltemperature,0.048,16237)*max_concentration/catc_avogadro*exp(-(16237+1684)/ *walltemperature); 
//reaction 18,Dissociative Adsorption of NO,NO+2s-> Ns+Os
rate[17]=1.*1/2*catc_boltzconst/catc_planckconst/catc_vibrationalpartitionfunction(*walltemperature,0.048,16237)*catc_equilibriumconstantNO(*walltemperature)*catc_equilibriumconstantadsorption(*walltemperature,28866.,1714.,molarmass[3])*catc_equilibriumconstantadsorption(*walltemperature,30068.,1714.,molarmass[4])*max_concentration/catc_avogadro*exp(-(16237+1864)/ *walltemperature)/max_concentration/max_concentration;;

//reaction 19,LH of N,2N+2s->2N+2s
rate[18]=1.*catc_boltzconst/catc_planckconst/catc_vibrationalpartitionfunction(*walltemperature,0.048,16237)*max_concentration/catc_avogadro*exp(-(16237+1684)/ *walltemperature); ;
//reaction 20,Dissociative Adsorption of N2,N2+2s->2N
rate[19]=1.*2/2*catc_boltzconst/catc_planckconst/catc_vibrationalpartitionfunction(*walltemperature,0.048,16237)*catc_equilibriumconstantN(*walltemperature)*pow(catc_equilibriumconstantadsorption(*walltemperature,30068.,1714.,molarmass[4]),2)*max_concentration/catc_avogadro*exp(-(16237+1864)/ *walltemperature)/max_concentration/max_concentration;;;

/*
//reaction 21,mediated Dissocitation of O2,O2->2O
rate[20]=0;
//reaction 22,mediated Rekombination of O,2O->O2
rate[21]=0;
 //reaction 23,mediated Dissociation of N2,N2->2N
rate[22]=0;
//reaction 24,mediated Recombination of N,2N-> 2N
rate[23]=0;
//reaction 25,mediated Dissociation of NO,NO->N+O
rate[24]=0;
//reaction 26,mediated Combination of NO,N+O->NO
rate[25]=0;
*/

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




/*
source for next three subroutines:
\bibitem{b:mf}  Fertig, M.,
"Modellierung reaktiver Prozesse auf Siliziumkarbid-Oberfl\"achen on verd\"unten Nichtgleichgewichts-Luftstr\"omungen,"
PhD thesis, University of Stuttgart, Germany, November 2005
status:
*/

double catc_equilibriumconstantadsorption(double temperature,double characteristicvibrationaltemperature,double characteristictemperature,double molarmass){	
	return	catc_avogadro/max_concentration* catc_vibrationalpartitionfunction(temperature,characteristicvibrationaltemperature,characteristictemperature )/pow(catc_translationalpartitionfunction(temperature, molarmass),3)*exp(characteristictemperature/temperature);
}
double catc_equlibriumconstantN(double temperature){
	return catc_equilibriumconstant(12.622585,1.325,-98.56,-17.4,8.,temperature);
}
double catc_equlibriumconstantO(double temperature){
	return catc_equilibriumconstant(18.945465,-0.988,-61.81,-2.3,-1,temperature);
}
double catc_equlibriumconstantNO(double temperature){
	return catc_equilibriumconstant(13.474639,0.492,-67.61,-9.1,4.,temperature);
}
double catc_equilibriumconstant(double p1,double p2,double p3,double p4,double p5,double temperature){
	return exp(p1+p2*log(temperature/1000)+p3*1000/temperature+p4*(1000/temperature)*(1000/temperature)+p5*(1000/temperature)*(1000/temperature)*(1000/temperature));
}
double catc_vibrationalpartitionfunction(double walltemperature,double characteristicvibrationaltemperature,double characteristictemperature){
	return (1-exp(-characteristictemperature/walltemperature))/(1-exp(-characteristicvibrationaltemperature/walltemperature));
}
double catc_translationalpartitionfunction(double temperature, double molarmass){
	return sqrt(2*catc_pi*molarmass*catc_Runiv*temperature/catc_planckconst*catc_planckconst);
}
double catc_reducedmass(double mass1,double mass2){
	return mass1*mass2/(mass1+mass2);
}
