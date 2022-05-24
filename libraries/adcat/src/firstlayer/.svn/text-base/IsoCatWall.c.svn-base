//void catc_invfickdiffusion_(int *nb_species,double moleproduction[*nb_species],double concentration[*nb_species],double D[*nb_species][*nb_species],double speciesmolegradients[*nb_species]){

void catc_invfickdiffusion_(int *nb_species,int diffusiontype,int cattype, double concentrationwall,double speciesconcentrationwallnew[*nb_species],double speciesconcentrationgas[*nb_species]){
/*
input
-number of species
-species mole fraction
-mole priduction of catalysis
output
-species gradients

to do:
*/
int i;
double speciesconcentrationupdate[*nb_species],speciesconcentrationwallold[*nb_species];

while((speciesconcentrationwallnew-speciesconcentrationwallold)/concentrationwall>0.1){

//	catc_ficklaw();
//	catc_ficklawjac();
//	catc_solve(*nb_species,speciesconcentrationupdate);

	for(i=0;i<*nb_species;i++){
		speciesconcentrationwallnew[i]=speciesconcentrationwallold[i]-speciesconcentrationupdate[i];

	}

	switch(diffusiontype){
	case 1:
	; 	
//		catc_ficklaw();
		

//		catc_ficklawjac();
	}	

}


}

	

void catc_solve(int nb_species,double concentrationupdate[nb_species]){
}
void catc_ficklaw(){
}
void catc_ficklawjac(){
}

/*
	int i,j;
	double sum;
	double efficientdiffusioncoefficient[*nb_species];
	double mixtureconcentration=0;
	
//0.	
	for(i=0;i<*nb_species;i++){
		mixtureconcentration=mixtureconcentration+concentration[i];
	}

//1.calculate efficient diffusion coefficients
	for(i=0;i<*nb_species;i++){
        	sum = 0.0;
         	for(j=*nb_species-1;j>=0;j--){
            		if(j!=i){ sum=sum+concentration[j]/mixtureconcentration/D[i][j];}
         	}
         efficientdiffusioncoefficient[i] = (1.0-concentration[i])/mixtureconcentration/sum;
      	}	


//2. calculate species gradient

	for(i=0;i<*nb_species;i++){
		speciesmolegradients[i]=moleproduction[i]/efficientdiffusioncoefficient[i];
	
	}

*/
