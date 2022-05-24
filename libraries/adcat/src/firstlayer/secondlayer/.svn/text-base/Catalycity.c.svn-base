#include <stdio.h>
void catc_catalycity_(int *cattype,int *nb_species,double *catalycity, double concentration[*nb_species],double  molarmass[*nb_species],double *walltemperature,double wallmoleproduction[*nb_species]){
/*------------------------
to do:
use only pointers
--------------------------*/
	double reactantconsumption;
	int i;
//	printf("adcat %e", *catalycity);
//	getchar();
	if(*nb_species!=5) {printf("*catalycity model valid only for 5 species gas");fflush(stdout);getchar();}
	
	for(i=0;i<*nb_species;i++){wallmoleproduction[i]=0.;}

	switch (*cattype) {
//air5
				 
		case 1:
		 	// Oxygen recombination only, pegase order
			// tested			
			catc_catalycitymoleproduction_(concentration[3],molarmass[3],*walltemperature,*catalycity,&reactantconsumption);
			wallmoleproduction[0]=-reactantconsumption/2;
			wallmoleproduction[3]=reactantconsumption;
			break;
		case 2:
		 	// Nitrogen recombination only, pegase order
			// tested			
			catc_catalycitymoleproduction_(concentration[4],molarmass[4],*walltemperature,*catalycity,&reactantconsumption);
			wallmoleproduction[1]=-reactantconsumption/2;
			wallmoleproduction[4]=reactantconsumption;
			break;
		case 3:
		 	// Oxygen and Nitrogen recombination only using equal *catalycity, pegase order
			// tested
			catc_catalycitymoleproduction_(concentration[3],molarmass[3],*walltemperature,*catalycity,&reactantconsumption);
			wallmoleproduction[0]=-reactantconsumption/2;
			wallmoleproduction[3]=reactantconsumption;
			
			catc_catalycitymoleproduction_(concentration[4],molarmass[4],*walltemperature,*catalycity,&reactantconsumption);
			wallmoleproduction[1]=-reactantconsumption/2;
			wallmoleproduction[4]=reactantconsumption;

			break;
//co2
		case 4:
		 	// co+o->co2
			if(concentration[3]>concentration[2]){
				catc_catalycitymoleproduction_(concentration[2],molarmass[2],*walltemperature,*catalycity,&reactantconsumption);
			}
			else{
				*catalycity=1.;
				catc_catalycitymoleproduction_(concentration[3],molarmass[3],*walltemperature,*catalycity,&reactantconsumption);
			}
			wallmoleproduction[0]=-reactantconsumption;
			wallmoleproduction[2]=reactantconsumption;
			wallmoleproduction[3]=reactantconsumption;
					
			break;
		case 5:
		 	// 2o->O2
			catc_catalycitymoleproduction_(concentration[3],molarmass[3],*walltemperature,*catalycity,&reactantconsumption);
			wallmoleproduction[1]=-reactantconsumption/2;
			wallmoleproduction[3]=reactantconsumption;
			break;
		case 6:
		 	// co+o->co2;2o->O2
			
			if(concentration[3]>concentration[2]){
				printf("-");
				catc_catalycitymoleproduction_(concentration[2],molarmass[2],*walltemperature,*catalycity,&reactantconsumption);
				wallmoleproduction[0]=-reactantconsumption;
				wallmoleproduction[2]=reactantconsumption;
				wallmoleproduction[3]=reactantconsumption;
				catc_catalycitymoleproduction_((concentration[3]-concentration[2]**catalycity),molarmass[2],*walltemperature,*catalycity,&reactantconsumption);
				wallmoleproduction[1]=-reactantconsumption/2;
				wallmoleproduction[3]=wallmoleproduction[3]+reactantconsumption;
				}
			else{
				printf(".");
				*catalycity=1.;
				catc_catalycitymoleproduction_(concentration[3],molarmass[3],*walltemperature,*catalycity,&reactantconsumption);
				wallmoleproduction[0]=-reactantconsumption;
				wallmoleproduction[2]=reactantconsumption;
				wallmoleproduction[3]=reactantconsumption;
				}

			break;
case 7:
		 	// o2->2O;co+o->co2;2o->O2
		 	*catalycity=1;
			catc_catalycitymoleproduction_(concentration[1],molarmass[1],*walltemperature,*catalycity,&reactantconsumption);
				wallmoleproduction[1]=reactantconsumption;
				wallmoleproduction[3]=-reactantconsumption*2;
			
			if((concentration[3]+2*concentration[1]**catalycity)>concentration[2]){
				printf("-");
				catc_catalycitymoleproduction_(concentration[2],molarmass[2],*walltemperature,*catalycity,&reactantconsumption);
				wallmoleproduction[0]=-reactantconsumption;
				wallmoleproduction[2]=reactantconsumption;
				wallmoleproduction[3]=wallmoleproduction[3]+reactantconsumption;
				catc_catalycitymoleproduction_((concentration[3]+2*concentration[1]**catalycity-concentration[2]*2**catalycity/(2-*catalycity)),molarmass[2],*walltemperature,*catalycity,&reactantconsumption);
				wallmoleproduction[1]=-reactantconsumption/2;
				wallmoleproduction[3]=wallmoleproduction[3]+reactantconsumption;
				}
			else{
				printf(".");
				*catalycity=1.;
				catc_catalycitymoleproduction_(concentration[3]+2*concentration[1]**catalycity,molarmass[3],*walltemperature,*catalycity,&reactantconsumption);
				wallmoleproduction[0]=-reactantconsumption;
				wallmoleproduction[2]=reactantconsumption;
				wallmoleproduction[3]=reactantconsumption;
				}


			break;
			
			
		 default:
			printf("\nadcat: boundary condition does not exists.\n");fflush(stdout);
			getchar();
			break;
		}



}
