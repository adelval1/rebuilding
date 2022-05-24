#include <physicalconstants.h>
#include<math.h>
void catc_catalycitymoleproduction_(double  concentrationreactant,double molarmassreactant,double Twall,double catalycity,double *reactantconsumption){
/*************************************************************************************
name: Subroutine scottcat_1
author: Jan Thoemel
date: October 2006
remark: computes mole production rate
input:
	concentrationreactant-concentration of reactant in [mol/m^3]
	molarmassreactant-molecular mass,[kg/particle]
	Twall-temperature of the wall,[K]
	catalycity-recombination coefficient,[1]
output:
	*reactantconsumption-moleproduction at the wall, [mol/m^2]
status:
	tested	
*************************************************************************************/
	catalycity=2*catalycity/(2-catalycity);
	*reactantconsumption=-catalycity*concentrationreactant*sqrt(catc_Runiv*Twall/(2.*catc_pi*molarmassreactant));

}
