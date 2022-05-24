subroutine PROPFROZ(T,yout,h)

!*****************************************************************************
!*****************************************************************************
!Inputs: T,P,nspecies
!Outputs: h,y,rho
!*****************************************************************************
!* Given T,P  ---------->> enthalpy (h) , mass fraction (y), density (rho)
!  for thermochemical equilibrium conditions
!*****************************************************************************
!*****************************************************************************

!*****************************************************************************
!T: Thermal Equilibrium Temperature
!Plocal: Chamber Pressure
!h: Enthalpy of the mixture (=trans+internal+formation)
!y: mass fractions of the species
!rho: Density of the mixture
!*****************************************************************************
!x: Mole fractions of the species
!enth(1:8): Enthalpy of the mixture (components)
!*****************************************************************************

USE interf_thermo
USE interf_traco
USE global_thermo

USE FLAGS_par
USE OPTIONS_par

implicit none

real*8 yout(1:nspecies)
real*8 T,h,rho

real*8 TNEQ(1:4),eps
integer :: flg_log,flg_neq,flg_guess,niter
real*8 zo

integer mode
real*8 enth(1:8)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!MOL FRACTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

TNEQ=T	      !equilibrium computation
flg_guess=0   !no initial guess for the mole fractions
eps=1.d-12    !tolerance in the iteration
flg_log=0     !printing out a convergence history in the unit 66
niter=1000	      !iteration up to convergence

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!mode=1 !1: to obtain the enthalpy per unit mole
mode=2  !2: to obtain the enthalpy per unit mass

! x/y: depending on per unit mole/per unit mass'

 call MIXTURE_ENTHALPY(T,TNEQ,yout(1:nspecies),flg_anha,mode,flg_neq,&
&   flg_stop,flg_termo,enth)
h=enth(1)

end subroutine PROPFROZ
