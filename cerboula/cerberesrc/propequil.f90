subroutine PROPEQUIL(T,h,y,rho)

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
!molarmass_rebuild: Dummy argument. Molar mass of the mixture.
!*****************************************************************************

USE interf_thermo
USE interf_traco
USE global_thermo

USE FLAGS_par
USE OPTIONS_par
USE EXPERIMENTS_par

implicit none

real*8, allocatable:: x_guess(:),x(:)

real*8 y(1:nspecies)
real*8 T,h,rho

real*8 TNEQ(1:4),eps
integer :: flg_log,flg_neq,flg_guess,niter
real*8 zo

integer mode
real*8 enth(1:8),molarmass_rebuild

allocate(x_guess(1:nspecies))
allocate(x(1:nspecies))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!MOL FRACTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

TNEQ=T	      !equilibrium computation
flg_guess=0   !no initial guess for the mole fractions
x_guess=1     !
eps=1.d-12    !tolerance in the iteration
flg_log=0     !printing out a convergence history in the unit 66
niter=1000    !iteration up to convergence


call MOLEFRAC(x(1:nspecies),Pstatic,T,zo,TNEQ,eps,flg_log,flg_anha,flg_neq,&
&    flg_stop,flg_termo,flg_guess,x_guess(1:nspecies),niter)


mode=1	!1: to use mole fraction input
molarmass_rebuild=mixture_molarmass(x(1:nspecies),mode,flg_stop) !molarmass

call MOLEFRAC_TO_MASSFRAC(x(1:nspecies),molarmass_rebuild,y(1:nspecies))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!mode=1 !1: to obtain the enthalpy per unit mole
mode=2  !2: to obtain the enthalpy per unit mass

! x/y: depending on per unit mole/per unit mass'

call MIXTURE_ENTHALPY(T,TNEQ,y(1:nspecies),flg_anha,mode,flg_neq,&
&   flg_stop,flg_termo,enth)
h=enth(1)

rho=mixture_density(Pstatic,T,x(1:nspecies),flg_stop)

deallocate(x)
deallocate(x_guess)

end subroutine PROPEQUIL
