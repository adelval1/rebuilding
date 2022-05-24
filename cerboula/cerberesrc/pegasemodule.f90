subroutine PEGASEMODULE(T)

USE EXPERIMENTS_par
USE FLAGS_par
USE HEATFLUX_var
USE OPTIONS_par

USE interf_thermo
USE interf_traco
USE global_thermo
USE global_equilibrium
USE global_traco

implicit none

real*8 T
real*8 TNEQ(1:4),eps
integer flg_guess,flg_log,flg_neq,niter
real*8 zo

real*8, allocatable :: x_guess(:),x(:),ndens(:)
real*8, allocatable :: cpr(:),cpv(:),cpe(:)
integer oper

real*8, allocatable :: h_spc_mole(:)
real*8  ent(1:8),enth(1:8)
integer k

integer mode
real*8 molarmass_rebuild

allocate(x_guess(1:nspecies))
allocate(x(1:nspecies))

allocate(ndens(1:nspecies))

allocate(cpr(1:nspecies))
allocate(cpv(1:nspecies))
allocate(cpe(1:nspecies))

allocate(h_spc_mole(1:nspecies))

TNEQ = T	!Thermal equillibrium computation
flg_guess = 0	!No initial guess for the mole fraction

x_guess = 1
eps = 1.d-12	!Tolerance in the iteration
flg_log = 0	!Printing out a convergence history in the unit 66
niter = 0	!Iteration up to convergence

call MOLEFRAC(x(1:nspecies),Pstatic,T,zo,TNEQ,eps,flg_log,flg_anha,flg_neq,&
& flg_stop,flg_termo,flg_guess,x_guess(1:nspecies),niter)

ndens = Pstatic/(kuniv*T)*x

cpr  = 0.
cpv  = 0.
cpe  = 0.
oper = 0

h_spc_mole = 0.
do k = 1,nspecies
  ent = 0.
  call SPECIES_ENTHALPY(k,TNEQ(1),TNEQ,flg_anha,flg_mode,flg_neq,flg_stop,&
  &flg_termo,ent(1:8))
  h_spc_mole(k) = ent(1)+ent(8)
enddo

call SET_TRACO(x+1.d-32,ndens,cpr,cpv,cpe,h_spc_mole,TNEQ,oper)

visco_rebuild = VISCOSITY()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mode = 1	!1 to use mole fraction input

molarmass_rebuild = MIXTURE_MOLARMASS(x(1:nspecies),mode,flg_stop)
call MOLEFRAC_TO_MASSFRAC(x(1:nspecies),molarmass_rebuild,concent_rebuild(1:nspecies))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

mode = 2	!2 to obtain the enthalpy per unit mass
call MIXTURE_ENTHALPY(T,TNEQ,concent_rebuild(1:nspecies),flg_anha,mode,&
&flg_neq,flg_stop,flg_termo,enth)

h_rebuild = enth(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

rho_rebuild = MIXTURE_DENSITY(Pstatic,T,x(1:nspecies),flg_stop)

deallocate(x_guess)
deallocate(x)

deallocate(cpr)
deallocate(cpv)
deallocate(cpe)

deallocate(ndens)
deallocate(h_spc_mole)

end subroutine PEGASEMODULE
