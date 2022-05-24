subroutine EQUILCONDITION(hex,pe,Teq,mue,concent,rhoe)

USE interf_thermo
USE interf_traco
USE global_thermo 

USE CONVERGENCE_ENTH_par
USE FLAGS_par
USE EXPERIMENTS_par
USE OPTIONS_par
USE INIT_DIST_par

implicit none

real*8 hex,pe,h2,Teq,mue,rhoe,rho
real*8 concent(1:nspecies)

real*8 TNEQ(1:4),eps
integer :: flg_log,flg_neq,flg_guess,niter
real*8 zo

integer mode
real*8 enth(1:8),ent(1:8)
real*8 molarmass_rebuild,enth_rebuild,residual
real*8 y(1:nspecies),T,T_1,T_2,h_1,h_2,deltaT,denth_dt

real * 8, allocatable :: cpr(:),cpv(:),cpe(:)
real * 8, allocatable :: x_guess(:),x(:),ndens(:)
real * 8, allocatable :: h_spc_mole(:)

integer icount
integer k,oper

allocate(cpr(1:nspecies))
allocate(cpv(1:nspecies))
allocate(cpe(1:nspecies))

allocate(x_guess(1:nspecies))
allocate(x(1:nspecies))

allocate(ndens(1:nspecies))
allocate(h_spc_mole(1:nspecies))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
T=T_0
call PROPEQUIL(T,enth_rebuild,y(1:nspecies),rho)

deltaT=0.
icount=0

residual=enth_rebuild-hex
do while((abs(residual)/hex).gt.enthalpyerror)
  icount=icount+1
  
  T=T+enthalpha*deltaT
  T_1=T
  call PROPEQUIL(T_1,h_1,y,rho)
  
  enth_rebuild=h_1
  residual=enth_rebuild-hex 
  
  T_2=T*(1.+enthepsilon)
  call PROPEQUIL(T_2,h_2,y,rho)
  
  denth_dt=(h_2-h_1)/T/enthepsilon
  deltaT=-residual/denth_dt
enddo

Teq=T

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

TNEQ=T	      !equilibrium computation
flg_guess=0   !no initial guess for the mole fractions
x_guess=1     !
eps=1.d-12     !tolerance in the iteration
flg_log=0     !printing out a convergence history in the unit 66
niter=0	      !iteration up to convergence

call MOLEFRAC(x(1:nspecies),pe,Teq,zo,TNEQ,eps,flg_log,flg_anha,flg_neq,&
&    flg_stop,flg_termo,flg_guess,x_guess(1:nspecies),niter)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ndens = Pstatic / (kuniv * T) * x

cpr=0.
cpv=0.
cpe=0.
oper=0

h_spc_mole=0.
DO k=1, nspecies
    ent = 0.
    CALL SPECIES_ENTHALPY (k,TNEQ(1),TNEQ,flg_anha,flg_mode,flg_neq,flg_stop,&
    &flg_termo,ent(1:8))
    h_spc_mole(k) = ent(1) + ent(8)
ENDDO  
call SET_TRACO(x+1.d-32,ndens,cpr,cpv,cpe,h_spc_mole,TNEQ,oper)

mue=viscosity()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

mode=1	!1: to use mole fraction input
molarmass_rebuild=mixture_molarmass(x,mode,flg_stop) !molarmass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call MOLEFRAC_TO_MASSFRAC(x(1:nspecies),molarmass_rebuild,y(1:nspecies))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

concent(1:nspecies)=y(1:nspecies)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rhoe=mixture_density(Pe,Teq,x(1:nspecies),flg_stop)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

deallocate(cpr)
deallocate(cpv)
deallocate(cpe)

deallocate(x_guess)
deallocate(x)

deallocate(ndens)
deallocate(h_spc_mole)

end subroutine EQUILCONDITION
