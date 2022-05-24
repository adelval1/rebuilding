subroutine ENTH_REBUILDING(gam)

!**********************************************************************
!**********************************************************************
!ENTHALPY REBUILDING PROCEDURE SUBROUTINE (given Kp and gamma)
!**********************************************************************
!READING OF THE DATA FROM ICP CODE
!READING OF THE DATA FROM THE EXPERIMENTS
!READING OF THE DATA FOR THE CONVERGENCE OF THE NEWTON LOOP
!READING OF THE INITIAL GUESS FOR THE TEMPERATURE AT THE OUTER EDGE
!
!    T_1 --------------------------> Qw_1
!                                             --->> d(Qw_num-Qw_exp_1)/dT
!    T_2=T_1 * (1 + epsilon) ------> Qw_2
!    NEWTON LOOP:  DeltaT=-(Qw_1-Qw_exp_1) / d(Qw_num-Qw_exp_1)/dT
!
!    T(n+1)=T(n)+alpha * DeltaT
!**********************************************************************
!**********************************************************************

USE interf_thermo
USE interf_traco
USE global_thermo
USE global_equilibrium
USE global_traco

USE CONVERGENCE_par
USE EXPERIMENTS_par
USE OUTEREDGE_var
USE HEATFLUX_var
USE REBUILDING_var
USE OPTIONS_par

implicit none

real*8 T_1,Qw_1,Qwc_1,Qwd_1,T_2,Qw_2,Qwc_2,Qwd_2,gam
real*8 deltaT,error,residual,heat_deriv
character*7 str


call CHEMAR(gam)

T_1=T_reb
deltaT=0d0

error=1000.d0
do while(error.gt.errorstop)
  count_iter=count_iter+1

  T_1=T_1+alpha*deltaT
  call HEATFLUX(T_1,Qw_1,Qwc_1,Qwd_1)
  residual=Qw_1-Qw_exp_1

  write(*,*) '                                       '
  write(*,*) 'ITERATIVE STEP NUMBER',' ',count_iter
  write(*,'(A2,3X,F15.8)') 'Te',T_1
  write(*,'(A2,3X,F20.8)') 'he',h_rebuild
  write(*,'(A2,3X,F20.8)') 'Vs',Vs_ext

  T_ext      =T_1
  rho_ext    =rho_rebuild
  mu_ext     =visco_rebuild
  h_ext      =h_rebuild
  concent_ext=concent_rebuild
  
  T_2=T_1*(1.d0+epsilon)
  call HEATFLUX(T_2,Qw_2,Qwc_2,Qwd_2)

  heat_deriv=(Qw_2-Qw_1)/T_1/epsilon
  deltaT=-residual/heat_deriv
  if (abs(deltaT)>200) then
     if (deltaT > 0) then
        deltaT = 200.d0
     elseif (deltaT < 0) then
        deltaT = -200.d0
     endif
  endif
  error=abs(residual)/Qw_exp_1

  write(*,*) ' '
  write(*,*) 'deltaT',deltaT
  write(*,*) 'ERROR',error

  call INITIALCOND
enddo

T_reb=T_1
Qwc_reb = Qwc_1
Qwd_reb = Qwd_1


END subroutine ENTH_REBUILDING
