subroutine KP_REBUILDING

USE REBUILDING_var
USE OUTEREDGE_var, ONLY : mu_ext, Vs_ext
USE EXPERIMENTS_par, ONLY : Rm
USE HEATFLUX_var, ONLY : rho_rebuild
USE CONVERGENCE_par

implicit none

real*8 Rey

if(count_iter.gt.0) then
  Rey =  rho_rebuild*Vs_ext*Rm/mu_ext
  Kp = 1.d0 + 6.d0/(Rey + 0.455*dsqrt(Rey))
endif

write(*,*) '                        '
write(*,*) 'Kp',Kp
write(*,*) '                        '

end subroutine KP_REBUILDING
