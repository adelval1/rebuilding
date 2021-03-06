subroutine HEATFLUX(T,Qw,Qwc,Qwd)

USE HEATFLUX_var
USE REBUILDING_var
USE EXPERIMENTS_par
USE ICP_par
USE OUTEREDGE_var

implicit none

real*8 T,Qw,Qwc,Qwd,dummy
integer I,system,jj


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!		CALL TO PEGASE LIBRARY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call PEGASEMODULE(T)

Vs_ext = sqrt(2d0*deltaPs_1/Kp/rho_rebuild)	!en supposant que rho_s = rho_reb
write(*,*) 'rho_rebuild=', rho_rebuild, 'Vs_ext=', Vs_ext
Ve_ext = Vs_ext*NDP_5_1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!		INPUTS FOR THE VKI BL CODE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call NUMERICAL_PARAMETER(T,rho_rebuild)
call BC_WALL(Tw)
call BC_OUT(rho_rebuild,visco_rebuild,h_rebuild,u1e_1,u1y_1,Ve_ext,delta_1,ve_1,concent_rebuild,Rm)

I = system("rm neboulaoutput/heatskin.dat")
I = system("../../cerboula/exe/neboula.exe")

open(93,file='neboulaoutput/heatskin.dat')
  read(93,*) dummy,dummy,dummy,dummy,dummy,dummy,dummy,Qw,dummy,dummy,Qwc,Qwd
close(93)

end subroutine HEATFLUX
