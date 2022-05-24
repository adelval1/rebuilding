subroutine write_sol

USE OUTEREDGE_var
USE HEATFLUX_var
USE EXPERIMENTS_par
USE OPTIONS_par
USE REBUILDING_var

IMPLICIT NONE

real * 8 heatwall_1,heatwall_cond_1,heatwall_diff_1
real * 8 T_1
INTEGER jj
character*80 dummytxt

external system

real*8 enth_wall,temp_wall

T_1 = T_ext

CALL HEATFLUX(T_1,heatwall_1,heatwall_cond_1,heatwall_diff_1)
  
OPEN(97,FILE='neboulaoutput/stagsol.dat',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',&
  &         ACTION='READ')  
   read(97,*) dummytxt
   read(97,*) temp_wall, enth_wall
CLOSE(97)    

OPEN(97,FILE='cerbereoutput/rebuilding.out',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',&
  &         ACTION='WRITE')

   WRITE(97,'(A2,3X,F20.10)') 'Te',T_ext
   WRITE(97,'(A2,3X,F20.10)') 'Pe',Pstatic
   WRITE(97,'(A2,3X,F20.10)') 'he',h_rebuild
   WRITE(97,'(A4,3X,F20.10)') 'rhoe',rho_rebuild
   WRITE(97,'(A3,3X,F20.10)') 'mue',visco_rebuild
   do jj=1,nspecies
      WRITE(97,'(A2,3X,F20.10)') 'ce',concent_rebuild(jj)
   enddo
   WRITE(97,'(A2,3X,F20.10)') 'Tw',temp_wall
   WRITE(97,'(A2,3X,F20.10)') 'hw',enth_wall
   WRITE(97,'(A2,3X,F20.10)') 'Vs',Vs_ext
   WRITE(97,'(A2,3X,F20.10)') 'Kp',Kp
   WRITE(97,'(A5,3X,F20.10)') 'gam_w',gamma
   if (choice.eq.2) then
      WRITE(97,'(A10,3X,F20.10)') 'gamma_smpl',result_gam
   endif
CLOSE(97)


END subroutine write_sol
