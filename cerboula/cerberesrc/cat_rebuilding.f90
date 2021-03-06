subroutine CAT_REBUILDING(flag_OK)

USE OUTEREDGE_var
USE REBUILDING_var
USE EXPERIMENTS_par
USE ICP_par
USE HEATFLUX_var

implicit none

integer i,nsteps,II,system,label,j
real*8 exponen_min,exponen_max,errorlim,error
real*8,allocatable:: exponen(:),gam(:)
real*8 Qw_num,Qw_experim,Boltz_ct
real*8 Ve_ext_2

real*8 T_wall,gam_init,tmp
integer flag_OK

nsteps = 40
allocate(exponen(1:nsteps))
allocate(gam(1:nsteps))


T_wall    = Tw_smpl
Boltz_ct = 5.67d-8
Qw_experim = Boltz_ct*emissivity_smpl*Tw_smpl**4.


call NUMERICAL_PARAMETER(T_ext,rho_ext)
call BC_WALL(T_wall)
Ve_ext_2 = Vs_ext*NDP_5_smpl! Vs_ext hasn't changed since the enthalpy is the same at the outer edge of the DYNAMIC bl.
call BC_OUT(rho_rebuild,visco_rebuild,h_rebuild,u1e_smpl,u1y_smpl,Ve_ext_2,delta_smpl,ve_smpl,concent_rebuild,Rm_smpl)

II = system("rm neboulaoutput/heatskin.dat")

 
exponen_min = -5.d0
exponen_max = 0.d0

errorlim = 1d-6
j = 0
j = j+1

gam_init = 10**exponen_max

write(*,*) '---------------------------------------------------------------'
write(*,*) '---------------------------------------------------------------'
write(*,*) 'CATALYCITY VALUE (between 0 and 1)',gam_init,j


call HEAT_IDENTIF(gam_init,Qw_num,j)

if (Qw_num.lt.Qw_experim) then
  write(*,*) '                                                    '
  write(*,*) '                                                    '
  write(*,*) 'THE EXPERIMENTAL Qw IS GREATER THAN THE MAX POSSIBLE'
  write(*,*) '                                                    '
  write(*,*) '                                                    '
  label = 30
  goto 30
endif


10 write(*,*) 'Starting new iteration'
do i = 2,nsteps
  j = j+1
  if (j.gt.nsteps.and.i.eq.nsteps) then
     exponen_max = exponen(i-1)
     goto 10
  endif

  tmp = exponen_max-(exponen_max-exponen_min)/float(nsteps-1)*dfloat(i-1)
  exponen(i) = tmp
  gam(i) = 10**exponen(i)
 
  write(*,*) '---------------------------------------------------------------'
  write(*,*) '---------------------------------------------------------------'
  write(*,*) 'CATALYCITY VALUE (between 0 and 1)',gam(i),j

  call HEAT_IDENTIF(gam(i),Qw_num,j)
  error = abs(Qw_num-Qw_experim)/Qw_experim

  write(*,*) ''
  write(*,*) 'heat', Qw_num, 'error',error

  if (error.le.errorlim) then
    label =20
    goto 20
  endif

  if (Qw_num.le.Qw_experim) then
    if (i.eq.2) then
       exponen_min = exponen(i)
       goto 10
    else
       exponen_max = exponen(i-1)
       exponen_min = exponen(i)
       goto 10
    endif
  endif

  if (i.eq.nsteps) then 
    write(*,*) '                                                    '
    write(*,*) '                                                    '
    write(*,*) 'THE EXPERIMENTAL Qw IS SMALLER THAN THE MIN POSSIBLE'
    write(*,*) '                                                    '
    write(*,*) '                                                    '
    label = 30
    goto 30
  endif

enddo

20 if (label.eq.20) then
  write(*,*) 'CONVERGED CASE'
endif

if (label.eq.20) then
  result_gam = gam(i)
  write(*,*) '                                                    '
  write(*,*) '                                                    '
  write(*,'(5x,a7,1x,F17.6)') 'Tw',Tw
  write(*,'(5x,a7,1x,F16.5)') 'Qw_experim',Qw_experim
  write(*,'(5x,a7,1x,F12.10)')'GAMMA_w',result_gam
  write(*,*) '                                                    '
  write(*,*) '                                                    '
  flag_OK = 1
endif

30 if (label.eq.30) then
   write(*,*) 'NO SOLUTION'
   flag_OK = 0
endif


deallocate(gam)
deallocate(exponen)

end subroutine CAT_REBUILDING
