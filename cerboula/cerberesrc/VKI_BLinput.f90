subroutine NUMERICAL_PARAMETER(T_1,rho)

implicit none

real*8 T_1,rho
character*80 dummy
integer icount,II,system

OPEN(95,FILE='neboulainput/numerical_parameter.in',status='OLD', &
& ACCESS='SEQUENTIAL',ACTION='READWRITE',form='FORMATTED')

OPEN(96,FILE='neboulainput/numerical_parameter.new',status='UNKNOWN', &
& ACCESS='SEQUENTIAL',ACTION='READWRITE',form='FORMATTED')

icount=1

do while(icount.le.45)
  read(95,FMT='(A80)') dummy
   
   if(icount.eq.18) then
     write(96,FMT='(F20.10)') T_1
   elseif(icount.eq.20) then
     write(96,FMT='(F20.10)') rho
   elseif(icount.eq.22) then
     write(96,FMT='(F20.10)') 5000000.
   else
     write(96,FMT='(A80)') dummy
   endif

  icount=icount+1
enddo

close(95)
close(96)

II=system("mv neboulainput/numerical_parameter.new neboulainput/numerical_parameter.in")


end subroutine NUMERICAL_PARAMETER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BC_OUT(rho,visco,enthal,NDP_u1e,NDP_u1y,V_e,NDP_delta,NDP_ve,conc,R_probe)

USE BL_par
USE EXPERIMENTS_par
USE OPTIONS_par

implicit none

real*8 rho,visco,enthal,conc(1:nspecies),V_torch
real*8 NDP_u1e,NDP_u1y,NDP_delta,NDP_ve,R_probe
real*8 delta2,duedx,duedy,V_e
integer II,system
integer nsize,icount
character*300 dummy

OPEN(95,FILE='neboulainput/bc_out.in',status='OLD', &
& ACCESS='SEQUENTIAL',ACTION='READWRITE',form='FORMATTED')

OPEN(96,FILE='neboulainput/bc_out.new',status='UNKNOWN', &
& ACCESS='SEQUENTIAL',ACTION='READWRITE',form='FORMATTED')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! V_torch: torch exit velocity                                           !!
!! V_s: Free stream velocity                                              !!
!! V_e: Velocity at the outer edge of the finite thickness boundary layer !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

delta2 = NDP_delta*R_probe
V_torch=V_e/NDP_ve
duedx = V_torch*NDP_u1e/R_probe         !V_torch is used for the NDP in ICP
duedy = NDP_u1y*V_torch/R_probe/R_probe

icount=1
do while(icount.le.10)
  read(95,FMT='(A300)') dummy
   
   if (icount.eq.2) then
     write(96,FMT='(F13.8,1X,F13.8,1X,F15.10,1X,F15.10)') &
     &0.,0.,rho,visco
   elseif(icount.eq.4) then
     write(96,FMT='(F10.5,1X,F20.10,1X,F20.10)') 0.,enthal,pstatic
   elseif(icount.eq.6) then
     write(96,FMT='(F20.10)') duedx
   elseif(icount.eq.8) then
     write(96,FMT='(F13.8,1X,F20.10,1X,F20.10)') -V_e,-duedy,delta2
   elseif(icount.eq.10) then
     write(96,*) conc(1:nspecies)
   else
     write(96,FMT='(A300)') dummy
   endif
  100  FORMAT ((E20.15,10X))

  icount=icount+1
enddo

close(95)
close(96)

II=system("mv neboulainput/bc_out.new neboulainput/bc_out.in")

end subroutine BC_OUT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BC_WALL(T_wall)

USE BL_par

implicit none

real*8 T_wall
character*80 dummy,dummy2,dummy3,dummy4,dummy5
integer icount,II,system

OPEN(95,FILE='neboulainput/bc_wall.in',status='OLD', &
& ACCESS='SEQUENTIAL',ACTION='READWRITE',form='FORMATTED')

OPEN(96,FILE='neboulainput/bc_wall.new',status='UNKNOWN', &
& ACCESS='SEQUENTIAL',ACTION='READWRITE',form='FORMATTED')

do icount=1,4
   if(icount.eq.2) then
     read(95,FMT='(A80)') dummy
     write(96,FMT='(F20.10)') T_wall
   elseif (icount.eq.4) then
      if (nreacwall.eq.1) then
         read(95,*) dummy
         write(96,*) dummy
      elseif (nreacwall.eq.2) then
         read(95,*) dummy,dummy2
         write(96,*) dummy,dummy2
      elseif (nreacwall.eq.3) then
         read(95,*) dummy,dummy2,dummy3
         write(96,*) dummy,dummy2,dummy3
      elseif (nreacwall.eq.4) then
         read(95,*) dummy,dummy2,dummy3,dummy4
         write(96,*) dummy,dummy2,dummy3,dummy4
      elseif (nreacwall.eq.5) then
         read(95,*) dummy,dummy2,dummy3,dummy4,dummy5
         write(96,*) dummy,dummy2,dummy3,dummy4,dummy5
      endif
   else
     read(95,FMT='(A80)') dummy
     write(96,FMT='(A80)') dummy
   endif
enddo

close(95)
close(96)

II=system("mv neboulainput/bc_wall.new neboulainput/bc_wall.in")

end subroutine BC_WALL

