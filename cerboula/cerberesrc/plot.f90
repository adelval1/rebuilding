subroutine PLOT(Twall,gam,Qw)

USE BL_par
USE OPTIONS_par

implicit none

real * 8 Twall,gam,Qw
INTEGER II,icount,system
character*80 dummy,dummyB,dummyC,dummyD,dummyE
real * 8 dummy2

OPEN(95,FILE='neboulainput/bc_wall.in',status='OLD', &
& ACCESS='SEQUENTIAL',ACTION='READWRITE',form='FORMATTED')

OPEN(96,FILE='neboulainput/bc_wall.new',status='UNKNOWN', &
& ACCESS='SEQUENTIAL',ACTION='READWRITE',form='FORMATTED')

   do icount=1,3
      read(95,FMT='(A80)') dummy
      if (icount.eq.2) then
	write(96,FMT='(F20.10)') Twall
      else
        write(96,FMT='(A80)') dummy
      endif
   enddo

   if (nreacwall.eq.1) then
      read(95,*) dummy
      write(96,*) dummy
   elseif (nreacwall.eq.2) then
      read(95,*) dummy,dummyB
      write(96,*) dummy,dummyB
   elseif (nreacwall.eq.3) then
      read(95,*) dummy,dummyB,dummyC
      write(96,*) dummy,dummyB,dummyC
   elseif (nreacwall.eq.4) then
      read(95,*) dummy,dummyB,dummyC,dummyD
      write(96,*) dummy,dummyB,dummyC,dummyD
   elseif (nreacwall.eq.5) then
      read(95,*) dummy,dummyB,dummyC,dummyD,dummyE
      write(96,*) dummy,dummyB,dummyC,dummyD,dummyE
   endif


close(95)
close(96)

II=system("mv neboulainput/bc_wall.new neboulainput/bc_wall.in")

call CHEMAR(gam)

II=system("rm neboulaoutput/heatskin.dat")
II=system("../../cerboula/exe/neboula.exe")

open(93,file='neboulaoutput/heatskin.dat')
   read(93,*) dummy2,dummy2,dummy2,dummy2,dummy2,dummy2,dummy2,Qw
close(93)

END subroutine PLOT
