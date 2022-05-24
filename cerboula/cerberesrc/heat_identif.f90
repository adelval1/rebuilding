subroutine HEAT_IDENTIF(gamm,Qw,j)

USE OPTIONS_par
USE BL_par
USE EXPERIMENTS_par

implicit none

real*8 gamm,Qw,initialWallTemperature
integer i,II,system,j
real*8 dummy



open(97,file='neboulainput/bc_wall.in',status='replace')

if (j<=3) then
      write(*,*) 'Wall temperature initial',j/3.*(Tw_smpl-350.)+350.
      write(97,*) 'Wall temperature initial'
      write(97,*) j/3.*(Tw_smpl-350.)+350.
else
      write(*,*) 'Wall temperature',Tw_smpl
      write(97,*) 'Wall temperature'
      write(97,*) Tw_smpl
endif

      write(97,*) 'Gamma vector'
      if (nreacwall.eq.1) then
         write(97,*) gamm
      elseif (nreacwall.eq.2) then
         write(97,*) gamm,gamm
      elseif (nreacwall.eq.3) then
         write(97,*) gamm,gamm,gamm
      elseif (nreacwall.eq.4) then
         write(97,*) gamm,gamm,gamm,gamm
      elseif (nreacwall.eq.5) then
         write(97,*) gamm,gamm,gamm,gamm,gamm
      endif
close(97)


II = system("../../cerboula/exe/neboula.exe")

open(93,file='neboulaoutput/heatskin.dat')
   read(93,*) dummy,dummy,dummy,dummy,dummy,dummy,dummy,Qw
close(93)
II = system("rm neboulaoutput/heatskin.dat")


call INITIALCOND

end subroutine HEAT_IDENTIF
