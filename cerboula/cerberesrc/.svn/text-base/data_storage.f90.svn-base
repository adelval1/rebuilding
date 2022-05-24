subroutine DATA_STORAGE(i,gam,T,h,heatwall_cond,heatwall_diff,rho,visco,Vs,barker)

USE NAMES_par

implicit none

integer step,j,i
real*8 dummy
real*8 catalycity,gam,T,h,heatwall,heatwall_cond,heatwall_diff,rho,visco,Vs,barker
character*80 filename

catalycity = gam

heatwall = heatwall_cond+heatwall_diff

filename = trim('cerbereoutput/')//trim(probe_name1)
OPEN(96,FILE=filename,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='READWRITE')
   if (i.eq.1) then
      write(96,*) catalycity,T,h,heatwall,heatwall_cond,heatwall_diff,rho,visco,Vs,barker
   else
      do j = 1,i-1
         read(96,*) dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy
      enddo
      write(96,*) catalycity,T,h,heatwall,heatwall_cond,heatwall_diff,rho,visco,Vs,barker
   endif
CLOSE(96)

end subroutine DATA_STORAGE
