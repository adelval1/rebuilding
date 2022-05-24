subroutine INITIALDISTRIB(Te)

!*****************************************************************************
!*****************************************************************************
!Inputs: irebuilding,iabacus,nspecies,Te
!Outputs: files 'bc_out.in' and 'bc_ini.in'
!*****************************************************************************
!* Reading of Tw, Pe
!* Guess of temperature distribution across the boundary layer
!* Enthalpy and concentrations for equilibrium conditions given by T, P
!* Writes down the distribution of temperature, enthalpy and concentrations 
!  across the boundary layer in the file 'bc_ini.in'  
!* Completes with blanks the file 'bc_ini.in'
!*****************************************************************************
!*****************************************************************************

USE BL_par
USE OPTIONS_par
USE EXPERIMENTS_par

implicit none

integer i

real * 8, allocatable:: Tdist(:),enthdistrib(:),concentdistrib(:,:),y(:),yout(:)
real*8 dummy,Te,h,rho

real*8 pi

allocate(Tdist(1:npoints))
allocate(enthdistrib(1:npoints))
allocate(concentdistrib(1:nspecies,1:npoints))
allocate(y(1:nspecies))
allocate(yout(1:nspecies))

pi=4d0*atan(1.d0)


open(29,file='neboulainput/bc_out.in')

   write(29,'(a30)') '   xgrid,  rbody,  rhoext,  muext'
   write(29,*) '*******************************'
   write(29,'(a30)') '   uext,  hext,  pext'
   write(29,*) '*******************************'
   write(29,'(a30)') '   duedx'
   write(29,*) 0.
   write(29,'(a30)') 'Data to build alphae and bl.l. thickness'
   write(29,*) '*******************************'

   write(29,'(a30)') 'Concent_ext (PEGASE ORDER)'

   write(29,*) '*******************************'
   100  FORMAT ((E14.8,5X))
close(29)

open(29,file='neboulainput/bc_ini.in')
   write(29,'(a28)') 'T, h at the stagnation point'


   if (neq_fr.eq.0) call PROPEQUIL(Te,h,yout,rho)

   do i=1,npoints
      Tdist(i)=Tw+(Te-Tw)*sin(.5d0*pi*float(i-1)/float(npoints-1)) !Guess for the distribution &
! & across the b.l. for the temperature

      if (neq_fr.eq.1) then
        call PROPEQUIL(Tdist(i),h,y,rho)
        concentdistrib(1:nspecies,i)=y(1:nspecies)
      elseif (neq_fr.eq.0) then
        call PROPFROZ(Tdist(i),yout,h)
        concentdistrib(1:nspecies,i)=yout(1:nspecies)
      endif

      enthdistrib(i)=h
   enddo

   do i=1,npoints
     write(29,'(2(E14.8,1X))') Tdist(i),enthdistrib(i)
   enddo

   write(29,'(a31)') 'Species at the stagnation point'

   do i=1,npoints
      write(29,*) concentdistrib(1:nspecies,i)
   enddo

   101  FORMAT ((E14.8,5X))

close(29)

deallocate(Tdist)
deallocate(enthdistrib)
deallocate(concentdistrib)
deallocate(y)
deallocate(yout)

end subroutine INITIALDISTRIB
