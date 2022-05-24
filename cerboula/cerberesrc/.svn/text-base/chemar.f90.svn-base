subroutine CHEMAR(gam)

USE BL_par
USE PEG_DIR_par
USE OPTIONS_par

implicit none

real*8 gam
integer jj,ii,kk,j
character*80 boundfile,brol
character*150 dummytext,dummydata
integer ll,system


jj=scan(mixname,'.',.false.)

boundfile=mixname(1:jj-1)//'.bc'
boundfile='../boundcond/'//trim(boundfile)

OPEN(96,FILE=boundfile,status='unknown', &
& ACCESS='SEQUENTIAL',ACTION='READWRITE')

OPEN(95,FILE='neboulainput/wall_chemistry.in',status='unknown', &
& ACCESS='SEQUENTIAL',ACTION='READWRITE')

write(95,'(a47)') 'Matrix for the boundary conditions in Scott form.'
write(95,'(a9)') 'Matrix nu'

read(96,'(a47)') dummytext !Reading of 'matrix for the boundary...'
read(96,'(a47)') dummytext !Reading of 'Matrix nu'

do ii=1,nreacwall
  read(96,'(a80)') dummydata
  write(95,'(a80)') dummydata
enddo

do kk=1,nreacwall
   read(96,'(a80)') dummytext
   write(95,'(a30,1X,I1)') 'Matrix mu of the reaction number',kk
   do ii=1,nspecies
     read(96,'(a80)') dummydata
     write(95,'(a80)') dummydata 
   enddo
enddo

CLOSE(95)
CLOSE(96)

!/////////////////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////////////////

OPEN(96,FILE='neboulainput/bc_wall.in',status='old', &
& ACCESS='SEQUENTIAL',ACTION='READWRITE')

OPEN(95,FILE='neboulainput/bc_wall.new',status='unknown', &
& ACCESS='SEQUENTIAL',ACTION='READWRITE')

do j=1,3
   read(96,'(A80)') brol
   if (j.eq.3) then
      write(95,'(a80)') brol
      if (nreacwall.eq.1) then
         write(95,*) gam
      elseif (nreacwall.eq.2) then
         write(95,*) gam,gam
      elseif (nreacwall.eq.3) then
         write(95,*) gam,gam,gam
      elseif (nreacwall.eq.4) then
         write(95,*) gam,gam,gam,gam
      elseif (nreacwall.eq.5) then
         write(95,*) gam,gam,gam,gam,gam
      endif
   else
      write(95,'(A80)') brol
   endif
enddo

CLOSE(95)
CLOSE(96)

ll = system("mv neboulainput/bc_wall.new neboulainput/bc_wall.in")

end subroutine CHEMAR
