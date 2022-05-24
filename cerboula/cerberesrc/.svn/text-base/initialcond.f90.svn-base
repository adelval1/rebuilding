subroutine INITIALCOND

USE BL_par
USE OPTIONS_par

integer I,II,system,jj
character * 300 textcopy
real * 8 , allocatable:: dummy(:)

allocate(dummy(1:nspecies))

open(89,file='neboulaoutput/stagsol.dat',status='old')
open(88,file='neboulainput/bc_out.in',status='old')
open(87,file='neboulainput/bc_out.new',status='unknown')

do I=1,9
   read(88,'(A300)') textcopy
   write(87,'(A300)') textcopy
enddo

read(88,*) dummy(1:nspecies)
write(87,*) dummy(1:nspecies)

100  format ((5X,E16.8))

close(88)
close(87)

open(88,file='neboulainput/bc_ini.in',status='old')
open(87,file='neboulainput/bc_ini.new',status='unknown')

READ(89,'(A300)') textcopy
WRITE(87,'(A300)') textcopy

do I=1,npoints
   read(89,*) (dummy(jj),jj=1,2)
   write(87,*) (dummy(jj),jj=1,2)
enddo

read(89,'(A300)') textcopy
write(87,'(A300)') textcopy

do I=1,npoints
   read(89,*) dummy(1:nspecies)
   write(87,*) dummy(1:nspecies)
enddo

101  format ((5X,E16.8))

close(87)
close(88)
close(89)

II = system("mv neboulainput/bc_out.new neboulainput/bc_out.in")
II = system("mv neboulainput/bc_ini.new neboulainput/bc_ini.in")

deallocate(dummy)

end subroutine INITIALCOND
