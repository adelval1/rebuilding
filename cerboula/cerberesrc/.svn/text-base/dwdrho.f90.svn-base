subroutine DWDRHO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine calls a pegase routine to compute the derivative of              !!
!! the species mass production with respect to the species partial density 'rho' !!
!! This subroutine requires 'flowfield.plt' provided by a NEBOULA computation    !!
!! Output: dW_drho_profile.dat in the cerbereoutput directory                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use FLAGS_par
use BL_par
use OPTIONS_par

use interf_chemco

implicit none
character*80 dummy
integer i,j,k,l,flg_neq
real*8,allocatable ::  x(:),T(:),rho(:),Tneq(:)
real*8,allocatable ::  mol_frac(:),part_dens(:)
real*8 brol
real*8,allocatable ::  dWdrho_profile(:),T_jac(:),rho_jacp(:,:),rho_jacm(:,:)

allocate(x(1:npoints))
allocate(T(1:npoints))
allocate(rho(1:npoints))
allocate(mol_frac(1:nspecies))
allocate(part_dens(1:nspecies))
allocate(Tneq(1:4))
allocate(T_jac(1:nspecies))
allocate(dWdrho_profile(1:nspecies))
allocate(rho_jacp(1:nspecies,1:nspecies))
allocate(rho_jacm(1:nspecies,1:nspecies))

call chemcodef
flg_neq = 0
OPEN(96,FILE="neboulaoutput/flowfield.plt",STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='READ')
   OPEN(97,FILE="cerbereoutput/dW_drho_profile.dat",STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')
   do i = 1,25
      read(96,*) dummy
   enddo
   do i = 1,npoints
      read(96,*) x(i),brol,brol,rho(i),T(i),brol,brol,brol,brol,brol,brol,brol,brol,mol_frac(1),mol_frac(2),mol_frac(3),mol_frac(4),mol_frac(5),mol_frac(6),mol_frac(7)
      do l = 1,nspecies
         part_dens(l) = mol_frac(l)*rho(i)
      enddo
      Tneq = T(i) !Equilibrium computation

      call jac_term_den(T(i),part_dens,rho_jacp,rho_jacm,T_jac,Tneq,flg_anha,flg_neq,flg_stop,flg_termo)

      write(97,*) "dW/drho at x = ",x(i)
      do j = 1,nspecies
         do k = 1,nspecies
            dWdrho_profile(k) = rho_jacm(j,k)+rho_jacp(j,k)
         enddo
         write(97,*) dWdrho_profile
      enddo
      write(97,*) ""
   enddo
   CLOSE(97)
CLOSE(96)


deallocate(x)
deallocate(T)
deallocate(rho)
deallocate(mol_frac)
deallocate(part_dens)
deallocate(Tneq)
deallocate(T_jac)
deallocate(dWdrho_profile)
deallocate(rho_jacp)
deallocate(rho_jacm)

end subroutine DWDRHO


