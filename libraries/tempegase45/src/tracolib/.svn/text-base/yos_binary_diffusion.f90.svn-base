SUBROUTINE diffusion_coefficients_multi(p,TNEQ,mm,d)
!
    USE global_thermo
    USE global_traco
    USE traco_ord_var, only: x_ord
    IMPLICIT NONE
!
    REAL(kind=8) p,TNEQ(4)
    REAL(kind=8) dij(nsp,nsp),d(nsp),d_ord(nsp),den,mm,mr
    INTEGER i,j
!
!   First, call a subroutine that provides the true array of diffusion
!   coefficients (the routine is at the bottom of this file)
!
    d_ord=0.0d0
    CALL set_binary_diffusion(p,TNEQ,dij)
!
!   Now compute the array of equivalent binary diffusion coefficients in
!   the multicomponent mixture
!
    DO i=1,nsp
        den=0.0d0
        DO j=1,nsp
            IF(i.EQ.j) CYCLE
            den=den+x_ord(j)/dij(i,j)
        END DO
        mr=MMOL_ord(i)/mm
        d_ord(i)= mr * ( 1.0d0-x_ord(i)*mr ) / den
!       d_ord(i)= ( 1.0d0-x_ord(i) ) / den
    END DO
!
    CALL Reorder_Data (nsp,index_array, d_ord, d)

END SUBROUTINE diffusion_coefficients_multi
!
!*********************************************************************
!
SUBROUTINE set_binary_diffusion(p,TNEQ,dij)
!
    USE global_thermo, only: nsp,KUNIV
    USE global_traco
    IMPLICIT NONE
    REAL(kind=8) :: p,TNEQ(4),dij(nsp,nsp),dijp(nsp,nsp)

    INTEGER :: i,j
    REAL(kind=8) :: T,Te
!
!   Initialize useful temperatures
!
    T=TNEQ(1)
    IF (nsp.EQ.nneut) THEN  !in this case the mixture is neutral
        Te=TNEQ(1)
    ELSE                    !in this case the mixture has electrons
        Te=TNEQ(4)
    END IF
!
!   Compute upper half of the matrix
!
!   Beware, the following lines are not at all trivial.
!   Electron temperature appears to be missing in electron-neutral collisions.
!   However, it IS already hidden in the reduction of collision integrals.
!   So DON'T change T back into Te ! Only witht the following lines does one
!   obtain the correct two-temperature binary diffusion coefficients.

    DO i=1,nsp-1
        DO j=i,nsp-1
            dij(i,j)=KUNIV*T / ( p*DELTA(1,i,j) ) ! Neutral-neutral collisions
        END DO
        dij(i,nsp)=KUNIV*T / ( p*DELTA(1,i,nsp) ) ! Electron-neutral collisions
    END DO
    dij(nsp,nsp)=KUNIV*Te / ( p*DELTA(1,i,nsp) ) ! Electron-electron collisions

!
!   Mirror the matrix
!
    DO i=1,nsp
        DO j=1,i-1
            dij(i,j)=dij(j,i)
        END DO
    END DO
!
END SUBROUTINE set_binary_diffusion

!******************************************************************************

SUBROUTINE diffusion_coefficients_binary(p,TNEQ,dij_peg)
   USE global_thermo, only: nsp
   USE global_traco, only: index_array
   IMPLICIT NONE
   REAL(kind=8) :: p,TNEQ(4),dij_peg(nsp,nsp)

   INTEGER :: i,j
   REAL(kind=8) :: dij_ord(nsp,nsp)

   dij_ord = 0.0

   CALL set_binary_diffusion(p,TNEQ,dij_ord)

   do i=1,nsp
      do j=1,nsp
         dij_peg(index_array(i),index_array(j)) = dij_ord(i,j)
      end do
   end do


END SUBROUTINE diffusion_coefficients_binary
