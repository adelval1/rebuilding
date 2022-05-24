SUBROUTINE Set_Hreact (nr,nsp, NU_ord, HINT_ord, HREACT)
!
!*********************************************************************
! Set total reaction enthalpies at nonzero temperatures
! (i.e. add internal energies).
!*********************************************************************
        IMPLICIT NONE
!*********************************************************************
!
      INTEGER :: i, j, nr, nsp
      INTEGER :: NU_ord (nr, nsp)
      REAL(kind=8)  :: HREACT (1:nr), HINT_ord(1:nsp)
!
!
!*********************************************************************

!  Set the reaction enthalpies.

      DO i = 1, nr
         HREACT (i) = 0.
         DO j = 1, nsp
            HREACT (i) = HREACT (i) + NU_ord (i, j) * HINT_ord (j)
         END DO
      END DO
!
END SUBROUTINE Set_Hreact
