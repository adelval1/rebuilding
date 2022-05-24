SUBROUTINE Order_Reactions (nr,nsp,index_array, NU, NU_ord)
!
!*********************************************************************
! Reorder reactions according to the traco-format (see routine
! init_particles.f90).
!*********************************************************************
      IMPLICIT NONE
!*********************************************************************
!
      INTEGER :: i, j, nr, nsp
      INTEGER :: index_array (1:nsp)
      INTEGER :: NU (nr, nsp+1), NU_ord (nr, nsp)
!
!
!*********************************************************************
!
!  Order the stoechiometric matrix.
!
      DO i = 1, nr
         DO j = 1, nsp
!
            NU_ord (i, j) = NU (i, index_array(j))
!
         END DO
      END DO
!
END SUBROUTINE Order_Reactions
