SUBROUTINE Order_Data (nsp,index_array, unordered_array, ordered_array)
!
!***************************************************************
! Reorder a generic array(1:nsp) in the traco-format (see
! routine init_particles.f90).
!***************************************************************
      IMPLICIT NONE
!***************************************************************
!
      INTEGER :: I, nsp
      REAL(kind=8) :: unordered_array (1:nsp)
      INTEGER :: index_array (1:nsp)
!
      REAL(kind=8) :: ordered_array (1:nsp)
!
!
! Order the data.
!
      DO I = 1, nsp
         ordered_array (I) = unordered_array (index_array(I))
      END DO
!
!***************************************************************
END SUBROUTINE Order_Data
!
SUBROUTINE Reorder_Data (nsp,index_array, ordered_array, unordered_array)
!
!***************************************************************
! Reorder a generic array(1:nsp) in the traco-format (see
! routine init_particles.f90).
!***************************************************************
      IMPLICIT NONE
!***************************************************************
!
      INTEGER :: I, nsp
      REAL(kind=8) :: unordered_array (1:nsp)
      INTEGER :: index_array (1:nsp)
!
      REAL(kind=8) :: ordered_array (1:nsp)
!
!
! Order the data.
!
      DO I = 1, nsp
         unordered_array (index_array(I)) = ordered_array (I)
      END DO
!
!***************************************************************
END SUBROUTINE Reorder_Data
!
