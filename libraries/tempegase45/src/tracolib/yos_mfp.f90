SUBROUTINE mfp_Yos (n_tot, mfp_array)
!
!******************************************************
! Computes an expression for the mean free path of all
! particles in the plasma.
! Mean free paths are stored in 'mfp_array'.
!******************************************************
      USE global_thermo, only: nsp
      USE global_traco
      USE traco_ord_var, only: x_ord
      IMPLICIT NONE
      REAL(kind=8) :: n_tot
      REAL(kind=8) :: mfp_array (1:nsp), mfp_ord(1:nsp)
      INTEGER I, J
!******************************************************
!
      DO I = 1, nsp
!
         mfp_ord (I) = 0.
!
         DO J = 1, nsp
!
            mfp_ord (I) = mfp_ord (I) + x_ord (J) * Omega (2, I, J) &
           & * Sqrt (1.+mmol_ord(I)/mmol_ord(J))
!
         END DO
!
         mfp_ord (I) = 1. / (n_tot*mfp_ord(I))
!
      END DO
!
      CALL Reorder_Data (nsp,index_array, mfp_ord, mfp_array)
!******************************************************
END SUBROUTINE mfp_yos
