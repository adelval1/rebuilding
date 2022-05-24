SUBROUTINE Set_a_av (nsp, x_ord, a_Yos, A, a_average, flg_viscosity)

!***********************************************************************
! Compute average off-diagonal matrix element to be used in Yos'
! mixture rule. See e.g. Gupta et al.
!***********************************************************************
      USE global_traco
      IMPLICIT NONE
!***********************************************************************
!
      INTEGER nsp
      REAL(kind=8) :: x_ord (1:nsp)
      REAL(kind=8) :: a_Yos (1:nsp, 1:nsp)
      REAL(kind=8) :: A (1:nsp)
      REAL(kind=8) :: a_average
      INTEGER :: flg_viscosity
!
      REAL(kind=8) :: expr (1:nsp, 1:nsp)
      REAL(kind=8) :: numerator
      REAL(kind=8) :: denominator
      INTEGER :: The_End
      INTEGER :: I, J
!
!***********************************************************************
!
! 1. For the calculation of the thermal conductivity, the contribution
!    of the electrons is calculated with a different formula. If there
!    electrons present we should only include the interactions which do
!    not include electrons in the average.
!    For thermal conductivity and if nelectrons=1, The_End must be nsp-1
!
      IF (flg_viscosity .EQ. 1) THEN
         The_End = nsp
      ELSE
         The_End = nsp - nelectrons
      END IF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 2. Defining some often used expressions, and perform initialisation
!
      expr = 0.
      DO I = 1, The_End
         IF (I .LT. The_End) THEN
            DO J = I + 1, The_End
               expr (I, J) = 2. * x_ord (I) * x_ord (J) * &
              & (1./A(I)-1./A(J)) ** 2
            END DO
         END IF
      END DO
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      numerator = 0.
      denominator = 0.
      DO I = 1, The_End
         IF (I .LT. The_End) THEN
            DO J = I + 1, The_End
               numerator = numerator + expr (I, J) * a_Yos (I, J)
               denominator = denominator + expr (I, J)
            END DO
         END IF
      END DO
      a_average = numerator / denominator
!
!***********************************************************************
END SUBROUTINE Set_a_av
