REAL(kind=8) FUNCTION viscosity()

!***********************************************************************
! Computes viscosity according to Yos' mixture rule.
! Off-diagonal elements taken into account in an approximate manner.
! Intermediate quantities such as the average off-diagonal element
! 'a_average' should be computed in advance.
! Approximation to result of Hirschfelder, Curtiss & Bird (first non-
! vanishing Sonine polynomial contribution). See e.g. Gupta et al.
!***********************************************************************
      USE global_thermo
      USE global_traco
      USE traco_ord_var, only: x_ord
      IMPLICIT NONE
!***********************************************************************


      REAL(kind=8) :: numerator
      INTEGER :: s1

!***********************************************************************

      numerator = 0.
      DO s1 = 1, nsp
         numerator = numerator + x_ord (s1) / (visc_factor(s1))
      END DO
      viscosity = numerator / (1-a_average_visc*(numerator))

!***********************************************************************
END FUNCTION viscosity
