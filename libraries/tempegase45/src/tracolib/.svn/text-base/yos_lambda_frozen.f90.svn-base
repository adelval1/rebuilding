SUBROUTINE Yos_lambda_frozen(k_trans, k_rot, k_vib, k_electronic)

!*****************************************************************************
! Translational thermal conductivity and internal (Eucken-type) thermal 
! conductivity according to Yos' mixture rule. Off-diagonal elements taken
! into account in an approximate manner. Some intermediate quantities
! have to be computed in advance (delta, a_average_lambda, lambda_factor, ...)
! Approximation to result of Hirschfelder, Curtiss & Bird (first non-
! vanishing Sonine polynomial contribution). See e.g. Gupta et al.
!*****************************************************************************
      USE global_thermo
      USE global_traco
      USE traco_ord_var
      IMPLICIT NONE
      REAL(kind=8) :: k_trans, k_rot, k_vib, k_electronic


      REAL(kind=8) :: numerator (1:nsp), numerator1
      INTEGER :: s1, s2


!
!*****************************************************************************
!
!  1. Translational
!
      numerator1 = 0.
      DO s1 = 1, nsp - nelectrons
         numerator1 = numerator1 + x_ord (s1) / (lambda_factor(s1))
      END DO
      k_trans = numerator1 / (1.0-a_average_lambda*(numerator1))
!
!  For use in internal thermal conductivity:
!
      numerator = 0.
      DO s1 = 1, nsp - nelectrons 
         DO s2 = 1, nsp 
            numerator (s1) = numerator (s1) + x_ord (s2) * &
           & delta (1, s1, s2)
         END DO
      END DO
!
!  2. Rotational & vibrational
!
      DO s1 = 1, nmolecn
         k_rot = k_rot + cpr_ord (s1) / RUNIV * x_ord (s1) / numerator (s1)
         k_vib = k_vib + cpv_ord (s1) / RUNIV * x_ord (s1) / numerator (s1)
      END DO

      IF (nions.ne.0) THEN 
         DO s1 = nneut + 1, nneut + nmoleci
            k_rot = k_rot + cpr_ord (s1) / RUNIV * x_ord (s1) / numerator (s1)
            k_vib = k_vib + cpv_ord (s1) / RUNIV * x_ord (s1) / numerator (s1)
         END DO
      ENDIF

      k_rot = k_rot * KUNIV
      k_vib = k_vib * KUNIV

!
!  3. Electronic (exitation of electrons)
!
      DO s1 = 1, nsp - nelectrons
         k_electronic = k_electronic + cpe_ord (s1) / RUNIV * x_ord (s1) &
        & / numerator (s1)
      END DO
      k_electronic = k_electronic * KUNIV
!
!*****************************************************************************
END SUBROUTINE Yos_lambda_frozen
