SUBROUTINE Set_Yos (nsp,x_ord)
!
!*******************************************************************
! Set intermediate quantities to be used in Yos' mixture rule
! for the frozen thermal conductivity and the viscosity:
!  - a_average_lambda, a_average_visc (approx. off-diag. elements)
!  - lambda_factor, visc_factor
! Intermediate quantities 'delta' to be computed in advance.
!*******************************************************************
      USE global_traco
      IMPLICIT NONE
!*******************************************************************
!
      INTEGER :: I, J, nsp
      REAL(kind=8) :: x_ord (1:nsp)
!
      REAL(kind=8) :: A (1:nsp)
      REAL(kind=8) :: B (1:nsp, 1:nsp)
      REAL(kind=8) :: a_array (1:nsp, 1:nsp)
      REAL(kind=8) :: C1
      REAL(kind=8),PARAMETER :: NA=6.0221367e23
      REAL(kind=8),PARAMETER :: KB=1.380658e-23  
!
!*******************************************************************
      A = 0.
      B = 0.
      a_array = 0.
!*******************************************************************
!
! Start with determining the factor appearing in the calculation
! of the viscosity
!
      DO I = 1, nsp
         A (I) = 0.
         DO J = 1, nsp
!
            A (I) = A (I) + x_ord (J) * (NA/mmol_ord(I)*delta(2, I, J)) 
!
            a_array (I, J) = (NA/(mmol_ord(I)+mmol_ord(J))) * &
           & (2*delta(1, I, J)-delta(2, I, J))
!
         END DO
      END DO
!
      CALL Set_a_av (nsp, x_ord, a_array, A, a_average_visc, 1)
!
      DO I = 1, nsp
         visc_factor (I) = A (I) + a_average_visc
      END DO
!
! Now determine the factor appearing in the expression for the
! thermal conductivity.
! If there are electrons present, their contribution should not be
! included, since it is calculated with a different formula. So only
! the contributions up to nsp-nelectrons must be calculated here.
! (nelectrons=1 if electrons are present, otherwise 0)
!
      DO I = 1, nsp - nelectrons
         A (I) = 0.
         DO J = 1, nsp - nelectrons
!
            C1 = (mmol_ord(I)+mmol_ord(J)) ** 2.
!
            A (I) = A (I) + (2./(15.*kb)) * (x_ord(J)/C1) * &
           & (8*mmol_ord(I)*mmol_ord(J)*delta(2, I, J)+(mmol_ord(I)-&
           & mmol_ord(J))*(9.*mmol_ord(I)-&
           & 15./2.*mmol_ord(J)+(18./5.)*BStar(I, &
           & J)*mmol_ord(J))*delta(1, I, J))
!
            a_array (I, J) = (2./(15.*kb)) * &
           & (mmol_ord(I)*mmol_ord(J)/C1) * &
           & ((33./2.-(18./5.)*BStar(I, &
           & J))*delta(1, I, J)-4.*delta(2, I, J))
!
         END DO
      END DO
!
      CALL Set_a_av (nsp, x_ord, a_array, A, a_average_lambda, 0)
!
      DO I = 1, nsp - nelectrons
         lambda_factor (I) = A (I) + a_average_lambda
      END DO
!
!*******************************************************************
END SUBROUTINE Set_Yos
