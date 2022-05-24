REAL(kind=8) FUNCTION electrical_conductivity (ne, T_e)
!
!*****************************************************
! Compute the electrical conductivity through Devoto's
! mixture rule. Two non-vanishing Sonine polynomial
! contributions are used. Use is made of the fact
! that the elctron mass is negligible w.r.t. to the
! heavy species masses; this enables us to considerably
! simplify the formula.
! Intermediate q-{mp} matrix needs to be computed in 
! advance.
!*****************************************************
      USE global_thermo
      USE global_traco
      IMPLICIT NONE
!*****************************************************
!
      REAL(kind=8) :: ne
      REAL(kind=8) :: me
      REAL(kind=8) :: T_e
!
      REAL(kind=8) :: sigma
!
      REAL(kind=8) :: m1
      REAL(kind=8) :: eps
!******************************************************

      me = 0.549d-6
      eps = .01
      m1 = me / NAUNIV

      SELECT CASE (norder)

! One Sonine  polynome -  first order.

      CASE (2)
      sigma = 3*(EUNIV*ne)**2/(2*KUNIV*T_e)*sqrt(2.*PIUNIV*KUNIV*T_e/m1)/qel(0,0)

! Two Sonine  polynomes - second order.

      CASE (3)
      sigma = 3 * (EUNIV*ne) ** 2 / (2*KUNIV*T_e) * Sqrt (2.*PIUNIV*KUNIV*T_e/m1) &
     & / (qel(0, 0)-qel(0, 1)**2/qel(1, 1))

! Three Sonine  polynomes - third order.

      CASE (4)
      sigma = 3*(EUNIV*ne)**2/(2*KUNIV*T_e)*                   &
             & sqrt(2.*PIUNIV*KUNIV*T_e/m1)*              &
             & (qel(1,1)*qel(2,2)-qel(1,2)**2)/(-qel(0,1)**2*qel(2,2) +  &
             & (2*qel(0,1)*qel(0,2)-qel(0,0)*qel(1,2))*qel(1,2) +        &
             & (qel(0,0)*qel(2,2)-qel(0,2)**2)*qel(1,1))

      END SELECT
     
      electrical_conductivity= sigma

!*****************************************************
END FUNCTION electrical_conductivity
