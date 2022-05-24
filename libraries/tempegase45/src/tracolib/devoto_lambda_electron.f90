REAL(kind=8) FUNCTION Devoto_lambda_electron (ne, T_e)

!***********************************************************************
! Compute the electron thermal conductivity through Devoto's
! mixture rule. Two non-vanishing Sonine polynomial
! contributions are used. Use is made of the fact
! that the electron mass is negligible w.r.t. to the
! heavy species masses; this enables us to considerably
! simplify the formula.
! Intermediate q[{mp} matrix needs to be computed in 
! advance.
!***********************************************************************
      USE global_thermo
      USE global_traco
      IMPLICIT NONE
!***********************************************************************

      REAL(kind=8) :: ne
      REAL(kind=8) :: me
      REAL(kind=8) :: T_e

      REAL(kind=8) :: C1, C
      REAL(kind=8) :: m1
      REAL(kind=8) :: lambda
      REAL(kind=8) :: eps

!***********************************************************************
      me=0.549d-6
      eps = .01

! Defining some often used expressions

      m1 = me / NAUNIV
      C1 = 2. * PIUNIV * KUNIV * T_e
      C = 75. * KUNIV * ne ** 2 * Sqrt (C1/m1) / 8.

! Calculate the contribution of the electrons to the thermal
! conductivity. This is the translational part. The part caused by
! exitation of electrons is not included, yet.

      SELECT CASE (norder)

! One non-vanishing Sonine-polynomial contribution.

      CASE (2)
          lambda = C / qel(1, 1)

! Two non-vanishing Sonine-polynomial contributions.

      CASE (3)
          lambda = C / (qel(1, 1)-qel(1, 2)**2/qel(2, 2))
    
! Three non-vanishing Sonine-polynomial contributions.

      CASE (4)
          lambda = C * (qel(2,2)*qel(3,3)-qel(2,3)**2)/(qel(1,1)*qel(2,2)*qel(3,3)-&
        & qel(1,1)*qel(2,3)**2-qel(1,2)**2*qel(3,3)+2*qel(1,2)*qel(1,3)*qel(2,3)-&
        & qel(1,3)**2*qel(2,2))

      END SELECT

      Devoto_lambda_electron = lambda

!***********************************************************************
END FUNCTION Devoto_lambda_electron
