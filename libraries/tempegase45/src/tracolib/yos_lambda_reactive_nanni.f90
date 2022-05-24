REAL(kind=8) FUNCTION Yos_lambda_reactive_nanni (T)

!***********************************************************
! Reactive thermal conductivity (i.e. diffusion of
! reaction enthalpies in an equilibrium mixture) 
! computed by Yos' mixture rule. Off-diagonal elements 
! taken into account in an approximate manner. 
! Intermediate quantities (delta) to be computed in advance.
! Approximation to result of Hirschfelder, Curtiss & Bird (first non-
! vanishing Sonine polynomial contribution). See e.g. Gupta et al.
!***********************************************************
      USE global_thermo
      USE global_traco
      USE traco_ord_var, only: x_ord
      IMPLICIT NONE
!***********************************************************

      INTEGER :: i, j, l
      REAL(kind=8) :: T
      REAL(kind=8) :: lambda, tempi, tempj
      REAL(kind=8) :: eps

!***********************************************************


!*********************************************************************
!
!     INTEGER :: i, j, nr, nsp
!     INTEGER :: NU_ord (nr, nsp)
!     REAL(kind=8)  :: HREACT (1:nr), HFOR(1:nsp)
!
!
!*********************************************************************

! write(*,*)"HFOR=",HFOR  ! ### NANNI ### 

!  Set the reaction enthalpies.

      DO i = 1, nr
         HREACT (i) = 0.
         DO j = 1, nsp
            HREACT (i) = HREACT (i) + NU_ord (i, j) * HFOR (j)
         END DO
      END DO
!
      
      eps = 1d-16

      lambda = 0.

      DO l = 1, nr

         tempi = 0.
         DO i = 1, nsp

            tempj = 0.
            DO j = 1, nsp
               tempj = tempj + (NU_ord(l, i)*x_ord(j)-NU_ord(l, j)&
              & *x_ord(i)) * delta (1, i, j)
            END DO

            tempi = tempi + NU_ord (l, i) / (x_ord(i)+eps) * tempj

         END DO

         lambda = lambda + (HREACT(l)/(RUNIV*T)) ** 2 / tempi

      END DO
      
      !write(*,*)"HREACT=",HREACT  ! ### NANNI ### 
      
      Yos_lambda_reactive_nanni = KUNIV * lambda

!***********************************************************
END FUNCTION Yos_lambda_reactive_nanni
!
