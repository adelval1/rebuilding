SUBROUTINE Set_q (nsp, T_e, n_ord,n)

!************************************************************************
! Set intermediate q-{mp} matrix, to be used in Devoto's formulas for the
! electrical conductivity and the electron thermal conductivity. Fill 
! enough q-{mp}'s to allow two non-vanishing sonine polynomial 
! contributions in both formulas.
!************************************************************************
      USE global_traco
      IMPLICIT NONE
!************************************************************************
      INTEGER nsp      
      REAL(kind=8) :: T_e
      REAL(kind=8) :: n_ord (1:nsp),n

      REAL(kind=8) :: Qbar (1:2, 1:7, 1:nsp),ne_lim
      REAL(kind=8) :: dummy
      INTEGER :: p, m, J, K

!************************************************************************

!  C_omega_el is the array which includes the various factors occurring in the
!  expressions for the q-mp (see the work of Devoto).
!  The meaning of the indices is:

!  1st : = 1 -> Q(1,s)
!        = 2 -> Q(2,s)
!  2nd : = m
!  3rd : = p
!  4th : = s

!  Qbar : contains all electron - particle collision integrals (also
!  higher order ones).

!************************************************************************

!  Fill Qbar with the already computed electron - particle collision
!  integrals; add moreover higher order electron - particle collision
!  integrals.

      ne_lim = n_ord(nsp)+1.d-32*n
      Qbar = 0.
      CALL Set_Qbar (nsp, T_e, ne_lim, Qbar)

!  Determine the contribution of the collisions of electrons with
!  other particles. In this the Omega(1,s) integrals are used.

      DO p = 0, norder - 1
         DO m = 0, p
            qel(m, p) = 0.
            DO J = 1, nsp - 1
               DO K = 1, m + p + 1
                  qel(m, p) = qel(m, p) + C_omega_el (1, m, p, K) * ne_lim * &
                 & n_ord (J) * Qbar (1, K, J)
               END DO
            END DO
            qel(m, p) = qel(m, p) / 128.
         END DO
      END DO

!  Determine the contribution of the collisions of electrons with
!  electrons. In this the Omega(2,s) integrals are used.

      DO p = 1, norder - 1
         DO m = 1, p
            dummy = 0.
            DO K = 2, m + p
               dummy = dummy + C_omega_el (2, m, p, K) * Qbar (2, K, nsp) * &
              & ne_lim * ne_lim
            END DO
            qel(m, p) = qel(m, p) + dummy * (Sqrt(2.)/128.)
         END DO
      END DO


END SUBROUTINE Set_q
