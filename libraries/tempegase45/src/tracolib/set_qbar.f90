SUBROUTINE Set_Qbar (nsp, T_e, ne, Qbar)

!***********************************************************************
! Compute a wide range of electron-particle collision integrals, to be
! used in the electrical conductivity and electron thermal conductivity
! calculations (See Devoto). Higher order charged-charged
! collision integrals taken from the work of Mason & Munn / Devoto.
!***********************************************************************
      USE global_traco
      IMPLICIT NONE
!***********************************************************************
      INTEGER nsp
      REAL(kind=8) :: T_e, TStare, lnTStare
      REAL(kind=8) :: ne
      REAL(kind=8) :: Qbar (1:2, 1:7, 1:nsp)
      REAL(kind=8) :: d0, b0, d1
      REAL(kind=8) :: t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12
      INTEGER :: I
      REAL(kind=8),PARAMETER :: PI=3.141592653589e0
      REAL(kind=8),PARAMETER :: KB=1.380658e-23
      REAL(kind=8),PARAMETER :: E=1.60217733e-19
      REAL(kind=8),PARAMETER :: EPSILON0=8.854187817e-12

!***********************************************************************

      DO I = 1, nneut

         Qbar (1, 1, I) = Omega (1, I, nsp)
         Qbar (1, 2, I) = CStar (nsp, I) * Qbar (1, 1, I)
         Qbar (1, 3, I) = 5. / 4. * Qbar (1, 2, I) - &
        & (BStar(nsp, I)/4.) * Qbar (1, 1, I)

!    Extrapolate in a hard-sphere manner for the remaining part.
!    (See Hirschfelder)

         Qbar (1, 4, I) = Qbar (1, 3, I)
         Qbar (1, 5, I) = Qbar (1, 4, I)
         Qbar (1, 6, I) = Qbar (1, 5, I)
         Qbar (1, 7, I) = Qbar (1, 6, I)

         Qbar (2, 2, I) = Omega (2, I, nsp)

         Qbar (2, 3, I) = Qbar (2, 2, I)
         Qbar (2, 4, I) = Qbar (2, 3, I)
         Qbar (2, 5, I) = Qbar (2, 4, I)
         Qbar (2, 6, I) = Qbar (2, 5, I)

      END DO

      d0 = Sqrt (epsilon0*kb*T_e/(2.*ne*e**2))! Debye-shielding distance
                                        ! ions and electrons

      d1 = ne ** (-1./3.)! Inter-electron distance

      d0 = Max (d0, d1)

      b0 = e ** 2 / (8.*pi*epsilon0*kb*T_e)! Average closest impact
                                        ! parameter (electons)

      d0 = min(d0, 20000.*b0)

      TStare = d0 / (2.*b0)
      Tstare= MAX (Tstare, 0.1d0)
      lnTStare = Log (TStare)

      I = nneut + 1

      t1 = Omega (1, I, nsp)
      t2 = CStar (nsp, I) * t1
      t3 = 5. / 4. * t2 - (BStar(nsp, I)/4.) * t1
      t4 = pi * d0 ** 2 * dexp (-2.3059+.3918*lnTStare-&
     & .0506*lnTStare**2+.005026*lnTStare**3-&
     & .0002161*lnTStare**4) / TStare ** 2
      t5 = pi * d0 ** 2 * dexp (-2.6202+.3688*lnTStare-&
     & .04610*lnTStare**2+.004482*lnTStare**3-&
     & .0001906*lnTStare**4) / TStare ** 2


!     Use estimate based on crude analytical expression

      t6=5./7.*t5
      t7=6./8.*t6

      t8 = Omega (2, I, nsp)
      t9 = max(.6747d0-.0632*lnTStare+.01228*lnTStare**2-&
     & .001151*lnTStare**3+4.001d-5*lnTStare**4, .5d0) * t8
      t10 = pi * d0 ** 2 * dexp (-1.5267+.5592*lnTStare-&
     & .1028*lnTStare**2+.01243*lnTStare**3-&
     & .0005828*lnTStare**4) / TStare ** 2


!     Use estimate based on crude analytical expression

      t11 = 4./6.*t10
      t12 = 5./7.*t11


      DO I = nneut + 1, nsp - 1

      Qbar (1, 1, I) = t1
      Qbar (1, 2, I) = t2
      Qbar (1, 3, I) = t3
      Qbar (1, 4, I) = t4
      Qbar (1, 5, I) = t5

      Qbar (1, 6, I) = t6
      Qbar (1, 7, I) = t7

      Qbar (2, 2, I) = t8
      Qbar (2, 3, I) = t9
      Qbar (2, 4, I) = t10

      Qbar (2, 5, I) = t11
      Qbar (2, 6, I) = t12

      END DO


      t1 = Omega (1, nsp, nsp)
      t2 = CStar (nsp, nsp) * t1
      t3 = 5. / 4. * t2 - (BStar(nsp, nsp)/4.) * t1
      t4 = pi * d0 ** 2 * dexp (-2.3059+.3918*lnTStare-&
     & .0506*lnTStare**2+.005026*lnTStare**3-&
     & .0002161*lnTStare**4) / TStare ** 2
      t5 = pi * d0 ** 2 * dexp (-2.8765+.5380*lnTStare-&
     & .08601*lnTStare**2+.008475*lnTStare**3-&
     & .0003343*lnTStare**4) / TStare ** 2

!     Use estimate based on crude analytical expression

      t6=5./7.*t5
      t7=6./8.*t6

      t8 = Omega (2, nsp, nsp)
      t9 = max(.6952d0-.05536*lnTStare+.006934*lnTStare**2-&
     & .0003886*lnTStare**3+7.949d-6*lnTStare**4, .5d0) * t8
      t10 = pi * d0 ** 2 * dexp (-1.7673+.6468*lnTStare-&
     & .09789*lnTStare**2+.008582*lnTStare**3-&
     & .0002997*lnTStare**4) / TStare ** 2


!     Use estimate based on crude analytical expression

      t11 = 4./6.*t10
      t12 = 5./7.*t11


      Qbar (1, 1, nsp) = t1
      Qbar (1, 2, nsp) = t2
      Qbar (1, 3, nsp) = t3
      Qbar (1, 4, nsp) = t4
      Qbar (1, 5, nsp) = t5

      Qbar (1, 6, nsp) = t6
      Qbar (1, 7, nsp) = t7

      Qbar (2, 2, nsp) = t8
      Qbar (2, 3, nsp) = t9
      Qbar (2, 4, nsp) = t10

      Qbar (2, 5, nsp) = t11
      Qbar (2, 6, nsp) = t12

!***********************************************************************
END SUBROUTINE Set_Qbar
