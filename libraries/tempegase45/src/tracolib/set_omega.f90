SUBROUTINE Set_Omega (n, ne, T, T_e)
!
!***********************************************************
! Set lower-order collision integrals omega-{11} and
! omega-{22}. Set also ratios BStar and CStar.
! Use fourth-order curve-fits (e.g. Gupta et al).
!***********************************************************
      USE global_thermo
      USE global_traco
      IMPLICIT NONE
      REAL(kind=8) :: n,ne,T,T_e

      REAL(kind=8) :: b0_ions, b0_electrons, d0, d1,ne_lim
      REAL(kind=8) :: e, lnT, lnTe, lnT_temp, kT, kTe, Tstar, Tstare, lnTstar, lnTstare
      REAL(kind=8) :: omega11_ion_ion, omega11_el_el, omega11_ion_el
      REAL(kind=8) :: omega22_ion_ion, omega22_el_el, omega22_ion_el
      REAL(kind=8) :: BStar_ion_ion, BStar_el_el, BStar_ion_el
      REAL(kind=8) :: CStar_ion_ion, CStar_el_el, CStar_ion_el
      INTEGER :: s1, s2
      INTEGER :: Omega_nr
      INTEGER :: temp1, temp2
!
!***********************************************************

!  Initialize

      Omega = 0.
      Bstar = 0.
      Cstar = 0.

!  Define some often used expressions

      lnT = Log (T)
      lnTe = Log (T_e)
      kT = KUNIV * T
      kTe = KUNIV * T_e
      e = EUNIV

      IF (nions.ne.0) THEN

!     d0 goes to infinity for zero electron density. We add a small quantity
!     to keep it finite

         ne_lim = ne+1.d-32*n
         d0 = Sqrt (EPSILONUNIV*kTe/(2.*ne_lim*e**2))
                                             ! Debye-shielding distance
                                             ! electrons and ions!

         d1 = ne_lim ** (-.333333333333)! Inter-electron distance

!     Debye shielding is effective only when interelectron distance is smaller
!     than the Debye length.

!         d0 = Max (d0, d1)   ! Like this d0 always Debye length. Fits better with
                              ! experiments

         b0_electrons = e*e/(8.*PIUNIV*EPSILONUNIV*kTe) ! Average closest impact 
                                                        ! parameter for collisions
                                                        ! involving an electron.

         b0_ions = e*e/(8.*PIUNIV*EPSILONUNIV*kT) ! Average closest impact parameter
                                                  ! for coll. involving only ions.

         d0 = MIN (d0, 10000.*(b0_electrons+b0_ions))

         Tstar = d0 / (2.*b0_ions)
         Tstare = d0 / (2.*b0_electrons)
         Tstar = MAX(Tstar, 0.1d0)
         Tstare= MAX(Tstare, 0.1d0)
         lnTstar = Log(Tstar)
         lnTstare = Log(Tstare)
!
!        Evaluate the collision integrals for charged-charged
!        interactions. Curve-fits based upon numerical results
!        of Mason, Munn, Smith / Devoto.
!
         omega11_el_el = PIUNIV * d0 ** 2 * dexp (-1.4041+&
              & .8004*lnTstare-.09307*lnTstare**2+&
              & .005406*lnTstare**3-.0001161*lnTstare**4) / &
              & Tstare ** 2

         omega11_ion_ion = PIUNIV * d0 ** 2 * dexp (-1.3513+&
              & 7.5847d-1*lnTstar-8.01279d-2*lnTstar**2+3.41942d-3*lnTstar**3) &
              & / Tstar ** 2

         omega11_ion_el = PIUNIV * d0 ** 2 * dexp (-.75917+&
              & .50283*lnTstare-4.33739d-2*lnTstare**2+1.67793d-3*lnTstare**3) &
              & / Tstare ** 2

         omega22_el_el = PIUNIV * d0 ** 2 * dexp (-1.1266+&
              & .7714*lnTstare-.09587*lnTstare**2+&
              & .005881*lnTstare**3-.0001312*lnTstare**4) / &
              & Tstare ** 2

         omega22_ion_ion = PIUNIV * d0 ** 2 * dexp (-1.1266+&
              & .7714*lnTstar-.09587*lnTstar**2+&
              & .005881*lnTstar**3-.0001312*lnTstar**4) / &
              & Tstar ** 2

         omega22_ion_el = PIUNIV * d0 ** 2 * dexp (-.8848+&
              & .71509*lnTstare-.0959*lnTstare**2+&
              & .006364*lnTstare**3-.00015072*lnTstare**4) / &
              & Tstare ** 2

         BStar_el_el = max( 1.4292 - .0869 * lnTstare + .007332 * &
              & lnTstare ** 2 - .0002626 * lnTstare ** 3 + 3.1204d-6 * &
              & lnTstare ** 4, 1.d0)

         BStar_ion_ion = max( 1.4292 - .0869 * lnTstar + .007332 * &
              & lnTstar ** 2 - .0002626 * lnTstar ** 3 + 3.1204d-6 * lnTstar ** 4, &
              & 1.d0 )

         BStar_ion_el = max( 1.3281 - .0691 * lnTstare + .007629 * &
              & lnTstare ** 2 - .0004144 * lnTstare ** 3 + 8.5491d-6 * &
              & lnTstare ** 4, 1.d0)

         CStar_el_el = max( .6008 - .05934 * lnTstare + .004867 * &
              & lnTstare ** 2 - .0001389 * lnTstare ** 3 + 4.8891d-7 * &
              & lnTstare ** 4, .333333d0)

         CStar_ion_ion = max( .6008 - .05934 * lnTstar + .004867 * &
              & lnTstar ** 2 - .0001389 * lnTstar ** 3 + 4.8891d-7 * lnTstar ** 4, &
              & .333333d0)

         CStar_ion_el = max( .6008 - .05934 * lnTstare + .004867 * &
              & lnTstare ** 2 - .0001389 * lnTstare ** 3 + 4.8891d-7 * &
              & lnTstare ** 4, .333333d0)

      END IF
!
!   Start with Omega-collision integrals.
!
      DO Omega_nr = 1, 2

!     Calculates the Omega11 and Omega22 for the interactions of a neutral species
!     with another species.
!     Curve-fits are used
         DO s1 = 1, nneut
            DO s2 = s1, nsp-1

               temp1 = neutral_index_array(s1)
               temp2 = index_array(s2)

               Omega(Omega_nr, s1, s2) = dexp (((fits(Omega_nr, temp1, &
              &    temp2, 1)*lnT+fits(Omega_nr, temp1, temp2, 2))*lnT &
              &   +fits(Omega_nr, temp1, temp2, 3))*lnT &
              &   +fits(Omega_nr, temp1, temp2, 4)) * 1.0d-20

            END DO

            s2=nsp
            lnT_temp=lnT
            IF(CHARGE(nsp).EQ.-1) lnT_temp=lnTe
                           
            temp1 = neutral_index_array(s1)
            temp2 = index_array(s2)

            Omega(Omega_nr, s1, s2) = dexp (((fits(Omega_nr, temp1, &
          &    temp2, 1)*lnT_temp+fits(Omega_nr, temp1, temp2, 2))*lnT_temp &
          &    +fits(Omega_nr, temp1, temp2, 3))*lnT_temp &
          &    +fits(Omega_nr, temp1, temp2, 4)) * 1.0d-20

         END DO

      END DO

      IF(nions .NE. 0) THEN

!     Omega integrals for ion-ion interactions

         DO s1 = nneut + 1, nneut + nions
            DO s2 = s1, nsp - 1

               Omega(1, s1, s2) = omega11_ion_ion
               Omega(2, s1, s2) = omega22_ion_ion

            END DO
         END DO

!     Omega integrals for ion-electron interactions

         DO s1 = nneut + 1, nneut + nions
            Omega(1, s1, nsp) = omega11_ion_el
            Omega(2, s1, nsp) = omega22_ion_el
         END DO
!
         Omega(1, nsp, nsp) = omega11_el_el
         Omega(2, nsp, nsp) = omega22_el_el
!
      END IF
!
!   Mirror the upper-triangular part of the matrix omega onto
!   the lower triangular-part.
!
      DO Omega_nr = 1, 2

         DO s1 = 2, nsp
            DO s2 = 1, s1 - 1
               Omega(Omega_nr, s1, s2) = Omega (Omega_nr, s2, s1)
            END DO
         END DO
!
      END DO

!
!   Proceed with Bstar and Cstar integrals.
!

!   Neutral with other species
      Omega_nr = 3
      DO s1 = 1, nneut
         DO s2 = s1, nsp-1

            temp1 = neutral_index_array(s1)
            temp2 = index_array(s2)

            Bstar(s1, s2) = dexp (((fits(Omega_nr, temp1, &
              & temp2, 1)*lnT+fits(Omega_nr, temp1, temp2, 2))*lnT &
              & +fits(Omega_nr, temp1, temp2, 3))*lnT &
              & +fits(Omega_nr, temp1, temp2, 4))*1.0d-20
!		write(*,*)'BSTAR',Bstar(s1, s2)!pietro
           ! For the time being, use hard-sphere extrapolation for
           ! Cstar where neutral particles are involved. In principle,
           ! all needed higher order collision integrals are available
           ! (Capitelli et al, AIAA98-2936 rep) and the actual values
           ! of Cstar could be easily added.
            Cstar(s1, s2) = 1.d0

         END DO

         s2=nsp
         lnT_temp=lnT
         IF(CHARGE(nsp).EQ.-1) lnT_temp=lnTe
                           
         temp1 = neutral_index_array(s1)
         temp2 = index_array(s2)
         Omega_nr = 3

         Bstar(s1, s2) = dexp (((fits(Omega_nr, temp1, &
          & temp2, 1)*lnT_temp+fits(Omega_nr, temp1, temp2, 2))*lnT_temp &
          & +fits(Omega_nr, temp1, temp2, 3))*lnT_temp &
          & +fits(Omega_nr, temp1, temp2, 4))* 1.0d-20

        ! For the time being, use hard-sphere extrapolation for
        ! Cstar where neutral particles are involved. In principle,
        ! all needed higher order collision integrals are available
        ! (Capitelli et al, AIAA98-2936 rep) and the actual values
        ! of Cstar could be easily added.
         Cstar(s1, s2) = 1.d0

      END DO

      IF(nions .NE. 0) THEN
!
!     Bstar and Cstar integrals. Ion-ion interactions
!
         DO s1 = nneut + 1, nneut + nions
            DO s2 = s1, nsp - 1

               Bstar(s1, s2) = Bstar_ion_ion
               Cstar(s1, s2) = Cstar_ion_ion

            END DO
         END DO

!     Ion-electron interactions
         DO s1 = nneut + 1, nneut + nions
               Bstar(s1, nsp) = Bstar_ion_el
               Cstar(s1, nsp) = Cstar_ion_el
         END DO

!     Electron-electron interactions
         Bstar(nsp, nsp) = Bstar_el_el
         Cstar(nsp, nsp) = Cstar_el_el

      END IF
!
!   Mirror the upper-triangular part of the matrix omega onto
!   the lower triangular-part.
!
      DO s1 = 2, nsp
         DO s2 = 1, s1 - 1
            Bstar(s1, s2) = Bstar (s2, s1)
            Cstar(s1, s2) = Cstar (s2, s1)
         END DO
      END DO

END SUBROUTINE Set_Omega
