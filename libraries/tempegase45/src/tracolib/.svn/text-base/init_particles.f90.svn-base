SUBROUTINE Init_Particles
!
!***********************************************************************
! Detects number of ions, electrons, molecules.
! The traco routines use their own format:

!      1. Neutral molecules in alphabetic order.   
!      2. Neutral atoms in alphabetic order.   
!      3. Ionic molecules in alphabetic order.   
!      4. Ionic atoms in alphabetic order.   
!      5. Electron

! One converts to this format using the 'index_array' (to be set in
! this routine). One stocks the order of the neutral particles in
! 'neutral_index_array'.
!***********************************************************************
      USE global_thermo
      USE global_traco
      IMPLICIT NONE
!***********************************************************************
!
      INTEGER :: nmolec1
      INTEGER :: nmolec2
      INTEGER :: ntemp
!
      INTEGER :: Nr_Species, i
      INTEGER :: neutral_array (1:nsp)
      INTEGER :: ion_array (1:nsp)
      INTEGER :: molec_array1 (1:nsp)
      INTEGER :: molec_array2 (1:nsp)

!***********************************************************************
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 1. Initialisation
!
      nions = 0
      nneut = 0
      nelectrons = 0
      nmolec1 = 0
      nmolec2 = 0
      index_array = 0
      neutral_index_array = 0
      neutral_array = 0
      ion_array = 0
      molec_array1 = 0
      molec_array2 = 0
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 2. Check the charge of all the species, count the number of ions
!    and neutrals and determine if there are electrons present.
!    Check the symmetry factor to determine whether any molecules are
!    present.
!
      DO Nr_Species = 1, nsp
!
         IF (Charge(Nr_Species) .GT. 0) THEN
            nions = nions + 1
            ion_array (nions) = Nr_Species
         ELSE IF (Charge(Nr_Species) .EQ.-1) THEN
            index_array (nsp) = Nr_Species
            nelectrons = 1
         ELSE
            nneut = nneut + 1
            neutral_array (nneut) = Nr_Species
         END IF
!
      END DO
!
      DO Nr_Species = 1, nneut
!
         IF (SIG(neutral_array(Nr_Species)) .GT. 0) THEN
            nmolec1 = nmolec1 + 1
            molec_array1 (nmolec1) = neutral_array (Nr_Species)
         END IF

      END DO


!  nmolecn stores number of neutral molecules

      nmolecn = nmolec1

      DO Nr_Species = 1, nneut
!
         IF (SIG(neutral_array(Nr_Species)) .EQ. 0) THEN
            nmolec1 = nmolec1 + 1
            molec_array1 (nmolec1) = neutral_array (Nr_Species)
         END IF
!
      END DO
!
      DO Nr_Species = 1, nions
!
         IF (SIG(ion_array(Nr_Species)) .GT. 0) THEN
            nmolec2 = nmolec2 + 1
            molec_array2 (nmolec2) = ion_array (Nr_Species)
         END IF
!
      END DO

!  nmoleci stores number of ionised molecules

      nmoleci = nmolec2

!
      DO Nr_Species = 1, nions
!
         IF (SIG(ion_array(Nr_Species)) .EQ. 0) THEN
            nmolec2 = nmolec2 + 1
            molec_array2 (nmolec2) = ion_array (Nr_Species)
         END IF
!
      END DO
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 3. Filling the index_array with the numbers of the neutral particles
!    followed by the numbers of the ions and as the last element (index
!    nsp) the number of the electron. These numbers are all the position
!    numbers of the particles as they occur in the input file.
!    Put molecular species before atomic species within the above ordering.
!
      DO Nr_Species = 1, nneut
         index_array (Nr_Species) = molec_array1 (Nr_Species)
      END DO
      DO Nr_Species = 1 + nneut, nions + nneut
         index_array (Nr_Species) = molec_array2 (Nr_Species-nneut)
      END DO
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 4. Keep track of the reordering put through within the set of
!    neutral particles.
!
      DO Nr_Species = 1, nneut
!
         ntemp = molec_array1 (Nr_Species)
         DO i = 1, nneut
            IF (neutral_array(i) .EQ. ntemp) neutral_index_array &
           & (Nr_Species) = i
         END DO
!
      END DO
!
!***********************************************************************
END SUBROUTINE Init_Particles
