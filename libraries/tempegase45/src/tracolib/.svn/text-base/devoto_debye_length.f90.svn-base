REAL(kind=8) FUNCTION Debye_Length_Devoto (ne, n, T_e)

!*********************************************************************
! Compute the Debye-length (for e.g. use by people designing 
! electrostatic probes) .
!*********************************************************************
      USE global_thermo
      IMPLICIT NONE
!*********************************************************************

      REAL(kind=8) :: ne,n,T_e

      REAL(kind=8) :: ne_lim
!*********************************************************************

! With this trick the electron number density is finite
  ne_lim = ne+1.d-32*n

! The next formula gives the Debye shielding if both the ions and the
! electrons contribute to the shielding

     Debye_Length_Devoto = sqrt(EPSILONUNIV*KUNIV*T_e/(2.*ne_lim*EUNIV**2))

! The next formula gives the Debye shielding if only the electrons
! make a contribution.

!      Debye_Length_Devoto = Sqrt (EPSILONUNIV*KUNIV*T_e/(ne_lim*EUNIV**2))

!*********************************************************************
END FUNCTION Debye_Length_Devoto
!
