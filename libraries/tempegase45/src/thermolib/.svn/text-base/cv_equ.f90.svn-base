! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE cv_equilibrium                                       //
! //                                                                  //
! //  Input:  mm          molar mass of the mixture (kg/mol)          //    
! //          dedt(8)     partial derivative de/dt at constant rho    //
! //          dmdt        partial derivative dm/dt at constant rho    //
! //          ener(8)     array of energies ALWAYS PER MASS           //
! //                                                                  //
! //  Flags:  mode        defines a per mole or per mass input (1/2)  //
! //                                                                  //
! //  Output: cv(8)       array of cv per mole or per mass (mode=1/2) //
! //              1           complete                                //
! //              2           translation                             //
! //              3           rotation (anha=0)                       //
! //              4           vibration (anha=0)                      //
! //              5           electronic (anha=0)                     //
! //              6           coupled rotation-vibration (anha=1)     //
! //              7           coupled internal (anha=2)               //
! //              8           formation energy                        //
! //                                                                  //
! //  This subroutine provides the specific heat at constant volume   //
! //  of the mixture in equilibrium. The output can be given per mole //
! //  or per mass depending on mode (ener is always passed per mass)  //
! //                                                                  //
! //  Benoit Bottin, 23/09/97                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE cv_equilibrium (mode,ener,mm,dedt,dmdt,cv)
!
      IMPLICIT NONE
      INTEGER i,mode
      REAL(kind=8) cv(8),ener(8),mm,dedt(8),dmdt
!
!     Case 1: mole fraction case
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (mode.EQ.1) THEN
          DO i=1,8
              cv(i)=mm*dedt(i)+ener(i)*dmdt
          END DO
!
!     Case 2: mass fraction case
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^
      ELSE
          DO i=1,8
              cv(i)=dedt(i)
          END DO
      END IF
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
