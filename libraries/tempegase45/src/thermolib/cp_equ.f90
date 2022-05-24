! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE cp_equilibrium                                       //
! //                                                                  //
! //  Input:  mm          molar mass of the mixture (kg/mol)          //    
! //          dhdt(8)     partial derivative dh/dt at constant rho    //
! //          dhdr(8)     partial derivative dh/drho at constant T    //
! //          dmdt        partial derivative dm/dt at constant rho    //
! //          dmdr        partial derivative dm/drho at constant T    //
! //          detdt(8)    partial derivative det/dt at constant rho   //
! //          detdr(8)    partial derivative det/drho at constant T   //
! //          enth(8)     array of enthalpies ALWAYS PER MASS         //
! //                                                                  //
! //  Flags:  mode        defines a per mole or per mass input (1/2)  //
! //                                                                  //
! //  Output: cp(8)       array of cp per mole or per mass (mode=1/2) //
! //              1           complete                                //
! //              2           translation                             //
! //              3           rotation (anha=0)                       //
! //              4           vibration (anha=0)                      //
! //              5           electronic (anha=0)                     //
! //              6           coupled rotation-vibration (anha=1)     //
! //              7           coupled internal (anha=2)               //
! //              8           formation energy                        //
! //                                                                  //
! //  This subroutine provides the specific heat at constant pressure //
! //  of the mixture in equilibrium. The output can be given per mole //
! //  or per mass depending on mode (enth is always passed per mass)  //
! //                                                                  //
! //  Benoit Bottin, 23/09/97                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE  cp_equilibrium (mode,enth,mm,dhdt,dhdr,dmdt,dmdr,&
     &dpdt,dpdr,cp)
!
      IMPLICIT NONE
      INTEGER i,mode
      REAL(kind=8) cp(8),enth(8),dhdt(8),dhdr(8),dmdt,dmdr,dpdt,dpdr
      REAL(kind=8) mm,var
!
      var=dpdt/dpdr
!
!     Case 1: mole fraction case
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (mode.EQ.1) THEN
          DO i=1,8
             cp(i)=mm*dhdt(i)+enth(i)*dmdt-var*(mm*dhdr(i)+enth(i)*dmdr)
          END DO
!
!     Case 2: mass fraction case
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^
      ELSE
          DO i=1,8
              cp(i)=dhdt(i)-dhdr(i)*var
          END DO
      END IF
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
