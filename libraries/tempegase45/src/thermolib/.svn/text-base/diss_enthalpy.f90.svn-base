! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTION dissoc_enthalpy                                        //
! //                                                                  //
! //  Input:  x           array of mole or mass fractions             //    
! //          zo          number of moles per mole of cold gas        //
! //                                                                  //
! //  Flags:  mode        defines a p,T or rho,T input                //
! //          method      defines the type of calculation method used //
! //              1           full use of all formation enthalpies    //
! //              2           limited use to dissociation only        //
! //              3           curve fit of Lewis for dissociated air  //
! //                                                                  //
! //  Output: hd          dissociation enthalpy in MJ/kg              //
! //                                                                  //
! //  This function provides the dissociation enthalpy of the mixture //
! //  to be used in Fay & Riddell's heat transfer law. The value is   //                       //
! //  returned per mole or per mass depending on the flag mode.       //
! //                                                                  //
! //  Benoit Bottin, 23/09/97                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION dissoc_enthalpy (x,zo,mode,method)
!
      USE global_thermo
      IMPLICIT NONE
      REAL(kind=8) x(:),zo,hd,mm
      INTEGER mode,method,i
!
      hd=0.0d0
      SELECT CASE (method)
      CASE (1,2)
          IF (mode.EQ.1) THEN
              DO i=1,nsp
                 IF ((method.EQ.2).AND.(CHARGE(i).NE.0)) THEN
                      CYCLE
                 END IF
                 hd=hd+x(i)*HFOR(i)
              END DO
          ELSE
              DO i=1,nsp
                 IF ((method.EQ.2).AND.(CHARGE(i).NE.0)) THEN
                      CYCLE
                 END IF
                 hd=hd+x(i)*HFOR(i)/MMOL(i)
              END DO
          END IF
      CASE (3)
          IF (zo.LE.1.21153d0) THEN
              hd=1.830923d8*(zo-1)
          ELSE
              hd=3.87295d7+3.512878d8*(zo-1.21153d0)
          END IF
          hd=hd*9.29d-2
          IF (mode.EQ.1) THEN
              mm=0.0d0
              DO i=1,nsp
                  mm=mm+x(i)*MMOL(i)
              END DO
              hd=hd*mm
          END IF
      END SELECT
!
      dissoc_enthalpy=hd
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
