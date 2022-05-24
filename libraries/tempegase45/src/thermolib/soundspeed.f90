! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE sound_speed_equilibrium                              //
! //                                                                  //
! //  Input:  enth(8)     enthalpy per unit mass                      //
! //          dpdt        partial derivative dp over dt at rho        //
! //          dpdr        partial derivative dp over drho at T        //
! //          detdt       partial derivative det over dt at rho       //
! //          detdr       partial derivative det over drho at T       //
! //                                                                  //
! //  Flags:  flg_stop    stops if errors in subroutines              //
! //                                                                  //
! //  Output: a           equilibrium speed of sound                  //
! //          kappa       equilibrium value of parameter kappa        //
! //          khi         equilibrium value of parameter khi          //
! //                                                                  //
! //  This subroutine provides the equilibrium speed of sound, kappa  //
! //  and khi of the mixture.                                         //
! //                                                                  //
! //  Benoit Bottin, 04/12/96. Updated Win95 Fortran 23/09/97.        //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE sound_speed_equilibrium (enth,dpdt,dpdr,detdt,detdr,&
     &a,kappa,khi)
!
      IMPLICIT NONE
      REAL(kind=8) enth(8),dpdt,dpdr,detdt(8),detdr(8),a,aa,kappa,khi
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
!
      pgname='pegaslib'
      kappa=dpdt/detdt(1)
      khi=dpdr-kappa*detdr(1)
      aa=khi+kappa*enth(1)
      IF (aa.LE.0.0d0) THEN
          localinfo(1:34)='           sound_speed_equilibrium'
          CALL PRINT_ERROR (40,pgname,localinfo,0)
          a=0.0d0
      ELSE
          a=DSQRT(aa)
      END IF
!
      RETURN
      END
! ////////////////////////////////////////////////////////////////////// 
