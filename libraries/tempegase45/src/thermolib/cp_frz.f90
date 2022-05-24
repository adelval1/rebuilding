! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE cp_frozen                                            //
! //                                                                  //
! //  Input:  t           equilibrium temperature (K)                 //
! //          tneq        nonequilibrium temperatures (K)             //
! //          x           array of mole or mass fractions (cf mode)   //
! //                                                                  //
! //  Flags:  flg_anha    anharmonicity corrections flag              //
! //          mode        defines a p,T or rho,T input                //
! //          flg_neq     if 1, considers multiple temperaturess      //
! //          flg_stop    stops if errors in subroutines              //
! //          flg_termo   statistical or reftable computation         //
! //                                                                  //
! //  Output: cp(8)        array of  specific heats                   //
! //              1           complete (except formation energy)      //
! //              2           translation                             //
! //              3           rotation (anha=0)                       //
! //              4           vibration (anha=0)                      //
! //              5           electronic (anha=0)                     //
! //              6           coupled rotation-vibration (anha=1)     //
! //              7           coupled internal (anha=2)               //
! //              8           formation energy/energy (not used)      //
! //                                                                  //
! //  This subroutine provides the specific heat at constant pressure //
! //  for the mixture of nsp according to the flags. The value is     //
! //  computed from the cp of individual species computed from        //
! //  statistical mechanics. It works for cp per mole or cp per mass  //
! //  provided that x is the mole fraction array or mass fraction     //
! //  array, depending on the case as defined by mode.                //
! //                                                                  //
! //  Benoit Bottin, 03/12/96 - modified Windows 95 5/8/97.           //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE cp_frozen (T,TNEQ,x,&
     &flg_anha,mode,flg_neq,flg_stop,flg_termo,cp)
!
      USE global_thermo
      IMPLICIT NONE
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      INTEGER flg_anha,mode,flg_neq,flg_stop,flg_termo,istop,i,j
      REAL(kind=8) T,TNEQ(4),cp(8),x(:),sent(8)
!
!     Checks of correct variable values and initialize
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      istop=0
      pgname='pegaslib'
      IF (flg_neq.EQ.0) THEN
          IF (T.LE.0) THEN
              WRITE (localinfo,101) T
101           FORMAT (11x,'Routine: cp_frozen. Value:',d15.8)
              CALL PRINT_ERROR (24,pgname,localinfo,0)
              istop=1
          ENDIF
      ELSE
          DO i=1,4
              IF (TNEQ(i).LE.0) THEN
                  WRITE (localinfo,102) i,TNEQ(i)
102   FORMAT (11x,'Routine: cp_frozen. Index:',i3,' Value:',d15.8)
                  CALL PRINT_ERROR (27,pgname,localinfo,0)
                  istop=1
              ENDIF
          END DO
      ENDIF
      DO i=1,8
          cp(i)=0.0d0
      END DO
!
!     Build cv array by recursive calls to species cv
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      DO i=1,nsp
          DO j=1,8
              sent(j)=0.0d0
          END DO
          CALL SPECIES_CP (i,T,TNEQ,&
     &flg_anha,mode,flg_neq,flg_stop,flg_termo,sent)
          DO j=1,8
              cp(j)=cp(j)+sent(j)*x(i)
          END DO
      END DO
!
!     Error handling
!     ^^^^^^^^^^^^^^
      IF ((istop.EQ.1).AND.(flg_stop.EQ.1)) THEN
          STOP 'Program ended due to PEGASE thermodynamic library error'
      ENDIF
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
