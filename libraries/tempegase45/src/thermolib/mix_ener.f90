! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE mixture_energy                                       //
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
! //  Output: ener(8)      array of  enthalpies                       //
! //              1           complete (except formation energy)      //
! //              2           translation                             //
! //              3           rotation (anha=0)                       //
! //              4           vibration (anha=0)                      //
! //              5           electronic (anha=0)                     //
! //              6           coupled rotation-vibration (anha=1)     //
! //              7           coupled internal (anha=2)               //
! //              8           formation energy/energy                 //
! //                                                                  //
! //  This subroutine provides the  energies of the mixture           //
! //  of nsp according to the flags specified. The array x passed     //
! //  must be the mole fractions if mode is 1 and the mass fractions  //
! //  if mode is 2.                                                   //
! //                                                                  //
! //  Benoit Bottin, 11/10/96. Modified and cleaned F90 28/7/97.      //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE mixture_energy (T,TNEQ,x,&
     &flg_anha,mode,flg_neq,flg_stop,flg_termo,ener)
!
      USE global_thermo
      IMPLICIT NONE
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      INTEGER flg_anha,mode,flg_neq,flg_stop,flg_termo,istop,i,j
      REAL(kind=8) T,TNEQ(4),ener(8),x(:),sent(8)
!
!     Checks of correct variable values and initialize
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      istop=0
      pgname='pegaslib'
      IF (flg_neq.EQ.0) THEN
          IF (T.LE.0) THEN
              WRITE (localinfo,101) T
101           FORMAT (11x,'Routine: mixture_energy. Value:',d15.8)
              CALL PRINT_ERROR (24,pgname,localinfo,0)
              istop=1
          ENDIF
      ELSE
          DO i=1,4
              IF (TNEQ(i).LE.0) THEN
                  WRITE (localinfo,102) i,TNEQ(i)
102   FORMAT (11x,'Routine: mixture_energy. Index:',i3,' Value:',d15.8)
                  CALL PRINT_ERROR (27,pgname,localinfo,0)
                  istop=1
              ENDIF
          END DO
      ENDIF
      DO i=1,8
          ener(i)=0.0d0
      END DO
!
!     Build energy array by recursive calls to species energies
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      DO i=1,nsp
          DO j=1,8
              sent(j)=0.0d0
          END DO
          CALL SPECIES_ENERGY (i,T,TNEQ,&
     &flg_anha,mode,flg_neq,flg_stop,flg_termo,sent)
          DO j=1,8
              ener(j)=ener(j)+sent(j)*x(i)
          END DO
      END DO
!
!     Including formation energy in the total for the mixture
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      ener(1)=ener(1)+ener(8)
!
!     Error handling
!     ^^^^^^^^^^^^^^
      IF ((istop.EQ.1).AND.(flg_stop.EQ.1)) THEN
          STOP 'Program ended due to PEGASE thermodynamic library error'
      ENDIF
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
