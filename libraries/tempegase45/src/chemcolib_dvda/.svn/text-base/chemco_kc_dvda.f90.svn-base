! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //  PEGASE 4. CHEMO LIBRARY                                         //
! //                                                                  //
! //  SUBROUTINE EQ_CNT_KC_DVDA                                       //
! //                                                                  //
! //  Inputs:  r           reaction number                       -    //
! //           T           temperature                           K    //
! //           T_NEQ       thermal nonequilibrium temperatures   K    //
! //                                                                  //
! //  Outputs: kc          array of equilibrium constants for    dep. //
! //                       concentrations                             //
! //                                                                  //
! //  Flags:   flg_anha    if 1, anharmonicity corrections            //
! //           flg_neq     if 1, thermal non-equilibrium              // 
! //           flg_stop    if 1, stops on subroutine errors           //
! //           flg_termo   mode of computation of equilibrium         // 
! //                       constants                                  //
! //                                                                  //
! //  This subroutine is called to compute the                        //
! //  equilibrium constants for concentration according to the        //
! //  kinetic model of reactions. It is just an excerpt of the        // 
! //  subroutine EQ_CONSTANTS of TERMO LIBRARY.                       //
! //                                                                  //
! //  J. P. Mellado, July 98                                          // 
! //  Modified by D. Vanden Abeele, April 2000 towards                //
! //  thermal non-equilibrium                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////

SUBROUTINE EQ_CNT_KC_DVDA(r, T, kc, flg_anha, flg_neq, flg_stop, flg_termo )

!   Global variables
!   ^^^^^^^^^^^^^^^^

    USE global_thermo
    USE global_chemco
    USE global_chemco_dvda

!   Local variables
!   ^^^^^^^^^^^^^^^

    IMPLICIT NONE
    REAL(kind=8) T, kc
    INTEGER r, flg_anha, flg_neq, flg_stop, flg_termo
    REAL(kind=8) T_NEQ(4)
    INTEGER s, nusum
    REAL(kind=8) muo(nsp), sum, po, pot(8)

!   Initialization
!   ^^^^^^^^^^^^^^

    po = 1.0d0
    muo = 0.
    pot = 0.
    T_NEQ = 0.d0

!   Compute the equilibrium constant
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    T_NEQ = T
    DO s = 1, nsp
       CALL SPECIES_GIBBS( s, po, T, T_NEQ, flg_anha, 1, flg_neq, flg_stop,&
   &                       flg_termo, pot )
       muo(s) = (pot(1)+pot(8))/RUNIV/T
    ENDDO

    sum = 0
    nusum = 0
    DO s = 1, nsp
        nusum = nusum + (CNE_P(r,s)-CNE_R(r,s))
        sum = sum + (CNE_P(r,s)-CNE_R(r,s))*muo(s)
    ENDDO
    kc = -sum - nusum*log(RUNIV*T)
    kc = exp( kc )
      
    RETURN
      
END SUBROUTINE EQ_CNT_KC_DVDA
