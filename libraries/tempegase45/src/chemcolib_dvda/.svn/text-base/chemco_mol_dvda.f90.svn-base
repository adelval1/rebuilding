! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //  PEGASE 4. CHEMO LIBRARY                                         //
! //                                                                  //
! //  SUBROUTINE PROD_TERM_MOL_DVDA                                   //
! //                                                                  //
! //  Inputs:  T_H, T_E    temperatures                      K        //
! //           p           pressure                          Pa       //
! //           X           mole fractions                    -        //     
! //                                                                  //
! //  Outputs: W           array with the production terms   kg/m^3/s //
! //                                                                  //
! //  Flags:   flg_anha    if 1, anharmonicity corrections            //
! //           flg_neq     if 1, thermal non-equilibrium              // 
! //           flg_stop    if 1, stops on subroutine errors           //
! //           flg_termo   mode of computation of equilibrium         // 
! //                       constants                                  //
! //                                                                  //
! //  This subroutine computes the production terms for each species  //
! //  once a kinetic model has been established by CHEMODEF.          //
! //                                                                  //
! //  J. P. Mellado, David Vanden Abeele (Dec 98)                     //
! //  Modified by David Vanden Abeele (July 99) to include thermal    //
! //  non-equilibrium                                                 //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////

SUBROUTINE PROD_TERM_MOL_DVDA ( T_H, T_E, p, X, W, We, &  
   &  flg_anha, flg_neq, flg_stop, flg_termo ) 
    
!   Global variables
!   ^^^^^^^^^^^^^^^^
  
    USE global_thermo
    USE global_chemco
    USE global_chemco_dvda

!   Local variables
!   ^^^^^^^^^^^^^^^

    IMPLICIT NONE

    REAL*8 T_H, T_E, p, X(nsp), W(nsp), We(nsp)
    REAL*8 TF, TB, N(nsp)
    REAL*8 r_prod, p_prod, kf, kfb, kb, kc

    INTEGER flg_anha, flg_neq, flg_stop, flg_termo
    INTEGER i, j, k

!   Initialization.
!   ^^^^^^^^^^^^^^^

    W = 0.d0; We = 0.d0

!   Set number densities
!   ^^^^^^^^^^^^^^^^^^^^

    N = 0.d0
    CALL NUMBER_DENS (p, T_H, T_E, X, N)

    DO j = 1, cne_nr

       TF   = T_E**QREAC(j,1)*T_H**(1.d0-QREAC(j,1))
       TB   = T_E**QREAC(j,1)*T_H**(1.d0-QREAC(j,1))

       kf = FRRC(j,1) * (TF**(FRRC(j,2))) * DEXP(-FRRC(j,3)/TF)
       kfb= FRRC(j,1) * (TB**(FRRC(j,2))) * DEXP(-FRRC(j,3)/TB)

       CALL EQ_CNT_KC_DVDA(j, TB, kc, flg_anha, flg_neq, flg_stop, flg_termo )
       kb = kfb/kc

       r_prod = 1.0d0
       p_prod = 1.0d0
       DO k = 1, nsp
          IF( CNE_R(j,k).NE.0 ) THEN
             r_prod = r_prod * ( N(k)**(CNE_R(j,k)) )
          ENDIF
          IF( CNE_P(j,k).NE.0 ) THEN
             p_prod = p_prod * ( N(k)**(CNE_P(j,k)) )
          ENDIF  
       ENDDO

       DO i = 1, nsp
          W(i) = W(i) + MMOL(i)*STO(j,i)*(kf*r_prod-kb*p_prod) 
       ENDDO

       ! Calculate contribution of electron-impact reactions separately
       IF ((QREAC(j,1).EQ.1.d0).AND.(QREAC(j,2).EQ.1.d0)) THEN
       DO i = 1, nsp
          We(i) = We(i) + MMOL(i)*STO(j,i)*(kf*r_prod-kb*p_prod) 
       ENDDO
       ENDIF

    ENDDO

END SUBROUTINE PROD_TERM_MOL_DVDA

SUBROUTINE NUMBER_DENS (p, T_H, T_E, X, N)

!   Global variables
!   ^^^^^^^^^^^^^^^^

    USE global_thermo

!   Local variables
!   ^^^^^^^^^^^^^^^

    IMPLICIT NONE
    REAL*8:: p, T_H, T_E, X(nsp)
    REAL*8:: N(nsp)
    REAL*8:: sum, N_tot, T
    INTEGER:: i

!   Take into account the effect of thermal NEQ !
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    N = 0.d0
    sum = 0.d0
    DO i = 1, nsp
       IF (charge(i).ge.0) THEN
          T = T_H
       ELSE
          T = T_E
       ENDIF
       sum = sum + X(i) * T
    ENDDO

    N_tot = p / (RUNIV * sum)
    N = N_tot * X

END SUBROUTINE NUMBER_DENS

