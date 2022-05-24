! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //  PEGASE 4. CHEMCO LIBRARY                                        //
! //                                                                  //
! //  SUBROUTINE PROD_TERM_MOL                                        //
! //                                                                  //
! //  Inputs:  T           temperature                       K        //
! //           p           pressure                          Pa       //
! //           X           molar fractions                   -        //     
! //           T_NEQ       thermal nonequilibrium                     //
! //                       temperatures                      K        //
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
! //  J. P. Mellado, David Vanden Abeele Dec 98                       //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////

SUBROUTINE PROD_TERM_MOL( T, p, X, W, T_NEQ, &  
   &  flg_anha, flg_neq, flg_stop, flg_termo ) 
    
!   Global variables
!   ^^^^^^^^^^^^^^^^
  
    USE global_thermo
    USE global_chemco

!   Local variables
!   ^^^^^^^^^^^^^^^

    IMPLICIT NONE

    REAL(kind=8) :: T, p, X(nsp), T_NEQ(4)
    REAL(kind=8) :: W(nsp)
    REAL(kind=8) :: r_prod, p_prod, kf, kb, kc(cne_nr), kx(cne_nr)

    INTEGER :: flg_anha, flg_neq, flg_stop, flg_termo
    INTEGER :: i, j

!   Initialization. The third block is to put Kc in cgs system
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    W = 0.0; kc = 0.0; kx = 0.0
    CALL EQ_CNT_KC( T, T_NEQ, kc, flg_anha, flg_neq, flg_stop, flg_termo )

!   Calculation (loop in the cne_nr reactions)
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    DO j = 1, cne_nr
 
!      Calculation of the rate coefficients in cgs
!      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

       kf = 1.0d0 * FRRC(j,1) * (T**(FRRC(j,2))) * DEXP(-FRRC(j,3)/T)
       kb = kf/kc(j)

       r_prod = 1.0d0
       p_prod = 1.0d0
       DO i = 1, nsp
          IF( CNE_R(j,i).NE.0 ) THEN
             r_prod = r_prod * ( (X(i)*p/RUNIV/T)**(CNE_R(j,i)) )
          ENDIF
          IF( CNE_P(j,i).NE.0 ) THEN
             p_prod = p_prod * ( (X(i)*p/RUNIV/T)**(CNE_P(j,i)) )
          ENDIF  
       ENDDO

!      Computing the values
!      ^^^^^^^^^^^^^^^^^^^^

       DO i = 1, nsp
          W(i) = W(i) + STO(j,i)*(kf*r_prod-kb*p_prod) 
       ENDDO

    ENDDO

!   Expressing the production terms in partial densities
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    DO i = 1, nsp
       W(i) = W(i)*MMOL(i)
       !WRITE(*,*)'production term for species ',i,' = ',W(i)!pietro
    ENDDO

    RETURN

END SUBROUTINE PROD_TERM_MOL

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //  PEGASE 4. CHEMCO LIBRARY                                        //
! //                                                                  //
! //  SUBROUTINE JAC_TERM_MOL                                         //
! //                                                                  //
! //  Inputs:  T           temperature                       K        //
! //           p           pressure                          Pa       //
! //           X           molar fractions                   -        //     
! //           T_NEQ       thermal nonequilibrium                     //
! //                       temperatures                      K        //
! //                                                                  //
! //  Outputs: JAC         jacobian                                   //
! //                                                                  //
! //  Flags:   flg_anha    if 1, anharmonicity corrections            //
! //           flg_neq     if 1, thermal non-equilibrium              // 
! //           flg_stop    if 1, stops on subroutine errors           //
! //           flg_termo   mode of computation of equilibrium         // 
! //                       constants                                  //
! //                                                                  //
! //  This subroutine computes the jacobian terms for each species    //
! //  once a kinetic model has been stablished by CHEMODEF. Besides,  //
! //  the jacobian is splitted up in the positive and negative terms. //
! //                                                                  //
! //  J. P. Mellado, David Vanden Abeele Dec 98                       //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////

SUBROUTINE JAC_TERM_MOL ( T, p, X, JAC, T_NEQ, &
   & flg_anha, flg_neq, flg_stop, flg_termo ) 
    
!   Global variables
!   ^^^^^^^^^^^^^^^^
  
    USE global_thermo
    USE global_chemco

!   Local variables
!   ^^^^^^^^^^^^^^^

    IMPLICIT NONE

    REAL(kind=8) :: T, p, X(nsp), T_NEQ(4), JAC(nsp,nsp+2)
    REAL(kind=8) :: dWdT(nsp), dWdp(nsp), dWdx(nsp,nsp) 
    REAL(kind=8) :: r_prod, p_prod, kf, kb, kc(cne_nr), kc1(cne_nr)
    REAL(kind=8) :: T1, dT, epsilon, dummy
    REAL(kind=8) :: dkcdT, dkbdT, dkfdT

    INTEGER :: flg_anha, flg_neq, flg_stop, flg_termo
    INTEGER :: i, j, k, r_sum, p_sum

!   Initialization. The second block is to the finite difference of
!   Kc (thermal eq. assumed). The third block is to put Kc in cgs system
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    dWdT  = 0.0d0   
    dWdp  = 0.0d0
    dWdx  = 0.0d0
    JAC = 0.0d0
    kc = 0.0d0; kc1 = 0.0d0
    CALL EQ_CNT_KC( T, T_NEQ, kc, flg_anha, flg_neq, flg_stop, flg_termo )

    epsilon = 1.0d-6
    dT = epsilon*T
    T1 = T + dT
    dT = T1 - T
    T_NEQ = T1
    CALL EQ_CNT_KC( T1, T_NEQ, kc1, flg_anha, flg_neq, flg_stop, flg_termo )
    T_NEQ = T

!   Calculation (loop in the cne_nr reactions)
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    DO j = 1, cne_nr
 
       kf = 1.0d0 * FRRC(j,1) * (T**(FRRC(j,2))) * DEXP(-FRRC(j,3)/T)
       kb = kf/kc(j)

       dkfdT = kf * (FRRC(j,2)+FRRC(j,3)/T)/T 
       dkcdT = (kc1(j)-kc(j)) / dT
       dkbdT = kb/kf * ( dkfdT - kb*dkcdT )

       r_prod = 1.0d0
       p_prod = 1.0d0
       DO i = 1, nsp
          IF( CNE_R(j,i).NE.0 ) THEN
             r_prod = r_prod * ( (X(i)*p/RUNIV/T)**(CNE_R(j,i)) )
          ENDIF
          IF( CNE_P(j,i).NE.0 ) THEN
             p_prod = p_prod * ( (X(i)*p/RUNIV/T)**(CNE_P(j,i)) )
          ENDIF  
       ENDDO

!      Computing the values
!      ^^^^^^^^^^^^^^^^^^^^

       DO i = 1, nsp
          r_sum = 0.0d0
          p_sum = 0.0d0
          DO k = 1, nsp
             IF( X(k).NE.0 ) THEN
               dWdX(i,k) = dWdX(i,k) + STO(j,i)*(kf*r_prod*CNE_R(j,k)-kb*p_prod*CNE_P(j,k))
             ENDIF
             r_sum = r_sum + CNE_R(j,k)
             p_sum = p_sum + CNE_P(j,k)
          ENDDO
          dummy = STO(j,i)*(kf*r_prod*r_sum-kb*p_prod*p_sum)
          dWdp(i) = dWdp(i) + dummy
          dWdT(i) = dWdT(i) - dummy/T + STO(j,i)*(dkfdT*r_prod-dkbdT*p_prod)
       ENDDO

    ENDDO

!   Final operations
!   ^^^^^^^^^^^^^^^^

    DO i = 1, nsp
       dWdp(i) = dWdp(i)*MMOL(i)/p
       dWdT(i) = dWdT(i)*MMOL(i)
       DO j = 1, nsp
          dWdX(i,j) = dWdX(i,j)*MMOL(i)/X(j)
       ENDDO
    ENDDO

!   Building of the jacobian
!   ^^^^^^^^^^^^^^^^^^^^^^^^

    DO i = 1, nsp
       JAC(i,1) = dWdT(i)
       JAC(i,2) = dWdp(i)
       DO k = 1, nsp
          JAC(i,k+2) = dWdx(i,k)
       ENDDO
    ENDDO

    RETURN

END SUBROUTINE JAC_TERM_MOL
