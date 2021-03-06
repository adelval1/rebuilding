! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //  PEGASE 4. CHEMCO LIBRARY                                        //
! //                                                                  //
! //  SUBROUTINE PROD_TERM_DEN                                        //
! //                                                                  //
! //  Inputs:  T           temperature                       K        //
! //           RHOSP       array with the partial                     //
! //                       densities of the species          kg/m^3   // 
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
! //  once a kinetic model has been stablished.                       //
! //                                                                   //
! //  J. P. Mellado,Paolo Barbante Dec 98                             //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////


SUBROUTINE PROD_TERM_DEN(T,RHOSP,W,T_NEQ,flg_anha,flg_neq,flg_stop,flg_termo)  
    
!   Global variables
!   ^^^^^^^^^^^^^^^^
  
    USE global_thermo
    USE global_chemco
    IMPLICIT NONE
    INTEGER :: flg_anha, flg_neq, flg_stop, flg_termo
    REAL(kind=8) :: T, RHOSP(:), W(:), T_NEQ(4)

!   Local variables
!   ^^^^^^^^^^^^^^^

    INTEGER :: i, j
    REAL(kind=8) :: r_prod, p_prod, kf, kb, kc(cne_nr)

!   Initialization.
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    W = 0.0d0

    CALL EQ_CNT_KC( T, T_NEQ, kc, flg_anha, flg_neq, flg_stop, flg_termo )


!   Calculation (loop in the cne_nr reactions)
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    DO j = 1, cne_nr
 
!      Calculation of the rate coefficients in cgs
!      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       kf = FRRC(j,1) * (T**(FRRC(j,2))) * exp(-FRRC(j,3)/T)
       kb = kf/kc(j)

       r_prod = 1.0d0
       p_prod = 1.0d0
       DO i = 1, nsp
          IF( CNE_R(j,i).NE.0 ) THEN
             r_prod = r_prod * ( (RHOSP(i)/MMOL(i))**(CNE_R(j,i)) )
          END IF
          IF( CNE_P(j,i).NE.0 ) THEN
             p_prod = p_prod * ( (RHOSP(i)/MMOL(i))**(CNE_P(j,i)) )
          END IF  
       END DO


!      Computing values
!      ^^^^^^^^^^^^^^^^

       DO i = 1, nsp
          W(i) = W(i) + STO(j,i)*(kf*r_prod-kb*p_prod) 
       END DO

    END DO

!   Expressing the production terms in partial densities form
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    DO i = 1, nsp
       W(i) = W(i)*MMOL(i)
    END DO
    

END SUBROUTINE PROD_TERM_DEN 

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //  PEGASE 4. CHEMCO LIBRARY                                         //
! //                                                                  //
! //  SUBROUTINE JAC_TERM_DEN                                         //
! //                                                                  //
! //  Inputs:  T           temperature                       K        //
! //           RHOSP       array with the partial                     //
! //                       densities of the species          kg/m^3   // 
! //           T_NEQ       thermal nonequilibrium                     // 
! //                       temperatures                      K        //
! //                                                                  //
! //  Outputs: RHO_JACP    positive jacobian wrt density              //
! //           RHO_JACM    negative jacobian wrt density              //
! //           T_JAC       jacobian wrt temperature                   //
! //                                                                  //
! //  Flags:   flg_anha    if 1, anharmonicity corrections            //
! //           flg_neq     if 1, thermal non-equilibrium              // 
! //           flg_stop    if 1, stops on subroutine errors           //
! //           flg_termo   mode of computation of equilibrium         // 
! //                       constants                                  //
! //                                                                  //
! //  This subroutine computes the jacobian terms for each species    //
! //  once a kinetic model has been stablished by CHEMODEF. Besides,  //
! //  the jacobian is  splitted up in the positive and negative terms.//
! //                                                                  //
! //  J. P. Mellado,Paolo Barbante Dec. 98                            //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////


    SUBROUTINE JAC_TERM_DEN(T,RHOSP,RHO_JACP,RHO_JACM,T_JAC,T_NEQ,   &  
   &                        flg_anha,flg_neq,flg_stop,flg_termo ) 
    
!   Global variables
!   ^^^^^^^^^^^^^^^^
  
    USE global_thermo
    USE global_chemco
    IMPLICIT NONE
    INTEGER :: flg_anha, flg_neq, flg_stop, flg_termo
    REAL(kind=8) :: T, RHOSP(:), T_NEQ(4)
    REAL(kind=8) :: RHO_JACP(:,:), RHO_JACM(:,:), T_JAC(:)

!   Local variables
!   ^^^^^^^^^^^^^^^

    INTEGER i, j, k
    REAL(kind=8) :: dWdT(nsp), dWdrhop(nsp,nsp), dWdrhom(nsp,nsp)
    REAL(kind=8) :: r_prod, p_prod, kf, kb, kc(cne_nr)
    REAL(kind=8) T1, dT, epsilon, kc1(cne_nr)
    REAL(kind=8) dkcdT, dkbdT, dkfdT

!   Initialization.
!   The second one is to use a finite-difference approximation
!   for the derivative of Kc (thermal eq. assumed)
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    dWdT   = 0.0d0   
    dWdrhop = 0.0d0
    dWdrhom = 0.0d0
    RHO_JACP  = 0.0d0
    RHO_JACM  = 0.0d0
    T_JAC  = 0.0d0


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
 
!      Calculation of the rate coefficients in cgs (now SI)
!      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

       kf = FRRC(j,1) * (T**(FRRC(j,2))) * exp(-FRRC(j,3)/T)
       kb = kf/kc(j)
       dkfdT = kf * (FRRC(j,2)+FRRC(j,3)/T)/T 
       dkcdT = (kc1(j)-kc(j)) / dT
       dkbdT = ( dkfdT - kb*dkcdT )/kc(j)

       r_prod = 1.0d0
       p_prod = 1.0d0
       DO i = 1, nsp
          IF( CNE_R(j,i).NE.0 ) THEN
             r_prod = r_prod * ( (RHOSP(i)/MMOL(i))**(CNE_R(j,i)) )
          END IF
          IF( CNE_P(j,i).NE.0 ) THEN
             p_prod = p_prod * ( (RHOSP(i)/MMOL(i))**(CNE_P(j,i)) )
          END IF  
       END DO

!      Computing values
!      ^^^^^^^^^^^^^^^^


       DO i = 1, nsp
          dWdT(i) = dWdT(i) + STO(j,i)*(dkfdT*r_prod-dkbdT*p_prod)
          DO k = 1, nsp
             IF( RHOSP(k).NE.0 ) THEN
                if(STO(j,i)>=0) then
                   dWdrhop(i,k) = dWdrhop(i,k) +   &
&                     STO(j,i)*kf*r_prod*CNE_R(j,k)/RHOSP(k)
                   dWdrhom(i,k) = dWdrhom(i,k) +   &
&                     STO(j,i)*(-kb*p_prod*CNE_P(j,k))/RHOSP(k)
                else
                   dWdrhop(i,k) = dWdrhop(i,k) +   &
&                     STO(j,i)*(-kb*p_prod*CNE_P(j,k))/RHOSP(k)
                   dWdrhom(i,k) = dWdrhom(i,k) +   &
&                     STO(j,i)*kf*r_prod*CNE_R(j,k)/RHOSP(k)
                end if
             END IF
          END DO
       END DO

    END DO

!   Final operations
!   ^^^^^^^^^^^^^^^^

    DO i = 1, nsp
       dWdT(i) = dWdT(i)*MMOL(i)
       DO j = 1, nsp
          dWdrhop(i,j) = dWdrhop(i,j)*MMOL(i)
          dWdrhom(i,j) = dWdrhom(i,j)*MMOL(i)
       END DO
    END DO

!   Writing in the jacobians

    DO i=1,nsp
       T_JAC(i) = dWdT(i)
       DO j=1,nsp
          RHO_JACP(i,j) = dWdrhop(i,j)
          RHO_JACM(i,j) = dWdrhom(i,j)
       END DO
    END DO


END SUBROUTINE JAC_TERM_DEN 
