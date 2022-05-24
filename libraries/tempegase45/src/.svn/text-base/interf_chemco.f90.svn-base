! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO    L I B R A R Y             //
! //                                                                  //
! //  MODULE interf_thermo                                            //
! //                                                                  //
! //  Initialization/deinitilization of chemical kinetic data         //
! //                                                                  //
! //  J. P. Mellado, July 98                                          //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////

MODULE interf_chemco

INTERFACE

    SUBROUTINE chemcodef 
    END SUBROUTINE chemcodef

    SUBROUTINE chemcostop
    END SUBROUTINE chemcostop

    SUBROUTINE prod_term_den( T, RHOSP, W, TNEQ,&
   &                           flg_anha, flg_neq, flg_stop, flg_termo )
       REAL(kind=8) T, RHOSP(:), W(:), TNEQ(4)
       INTEGER flg_anha, flg_neq, flg_stop, flg_termo
    END SUBROUTINE

    SUBROUTINE jac_term_den( T, RHOSP, RHO_JACP, RHO_JACM, T_JAC, TNEQ,&
   &                           flg_anha, flg_neq, flg_stop, flg_termo )
       REAL(kind=8) T, RHOSP(:), RHO_JACP(:,:), RHO_JACM(:,:), T_JAC(:), TNEQ(4)
       INTEGER flg_anha, flg_neq, flg_stop, flg_termo
    END SUBROUTINE

    SUBROUTINE prod_term_mol( T, p, X, W, TNEQ,&
   &                          flg_anha, flg_neq, flg_stop, flg_termo )
       REAL(kind=8) T, p, X(:), W(:), TNEQ(4)
       INTEGER flg_anha, flg_neq, flg_stop, flg_termo
    END SUBROUTINE

    SUBROUTINE jac_term_mol( T, p, X, JAC, TNEQ,&
   &                          flg_anha, flg_neq, flg_stop, flg_termo )
       REAL(kind=8) T, p, X(:), JAC(:,:), TNEQ(4)
       INTEGER flg_anha, flg_neq, flg_stop, flg_termo
    END SUBROUTINE

END INTERFACE

END MODULE
