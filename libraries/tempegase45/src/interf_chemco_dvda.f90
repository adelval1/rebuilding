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

MODULE interf_chemco_dvda

INTERFACE

    SUBROUTINE chemcodef_dvda 
    END SUBROUTINE chemcodef_dvda

    SUBROUTINE chemcostop_dvda
    END SUBROUTINE chemcostop_dvda

    SUBROUTINE prod_term_mol_dvda ( T_H, T_E, p, X, W, We, &
    &  flg_anha, flg_neq, flg_stop, flg_termo )
       USE global_thermo
       REAL*8 T_H, T_E, p, X(nsp), W(nsp), We(nsp)
       INTEGER flg_anha, flg_neq, flg_stop, flg_termo
    END SUBROUTINE prod_term_mol_dvda

END INTERFACE

END MODULE
