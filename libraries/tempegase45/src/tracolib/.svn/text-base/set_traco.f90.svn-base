SUBROUTINE Set_traco(x,n,cpr,cpv,cpe,HSPC,T_NEQ,oper)
!
! ***************************************************************
! The subroutine contains all the routines required first to 
! initialize all arrays needed by the traco engine. This initialization
! is to be done again if the composition or temperature changes, so it
! cannot be done once and for all !!
!
!   Intent IN: x,n,cpr,cpv,cpe,T_NEQ,flg_NEQ,oper
! ***************************************************************
    USE global_traco
    USE global_thermo, only: nsp,nr
    USE traco_ord_var
    IMPLICIT NONE
    INTEGER :: oper
    REAL(kind=8) :: x(1:nsp),n(1:nsp),T_NEQ(4)    
    REAL(kind=8) :: HSPC(1:nsp),cpr(1:nsp),cpv(1:nsp),cpe(1:nsp)


    INTEGER :: i
    REAL(kind=8) :: HSPC_ord(1:nsp),n_tot
!

    n_tot = 0.0
    DO i=1,nsp
       n_tot = n_tot + n(i)
    END DO

    x_ord=0.
    CALL Order_Data(nsp,index_array,x(1:nsp),x_ord(1:nsp))
    n_ord=0.
    CALL Order_Data(nsp,index_array,n(1:nsp),n_ord(1:nsp))
    HSPC_ord=0.
    CALL Order_Data(nsp,index_array,HSPC(1:nsp),HSPC_ord(1:nsp))
!
!   If OPER=0, then cp - cp(trans) can be passed in cpe_ord
!   cpr_ord then contains cp total (per mole). After cpe_ord
!   is computed then cpr_ord and cpv_ord are set to zero (case of
!   curve fits, valid in equilibrium only)
!
    cpr_ord=0.
    CALL Order_Data (nsp,index_array,cpr(1:nsp),cpr_ord(1:nsp))
!
    IF (oper.EQ.0) THEN
        cpe_ord=cpr_ord-20.786275d0
        cpr_ord=0.0d0
        cpv_ord=0.0d0
    ELSE
        cpv_ord=0.
        CALL Order_Data (nsp,index_array,cpv(1:nsp),cpv_ord(1:nsp))
        cpe_ord=0.
        CALL Order_Data (nsp,index_array,cpe(1:nsp),cpe_ord(1:nsp))
    END IF
!
    omega=0.; Bstar=0.; Cstar=0.
    CALL Set_Omega(n_tot,n_ord(nsp),T_NEQ(1),T_NEQ(4))
!
    delta=0.
    CALL Set_Delta(nsp,T_NEQ(1),T_NEQ(4))
!
    a_average_lambda=0.; a_average_visc=0.
    lambda_factor=0.; visc_factor=0.
    CALL Set_Yos(nsp,x_ord)
!
    HREACT=0.
    CALL Set_Hreact(nr,nsp, NU_ord, HSPC_ord, HREACT)
!
    IF (nions.NE.0) THEN
        qel=0.
        CALL Set_q (nsp,T_NEQ(4),n_ord,n_tot)
    END IF
!

END SUBROUTINE Set_traco
