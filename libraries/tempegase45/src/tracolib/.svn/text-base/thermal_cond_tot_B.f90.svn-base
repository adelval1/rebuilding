REAL(kind=8) FUNCTION thermal_cond_tot_b(T_NEQ)

!*****************************************************************
! Subroutine added by B. Bottin to get one single call computing the
! total thermal conductivity contributions. The individual contributions
! can be accessed to by individual calls to the functions below !
!*****************************************************************
    USE global_thermo, only: nsp
    USE global_traco, only: nions
    USE traco_ord_var, only: n_ord
    IMPLICIT NONE
    REAL(kind=8) :: T_NEQ(1:4)

    REAL(kind=8) :: Devoto_lambda_electron,Brokaw_lambda_reactive_nanni 
    REAL(kind=8) :: k_trans,k_rot,k_vib,k_electronic,k_electron,k_reactive

    k_trans=0.
    k_rot=0.
    k_vib=0.
    k_electron=0.
    k_electronic=0.
    k_reactive=0.

   ! write(*,*)"thermal_cond_tot_nanni started - B&B"   

    CALL Yos_lambda_frozen(k_trans,k_rot,k_vib,k_electronic)


    IF (nions.NE.0) THEN
        k_electron = Devoto_lambda_electron(n_ord(nsp),T_NEQ(4))
    ELSE
        k_electron = 0.
    END IF
    
    k_reactive = Brokaw_lambda_reactive_nanni (T_NEQ(1))

  !  write(*,*)"k_reactive",k_reactive    

    thermal_cond_tot_b = k_trans + k_rot + k_vib + k_electronic &
   & + k_electron + k_reactive


END FUNCTION thermal_cond_tot_b
