MODULE interf_traco
!
INTERFACE
!
!   Definition and allocation of traco variables when a mixture
!   is set for the calculation, and deallocation after use.
!   -----------------------------------------------------------
    SUBROUTINE tracodef (n)
        INTEGER n
    END SUBROUTINE tracodef
!
    SUBROUTINE tracostop
    END SUBROUTINE tracostop
!
!   Initialization of all values required by the transport properties 
!   functions and subroutines, i.e. collision integral stuff and the
!   like.
!   -----------------------------------------------------------------
    SUBROUTINE Set_traco(x,n,cpr,cpv,cpe,HSPC,T_NEQ,oper)
        USE global_thermo, only: nsp,nr
        INTEGER :: oper
        REAL(kind=8) :: x(1:nsp),n(1:nsp),T_NEQ(4)
        REAL(kind=8) :: HSPC(1:nsp),cpr(1:nsp),cpv(1:nsp),cpe(1:nsp)
    END SUBROUTINE Set_traco
!
!   Transport properties calculation subroutines and functions
!   ----------------------------------------------------------
    REAL(kind=8) FUNCTION viscosity()

    END FUNCTION viscosity
!
    REAL(kind=8) FUNCTION thermal_cond_frz(T_NEQ)
        REAL(kind=8) :: T_NEQ(1:4)
    END FUNCTION thermal_cond_frz 
!
    REAL(kind=8) FUNCTION thermal_cond_tot(T_NEQ)
        REAL(kind=8) :: T_NEQ(1:4)
    END FUNCTION thermal_cond_tot 

    REAL(kind=8) FUNCTION electrical_conductivity (ne, T_e)
          REAL(kind=8) :: ne,T_e
    END FUNCTION electrical_conductivity
!
    REAL(kind=8) FUNCTION Debye_Length_Devoto (ne,n,T_e)
        REAL(kind=8) :: ne,n,T_e
    END FUNCTION Debye_Length_Devoto
!
    REAL(kind=8) FUNCTION Devoto_taueh (T_h, T_e)
        REAL(kind=8) :: T_h, T_e
    END FUNCTION Devoto_taueh

    REAL(kind=8) FUNCTION Millikan_tauvh (p, T_h)
        REAL(kind=8) :: p, T_h
    END FUNCTION Millikan_tauvh
!
    SUBROUTINE mfp_Yos (n_tot, mfp_array)
        USE global_thermo
        REAL(kind=8) :: n_tot, mfp_array (1:nsp)
    END SUBROUTINE mfp_Yos
!
    REAL(kind=8) FUNCTION PRANDTL_NUMBER (mu,cp,k)
        REAL(kind=8) mu,cp,k
    END FUNCTION PRANDTL_NUMBER
!
    REAL(kind=8) FUNCTION LEWIS_NUMBER (k,rho,cp,d)
        REAL(kind=8) k,rho,cp,d
    END FUNCTION LEWIS_NUMBER
!
    SUBROUTINE diffusion_coefficients_multi(p,TNEQ,mm,d)
        USE global_thermo
        REAL(kind=8) p,TNEQ(4),mm,d(nsp)
    END SUBROUTINE diffusion_coefficients_multi

    SUBROUTINE Dvda_stefmax (p,x,z,T,d,diff_ord,J_old,J_new,E_amb, &
      & tresh,n_it,l_ref,resit,resmass,rescharge,k)
        USE global_thermo
        REAL(kind=8) p, x (nsp), z (nsp), T (nsp), d (nsp)
        REAL(kind=8) diff_ord(nsp,nsp), J_old(nsp), J_new(nsp), E_amb
        REAL(kind=8) tresh, l_ref, resit, resmass, rescharge
        INTEGER   :: n_it, k
    END SUBROUTINE Dvda_stefmax
!
!   Derived functions for the individual calculation of intermediate
!   components
!   --------------------------------------------------------------
    SUBROUTINE Yos_lambda_frozen (k_trans, k_rot, k_vib, k_electronic)
        REAL(kind=8) :: k_trans, k_rot, k_vib, k_electronic
    END SUBROUTINE Yos_lambda_frozen 
!
    REAL(kind=8) FUNCTION Yos_lambda_reactive (T)
        REAL(kind=8) :: T
    END FUNCTION Yos_lambda_reactive
!
    REAL(kind=8) FUNCTION Brokaw_lambda_reactive (T)
        REAL(kind=8) :: T
    END FUNCTION Brokaw_lambda_reactive
!
    REAL(kind=8) FUNCTION Devoto_lambda_electron (ne, T_e)
        REAL(kind=8) :: ne,T_e
    END FUNCTION Devoto_lambda_electron
!
    SUBROUTINE set_binary_diffusion(p,TNEQ,dij)
        USE global_thermo, only: nsp,KUNIV
        REAL(kind=8) :: p,TNEQ(4),dij(nsp,nsp)
    END SUBROUTINE set_binary_diffusion

    SUBROUTINE diffusion_coefficients_binary(p,TNEQ,dij_peg)
        USE global_thermo, only: nsp
        REAL(kind=8) :: p,TNEQ(4),dij_peg(nsp,nsp)
    END SUBROUTINE diffusion_coefficients_binary

!
!   Ordering and reordering of the species
!   --------------------------------------
    SUBROUTINE Order_Data (nsp,index_array, unordered_array, ordered_array)
        INTEGER :: nsp,index_array (1:nsp)
        REAL(kind=8) :: unordered_array (1:nsp),ordered_array (1:nsp)
    END SUBROUTINE Order_Data
!
    SUBROUTINE Reorder_Data (nsp,index_array, ordered_array, unordered_array)
        INTEGER :: nsp,index_array (1:nsp)
        REAL(kind=8) :: unordered_array (1:nsp),ordered_array (1:nsp)
    END SUBROUTINE Reorder_Data
!
END INTERFACE

END MODULE interf_traco
!
