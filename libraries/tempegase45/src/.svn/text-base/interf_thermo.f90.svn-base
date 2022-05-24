! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO    L I B R A R Y             //
! //                                                                  //
! //  MODULE interf_thermo                                            //
! //                                                                  //
! // Definition of the function and procedures interfaces accessible  //
! // from the thermodynamic library of Pegase.                        //
! //                                                                  //
! // Benoit Bottin, 30/09/97                                          //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
MODULE interf_thermo
!
INTERFACE
!
!   initialization/deinitialization of thermodynamic data
!   =====================================================
    SUBROUTINE termodef (thermo_dir,traco_dir,chemo_dir,mixname,flg_anha,&
       &flg_oper,flg_termo,flg_traco,flg_stop,maxlvl)
        INTEGER flg_anha,flg_oper,flg_termo,flg_traco,flg_stop,maxlvl
        CHARACTER*80 thermo_dir,traco_dir,chemo_dir,mixname
    END SUBROUTINE termodef
!
    SUBROUTINE termostop
    END SUBROUTINE termostop
!
!   Conversion procedures
!   =====================
    REAL(kind=8) FUNCTION mol_to_mass (val,isp,mm,flg_stop)
        REAL(kind=8) val,mm
        INTEGER isp,flg_stop
    END FUNCTION mol_to_mass
!
    REAL(kind=8) FUNCTION mass_to_mol (val,isp,mm,flg_stop)
        REAL(kind=8) val,mm
        INTEGER isp,flg_stop
    END FUNCTION mass_to_mol
!
    SUBROUTINE molefrac_to_massfrac (x,mm,y)
        REAL(kind=8) x(:),y(:),mm
    END SUBROUTINE molefrac_to_massfrac
!
    SUBROUTINE massfrac_to_molefrac (y,mm,x)
        REAL(kind=8) x(:),y(:),mm
    END SUBROUTINE massfrac_to_molefrac
!
    REAL(kind=8) FUNCTION number_density (p,T)
        REAL(kind=8) p,T
    END FUNCTION number_density
!
    REAL(kind=8) FUNCTION concentration (p,T)
        REAL(kind=8) p,T
    END FUNCTION concentration 
!
    REAL(kind=8) FUNCTION species_density (isp,p,T)
        REAL(kind=8) p,T
        INTEGER isp
    END FUNCTION species_density
!
    REAL(kind=8) FUNCTION species_pressure (isp,rho,T)
        REAL(kind=8) rho,T
        INTEGER isp
    END FUNCTION species_pressure
!
!   Perfect gas equation of state and its parameters
!   ================================================
    REAL(kind=8) FUNCTION mixture_pressure (rho,T,y,flg_stop)
        REAL(kind=8) rho,T,y(:)
        INTEGER flg_stop
    END FUNCTION mixture_pressure
!
    REAL(kind=8) FUNCTION mixture_density (p,T,x,flg_stop)
        REAL(kind=8) p,T,x(:)
        INTEGER flg_stop
    END FUNCTION mixture_density
!
    REAL(kind=8) FUNCTION mixture_molarmass (x,mode,flg_stop)
        REAL(kind=8) x(:)
        INTEGER flg_stop,mode
    END FUNCTION mixture_molarmass
!
    REAL(kind=8) FUNCTION mixture_rspecific (x,mode,flg_stop)
        REAL(kind=8) x(:)
        INTEGER flg_stop,mode
    END FUNCTION mixture_rspecific
!
    REAL(kind=8) FUNCTION mixture_moleratio (x,mode,flg_stop)
        REAL(kind=8) x(:)
        INTEGER flg_stop,mode
    END FUNCTION mixture_moleratio
!
!   Energy and enthalpy
!   ===================
    SUBROUTINE species_energy (isp,T,TNEQ,&
       &flg_anha,mode,flg_neq,flg_stop,flg_termo,ener)
        INTEGER isp,flg_anha,mode,flg_neq,flg_stop,flg_termo
        REAL(kind=8) T,TNEQ(4),ener(8)
    END SUBROUTINE species_energy
!
    SUBROUTINE mixture_energy (T,TNEQ,x,&
       &flg_anha,mode,flg_neq,flg_stop,flg_termo,ener)
        INTEGER flg_anha,mode,flg_neq,flg_stop,flg_termo
        REAL(kind=8) T,TNEQ(4),ener(8),x(:)
    END SUBROUTINE mixture_energy
!
    SUBROUTINE species_enthalpy (isp,T,TNEQ,&
       &flg_anha,mode,flg_neq,flg_stop,flg_termo,enth)
        INTEGER isp,flg_anha,mode,flg_neq,flg_stop,flg_termo
        REAL(kind=8) T,TNEQ(4),enth(8)
    END SUBROUTINE species_enthalpy
!
    SUBROUTINE mixture_enthalpy (T,TNEQ,x,&
       &flg_anha,mode,flg_neq,flg_stop,flg_termo,enth)
        INTEGER flg_anha,mode,flg_neq,flg_stop,flg_termo
        REAL(kind=8) T,TNEQ(4),enth(8),x(:)
    END SUBROUTINE mixture_enthalpy
!
!   Entropy and free energies
!   =========================
    SUBROUTINE species_entropy (isp,var,T,TNEQ,&
       &flg_anha,mode,flg_neq,flg_stop,flg_termo,entro)
        INTEGER isp,flg_anha,mode,flg_neq,flg_stop,flg_termo
        REAL(kind=8) var,T,TNEQ(4),entro(8)
    END SUBROUTINE species_entropy
!
    SUBROUTINE mixture_entropy (var,T,TNEQ,x,&
       &flg_anha,mode,flg_neq,flg_stop,flg_termo,entro)
        INTEGER flg_anha,mode,flg_neq,flg_stop,flg_termo
        REAL(kind=8) var,T,TNEQ(4),entro(8),x(:)
    END SUBROUTINE mixture_entropy 
!
    SUBROUTINE species_gibbs (isp,var,T,TNEQ,&
       &flg_anha,mode,flg_neq,flg_stop,flg_termo,gibb)
        INTEGER isp,flg_anha,mode,flg_neq,flg_stop,flg_termo
        REAL(kind=8) var,T,TNEQ(4),gibb(8)
    END SUBROUTINE species_gibbs
!
    SUBROUTINE mixture_gibbs (var,T,TNEQ,x,&
       &flg_anha,mode,flg_neq,flg_stop,flg_termo,gibbs)
        INTEGER flg_anha,mode,flg_neq,flg_stop,flg_termo
        REAL(kind=8) var,T,TNEQ(4),gibbs(8),x(:)
    END SUBROUTINE mixture_gibbs
!
      SUBROUTINE species_helmholtz (isp,var,T,TNEQ,&
     &flg_anha,mode,flg_neq,flg_stop,flg_termo,helm)
        INTEGER isp,flg_anha,mode,flg_neq,flg_stop,flg_termo
        REAL(kind=8) var,T,TNEQ(4),helm(8)
    END SUBROUTINE species_helmholtz
!
    SUBROUTINE mixture_helmholtz (var,T,TNEQ,x,&
       &flg_anha,mode,flg_neq,flg_stop,flg_termo,helm)
        INTEGER flg_anha,mode,flg_neq,flg_stop,flg_termo
        REAL(kind=8) var,T,TNEQ(4),helm(8),x(:)
    END SUBROUTINE mixture_helmholtz
!
!   Prerequisite finite difference "sensitivities" procedures
!   =========================================================
    SUBROUTINE sensitivities_r (p,rho,T,TNEQ,ener,enth,mm,y,findiff, &
       &epsilon,flg_anha,flg_neq,flg_stop,flg_termo,flg_log,niter,   &
       &dmdr,dpdr,dedr,dhdr,detdr)                               
        INTEGER flg_anha,flg_neq,flg_stop,flg_termo,flg_log
        REAL(kind=8) p,rho,T,TNEQ(4),ener(8),enth(8),mm,y(:),epsilon,findiff
        REAL(kind=8) dmdr,dpdr,dedr(8),dhdr(8),detdr(8)
    END SUBROUTINE sensitivities_r
!
    SUBROUTINE sensitivities_t (p,rho,T,TNEQ,ener,enth,mm,y,findiff, &
       &epsilon,flg_anha,flg_neq,flg_stop,flg_termo,flg_log,niter,   &
       &dmdt,dpdt,dedt,dhdt,detdt)
        INTEGER flg_anha,flg_neq,flg_stop,flg_termo,flg_log
        REAL(kind=8) p,rho,T,TNEQ(4),ener(8),enth(8),mm,y(:),epsilon,findiff
        REAL(kind=8) dmdt,dpdt,dedt(8),dhdt(8),detdt(8)
    END SUBROUTINE sensitivities_t
!
!   Specific heat at constant volume
!   ================================
    SUBROUTINE species_cv (isp,T,TNEQ,&
       &flg_anha,mode,flg_neq,flg_stop,flg_termo,cv)
        INTEGER :: isp,flg_anha,mode,flg_neq,flg_stop,flg_termo
        REAL(kind=8)  :: T,TNEQ(4),cv(8)
    END SUBROUTINE species_cv
!
    SUBROUTINE cv_frozen (T,TNEQ,x,&
       &flg_anha,mode,flg_neq,flg_stop,flg_termo,cv)
        INTEGER flg_anha,mode,flg_neq,flg_stop,flg_termo
        REAL(kind=8) T,TNEQ(4),cv(8),x(:)
    END SUBROUTINE cv_frozen
!
    SUBROUTINE cv_equilibrium (mode,ener,mm,dedt,dmdt,cv)
        INTEGER mode
        REAL(kind=8) cv(8),ener(8),mm,dedt(8),dmdt
    END SUBROUTINE cv_equilibrium
!
!   Specific heat at constant pressure
!   ==================================
    SUBROUTINE species_cp (isp,T,TNEQ,&
       &flg_anha,mode,flg_neq,flg_stop,flg_termo,cp)
        INTEGER :: isp,flg_anha,mode,flg_neq,flg_stop,flg_termo
        REAL(kind=8)  :: T,TNEQ(4),cp(8)
    END SUBROUTINE species_cp
!
    SUBROUTINE cp_frozen (T,TNEQ,x,&
       &flg_anha,mode,flg_neq,flg_stop,flg_termo,cp)
        INTEGER flg_anha,mode,flg_neq,flg_stop,flg_termo
        REAL(kind=8) T,TNEQ(4),cp(8),x(:)
    END SUBROUTINE cp_frozen
!
    SUBROUTINE cp_equilibrium (mode,enth,mm,dhdt,dhdr,dmdt,dmdr,&
       &dpdt,dpdr,cp)
        INTEGER :: mode
        REAL(kind=8) :: cp(8),enth(8),dhdt(8),dhdr(8),dmdt,dmdr,dpdt,dpdr,mm
    END SUBROUTINE cp_equilibrium
!
!   Equilibrium composition procedures
!   ==================================
    SUBROUTINE eq_constants (T,TNEQ,p,rho,k,&
       &flg_anha,flg_neq,flg_stop,flg_termo)
        INTEGER flg_anha,flg_neq,flg_stop,flg_termo
        REAL(kind=8) T,TNEQ(4),p,rho,k(:,:)
    END SUBROUTINE eq_constants
!
    SUBROUTINE massfrac (y,rho,T,TNEQ,eps,flg_log,&
       &flg_anha,flg_neq,flg_stop,flg_termo,flg_guess,y_guess,niter)
        INTEGER flg_log,flg_stop,flg_anha,flg_neq,flg_termo,flg_guess,niter
        REAL(kind=8) y(:),rho,T,TNEQ(4),y_guess(:),eps
    END SUBROUTINE massfrac
!
    SUBROUTINE molefrac (x,p,T,zo,TNEQ,eps,flg_log,&
       &flg_anha,flg_neq,flg_stop,flg_termo,flg_guess,x_guess,niter)
        INTEGER flg_log,flg_stop,flg_anha,flg_neq,flg_termo,flg_guess,niter
        REAL(kind=8) x(:),p,T,TNEQ(4),zo,x_guess(:),eps
    END SUBROUTINE molefrac
!
!   Miscellaneous procedures
!   ========================
    SUBROUTINE sound_speed_equilibrium (enth,dpdt,dpdr,detdt,detdr,&
       &a,kappa,khi)
        REAL(kind=8) enth(8),dpdt,dpdr,detdt(8),detdr(8),a,kappa,khi
    END SUBROUTINE sound_speed_equilibrium 
!
    REAL(kind=8) FUNCTION dissoc_enthalpy (x,zo,mode,method)
        REAL(kind=8) x(:),zo
        INTEGER mode,method
    END FUNCTION dissoc_enthalpy
!
    SUBROUTINE partition_function (isp,var,T,TNEQ,flg_anha,mode,&
       &flg_neq,flg_stop,flg_termo,parfn,dparfn,ddparfn,dlnpf,ddlnpf)
        INTEGER isp,flg_anha,mode,flg_neq,flg_stop,flg_termo
        REAL(kind=8) var,T,TNEQ(4)
        REAL(kind=8) parfn(8),dparfn(8),ddparfn(8),dlnpf(8),ddlnpf(8)
    END SUBROUTINE partition_function
!
!   Gupta curve fits procedures
!   ===========================
    SUBROUTINE load_gupta_fits (traco_dir,flg_stop)
        CHARACTER*80 :: traco_dir
        INTEGER :: flg_stop
    END SUBROUTINE load_gupta_fits
!
    REAL(kind=8) FUNCTION gupta_visco(p,T)                                          
        REAL(kind=8) p,T
    END FUNCTION gupta_visco
!
    REAL(kind=8) FUNCTION gupta_k(p,T)                                          
        REAL(kind=8) p,T
    END FUNCTION gupta_k
!
    REAL(kind=8) FUNCTION gupta_cp(p,T)                                          
        REAL(kind=8) p,T
    END FUNCTION gupta_cp
!
    REAL(kind=8) FUNCTION gupta_prandtl(p,T)                                          
        REAL(kind=8) p,T
    END FUNCTION gupta_prandtl
!
END INTERFACE
END MODULE
