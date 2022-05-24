!***********************************************************************
module OPTIONS_par

integer choice,nspecies
integer Kp_choice,neq_fr

end module OPTIONS_par
!***********************************************************************

!***********************************************************************
module PEG_DIR_par

character*80 thermo_dir,traco_dir,mixname,chem_dir

end module PEG_DIR_par
!***********************************************************************

!**********************************************************************
module FLAGS_par

integer :: flg_mode,flg_anha,flg_oper,flg_termo,flg_traco,flg_stop
integer :: sonine,nspec,maxlvl

end module FLAGS_par
!**********************************************************************

!**********************************************************************
module REBUILDING_var

real*8 Kp,T_reb,gamma,Qwc_reb,Qwd_reb,result_gam
real*8 error_Kp,stop_Kp,Kp_old
real*8, allocatable :: gamma_vector(:)
integer N_gamma

end module REBUILDING_var
!**********************************************************************

!**********************************************************************
module CONVERGENCE_par

integer count_iter
real*8 errorstop,alpha,epsilon

end module CONVERGENCE_par
!**********************************************************************

!**********************************************************************
module EXPERIMENTS_par

real*8 Tw,Pstatic
real*8 deltaPs_1,Qw_exp_1
real*8 emissivity_smpl,Tw_smpl
real*8 Rm,Rm_2,Rm_3,Rm_smpl

end module EXPERIMENTS_par
!**********************************************************************

!**********************************************************************
module OUTEREDGE_var

real*8 T_ext,Vs_ext,rho_ext,mu_ext,h_ext,Ve_ext
real*8, allocatable:: concent_ext(:)

end module OUTEREDGE_var
!**********************************************************************

!**********************************************************************
Module HEATFLUX_var

real*8 P,rho_rebuild,visco_rebuild,h_rebuild,flux
real*8,allocatable:: concent_rebuild(:)

end module HEATFLUX_var
!**********************************************************************

!**********************************************************************
module BL_par

integer npoints,nreacwall

end module BL_par
!**********************************************************************

!**********************************************************************
module ICP_par

real*8 delta_1,u1e_1,ve_1,u1y_1,NDP_5_1
real*8 u1e_2,u1e_3
real*8 delta_smpl,u1e_smpl,ve_smpl,u1y_smpl,NDP_5_smpl

end module ICP_par
!**********************************************************************

!**********************************************************************
module CONVERGENCE_ENTH_par

real*8 enthalpyerror,enthepsilon,enthalpha

end module CONVERGENCE_ENTH_par
!**********************************************************************

!**********************************************************************
module NAMES_par

character*30 probe_name1,probe_name2,probe_name3
end module NAMES_par
!**********************************************************************

!**********************************************************************
module ABACUS_par

real*8 Tw_min,Tw_max
integer ntw,ngamma
real*8, allocatable:: gammaarray(:)

end module
!**********************************************************************

!**********************************************************************
module INIT_DIST_par

real*8 he,T_0

end module INIT_DIST_par
