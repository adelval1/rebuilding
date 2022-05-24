subroutine c_mixture_enthalpy(t,tneq,c,flg_anha,mode,flg_neq,flg_stop,flg_termo,mixh)

   use interf_thermo
   use global_thermo, only: nsp
   implicit none
   integer :: flg_anha,mode,flg_neq,flg_stop,flg_termo
   real(kind=8) :: t,tneq(4),mixh(8),c(nsp)

   call mixture_enthalpy(t,tneq,c,flg_anha,mode,flg_neq,flg_stop,flg_termo,mixh)

end subroutine c_mixture_enthalpy




subroutine c_mixture_energy(t,tneq,c,flg_anha,mode,flg_neq,flg_stop,flg_termo,mixe)

   use interf_thermo
   use global_thermo, only: nsp
   implicit none
   integer :: flg_anha,mode,flg_neq,flg_stop,flg_termo
   real(kind=8) :: t,tneq(4),mixe(8),c(nsp)

   call mixture_energy(t,tneq,c,flg_anha,mode,flg_neq,flg_stop,flg_termo,mixe)

end subroutine c_mixture_energy




subroutine c_cp_frozen(t,tneq,ctmp,flg_anha,mode,flg_neq,flg_stop,flg_termo,cp_m)

   use interf_thermo
   use global_thermo, only: nsp
   implicit none
   integer :: flg_anha,mode,flg_neq,flg_stop,flg_termo
   real(kind=8) :: t,tneq(4),cp_m(8),ctmp(nsp)

   call cp_frozen(t,tneq,ctmp,flg_anha,mode,flg_neq,flg_stop,flg_termo,cp_m)

end subroutine c_cp_frozen




subroutine c_massfrac(Y,p,T,Tneq,eps,flg_log,flg_anha,flg_neq,flg_stop,flg_termo,   &
&                     flg_guess,x_guess,niter)

   use interf_thermo
   use global_thermo, only: nsp
   implicit none
   integer :: flg_log,flg_anha,flg_neq,flg_stop,flg_termo,flg_guess,niter
   real(kind=8) :: Y(nsp),p,T,Tneq(4),eps,x_guess(nsp)

   real(kind=8) :: X(nsp),zo,mm

   call molefrac(X,p,T,zo,Tneq,eps,flg_log,flg_anha,flg_neq,flg_stop,flg_termo,   &
&                flg_guess,x_guess,niter)
   !flg_guess=0!pietro
   mm = mixture_molarmass(X,1,flg_stop)

   call molefrac_to_massfrac(X,mm,Y)


end subroutine c_massfrac




real(kind=8) function c_mixture_rspecific(Y,mode,flg_stop)
   use interf_thermo
   use global_thermo, only: nsp
   implicit none
   integer :: mode,flg_stop
   real(kind=8) :: Y(nsp)

   c_mixture_rspecific = mixture_rspecific(Y,mode,flg_stop)

end function c_mixture_rspecific



!#####################################################################
!             Added by NANNI 
!#####################################################################
!

!subroutine c_thermal_cond_tot_nanni(k_t,tneq,h_react)
!
!       real(kind=8) :: k_t,tneq(4),h_react(2)
!
!        call thermal_cond_tot_nanni(k_t,tneq,h_react)
!	
!end subroutine c_thermal_cond_tot_nanni
!



subroutine c_massfrac_eq(Y,p,T,Tneq,eps,flg_log,flg_anha,flg_neq,flg_stop,flg_termo,   &
&                     flg_guess,x_guess,niter,Xc)

   use interf_thermo
   use global_thermo, only: nsp
   implicit none
   integer :: flg_log,flg_anha,flg_neq,flg_stop,flg_termo,flg_guess,niter
   real(kind=8) :: Y(nsp),p,T,Tneq(4),eps,x_guess(nsp),Xc(2)

   real(kind=8) :: X(nsp),zo,mm

   !write(*,*)'Inside c_massfrac_eq', p, T 
   
   flg_guess=0
   
   call molefrac_mod(X,p,T,Tneq,eps,flg_log,flg_anha,flg_neq,flg_stop,flg_termo,   &
&                flg_guess,x_guess,niter,Xc)

 ! pause
  
   mm = mixture_molarmass(X,1,flg_stop)

   call molefrac_to_massfrac(X,mm,Y)

  
end subroutine c_massfrac_eq



!------------------------------------------------------------------------------

! This subroutine is like the previous one except that here the added variable
! zo is taken out as an index of the numerical errors committed in evaluating
! the Xc quantities.







subroutine c_massfrac_eq_err(Y,p,T,Tneq,eps,flg_log,flg_anha,flg_neq,flg_stop,flg_termo,   &
&                     flg_guess,x_guess,niter,Xc,var_add_now)

   use interf_thermo
   use global_thermo, only: nsp
   implicit none
   integer :: flg_log,flg_anha,flg_neq,flg_stop,flg_termo,flg_guess,niter
   real(kind=8) :: Y(nsp),p,T,Tneq(4),eps,x_guess(nsp),Xc(2)

   real(kind=8) :: X(nsp),mm,zo,var_add_now

   write(*,*)'Inside c_massfrac_eq', Xc
   
   flg_guess=0
   
   call molefrac_mod_constr(X,p,T,zo,Tneq,eps,flg_log,flg_anha,flg_neq,flg_stop,flg_termo,   &
&                flg_guess,x_guess,niter,Xc)

   var_add_now=zo
   
   write(*,*)'zo = ',var_add_now
   
   !pause
   
  
   mm = mixture_molarmass(X,1,flg_stop)

   call molefrac_to_massfrac(X,mm,Y)

  
end subroutine c_massfrac_eq_err


!#######################################################################



!real(kind=8) function c_mix_den_mass(p,T,y,flg_stop)

!   use global_thermo, only: nsp
!   use interf_thermo
!   implicit none
!   integer :: flg_stop
!   real(kind=8) :: p,T,y(1:nsp)

!   c_mix_den_mass = mix_den_mass(p,T,y,flg_stop)
!   c_mix_den_mass = mixture_density(p,T,y,flg_stop)


!end function c_mix_den_mass
