subroutine diff_coeff_f(Dij_vec,tneq,p,nspec)

   use interf_traco
   implicit none
   integer :: nspec
   real(kind=8) :: Dij_vec(1:nspec*nspec),tneq(1:4),p

   integer :: i,j
   real(kind=8) :: Dij(nspec,nspec)

   Dij = 0

   call diffusion_coefficients_binary(p,tneq,Dij)

   do i=1,nspec
      do j=1,nspec
         Dij_vec((i-1)*nspec+j) = Dij(i,j)
      end do
   end do

end subroutine diff_coeff_f


subroutine set_1(flg_anha,flg_oper,flg_termo,flg_traco,flg_stop,sonine,   &
       &         R,pi,kb,Na,nspec)

   use interf_thermo
   use interf_traco
   use global_thermo, only: nsp,RUNIV,PIUNIV,NAUNIV,KUNIV
   implicit none
   real(kind=8) :: R,pi,kb,Na
   integer :: flg_anha,flg_oper,flg_termo,flg_traco,flg_stop,sonine,nspec
   integer :: maxlvl

   character*80 thermo_dir,traco_dir,mixname,chem_dir
   integer :: flg_oper2


   flg_oper2=1

   open(unit=20,file='neboulainput/peg_dir.in',status='unknown',err=100)
      read(20,*) maxlvl
      read(20,*) sonine
      read(20,'(a)') thermo_dir
      read(20,'(a)') traco_dir
      read(20,'(a)') chem_dir
      read(20,'(a)') mixname
   close(20)

   call termodef(thermo_dir,traco_dir,chem_dir,mixname,flg_anha,flg_oper2,flg_termo, &
        &        flg_traco,flg_stop,maxlvl)
   call tracodef(sonine)
   call chemcodef

   nspec=nsp
   R=RUNIV
   pi=PIUNIV
   kb=KUNIV
   Na=NAUNIV

   return
100 write(*,*) 'Error in opening file blinput/peg_dir.in'
   stop
end subroutine set_1


subroutine set_2(Mw,Mw_scaled,c_charge,chflag)

   use global_thermo, only: nsp,MMOL,CHARGE,NAUNIV,EUNIV
   implicit none
   integer :: chflag,i
   real(kind=8) :: Mw(1:nsp),Mw_scaled(1:nsp),c_charge(1:nsp)
   real(kind=8) :: maxi

   chflag=0
   do i=1,nsp
      Mw(i)=MMOL(i)
      if(CHARGE(i).ne.0) chflag=1
   end do

!  Computes the vector of charges (C/kg)
   c_charge=0.0d0
   if(chflag.eq.1) then
      do i=1,nsp
         c_charge(i) = EUNIV/(MMOL(i)/NAUNIV)*dble(CHARGE(i))
      end do
   end if

! Computes scaled molecular weight

   maxi = 0.0d0
   do i=1,nsp
      if(MMOL(i) >= maxi) maxi=MMOL(i)
   end do

   do i=1,nsp
      Mw_scaled(i) = MMOL(i)/maxi
   end do

end subroutine set_2



subroutine f_chemistry(T,Tneq,rhosp,W,dWdrp,dWdrm,dWdT,flg_anha,flg_neq,flg_stop,flg_termo)

   use global_thermo, only:nsp
   use interf_chemco
   implicit none
   integer :: flg_anha,flg_neq,flg_stop,flg_termo
   real(kind=8) :: T,Tneq(4),rhosp(1:nsp),W(1:nsp)
   real(kind=8) :: dWdrp(1:nsp*nsp),dWdrm(1:nsp*nsp),dWdT(1:nsp)

   integer :: i,j
   real(kind=8) :: rho_jacp(nsp,nsp),rho_jacm(nsp,nsp),T_jac(nsp)


   call prod_term_den(T,rhosp,W,Tneq,flg_anha,flg_neq,flg_stop,flg_termo)

   call jac_term_den(T,rhosp,rho_jacp,rho_jacm,T_jac,Tneq,flg_anha,flg_neq,flg_stop,flg_termo)

   do i=1,nsp
      dWdT(i)=T_jac(i)
      do j=1,nsp
         dWdrp((i-1)*nsp+j) = rho_jacp(i,j)
         dWdrm((i-1)*nsp+j) = rho_jacm(i,j)
      end do
   end do

end subroutine f_chemistry


