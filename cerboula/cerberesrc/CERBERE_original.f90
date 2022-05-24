program CERBERE

USE OPTIONS_par
USE PEG_DIR_par
USE FLAGS_par
USE REBUILDING_var
USE CONVERGENCE_par
USE EXPERIMENTS_par
USE OUTEREDGE_var

implicit none

integer carnaval,i
character*80 mixname2
real*8 gam
integer flag_OK

call INTERF
call INPUT_READING

mixname2 = mixname

call TERMODEF(thermo_dir,traco_dir,chem_dir,mixname2,flg_anha,flg_oper,&
   &flg_termo,flg_traco,flg_stop,maxlvl)
call TRACODEF(sonine)

count_iter = 0     !To count the total number of iteration
carnaval = 0       !To know if we are in the 1st iteration of the main loop (LHTS loop)

error_Kp = 1000

if (choice.eq.0) then
   if (Kp_choice.eq.0) then
      call ENTH_REBUILDING(gamma)
   elseif (Kp_choice.eq.1) then
      do while (error_Kp.gt.stop_Kp)
        Kp_old = Kp
        call KP_REBUILDING
        call ENTH_REBUILDING(gamma)
        if (carnaval.gt.0) then
          error_Kp = abs(Kp-Kp_old)
          write(*,*) 'BARKER ERROR', error_Kp
        endif
        carnaval = carnaval+1
      enddo
   endif
   call WRITE_SOL
elseif (choice.eq.1) then
   if (Kp_choice.eq.0) then
      do i = 1,N_gamma
         call ENTH_REBUILDING(gamma_vector(i))
         call DATA_STORAGE(i,gamma_vector(i),T_reb,h_ext,Qwc_reb,Qwd_reb,rho_ext,mu_ext,Vs_ext,Kp)
      enddo
   elseif (Kp_choice.eq.1) then
     do i = 1,N_gamma
         do while (error_Kp.gt.stop_Kp)
           Kp_old = Kp
           call KP_REBUILDING
           call ENTH_REBUILDING(gamma_vector(i))
           if (carnaval.gt.1) then
             error_Kp = abs(Kp-Kp_old)
             write(*,*) 'BARKER ERROR', error_Kp
           endif
           carnaval = carnaval + 1
         enddo
         call DATA_STORAGE(i,gamma_vector(i),T_reb,h_ext,Qwc_reb,Qwd_reb,rho_ext,mu_ext,Vs_ext,Kp)
      enddo
   endif
elseif (choice.eq.2) then
   if (Kp_choice.eq.0) then
        call ENTH_REBUILDING(gamma)
   elseif (Kp_choice.eq.1) then
      do while(error_Kp.gt.stop_Kp)
        Kp_old = Kp
        call KP_REBUILDING
   
        call ENTH_REBUILDING(gamma)
        if (carnaval.gt.0) then
          error_Kp = abs(Kp-Kp_old)
          write(*,*) 'BARKER ERROR', error_Kp
        endif
        carnaval = carnaval+1
      enddo
   endif
   call CAT_REBUILDING(flag_OK)
   if (flag_OK.eq.1) then
      call WRITE_SOL
   endif
elseif (choice.eq.3) then
   call SET_NEBOULA_INPUT
   call ABACUS
elseif (choice.eq.4) then
   call INIT_DIST
elseif (choice.eq.5) then
   call RICHARDSON_EXTRAPOLATION
elseif (choice.eq.6) then
   call DWDRHO
endif

write(*,*) 'END OF CERBERE'

call TERMOSTOP
call TRACOSTOP

end program CERBERE
