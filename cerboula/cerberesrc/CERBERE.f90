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

call INPUT_READING

mixname2 = mixname

call TERMODEF(thermo_dir,traco_dir,chem_dir,mixname2,flg_anha,flg_oper,&
   &flg_termo,flg_traco,flg_stop,maxlvl)
call TRACODEF(sonine)

count_iter = 0     !To count the total number of iteration
carnaval = 0       !To know if we are in the 1st iteration of the main loop (LHTS loop)

error_Kp = 1000

call SET_NEBOULA_INPUT

call TERMOSTOP
call TRACOSTOP

end program CERBERE
