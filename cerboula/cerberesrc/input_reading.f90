subroutine INPUT_READING

USE PEG_DIR_par
USE FLAGS_par
USE REBUILDING_var
USE CONVERGENCE_par
USE EXPERIMENTS_par
USE BL_par
USE ICP_par
USE OPTIONS_par
USE HEATFLUX_var
USE OUTEREDGE_var
USE NAMES_par
USE ABACUS_par
USE INIT_DIST_par
USE CONVERGENCE_ENTH_par

implicit none

integer i
character*80 dummy,mixfile
real*8 anspecies,expon_gamma_min,expon_gamma_max,dummy2



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	READING OF THE FILE PEG_DIR.DAT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

OPEN(20,FILE='neboulainput/peg_dir.in',STATUS='OLD',ACCESS='SEQUENTIAL',&
& ACTION='READ')
  read(20,*) maxlvl
  read(20,*) sonine
  read(20,'(a)') thermo_dir
  read(20,'(a)') traco_dir
  read(20,'(a)') chem_dir
  read(20,'(a)') mixname
close(20)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	READING OF THE FLAGS FOR PEGASE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

OPEN(20,FILE='cerbereinput/flags.in',STATUS='OLD',ACCESS='SEQUENTIAL',&
& ACTION='READ')
  read(20,*) flg_anha 
  read(20,*) flg_oper 
  read(20,*) flg_termo 
  read(20,*) flg_traco 
  read(20,*) flg_stop 
  read(20,*) flg_mode
close(20)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	READING OF THE BL PARAMETERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

OPEN(20,FILE='neboulainput/bl_mode.in',STATUS='OLD',ACCESS='SEQUENTIAL',&
& ACTION='READ')
  read(20,*) dummy
  read(20,*) dummy
  read(20,*) dummy
  read(20,*) dummy
  read(20,*) dummy
  read(20,*) dummy
  read(20,*) dummy
  read(20,*) npoints
  read(20,*) dummy
  read(20,*) dummy
  read(20,*) dummy
  read(20,*) dummy
  read(20,*) dummy
  read(20,*) nreacwall
CLOSE(20)

mixfile = trim(thermo_dir)//trim(mixname)

OPEN(20,FILE=mixfile,STATUS='OLD',ACCESS='SEQUENTIAL',ACTION='READ')
  read(20,*) dummy
  read(20,*) anspecies
CLOSE(20)

nspecies = int(anspecies)
allocate(concent_rebuild(1:nspecies))
allocate(concent_ext(1:nspecies))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	READING OF THE CONVERGENCE PARAMETERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
OPEN(20,FILE='cerbereinput/convergence.in',STATUS='OLD',ACCESS='SEQUENTIAL',&
& ACTION='READ')
  read(20,*) epsilon
  read(20,*) alpha
  read(20,*) errorstop
  read(20,*) stop_Kp
CLOSE(20)


   OPEN(24,FILE='cerbereinput/input.in',STATUS='OLD',ACCESS='SEQUENTIAL',&
   & ACTION='READ')
      read(24,*) Kp
      read(24,*) he
      read(24,*) Pstatic
      read(24,*) Vs_ext
      read(24,*) Rm,Tw
      read(24,*) delta_1,u1e_1,u1y_1,ve_1,NDP_5_1
      read(24,*) enthalpyerror
      read(24,*) enthepsilon
      read(24,*) enthalpha
      read(24,*) T_0
      read(24,*) Tw_min,Tw_max
   close(24)
   
end subroutine INPUT_READING
