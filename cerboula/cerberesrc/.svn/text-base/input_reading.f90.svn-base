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



if (choice.eq.0) then
   OPEN(21,FILE='cerbereinput/rebuilding.in',STATUS='OLD',ACCESS='SEQUENTIAL',&
   & ACTION='READ')
      read(21,*) Kp
      read(21,*) gamma
      read(21,*) Tw
      read(21,*) deltaPs_1
      read(21,*) Pstatic
      read(21,*) Rm,Qw_exp_1
      read(21,*) delta_1,u1e_1,u1y_1,ve_1,NDP_5_1
   close(21)
elseif (choice.eq.1) then
   OPEN(22,FILE='cerbereinput/enthalpy_S-curve.in',STATUS='OLD',ACCESS='SEQUENTIAL',&
   & ACTION='READ')
      read(22,*) expon_gamma_min
      read(22,*) expon_gamma_max
      read(22,*) Kp
      read(22,*) N_gamma
      allocate(gamma_vector(1:N_gamma))
      do i = 1,N_gamma
         gamma_vector(i) = expon_gamma_max - (expon_gamma_max-expon_gamma_min)/float(N_gamma-1)*dfloat(i-1)
         gamma_vector(i) = 10**gamma_vector(i)
      enddo
      read(22,*) Tw
      read(22,*) deltaPs_1
      read(22,*) Pstatic
      read(22,*) Rm,Qw_exp_1
      read(22,*) delta_1,u1e_1,u1y_1,ve_1,NDP_5_1
      read(22,*) probe_name1
   close(22)
elseif (choice.eq.2) then
   OPEN(23,FILE='cerbereinput/identification.in',STATUS='OLD',ACCESS='SEQUENTIAL',&
   & ACTION='READ')
      read(23,*) Kp
      read(23,*) gamma
      read(23,*) deltaPs_1
      read(23,*) Pstatic
      read(23,*) Rm,Qw_exp_1,Tw
      read(23,*) Rm_smpl,emissivity_smpl,Tw_smpl
      read(23,*) delta_1,u1e_1,u1y_1,ve_1,NDP_5_1
      read(23,*) delta_smpl,u1e_smpl,u1y_smpl,ve_smpl,NDP_5_smpl
   close(23)
elseif (choice.eq.3) then
   OPEN(24,FILE='cerbereinput/abacus.in',STATUS='OLD',ACCESS='SEQUENTIAL',&
   & ACTION='READ')
      read(24,*) Kp
      read(24,*) he
      read(24,*) Pstatic
      read(24,*) Vs_ext
      read(24,*) Rm,Qw_exp_1,Tw
      read(24,*) delta_1,u1e_1,u1y_1,ve_1,NDP_5_1
      read(24,*) enthalpyerror
      read(24,*) enthepsilon
      read(24,*) enthalpha
      read(24,*) T_0
      read(24,*) Tw_min,Tw_max
      read(24,*) ntw
      read(24,*) ngamma
      ngamma = int(ngamma)
      allocate(gammaarray(1:ngamma))
      do i=1,ngamma
         read(24,*) gammaarray(i)
      enddo
   close(24)
elseif (choice.eq.4) then
   OPEN(25,FILE='cerbereinput/initial_distribution.in',STATUS='OLD',ACCESS='SEQUENTIAL',&
   & ACTION='READ')
      read(25,*) Tw
      read(25,*) he
      read(25,*) Pstatic
      read(25,*) Vs_ext
      read(25,*) delta_1,u1e_1,u1y_1,ve_1,NDP_5_1
      read(25,*) Kp
      read(25,*) Rm
      read(25,*) enthalpyerror
      read(25,*) enthepsilon
      read(25,*) enthalpha
      read(25,*) T_0
      read(25,*) neq_fr
   close(25)
elseif (choice.eq.5) then
   OPEN(26,FILE='cerbereinput/equilibrium_extrapolation.in',STATUS='OLD',ACCESS='SEQUENTIAL',&
   & ACTION='READ')
      read(26,*) N_gamma
      read(26,*) probe_name1
      read(26,*) probe_name2
      read(26,*) probe_name3
      read(26,*) u1e_1,Rm
      read(26,*) u1e_2,Rm_2
      read(26,*) u1e_3,Rm_3
   close(26)
endif


end subroutine INPUT_READING
