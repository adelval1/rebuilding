! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  MAIN    L I B R A R Y              //
! //                                                                  //
! //  MODULE global_pegase                                            //
! //                                                                  //
! // Definition of the global variables used in the PEGASE program    //
! // and which are not to be reused with the library only.            //
! //                                                                  //
! //  Benoit Bottin, 30/09/97                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
MODULE global_pegase
!
!   Old common block /custom_file/ for the custom output file
!   =========================================================
!   custmax     number of output items requested
!   custom      the output array of items requested:
!                   1:99 the number of the output item
!                   1:4  the reference of the output item
!                           1   the output array number
!                           2   first index in the output array
!                           3   second index in the output array
!                           4   third index in the output array
!   label       the title of the column for the output item
!   ------------------------------------------------------------
    INTEGER,SAVE :: custmax
    REAL(kind=8), SAVE :: custom(99,4)
    CHARACTER*14, SAVE :: label(99)
!
!   Old common block /filesopen/ for which files are open
!   =====================================================
!   fop_xxx     valued 1 if the file is open, 0 if closed
!                     -1 if the file had been opened in the run but
!                      has been closed in the meantime
!   ---------------------------------------------------------------    
    INTEGER,SAVE :: fop_out,fop_dat,fop_da2,fop_cust,fop_log
!
!   Old common block /ostype/ defining the OS system
!   ================================================
!   flg_os      1 if PC, 0 if UNIX
!   ------------------------------
    INTEGER,SAVE :: flg_os
!
!   Pegase program storage variable arrays
!   ======================================
!   array1      dimension 10,nsp
!   array2      dimension 16,nsp,8
!   array3      dimension 22
!   array4      dimension 18,8
!   array5      dimension nr,4
!   array6      dimension 6,nsp (nsp min 7 in this case)
!   userarray   dimension nusr
!   nusr        maximum number of user-defined variables in user array
!   ------------------------------
    REAL(kind=8),ALLOCATABLE,SAVE :: array1(:,:)
    REAL(kind=8),ALLOCATABLE,SAVE :: array2(:,:,:)
    REAL(kind=8),SAVE             :: array3(24)
    REAL(kind=8),SAVE             :: array4(18,8)
    REAL(kind=8),ALLOCATABLE,SAVE :: array5(:,:)
    REAL(kind=8),ALLOCATABLE,SAVE :: array6(:,:)
    REAL(kind=8),ALLOCATABLE,SAVE :: userarray(:)
    INTEGER,SAVE :: nusr
!
!   Pegase nonequilibrium/frozen input composition
!   ==============================================
!   X_NEQ       allocated in the input file read sequence
!   -----------------------------------------------------
    REAL(kind=8),ALLOCATABLE,SAVE :: X_NEQ(:)
!
END MODULE global_pegase

MODULE input_pegase
!
!   Pegase directory structure returned from the control file
!   =========================================================
!   in_dir      directory of input files
!   out_dir     directory of output files
!   thermo_dir  directory of thermodynamics data files
!   traco_dir   directory of transport coefficients data files
!   chemco_dir  directory of chemical kinetics data files
!   ----------------------------------------------------------
      CHARACTER*80 in_dir,out_dir,thermo_dir,chemco_dir,traco_dir
!
!   Variables defined in the input file
!   ===================================
!     a) strings:
!         mixname     name if the mixture file (including path)
!         casename    name of the test case being considered
!         formatname  name of the custom format description file
!     b) flags (integers):
!         flg_neq     set to 1 = thermal nonequilibrium operation
!         flg_oper    set to 0 = single case computation
!                     set to 1 = equilibrium calculation
!                     set to 2 = frozen flow calculation
!         flg_scan    set to 0 = single input point
!                     set to 1 = do loop on the second input
!                     set to 2 = do loops on both inputs
!         flg_mode    set to 1 = pressure-temperature input
!                     set to 2 = density-temperature input
!         flg_anha    set to 0 = no anharmonicity correction
!                     set to 1 = simple correction (constant parameters)
!                     set to 2 = full correction (variable parameters)
!         flg_out     set to 1 = output file is required
!         flg_dat     set to 1 = tabulated output file is required
!         flg_da2     set to 1 = second tabulated output file is required
!         flg_cust    set to 1 = custom output file is required
!         flg_tec     set to 1 = output of .dat is in TECPLOT format
!         flg_termo   set to 0 = use of statistical mechanics
!                     set to 1 = use of reference tables (Gurvich-style)
!                     set to 2 = use Srinivasan curve-fits for air
!         flg_traco   set to 0 = use of kinetic theory
!                     set to 1 = use of Gupta curve-fits for equilibrium air
!         flg_stop    set to 1 = stops programs in case of errors in routines
!         flg_more    set to 1 = another computation case follows this one
!         flg_user    set to 1 = use the user_compute routines linked to code
!     c) input values (real*8):
!         var1_from   value of var1 input or lower bound if scanning is on
!         var1_to     upper bound of var1 scan if flg_scan = 2
!         var1_step   step for scanning var1 if flg_scan = 2
!         var2_from   value of var2 input or lower bound if scanning is on
!         var2_to     upper bound on var2 if flg_scan = 1 or 2
!         var2_step   step for scanning var2 if flg_scan = 1 or 2
!         epsilon     finite differences tolerance value
!     d) nonequilibrium and frozen flow arrays (REAL(kind=8)):
!         t_neq(4)    trans, rot, vib and elec temperatures
!   --------------------------------------------------------------------------
      CHARACTER*80 mixname,casename,formatname
      INTEGER flg_neq,flg_oper,flg_scan,flg_mode,flg_anha,flg_out
      INTEGER flg_dat,flg_da2,flg_cust,flg_termo,flg_traco,flg_stop
      INTEGER flg_more,flg_user,flg_log,flg_tec
      REAL(kind=8) var1_from,var1_to,var1_step,var2_from,var2_to,var2_step
      REAL(kind=8) epsilon,findiff,t_neq(4)
!
END MODULE input_pegase
