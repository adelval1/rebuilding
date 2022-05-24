! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO    L I B R A R Y             //
! //                                                                  //
! //  MODULE global_thermo                                            //
! //                                                                  //
! // Definition of the global variables used in the PEGASE library    //
! // and which can be used either in pegase or in other codes.        //
! //                                                                  //
! // Benoit Bottin, 30/09/97                                          //
! //                                                                  //  
! // The code is modified by taking into account multi vibrational    //
! // modes.This is done by converting the 2-D array TV to a 3-D array //
! // the third dimension being the number of vibrational mode         //
! //                                                                  //
! // Alp Marangoz,  26/06/2001                                        //
! //                                                                  // 
! //////////////////////////////////////////////////////////////////////
MODULE global_thermo
    IMPLICIT NONE
!
!   Old parameters for the sizing of the arrays (old size.inc)
!   ==========================================================
!   maxlvl      maximum number of electronic levels
!   maxsp       maximum number of species + 1
!   maxc        maximum number of basic nuclei + 1
!   maxr        maximum number of chemical reactions
!   maxdata     maximum number of data points for thermodynamic data
!   ----------------------------------------------------------------
    INTEGER,SAVE :: nlvl,nsp,nc,nr,ndata,mxnvm
    INTEGER,ALLOCATABLE,SAVE :: nvibmode(:)                 

!   Old common block /ini_composition/ for the initial composition
!   ==============================================================
!   xini        array of initial compositions
!   xcons       array of RHS terms for the system in terms of mole
!   ycons       array of RHS terms for the system in terms of mass
!   The 3 arrays have maximum size maxc
!   --------------------------------------------------------------
    REAL(kind=8),ALLOCATABLE,SAVE :: XINI(:),XCONS(:),YCONS(:)
!
!   Old common block /thermo_properties/ for thermodynamic data
!   ===========================================================
!   TR          characteristic rotation temperature     (nlvl,nsp)
!   TV          characteristic vibration temperature    (nvibmode,nlvl,nsp)
!   TE          characteristic electronic temperatures  (nlvl,nsp)
!   GE          characteristic electronic degeneracies  (nlvl,nsp)
!   HFOR        formation enthalpy of the species       (nsp)
!   ST          Sackur-Tetrode constant                 (nsp)
!   SIG         symmetry factor                         (nsp)
!   NELVL       actual number of electronic levels used (nsp)
!   --------------------------------------------------------------
    REAL(kind=8),ALLOCATABLE,SAVE :: TR(:,:),TV(:,:,:),TE(:,:),GE(:,:)
    REAL(kind=8),ALLOCATABLE,SAVE :: HFOR(:),ST(:),SIG(:)
    INTEGER,ALLOCATABLE,SAVE :: NELVL(:)
!
!   Old common block /mass_properties/ for thermodynamic data
!   =========================================================
!   RSP         specific gas constant of the species    (nsp)
!   MMOL        molar mass of the species               (nsp)
!   --------------------------------------------------------------
    REAL(kind=8),ALLOCATABLE,SAVE :: RSP(:),MMOL(:)                  
!
!   Old common block /ion_charge/ for thermodynamic data
!   =========================================================
!   CHARGE      charge of the species                   (nsp)
!   --------------------------------------------------------------
    INTEGER,ALLOCATABLE,SAVE :: CHARGE(:)
!
!   Old common block /reference_table/ for thermodynamic data
!   =========================================================
!   TABLE_NAME  name of the reference table file name   (nsp)
!   --------------------------------------------------------------
    CHARACTER*80,ALLOCATABLE,SAVE :: TABLE_NAME(:)
!
!   Old common block /spectro_data/ for anharmonicity corrections
!   =============================================================
!   A_E         anharmonicity constant alpha            (nlvl,nsp)
!   B_E         rotation constant B                     (nlvl,nsp)
!   D_str       anharmonicity constant D                (nlvl,nsp)
!   D_E         dissociation/ionization energy          (nlvl,nsp)
!   W_E         natural vibration frequency omega       (nlvl,nsp)
!   XW_E        product of omega and anh. ct. x         (nlvl,nsp)
!   --------------------------------------------------------------
    REAL(kind=8),ALLOCATABLE,SAVE :: A_E(:,:),B_E(:,:),D_E(:,:)
    REAL(kind=8),ALLOCATABLE,SAVE :: D_str(:,:),W_E(:,:),XW_E(:,:)
!
!   Old common block /anharmonicity/ for anharmonicity corrections
!   ==============================================================
!   D_A         useful anharmonicity constant D         (nlvl,nsp)
!   G_A         useful anharmonicity constant gamma     (nlvl,nsp)
!   X_A         useful anharmonicity constant x         (nlvl,nsp)
!   --------------------------------------------------------------
    REAL(kind=8),ALLOCATABLE,SAVE :: D_A(:,:),G_A(:,:),X_A(:,:)
!
!   Old common block /stoechio/ for the stoechiometric matrix
!   ==============================================================
!   NU          stoechiometric matrix                   (nr,nsp+1)
!   --------------------------------------------------------------
    INTEGER,ALLOCATABLE,SAVE :: NU(:,:)
!
!   Old common block /conserv/ for the conservation of nuclei matrix
!   ================================================================
!   KHI         conservation matrix                     (nc+1,nsp+1)
!   ----------------------------------------------------------------
    REAL(kind=8),ALLOCATABLE,SAVE :: KHI(:,:)
!
!   Old common block /species_name/ for species names
!   =================================================
!   spc_label   label tag string for each species
!   spc_symbol  species symbols in order to read the reactions
!   ---------------------------------------------    
    CHARACTER*20,ALLOCATABLE,SAVE :: SPC_LABEL(:)
    CHARACTER*20,ALLOCATABLE,SAVE :: SPC_SYMBOL(:)
!
!   Old common block /reference_table/ for thermodynamic input
!   ==========================================================
!   reftable    reference table data storage - dim nsp,4,ndata
!                   1st index: species number
!                   2nd index:  1 = temperature
!                               2 = sensible enthalpy
!                               3 = entropy
!                               4 = cp specific heat
!                   3rd index: data items
!   refpress    reference pressure for the table
!   idata       the actual number of data points (actual range
!                   of the 3rd index of reftable)
!   ----------------------------------------------------------
    REAL(kind=8),ALLOCATABLE,SAVE :: REFTABLE (:,:,:), REFPRESS (:)
    INTEGER,ALLOCATABLE,SAVE :: IDATA (:)
!
!   Old common block /gupta/ for gupta transport curve fits
!   =======================================================
!   VMAT        viscosity matrix
!   KMAT        thermal conductivity matrix
!   PRMAT       Prandtl number matrix
!   CPMAT       Cp specific heat matrix
!   -------------------------------------------------------
      REAL(kind=8),SAVE :: VMAT(37,8),KMAT(37,8),PRMAT(37,8),CPMAT(37,8)
!
!   Physical parameters
!   ================================================================
!   RUNIV       universal gas constant
!   PIUNIV      trigonometric constant
!   NAUNIV      Avogadro's number
!   KUNIV       Boltzmann's constant
!   HUNIV       Planck's constant
!   ----------------------------------------------------------------
    REAL(kind=8),PARAMETER :: RUNIV=8.31451d0             ! J/(mol.K)
    REAL(kind=8),PARAMETER :: PIUNIV=3.141592653589d0     ! [1]
    REAL(kind=8),PARAMETER :: NAUNIV=6.0221367d23         ! part/mol
    REAL(kind=8),PARAMETER :: KUNIV=1.380658d-23          ! J/K
    REAL(kind=8),PARAMETER :: HUNIV=6.6260755d-34         ! J.s
    REAL(kind=8),PARAMETER :: EUNIV=1.60217733d-19        ! C
    REAL(kind=8),PARAMETER :: EPSILONUNIV=8.854187817d-12 ! F/m
!
!   Transport properties collision integrals file name
!   ==================================================
    CHARACTER*80,SAVE :: traconame
!
!   Chemical kinetics data file
!   ===========================
    CHARACTER*80,SAVE :: chemconame

END MODULE
