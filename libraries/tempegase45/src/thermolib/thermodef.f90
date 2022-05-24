! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  SETUP   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE termodef                                             //
! //                                                                  //
! //  Input:  thermo_dir  directory of the mixture and species files  //    
! //          mixname     name of the mixture file                    //
! //                                                                  //
! //  Outs:   nsp         actual number of species used               // 
! //          nr          actual number of reaction equations used    //
! //          nc          actual number of conservation equations     //
! //          mixname     the name of the mixture file returns with   //
! //                      the full name including directory           //
! //                                                                  //
! //  Flags:  flg_stop    forces to stop on library error             //
! //          flg_termo   passed to read_species_file directly        //
! //          flg_oper    used to check if it is 1 (equ. calculation) //
! //          flg_anha    passed to read_species_file directly        //
! //                                                                  //
! //  The subroutine termodef is simply used to perform the retrieval //
! //  of all the relevant information by use of a single subroutine   //
! //  call. The mixture file is opened, read and closed and all       //
! //  mixture parameters copied in the suitable commons in the        //
! //  child subroutine read_mixture_file. During the execution of     //
! //  this routine, each species file found is opened, read and then  //
! //  closed. All the species parameters are in turn copied in the    //
! //  common blocks in the routine read_species_file.                 //
! //                                                                  //
! //  The only new parameters returned by the routine are the actual  //
! //  sizes of the problem: nsp, nr and nc, defining the number of    //
! //  species, reactions and conservation equations respectively.     //
! //                                                                  //
! //  Benoit Bottin, 08/10/96. Rewritten Win95 24/7/97.               //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////

      SUBROUTINE thermodef (thermo_dir,traco_dir,chemco_dir,mixname,flg_anha,&
     &flg_oper,flg_termo,flg_traco,flg_stop,maxlvl)

!     Dummy routine to enable users to call the more logical routine
!     `thermodef' instead of `termodef'.

      USE global_thermo
      USE global_equilibrium
      IMPLICIT NONE
!
!     Definition of subroutine input/output variables
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      INTEGER flg_anha,flg_oper,flg_termo,flg_traco,flg_stop,maxlvl
      INTEGER ierr,FILE_EXISTS,FILE_OPEN,FILE_CLOSE
      CHARACTER*80 thermo_dir,traco_dir,chemco_dir,mixname

      CALL termodef (thermo_dir,traco_dir,chemco_dir,mixname,flg_anha,&
     &flg_oper,flg_termo,flg_traco,flg_stop,maxlvl)

      END SUBROUTINE thermodef

      SUBROUTINE termodef (thermo_dir,traco_dir,chemco_dir,mixname,flg_anha,&
     &flg_oper,flg_termo,flg_traco,flg_stop,maxlvl)

      USE global_thermo
      USE global_equilibrium
      IMPLICIT NONE
!
!     Definition of subroutine input/output variables
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      INTEGER flg_anha,flg_oper,flg_termo,flg_traco,flg_stop,maxlvl
      INTEGER ierr,FILE_EXISTS,FILE_OPEN,FILE_CLOSE
      CHARACTER*80 thermo_dir,traco_dir,chemco_dir,mixname
!
!     Definition of variables intrinsic to the subroutine
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      CHARACTER*80 dummyname
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
!
!     Initialization
!     ^^^^^^^^^^^^^^
      pgname='pegaslib'
      localinfo(12:70)=mixname(1:59)
      CALL FILL_WITH_BLANKS (dummyname)
      CALL CONCATENATE (thermo_dir,mixname,dummyname,ierr)
      mixname=dummyname
      nlvl = maxlvl
!
!     Open the mixture file
!     ^^^^^^^^^^^^^^^^^^^^^
      ierr=FILE_EXISTS (mixname)
      IF (ierr.EQ.0) THEN
          CALL PRINT_ERROR (8,pgname,localinfo,1)
      ELSEIF (ierr.EQ.-1) THEN
          CALL PRINT_ERROR (9,pgname,localinfo,1)
      ENDIF
      ierr=FILE_OPEN (mixname,12,1)
      IF (ierr.NE.1) THEN
          CALL PRINT_ERROR (10,pgname,localinfo,1)
      ENDIF
!
!     Analyze the mixture file
!     ^^^^^^^^^^^^^^^^^^^^^^^^
      CALL READ_MIXTURE_FILE (mixname,thermo_dir,traco_dir,chemco_dir,&
     &flg_anha,flg_oper,flg_termo,flg_traco,flg_stop)
!
!     Close the mixture file
!     ^^^^^^^^^^^^^^^^^^^^^^
      ierr=FILE_CLOSE (12)
      IF (ierr.EQ.-1) THEN
          CALL PRINT_ERROR (13,pgname,localinfo,1)
      ENDIF
!
!     If an equilibrium computation is required, compute fixedjac
!     and allocate all arrays required by equilibrium computations
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (flg_oper.EQ.1) THEN
      ALLOCATE (jac11(1:nr,1:nr),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of jac11            '
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE thermodynamic library error'
      END IF
      ALLOCATE(jac12(1:nr,1:nc+1),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of jac12            '
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE thermodynamic library error'
      END IF
          ALLOCATE (jaclow(1:nc+1,1:nsp+1),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of jaclow           '
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE thermodynamic library error'
      END IF
      ALLOCATE(jacmn(1:nc+1,1:nc+1),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of jacmn            '
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE thermodynamic library error'
      END IF
          ALLOCATE (f1(1:nr),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of f1               '
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE thermodynamic library error'
      END IF
      ALLOCATE(f2(1:nc+1),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of f2               '
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE thermodynamic library error'
      END IF
      ALLOCATE(delv(1:nsp+1),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of delv             '
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE thermodynamic library error'
      END IF
      ALLOCATE (indmn(1:nc+1),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of indmn            '
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE thermodynamic library error'
      END IF
      ALLOCATE(indsd(1:nr),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of indsd            '
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE thermodynamic library error'
      END IF
      ALLOCATE(luindx(1:nc+1),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of luindx           '
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE thermodynamic library error'
      END IF
      ALLOCATE (v(1:nsp+1),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of v                '
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE thermodynamic library error'
      END IF
      ALLOCATE (truevar(1:nsp+1),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of truevar          '
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE thermodynamic library error'
      END IF
      ALLOCATE(sys_keq(1:nr),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of sys_keq          '
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE thermodynamic library error'
      END IF
          jac11=0.0;jac12=0.0
          jaclow=0.0;jacmn=0.0
          f1=0.0;f2=0.0;indsd=0.0;indmn=0.0
          v=0.0;delv=0.0;sys_keq=0.0
          luindx=0;truevar=0.0
          CALL FIXEDJAC
      END IF
!     ------------------------------------------------------------------
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  SETUP   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE read_mixture_file                                    //
! //                                                                  //
! //  Input:  thermo_dir  directory of the mixture and species files  //    
! //          traco_dir   directory of the collision integrals files  //
! //                                                                  //
! //  Outs:   nsp         actual number of species used               // 
! //          nr          actual number of reaction equations used    //
! //          nc          actual number of conservation equations     //
! //          traconame   full name of collision integrals file       //
! //          chemconame   full name of chemical kinetics file         //            
! //                                                                  //
! //  Flags:  flg_stop    forces to stop on library error             //
! //          flg_termo   passed to read_species_file directly        //
! //          flg_traco   passed for consistency check purposes       //
! //          flg_oper    used to check if it is 1 (equ. calculation) //
! //          flg_anha    passed to read_species_file directly        //
! //                                                                  //
! //  Benoit Bottin, 21/10/96. Rewritten Win95 24/7/97.               //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////

      SUBROUTINE read_mixture_file (mixname,thermo_dir,traco_dir,chemco_dir,&
     &flg_anha,flg_oper,flg_termo,flg_traco,flg_stop)
      USE global_thermo
      IMPLICIT NONE
!     ------------------------------------------------------------------
!     PART 1: DEFINITION OF VARIABLES
!     ------------------------------------------------------------------
!
!     Definition of subroutine input/output variables
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      
      CHARACTER*80 thermo_dir,traco_dir, chemco_dir
      INTEGER flg_anha,flg_oper,flg_termo,flg_traco,flg_stop
!
!     Definition of variables intrinsic to the subroutine
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      
      INTEGER i,j,istop,ok,isp,ierr,pinar
      INTEGER,ALLOCATABLE :: K(:,:)
      REAL(kind=8) m
      REAL(kind=8) val
      CHARACTER*4 com
      CHARACTER*80 filename,spcname,fullcom,mixname
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
!     ------------------------------------------------------------------
!     END PART 1
!
!     PART 2: READING THE MIXTURE FILE 
!     ------------------------------------------------------------------
!     Initialization
!     ^^^^^^^^^^^^^^
   
      pgname='pegaslib'
      CALL fill_with_blanks (spcname)
      CALL fill_with_blanks (filename)
      CALL fill_with_blanks (traconame)
      CALL fill_with_blanks (chemconame)

      com='----'
      istop=0
      nsp=0
      nr=0
      nc=0
      m=0.0d0
      isp=0
      mxnvm=0
!pietro      write(*,*) 'reading mixture file'
! 
!     Interpretation of data file
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^
! Additonal loop for reading the max no of vib modes
      DO WHILE (com.NE.'stop')  
        CALL FILL_WITH_BLANKS(fullcom)
           com='    '
          READ (12,1000,err=901) fullcom
          CALL LCASE(fullcom(1:4),com,ierr)
          IF (com.EQ.'mnvm') THEN
             READ (12,1001,err=901) val
             mxnvm=INT(val)   
          ENDIF 
          localinfo(12:70)=mixname(1:59)
          CALL FILL_WITH_BLANKS (localinfo)
       END DO
       REWIND(12) 
       com='----'
!       
      DO WHILE (com.NE.'stop')
          CALL FILL_WITH_BLANKS(fullcom)
          com='    '
          READ (12,1000,err=901) fullcom
          CALL LCASE(fullcom(1:4),com,ierr)
          IF (com.EQ.'trac') THEN
!pietro	  write(*,*) 'reading trac'
              READ (12,1000,err=901) fullcom
              i=INDEX(fullcom,' ')
              IF ((i.EQ.0).AND.(flg_traco.NE.0)) THEN
                  localinfo(12:70)=mixname(1:59)
                  CALL PRINT_ERROR (2,pgname,localinfo,0)
                  flg_traco=0
              ELSE
                  CALL CONCATENATE (traco_dir,fullcom,traconame,ierr)
              ENDIF
!pietro	      write(*,*) 'no problem traco' 
          ELSEIF (com.EQ.'chem') THEN
!pietro          write(*,*) 'reading chem'
	      READ (12,1000,err=901) fullcom
              i=INDEX(fullcom,' ')
              IF (i.EQ.0) THEN
                  localinfo(12:70)=mixname(1:59)
                  CALL PRINT_ERROR (54,pgname,localinfo,0)
              ELSE
                  CALL CONCATENATE (chemco_dir,fullcom,chemconame,ierr)
              ENDIF
!pietro	      write(*,*) 'no problem chem'
          ELSEIF (com.EQ.'nsp ') THEN
!pietro	  write(*,*) 'reading nsp'
              READ (12,1001,err=901) val
              nsp=INT(val)
              CALL SIZE_ARRAYS (1)
              DO i=1,nsp
                  isp=i
                  READ (12,1000,err=901) spcname
                  j = LEN_TRIM(spcname)
                  SPC_SYMBOL(i)(1:j) = spcname(1:j)
                  CALL LCASE(spcname,filename,ierr)
                  filename(j+1:j+4) = '.spc'
                  CALL CONCATENATE (thermo_dir,filename,spcname,ierr)
                  CALL READ_SPECIES_FILE (isp,flg_anha,flg_termo,&
                 &flg_stop,thermo_dir,spcname)
                  CALL fill_with_blanks (filename)
                  CALL fill_with_blanks (spcname)
              END DO
!pietro	      write(*,*) 'no problem nsp'
          ELSEIF (com.EQ.'nc  ') THEN
!pietro	  write(*,*) 'reading nc'
              READ (12,1001,err=901) val
              nc=INT(val)
!pietro	   write(*,*) 'nc= ',nc
              CALL SIZE_ARRAYS (3)
!pietro	   write(*,*) 'size arrays called'
              ALLOCATE (K(1:nc,1:nsp))
!pietro	      write(*,*) 'K is allocated'
              K=0
              DO i=1,nc
!pietro	      write(*,*) 'just before the reactions reading',i
!pietro	      write(*,*) 'nsp',nsp,'nc',nc,i
                  READ (12,*,err=901) K(i,1:nsp)
!                  READ (12,1002,err=901) (K(i,j),j=1,nsp)
!pietro		  write(*,*) 'rection is read',i
              END DO
!pietro	      write(*,*) 'no problem nc'
          ELSEIF (com.EQ.'nr  ') THEN
!pietro	  write(*,*) 'reading nr'
              READ (12,1001,err=901) val
              nr=INT(val)
              CALL SIZE_ARRAYS (2)
              DO i=1,nr
!pietro	        write(*,*) 'Just befor ethe reading of nr reactions',i
                  READ (12,*) (NU(i,j),j=1,nsp)
!pietro		write(*,*) 'reaction is read',i
              END DO
!pietro	      write(*,*) 'no problem nr'
          ELSEIF (com.EQ.'xini') THEN
!pietro	  write(*,*) 'reading xini'
              IF (.NOT.ALLOCATED(XINI)) THEN
                  CALL PRINT_ERROR (44,pgname,localinfo,1)
                  STOP 'Errors encountered in PEGASE-LIB mixture file analysis'
              ELSE
              DO i=1,nsp
                  READ (12,1001,err=901) val
                  XINI(i)=DBLE(val)
              END DO
              END IF
!pietro	      write(*,*) 'no problem xini'
          ENDIF
          CALL FILL_WITH_BLANKS (localinfo)
      END DO 
!pietro      write(*,*) 'mixture file is read' 
!     ------------------------------------------------------------------
!     END PART 2
!
!     PART 3: CONSISTENCY CHECKS AND FINAL DATA CONVERSION/COMPUTATION
!     ------------------------------------------------------------------
!     Number of species is zero
!     ^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (nsp.EQ.0) THEN
          localinfo(11:70)=mixname(1:60)
          CALL PRINT_ERROR (6,pgname,localinfo,0)
          istop=1
      ENDIF
!
!     Number of equations has to be nonzero for equilibrium calculations
!     and the left block of the stoechiometric matrix has to be diagonal
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (flg_oper.EQ.1) THEN
          ok=1
          IF (nr.EQ.0) THEN
              localinfo(12:70)=mixname(1:59)
              CALL PRINT_ERROR (7,pgname,localinfo,0)
              istop=1
          ENDIF
          IF (nc.EQ.0) THEN
              localinfo(12:70)=mixname(1:59)
              CALL PRINT_ERROR (14,pgname,localinfo,0)
              istop=1
          ENDIF
          DO i=1,nr
              DO j=1,nr
                  IF (NU(i,j).NE.0) THEN
                      IF (i.NE.j) THEN
                          ok=0
                      ENDIF
                  ENDIF
              END DO
          END DO
          IF (ok.EQ.0) THEN
              localinfo(12:70)=mixname(1:59)
              CALL PRINT_ERROR (15,pgname,localinfo,0)
              istop=1
          ENDIF
      ENDIF
!
!     Conversion of KHI to REAL(kind=8) (only for nc=nsp-nr!)
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (flg_oper.EQ.1) THEN
          IF ((nr+nc).NE.nsp) THEN
              WRITE(localinfo,110) (nr+nc),nsp
110           FORMAT (11x,'Sum equal to:',i3,' Nsp equal to:',i3)
              CALL PRINT_ERROR (16,pgname,localinfo,0)
              istop=1
          ENDIF
          DO i=1,nsp
              DO j=1,nc
                  KHI(j,i)=DBLE(K(j,i))
              END DO
          END DO
!
!     Preparation of XCONS and YCONS vectors (only for nc=nsp-nr!)
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          DO i=1,nsp
              m=m+XINI(i)*MMOL(i)
              DO j=1,nc
                  XCONS(j)=XCONS(j)-KHI(j,i)*XINI(i)
              END DO
          END DO
          DO i=1,nc
              YCONS(i)=XCONS(i)/m
          END DO
!
!     Addition of X-formalism data "outside" the normal frame
!     (again, only for nc=nsp-nr!)
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          XCONS(nc+1)=-1.0d0
          DO i=1,nsp
              KHI(nc+1,i)=1.0d0
          END DO
      ENDIF
      IF ((istop.EQ.1).AND.(flg_stop.EQ.1)) THEN
          STOP 'Errors encountered in PEGASE-LIB mixture file analysis'
      ENDIF
      IF (ALLOCATED(K)) THEN
          DEALLOCATE (K)
      END IF
!pietro      write(*,*) 'mixture file is ended',mxnvm 
      
      RETURN
!     ------------------------------------------------------------------
!     END PART 3
!     ------------------------------------------------------------------
901   localinfo(11:70)=mixname(1:60)
!pietro      write(*,*) 'nsp=',nsp,'i=',i,'j=',j
      CALL PRINT_ERROR (1,pgname,localinfo,1)
!     ------------------------------------------------------------------
      
1000      FORMAT (a80)
1001      FORMAT (f12.6)
1002      FORMAT (13i6)
!     ------------------------------------------------------------------
      END     
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  SETUP   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE read_species_file                                    //
! //                                                                  //
! //  Input:  thermo_dir  directory of the mixture and species files  //    
! //                                                                  //
! //  Outs:   nsp         actual number of species used               // 
! //          nr          actual number of reaction equations used    //
! //          nc          actual number of conservation equations     //
! //                                                                  //
! //  Flags:  flg_stop    forces to stop on library error             //
! //          flg_oper    used to check if it is 1 (equ. calculation) //
! //                                                                  //
! //  The species file is read using a line input instruction (len 80 //
! //  string), with subsequent values on the next lines.              //
! //  Lines not containing commands are ignored.                      //
! //  n_lvls is used to loop on the number of electronic levels.      //
! //  The anharmonicity corrections and Sackur-Tetrode constants are  //
! //  computed within this routine as well.                           //
! //                                                                  //
! //  Benoit Bottin, 08/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE read_species_file (isp,flg_anha,flg_termo,&
     &flg_stop,thermo_dir,spcname)
      USE global_thermo
      IMPLICIT NONE
     
!     ------------------------------------------------------------------
!     PART 1: DEFINITION OF VARIABLES
!     ------------------------------------------------------------------
!
!     Definition of input/output variables
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      INTEGER isp,flg_anha,flg_termo,flg_stop
      CHARACTER*80 thermo_dir,spcname
!
!     Definition of variables intrinsic to the routine
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      
      CHARACTER*80 fullcom,dummy,tablename
      CHARACTER*4 com
      CHARACTER*8 pgname
      CHARACTER*70 localinfo
      REAL(kind=8) val(7)
      INTEGER i,j,n_lvls,istop,FILE_EXISTS,FILE_OPEN,FILE_CLOSE,ierr
!
!pietro      write(*,*) 'start  read_spec',mxnvm
      pgname='pegaslib'
      localinfo(12:70)=spcname(1:59)
!     ------------------------------------------------------------------
!     END PART 1
!
!     PART 2: OPENING THE SPECIES FILE
!     ------------------------------------------------------------------
!     Detection of existence
!     ^^^^^^^^^^^^^^^^^^^^^^
      ierr=FILE_EXISTS (spcname)
      IF (ierr.EQ.0) THEN
          CALL PRINT_ERROR (8,pgname,localinfo,1)
      ELSEIF (ierr.EQ.-1) THEN
          CALL PRINT_ERROR (9,pgname,localinfo,1)
      ENDIF
!
!     Opening the species file
!     ^^^^^^^^^^^^^^^^^^^^^^^^
      ierr=FILE_OPEN (spcname,13,1)
      IF (ierr.NE.1) THEN
          CALL PRINT_ERROR (10,pgname,localinfo,1)
      ENDIF
!     ------------------------------------------------------------------
!     END PART 2
!
!     PART 3: READING THE SPECIES FILE
!     ------------------------------------------------------------------
!     Initialization
!     ^^^^^^^^^^^^^^
      CALL FILL_WITH_BLANKS (tablename)
      CALL FILL_WITH_BLANKS (dummy)
      n_lvls=0
      com='----'
!
!     Main loop on all instructions contained in the species file
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!pietro      write(*,*) 'start reading species files' , isp
      DO WHILE (com.NE.'stop')
          CALL FILL_WITH_BLANKS(fullcom)
          CALL FILL_WITH_BLANKS(dummy)
          CALL FILL_WITH_BLANKS(tablename)
          com='    '
!
!     Reading the species file
!     ^^^^^^^^^^^^^^^^^^^^^^^^
          READ (13,1001,err=901) fullcom
          CALL LCASE(fullcom(1:4),com,ierr)
          IF (com.EQ.'name') THEN
              i=INDEX(fullcom,' ')
              IF (i.EQ.0) THEN
                  CALL PRINT_ERROR (17,pgname,localinfo,0)
              ELSE
                  SPC_LABEL(isp)=fullcom(i+1:i+20)
              ENDIF
!	      write(*,*) 'no problem with name ' ,SPC_LABEL(isp)
          ELSEIF (com.EQ.'file') THEN
              i=INDEX(fullcom,' ')
              IF ((i.EQ.0).AND.(flg_termo.EQ.1)) THEN
                  CALL PRINT_ERROR (18,pgname,localinfo,0)
                  istop=1
              ELSEIF (flg_termo.EQ.1) THEN
                  dummy(1:80-i)=fullcom(i+1:80)
                  CALL CONCATENATE (thermo_dir,dummy,tablename,ierr)
                  CALL READ_REFERENCE_FILE (tablename,isp)
              ENDIF
!	       write(*,*) ' no problem with file ' , SPC_LABEL(isp)
          ELSEIF (com.EQ.'trot') THEN
              READ (13,1002,err=901) val(1)
              TR(1,isp)=DBLE(val(1))
!	      write(*,*) ' no problem with trot ' , SPC_LABEL(isp)
          ELSEIF (com.EQ.'tvib') THEN
!  This part of the code is modified for nvibmode(isp) vibrational modes     
!pietro              write(*,*) 'Coming to the TV'
              READ (13,1002,err=901) val(1)
              nvibmode(isp)=INT(val(1))
                  if (flg_anha.ne.0.and.nvibmode(isp).ne.1) THEN
	 		  CALL PRINT_ERROR (56,pgname,localinfo,1) 
!	               istop=1
                  ENDIF
!pietro	      write(*,*) 'nvibmode' , nvibmode(isp),SPC_LABEL(isp)
              DO i=1, nvibmode(isp) 
                READ (13,1002,err=901) val(1)
!pietro		write(*,*) i,val(1),TV(i,1,isp)
                TV(i,1,isp)=DBLE(val(1))
!pietro		write(*,*) i, 'OK',TV(i,1,isp)
              END DO 
              DO j=1,nsp
                 IF (nvibmode(j).GT.mxnvm) THEN
                   localinfo(12:51)='Context: allocation of nvibmode         '
                   CALL PRINT_ERROR (55,pgname,localinfo,1)
                 END IF
              END DO
!pietro	      write(*,*) ' no problem with  tvib ' , SPC_LABEL(isp)
          ELSEIF (com.EQ.'hfor') THEN
              READ (13,1002,err=901) val(1)
              HFOR(isp)=DBLE(val(1))
!	      write(*,*) ' no problem with hfor ' , SPC_LABEL(isp)
          ELSEIF (com.EQ.'sig ') THEN
              READ (13,1002,err=901) val(1)
              SIG(isp)=DBLE(val(1))
!	      write(*,*) ' no problem with sig ' , SPC_LABEL(isp)
          ELSEIF (com.EQ.'mmol') THEN
              READ (13,1002,err=901) val(1)
              MMOL(isp)=DBLE(val(1))
!	      write(*,*) ' no problem with mmol ' , isp
          ELSEIF (com.EQ.'lvls') THEN
              READ (13,1002,err=901) val(1)
              n_lvls=INT(val(1))
              IF (n_lvls.GT.nlvl) THEN
                  CALL PRINT_ERROR (19,pgname,localinfo,0)
                  n_lvls=nlvl
              ENDIF              
              IF (n_lvls.EQ.0) THEN
                  CALL PRINT_ERROR (20,pgname,localinfo,0)
                  n_lvls=1
              ENDIF
              NELVL(isp)=n_lvls
              READ (13,1001) fullcom
              DO i=1,n_lvls
 !                 write(*,*) 'start'
                  READ (13,1003) (val(j),j=1,7)
 !                 write(*,*) 'here' 
                  GE(i,isp)=DBLE(val(1))
 !                 write(*,*) 'GE'
                  TE(i,isp)=DBLE(val(2))
 !                 write(*,*) 'TE'
                  A_E(i,isp)=DBLE(val(3))
 !                 write(*,*) 'A_E'
                  B_E(i,isp)=DBLE(val(4))
 !                 write(*,*) 'B_E'
                  D_E(i,isp)=DBLE(val(5))
 !                 write(*,*) 'D_E'
                  W_E(i,isp)=DBLE(val(6))
 !                 write(*,*) 'W_E'
                  XW_E(i,isp)=DBLE(val(7))
 !                 write(*,*) 'no error',i 
                  IF (SIG(isp).NE.0.0) THEN
                      D_str(i,isp)=4.0d0*B_E(i,isp)**3/W_E(i,isp)**2
                  ENDIF
 !                 write(*,*) i
              END DO
          ELSEIF (com.EQ.'chrg') THEN
              READ (13,1002,err=901) val(1)
              CHARGE(isp)=INT(val(1))
!	      write(*,*) ' no problem with chrg ' , SPC_LABEL(isp)
          ENDIF
      END DO
 !       write(*,*) 'Alp' 
 !       write(*,*) nvibmode(isp)
 !       write(*,*) mxnvm 
 !     DO i=1, nvibmode(isp)
 !       write(*,*) TV(i,1,isp)
 !     END DO
      
!     ------------------------------------------------------------------
!     END PART 3
!
!     PART 4: OPERATIONS ON THE INPUT DATA
!     ------------------------------------------------------------------
!     Calculation of Sackur-Tetrode constant for entropy/potential
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      ST(isp)=((2.0d0*PIUNIV*MMOL(isp)/NAUNIV/HUNIV/HUNIV)**1.5d0)
      ST(isp)=ST(isp)*(KUNIV**2.5d0)
!          
!     Calculation of specific gas constant
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (MMOL(isp).EQ.0.0d0) THEN
          CALL PRINT_ERROR (21,pgname,localinfo,0)
          istop=1
      ELSE
          RSP(isp)=RUNIV/MMOL(isp)
      ENDIF
!
!     Calculation of simplified anharmonicity correction terms
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (flg_anha.GE.2) THEN
          j=n_lvls
      ELSEIF (flg_anha.NE.0) THEN
          j=1
      ELSE
          j=0
      ENDIF
      DO i=1,j
!         Corrected rotation temperature
          TR(i,isp)=1.43876866d0*(B_E(i,isp)-A_E(i,isp)/2.0d0)
!         Corrected vibration temperature
          TV(1,i,isp)=1.43876866d0*(W_E(i,isp)-2.0d0*XW_E(i,isp))
!         Anharmonicity constants
          IF (B_E(i,isp).GT.0.0d0) THEN
          D_A(i,isp)=A_E(i,isp)/(B_E(i,isp)-0.5d0*A_E(i,isp))
          G_A(i,isp)=2.0d0*D_str(i,isp)/(B_E(i,isp)-0.5d0*A_E(i,isp))
          X_A(i,isp)=XW_E(i,isp)/W_E(i,isp)
          ENDIF
      END DO
      IF (flg_anha.EQ.-1) THEN
          flg_anha=0
      END IF
!
!     Conversion, of electronic energy levels into characteristic temps
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      DO i=1,n_lvls
          TE(i,isp)=1.43876866d0*TE(i,isp)
      END DO
!     ------------------------------------------------------------------
!     END PART 4
!     
!     PART 5: ERROR HANDLING AND FILE CLOSE
!     ------------------------------------------------------------------
      ierr=FILE_CLOSE (13)
      IF (ierr.EQ.-1) THEN
          CALL PRINT_ERROR (13,pgname,localinfo,1)
      ENDIF
      IF ((istop.EQ.1).AND.(flg_stop.EQ.1)) THEN
          STOP 'Errors encountered in PEGASE-LIB species file analysis'
      ELSE
          RETURN
      ENDIF
!     ------------------------------------------------------------------
901   CALL PRINT_ERROR (22,pgname,localinfo,1)
      STOP 'Errors encountered in PEGASE-LIB species file analysis'
!     ------------------------------------------------------------------
1001  FORMAT (a80)
1002  FORMAT (f12.6)              
1003  FORMAT (7f12.6)
!     ------------------------------------------------------------------
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  SETUP   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE read_reference_file                                  //
! //                                                                  //
! //  Input:  thermo_dir  directory of the mixture and species files  //    
! //          traco_dir   directory of the collision integrals files  //
! //                                                                  //
! //  Outs:   nsp         actual number of species used               // 
! //          nr          actual number of reaction equations used    //
! //          nc          actual number of conservation equations     //
! //          traconame   full name of collision integrals file       //
! //                                                                  //
! //  Flags:  flg_stop    forces to stop on library error             //
! //          flg_termo   passed to read_species_file directly        //
! //          flg_traco   passed for consistency check purposes       //
! //          flg_oper    used to check if it is 1 (equ. calculation) //
! //          flg_anha    passed to read_species_file directly        //
! //                                                                  //
! //  Benoit Bottin, 21/10/96. Rewritten Win95 24/7/97.               //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE read_reference_file (tablename,isp)
      USE global_thermo
      IMPLICIT NONE
!
      CHARACTER*80 tablename
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      INTEGER isp,i,j,ierr
      INTEGER FILE_EXISTS,FILE_OPEN,FILE_CLOSE
!
!     Initialization
!     ^^^^^^^^^^^^^^
      pgname='pegaslib'
      localinfo(12:70)=tablename(1:59)
!
!     Open the reference file
!     ^^^^^^^^^^^^^^^^^^^^^
      ierr=FILE_EXISTS (tablename)
      IF (ierr.EQ.0) THEN
          CALL PRINT_ERROR (8,pgname,localinfo,1)
      ELSEIF (ierr.EQ.-1) THEN
          CALL PRINT_ERROR (9,pgname,localinfo,1)
      ENDIF
      ierr=FILE_OPEN (tablename,19,1)
      IF (ierr.NE.1) THEN
          CALL PRINT_ERROR (10,pgname,localinfo,1)
      ENDIF   
!
!     Analyze the reference table file
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      READ (19,1004,err=903) REFPRESS(isp)
      READ (19,1005,err=903) IDATA(isp)
      IF (IDATA(isp).GT.ndata) THEN
          WRITE (localinfo,100) idata,ndata,isp
100       FORMAT (11x,'Requested:',i4,' Available:',i4,'Species:',i4)
          CALL PRINT_ERROR (42,pgname,localinfo,1)
      END IF
      DO i=1,IDATA(isp)
          READ (19,1006,err=903) (REFTABLE(isp,j,i),j=1,4)
      END DO
1004  FORMAT (f16.6)
1005  FORMAT (i6)
1006  FORMAT (4f12.6)
!
!     Close the reference file
!     ^^^^^^^^^^^^^^^^^^^^^^^^
      ierr=FILE_CLOSE (19)
      IF (ierr.EQ.-1) THEN
         CALL PRINT_ERROR (13,pgname,localinfo,0)
         STOP 'Errors encountered in PEGASE-LIB reference file analysis'
      ENDIF
!     ------------------------------------------------------------------
      RETURN
903   CALL PRINT_ERROR (23,pgname,localinfo,0)
      STOP 'Errors encountered in PEGASE-LIB reference file analysis'

      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  SETUP   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE size_arrays                                          //
! //                                                                  //
! //  Input:  imode       operation required from this subroutine     //    
! //                                                                  //
! //  The subroutine is used to allocate the dynamic arrays used by   //
! //  the thermodynamic library PEGASELIB.                            //
! //                                                                  //
! //  Benoit Bottin, 30/09/97.                                        //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE size_arrays (imode)
!
      USE global_thermo
      IMPLICIT NONE
	  INTEGER imode,i,ierr,j,k
      CHARACTER*8 pgname
      CHARACTER*70 localinfo
!
      pgname='pegaslib'
      IF (nsp.LE.0) THEN
          WRITE (localinfo,100) nsp
100       FORMAT (11x,i5)
          CALL PRINT_ERROR (3,pgname,localinfo,1)
      END IF
      SELECT CASE (imode)
!
!     Case 1 is the sizing of arrays depending on nsp alone
      CASE (1)
          IF (nlvl.LE.0) THEN
              WRITE (localinfo,100) nlvl
              CALL PRINT_ERROR (41,pgname,localinfo,1)
	      END IF
          ALLOCATE(TR(1:nlvl,1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of TR               '
             CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF     
          ALLOCATE (TV(1:mxnvm,1:nlvl,1:nsp),STAT=ierr)
!pietro	  write(*,*) 'TV is allocated'
!pietro	  write(*,*) mxnvm
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of TV               '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF         
          ALLOCATE (TE(1:nlvl,1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of TE               '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          ALLOCATE(GE(1:nlvl,1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of GE               '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          ALLOCATE (A_E(1:nlvl,1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of A_E              '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          ALLOCATE (B_E(1:nlvl,1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of B_E              '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          ALLOCATE(D_str(1:nlvl,1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of D_str            '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          ALLOCATE(D_E(1:nlvl,1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of D_E              '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          ALLOCATE (W_E(1:nlvl,1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of W_E              '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          ALLOCATE (XW_E(1:nlvl,1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of XW_E             '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          ALLOCATE (D_A(1:nlvl,1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of D_A              '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          ALLOCATE(G_A(1:nlvl,1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of G_A              '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          ALLOCATE (X_A(1:nlvl,1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of X_A              '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
	    ALLOCATE (HFOR(1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of HFOR             '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          ALLOCATE (ST(1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of ST               '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
! This part is added 
          ALLOCATE (nvibmode(1:nsp),STAT=ierr)
          nvibmode=0
          IF (ierr.NE.0) THEN
              localinfo(12:51)='Context: allocation of nvibmode        '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
! End of modification      
          ALLOCATE (SIG(1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of SIG              '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          ALLOCATE (NELVL(1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of NELVL            '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
	      ALLOCATE (RSP(1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of RSP              '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          ALLOCATE (MMOL(1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of MMOL             '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          ALLOCATE (CHARGE(1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of CHARGE           '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          ALLOCATE (XINI(1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of XINI             '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          ALLOCATE (TABLE_NAME(1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of TABLE_NAME       '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          ALLOCATE (SPC_LABEL(1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of SPC_LABEL        '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          ALLOCATE (SPC_SYMBOL(1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of SPC_SYMBOL       '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          ALLOCATE (REFTABLE(1:nsp,1:4,1:ndata),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of REFTABLE         '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          ALLOCATE (REFPRESS(1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of REFPRESS         '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          ALLOCATE (IDATA(1:nsp),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of IDATA            '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          TR=0.0;TE=0.0;GE=0.0;TV=0.0
          A_E=0.0;B_E=0.0;D_str=0.0;D_E=0.0;W_E=0.0
          XW_E=0.0;D_A=0.0;G_A=0.0;X_A=0.0
          HFOR=0.0;ST=0.0;SIG=0.0;NELVL=0
          RSP=0.0;MMOL=0.0;CHARGE=0;XINI=0.0
!pietro	  write(*,*) 'TV is set to zero'
!pietro	  write(*,*) mxnvm
	  DO i=1,nsp
              CALL FILL_WITH_BLANKS(TABLE_NAME(i))
              CALL FILL_WITH_BLANKS(SPC_LABEL(i))
              CALL FILL_WITH_BLANKS(SPC_SYMBOL(i))
          END DO
          REFTABLE=0.0;REFPRESS=0.0;IDATA=0
!
!     Case 2 is the sizing of arrays depending on nsp and nr
      CASE (2)
          IF (nr.LE.0) THEN
              WRITE (localinfo,100) nr
              CALL PRINT_ERROR (5,pgname,localinfo,1)
          END IF
	      ALLOCATE (NU(1:nr,1:nsp+1),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of NU               '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          NU=0
!
!     Case 3 is the sizing of arrays depending on nsp and nc
      CASE (3)
          IF (nc.LE.0) THEN
              WRITE (localinfo,100) nc
              CALL PRINT_ERROR (4,pgname,localinfo,1)
	      END IF
	      ALLOCATE (KHI(1:nc+1,1:nsp+1),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of KHI              '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          ALLOCATE (XCONS(1:nc+1),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of XCONS            '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF
          ALLOCATE (YCONS(1:nc+1),STAT=ierr)
          IF (ierr.NE.0) THEN
             localinfo(12:51)='Context: allocation of YCONS            '
              CALL PRINT_ERROR (43,pgname,localinfo,1)
          END IF

          KHI=0.0
          XCONS=0.0
          YCONS=0.0
      END SELECT
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  SETUP   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE termostop                                            //
! //                                                                  //
! //  The subroutine is used to deallocate the dynamic arrays used by //
! //  the thermodynamic library PEGASELIB.                            //
! //                                                                  //
! //  Benoit Bottin, 30/09/97.                                        //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE termostop
!
      USE global_thermo
      USE global_equilibrium
      IMPLICIT NONE
      DEALLOCATE(TR,TV,TE,GE,A_E,B_E,D_str,D_E,W_E,XW_E,D_A,G_A,X_A)
      DEALLOCATE (HFOR,ST,SIG,NELVL,RSP,MMOL,CHARGE,XINI,SPC_SYMBOL)
      DEALLOCATE (TABLE_NAME,SPC_LABEL,REFTABLE,REFPRESS,IDATA)
      IF (ALLOCATED(NU)) THEN
          DEALLOCATE (NU,KHI,XCONS,YCONS)
      END IF
      IF (ALLOCATED(jac11)) THEN
          DEALLOCATE (jac11,jac12,jaclow,jacmn)
          DEALLOCATE (f1,f2,indsd,indmn,delv,luindx)
          DEALLOCATE (v,truevar,sys_keq)
      END IF
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
