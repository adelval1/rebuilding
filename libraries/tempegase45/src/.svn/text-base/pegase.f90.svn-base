! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //                                                                  //
! //               P E G A S E   4.   M A I N   ! O D E               //
! //                                                                  //
! //                     M A I N      P R O G R A M                   //
! //                                                                  //
! //                                                                  //
! //                           Benoit Bottin.                         //
! //                   (pegase suite, thermodynamics)                 //
! //                                                                  //
! //                        David Vanden Abeele.                      //
! //                      (transport coefficients)                    //
! //                                                                  //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //  This is version 4.4                                             //
! //                                                                  //
! //  4.4     Finished October 15, 1997                               //
! //          Contains the full thermodynamics code operational and   //
! //          validated, including cp cv and speed of sound of the    //
! //          mixture in equilibrium. Also contains the transport     //
! //          properties part ("traco") in fully operational status.  //
! //                                                                  //
! //  First combined UNIX/DOS release version with full manual and    //
! //  documentation.                                                  //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //  Benoit Bottin, 15/10/97                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      PROGRAM pegase 
!
      USE global_pegase
      USE input_pegase
      USE global_thermo     
      USE interf_thermo
      USE global_chemco      
      USE interf_chemco
      IMPLICIT NONE 
!      
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
!     PART 1: DEFINITION OF VARIABLES
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
!
      CHARACTER*80 inputname,outname
      CHARACTER*70 localinfo
      CHARACTER*6 pgname
      INTEGER i,cntr,cntr2,cntrmax,ierr,niter,ntraco,anha_type
      REAL(kind=8) p,rho,t,v,test,zo
      REAL(kind=8),ALLOCATABLE :: x(:),z(:)
      INTEGER :: e1_fm,e1_to,e1_st,e2_fm,e2_to,e2_st,e1,e2,nstep
      REAL(kind=8) :: m1_fm,m1_to,m1_st,m2_fm,m2_to,m2_st,m1,m2,cnt
!
      INTERFACE
          SUBROUTINE fill_array1 (oper,mode,flg_neq,flg_stop,&
         &p,rho,v,T,TNEQ,x)
              REAL(kind=8) p,rho,v,T,TNEQ(4),x(:)
              INTEGER oper,mode,flg_neq,flg_stop
          END SUBROUTINE fill_array1
          SUBROUTINE fill_array2 (flg_anha,flg_neq,flg_stop,flg_termo,&
         &T,TNEQ,p,x,mode)
              REAL(kind=8) T,TNEQ(4),p,x(:)
              INTEGER mode,flg_anha,flg_neq,flg_stop,flg_termo
          END SUBROUTINE fill_array2
          SUBROUTINE fill_array_mix (mode,flg_anha,flg_neq,flg_log,&
         &flg_oper,flg_stop,flg_termo,p,rho,T,TNEQ,x,epsilon,findiff,niter,zo)
              REAL(kind=8) p,rho,T,TNEQ(4),x(:),epsilon,findiff,zo
              INTEGER mode,flg_anha,flg_neq,flg_oper,flg_stop
              INTEGER flg_termo,niter,flg_log
          END SUBROUTINE fill_array_mix
      END INTERFACE
!
!     initialization of strings
!     ^^^^^^^^^^^^^^^^^^^^^^^^^
      CALL FILL_WITH_BLANKS (inputname)
      CALL FILL_WITH_BLANKS (mixname)
      CALL FILL_WITH_BLANKS (casename)
      CALL FILL_WITH_BLANKS (in_dir)
      CALL FILL_WITH_BLANKS (out_dir)
      CALL FILL_WITH_BLANKS (thermo_dir)
      CALL FILL_WITH_BLANKS (traco_dir)
      CALL FILL_WITH_BLANKS (chemco_dir)
      CALL FILL_WITH_BLANKS (localinfo)
      pgname='pegase'
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
!     END PART 1
!
!     PART 2: USER INTERFACE
!
!     At the middle of the user interface, the main CONTINUE keyword will
!     be found. Because the input contains several computation cases, it
!     is necessary to open it outside of a loop and then to create a loop
!     of unknown number of recurrences, in which the first thing done is
!     to read the input file up to the last instruction of the test case,
!     which will define the flg_more flag.
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
!     Printing information screen
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      CALL presentation
!
!     Opening PEGASE.CTL for getting the directories
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      CALL get_input_param (in_dir,out_dir,thermo_dir,traco_dir,chemco_dir)
!
!     Opening the input file 
!     ^^^^^^^^^^^^^^^^^^^^^^
      CALL prompt_input_file (in_dir,inputname,outname)
      CALL open_input_file (inputname)
!
!     Setting output file open/close flags to close
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      fop_out=0
      fop_dat=0
      fop_da2=0
      fop_cust=0
!
!     The main CONTINUE iterating on all input cases
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      cntr=1
1     CONTINUE
      WRITE (*,3001) cntr
      WRITE (*,3002)
!
!     Reading the input file
!     ^^^^^^^^^^^^^^^^^^^^^^
      CALL explicit_input(mixname,casename,formatname,flg_neq,       &
     &flg_oper,flg_scan,flg_mode,flg_anha,flg_out,flg_dat,flg_da2,   &
     &flg_cust,flg_termo,flg_traco,flg_stop,flg_more,flg_user,t_neq, &
     &var1_from,var1_to,var1_step,var2_from,var2_to,                 &
     &var2_step,epsilon,findiff,ndata,nlvl,niter,ntraco,flg_log,flg_tec)
      IF (flg_cust.NE.0) THEN
          CALL define_custom (in_dir,formatname,flg_cust)
      END IF
      anha_type=flg_anha
!
!     Opening the output files as required
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (flg_log.EQ.1) THEN
          CALL open_output_file (66,out_dir,outname,'.log',flg_tec)
      ENDIF
      IF (flg_out.EQ.1) THEN
          CALL open_output_file (21,out_dir,outname,'.out',flg_tec)
      ENDIF
      IF (flg_dat.EQ.1) THEN
          CALL open_output_file (22,out_dir,outname,'.dat',flg_tec)
      ENDIF
      IF (flg_da2.EQ.1) THEN
          CALL open_output_file (23,out_dir,outname,'.da2',flg_tec)
      ENDIF
      IF (flg_cust.EQ.1) THEN
          CALL open_output_file (24,out_dir,outname,'.res',flg_tec)
      ENDIF
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
!     END PART 2
!
!     PART 3: THERMODYNAMIC DEFINITION
!
!     Test case by test case, the arrays containing the physical
!     and thermodynamic parameters of the mixture used are wiped clean
!     and filled with the needed information. 
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
!     Defining the thermodynamic and transport property data to be used
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      WRITE (*,3003)
      IF (fop_log.EQ.1) WRITE (66,3003)
!
      CALL termodef (thermo_dir,traco_dir,chemco_dir,mixname,flg_anha, &
     &flg_oper,flg_termo,flg_traco,flg_stop,nlvl)

      IF (flg_oper.NE.0) THEN
!!!!!!!!!! Olay BURDAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA ALp !!!!!!!!!!!!!!!!!!!!!!!!!
          CALL tracodef (ntraco)
          WRITE (*,3006)
	 
          CALL chemcodef 
	  
      END IF
!                                
!     Allocation of pegase arrays
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      ALLOCATE (array1(1:10,1:nsp),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of array1           '
          CALL PRINT_ERROR (37,pgname,localinfo,1)
          STOP 'PEGASE terminated due to an error in the main code'
      END IF
      ALLOCATE(array2(1:21,1:nsp,1:8),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of array2           '
          CALL PRINT_ERROR (37,pgname,localinfo,1)
          STOP 'PEGASE terminated due to an error in the main code'
      END IF
      i=MAX (nr,1)
      ALLOCATE (array5(1:i,1:4),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of array5           '
          CALL PRINT_ERROR (37,pgname,localinfo,1)
          STOP 'PEGASE terminated due to an error in the main code'
      END IF
      i= MAX (nsp,7)        
      ALLOCATE (array6(1:9,1:i),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of array6           '
          CALL PRINT_ERROR (37,pgname,localinfo,1)
          STOP 'PEGASE terminated due to an error in the main code'
      END IF
      ALLOCATE (x(1:nsp),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of x                '
          CALL PRINT_ERROR (37,pgname,localinfo,1)
          STOP 'PEGASE terminated due to an error in the main code'
      END IF
      ALLOCATE (z(1:nsp),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of z                '
          CALL PRINT_ERROR (37,pgname,localinfo,1)
          STOP 'PEGASE terminated due to an error in the main code'
      END IF
      array1=0.0
      array2=0.0
      array5=0.0
      array6=0.0
      x=0.0
      !z(1)=XINI(1)	!pietro
      !z(2)=XINI(2)	!pietro
      z(1:nsp)=XINI(1:nsp)+1e-10  !pietro
      IF (flg_user.NE.0) THEN
          ALLOCATE(userarray(1:nusr),STAT=ierr)
          IF (ierr.NE.0) THEN
              localinfo(12:51)='Context: allocation of userarray        '
              CALL PRINT_ERROR (37,pgname,localinfo,1)
              STOP 'PEGASE terminated due to an error in the main code'
          END IF
          userarray=0.0
      END IF
!
!     Loading the Gupta fits information
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF ((flg_traco.EQ.1).OR.(flg_termo.EQ.2)) THEN
          CALL load_gupta_fits (traco_dir,flg_stop)
      ENDIF
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
!     END PART 3
!
!     PART 4: THERMODYNAMIC CALCULATIONS
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
!
!     Calling the user_open routine
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      IF (flg_user.NE.0) THEN
      CALL user_open(in_dir,out_dir,thermo_dir,traco_dir,cntr,flg_user)
      ENDIF
!
!     Analysis of input variables and initialization if needed
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (var1_from.LE.0.0d0) THEN
          WRITE (localinfo,101) cntr,var1_from
101       FORMAT (11x,'Input case ',i3,' Value of variable 1:',d15.8)
          CALL PRINT_ERROR (23,pgname,localinfo,0)
          STOP 'PEGASE terminated due to an error in the main code'
      ENDIF
      IF (var2_from.LE.0.0d0) THEN
          WRITE (localinfo,102) cntr,var2_from
102       FORMAT (11x,'Input case ',i3,' Value of variable 2:',d15.8)
          CALL PRINT_ERROR (24,pgname,localinfo,0)
          STOP 'PEGASE terminated due to an error in the main code'
      ENDIF
      IF ((flg_scan.EQ.0).OR.(flg_scan.EQ.1)) THEN
          var2_to=var2_from
          var2_step=var2_from
      END IF
      IF ((flg_scan.EQ.0).OR.(flg_scan.EQ.2)) THEN
          var1_to=var1_from
          var1_step=var1_from
      END IF

      IF (flg_scan.GE.2) THEN
          IF (var2_from.EQ.var2_to) THEN
              CALL PRINT_ERROR (27,pgname,localinfo,0)
          ENDIF
          IF (var2_step.LE.0.0d0) THEN
              CALL PRINT_ERROR (28,pgname,localinfo,0)
              var2_step=var2_to
          ENDIF
      ENDIF
      IF ((flg_scan.EQ.1).OR.(flg_scan.GE.3)) THEN
          IF (var1_from.EQ.var1_to) THEN
              CALL PRINT_ERROR (25,pgname,localinfo,0)
          ENDIF
          IF (var1_step.LE.0.0d0) THEN
              CALL PRINT_ERROR (26,pgname,localinfo,0)
              var1_step=var1_to
          ENDIF
      ENDIF

!
!     Preparing loop bound values for the scanning
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      cntrmax=1
      IF (flg_scan.LE.2) THEN
          e1_fm=IDINT(DLOG10(var1_from))
          e1_to=e1_fm;e1_st=1
          m1_fm=var1_from/(10.0d0**e1_fm)
          m1_to=m1_fm;m1_st=m1_fm
          e2_fm=IDINT(DLOG10(var2_from))
          e2_to=e2_fm;e2_st=1
          m2_fm=var2_from/(10.0d0**e2_fm)
          m2_to=m2_fm;m2_st=m2_fm
      END IF
!
      IF ((flg_scan.EQ.1).OR.(flg_scan.EQ.3).OR.(flg_scan.EQ.5)) THEN
          e1_fm=IDINT(DLOG10(var1_from))
          e1_to=e1_fm;e1_st=e1_fm
          m1_fm=var1_from/(10.0d0**e1_fm)
          m1_to=var1_to/(10.0d0**e1_fm)
          m1_st=var1_step/(10.0d0**e1_fm)
          cntrmax=cntrmax * ( IDINT((m1_to-m1_fm)/m1_st) + 1 )
      END IF
! 
      IF ((flg_scan.EQ.2).OR.(flg_scan.EQ.3).OR.(flg_scan.EQ.4)) THEN
          e2_fm=IDINT(DLOG10(var2_from))
          e2_to=e2_fm;e2_st=e2_fm
          m2_fm=var2_from/(10.0d0**e2_fm)
          m2_to=var2_to/(10.0d0**e2_fm)
          m2_st=var2_step/(10.0d0**e2_fm)
          cntrmax=cntrmax * ( IDINT((m2_to-m2_fm)/m2_st) + 1 )
      END IF
!
      IF ((flg_scan.EQ.4).OR.(flg_scan.EQ.6)) THEN
          e1_fm=IDINT(DLOG10(var1_from))
          e1_to=IDINT(DLOG10(var1_to))
          e1_st=1
          m1_fm=var1_from/(10.0d0**e1_fm)
          m1_to=var1_from/(10.0d0**(e1_fm-1))
          m1_st=var1_step
          cnt=(m1_to-m1_fm)/m1_st
          IF (CEILING(cnt).EQ.FLOOR(cnt)) THEN
              m1_to=(cnt-1.0d0)*m1_st+m1_fm
          END IF
          nstep=(IDINT((m1_to-m1_fm)/m1_st) + 1) * (e1_to-e1_fm) + 1
          cntrmax=cntrmax * nstep
      END IF
!
      IF ((flg_scan.EQ.5 ).OR.(flg_scan.EQ.6)) THEN
          e2_fm=IDINT(DLOG10(var2_from))
          e2_to=IDINT(DLOG10(var2_to))
          e2_st=1
          m2_fm=var2_from/(10.0d0**e2_fm)
          m2_to=var2_from /(10.0d0**(e2_fm-1))
          m2_st=var2_step
          cnt=(m2_to-m2_fm)/m2_st
          IF (CEILING(cnt).EQ.FLOOR(cnt)) THEN
              m2_to=(cnt-1.0d0)*m2_st+m2_fm
          END IF
          nstep=(IDINT((m2_to-m2_fm)/m2_st + 1)) * (e2_to-e2_fm) + 1
          cntrmax=cntrmax * nstep
      END IF
      cntr2=0
!
!     The main nested DO loops scanning on variable 1
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      exploop1: DO e1=e1_fm,e1_to,e1_st
      mantloop1: DO m1=m1_fm,m1_to,m1_st
!
!     Computing variable 1 value depending on the type of scan performed
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      v=m1*10.0d0**e1
      IF (v.GT.(var1_to+epsilon)) THEN
          EXIT exploop1
      ENDIF
!
!     The main nested DO loops scanning on variable 2
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      exploop2: DO e2=e2_fm,e2_to,e2_st
      mantloop2: DO m2=m2_fm,m2_to,m2_st
!
!     Computing variable 2 value depending on the type of scan performed
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      t=m2*10.0d0**e2
      IF (t.GT.(var2_to+epsilon)) THEN
          EXIT exploop2
      ENDIF
!
      cntr2=cntr2+1
      WRITE (*,3004) cntr2,cntrmax
      IF (fop_log.EQ.1) WRITE (66,3004) cntr2,cntrmax
      IF (fop_log.EQ.1) WRITE (66,*)
!
!     Resetting problem size to one species in single-species case
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (flg_oper.EQ.0) THEN
          IF (nsp.NE.1) THEN
              WRITE (localinfo,103) cntr,nsp
103           FORMAT (11x,'Input case:',i3,' Number of species:',i3)
              CALL PRINT_ERROR (29,pgname,localinfo,0)
          ENDIF
          nsp=1
      ENDIF

!
!     Computing equilibrium composition in equilibrium multispecies
!     If mode = 1, x is a mole fraction array      
!     If mode = 2, x is a mass fraction array
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (flg_oper.EQ.1) THEN
          IF (flg_mode.EQ.1) THEN
              CALL MOLEFRAC(x,v,T,zo,T_NEQ,epsilon,flg_log, &
            & flg_anha,flg_neq,flg_stop,flg_termo,1,z,niter)
          ELSE
              CALL MASSFRAC(x,v,T,T_NEQ,epsilon,flg_log,    &
            & flg_anha,flg_neq,flg_stop,flg_termo,1,z,niter)
          ENDIF
      ENDIF
!
!     Storing composition in z array as initial guess for next step
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (T.LT.600) THEN !pietro -1-
      !write(*,'(e14.9)')T !pietro -2-
      !x(1)=XINI(1)	!pietro -3-
      !x(2)=XINI(2)	!pietro -4-
      x(1:nsp)=XINI(1:nsp)+1e-10  !pietro -5-
      ELSE
	z=x
      END IF
!	before instead of 1-2-3-4-5- there was only z=x
      
!
!     Setting the initial composition in the case of a nonequilibrium
!     or frozen calculation. Test checks that X_NEQ is not empty. If
!     it would be the case, then XINI from the mixture file is used.
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      test=0.0d0
      IF (flg_oper.EQ.2) THEN
          DO i=1,nsp
              x(i)=X_NEQ(i)
              test=test+X_NEQ(i)
          END DO
          IF (IDNINT(test).NE.1) THEN
              test=0.0d0
              DO i=1,nsp
                  x(i)=XINI(i)
                  test=test+XINI(i)
              END DO
              IF (IDNINT(test).NE.1) THEN
                  WRITE (localinfo,104) test
104               FORMAT (11x,'Total of mole/mass fractions is:',d15.8)
                  CALL PRINT_ERROR (30,pgname,localinfo,0)
              STOP 'PEGASE terminated due to an error in the main code'
              ENDIF
          ENDIF
      ENDIF
!
!     Computing and filling out array1 values
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      CALL fill_array1 (flg_oper,flg_mode,flg_neq,flg_stop,&
     &p,rho,v,T,t_neq,x)
!
!     Computing and filling out array2 values
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      CALL fill_array2 (flg_anha,flg_neq,flg_stop,flg_termo,&
     &T,t_neq,p,x,flg_mode)
!
!     Computing and filling out array3 and array4 values
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      CALL fill_array_mix (flg_mode,flg_anha,flg_neq,flg_log, &
     &flg_oper,flg_stop,flg_termo,p,rho,T,t_neq,x,epsilon,findiff,niter,zo)
!
!     Computing equilibrium constants if flg_oper=1
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (flg_oper.EQ.1) THEN
          CALL EQ_CONSTANTS (T,t_neq,p,rho,array5, &
     &    flg_anha,flg_neq,flg_stop,flg_termo)
      ENDIF

!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
!     END PART 4
!
!     PART 5: TRANSPORT PROPERTIES CALCULATIONS
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
!
      IF ((flg_traco.GE.0).AND.(flg_oper.NE.0)) THEN
          CALL fill_array6 (p,T,T_NEQ,flg_neq, &
        & flg_termo,flg_traco,flg_anha)
      END IF
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
!     END PART 5
!
!     PART 6: USER COMPUTATION AND PRINTOUT OF OUTPUT FILES
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
!
!     Allow the user to perform its own computations
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (flg_user.NE.0) THEN
          CALL user_compute 
      ENDIF
!
!     Printing the ouptut files
!     ^^^^^^^^^^^^^^^^^^^^^^^^^
      WRITE (*,3005)
      IF (fop_log.EQ.1) WRITE (66,3005)
!
      IF (flg_out.EQ.1) THEN
          CALL print_output(cntr,cntr2,cntrmax,casename,mixname, &
         &anha_type,flg_neq,flg_mode,flg_oper,flg_termo,flg_traco,flg_user)
      ENDIF
      IF (flg_dat.EQ.1) THEN
          CALL print_data(flg_mode,cntr,cntr2,casename,flg_tec,cntrmax)
      ENDIF
      IF (flg_da2.EQ.1) THEN
          CALL print_data2(flg_mode,cntr,cntr2,casename,flg_tec,cntrmax)
      ENDIF
      IF (flg_cust.EQ.1) THEN
          CALL print_custom (cntr,cntr2,casename,flg_tec,cntrmax)
      ENDIF
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
!     END PART 6
!
!     PART 7: CLOSING THE LOOP (and the input file)
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
!
!     End of the scanning loop for one test case
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      END DO mantloop2
      END DO exploop2
      END DO mantloop1
      END DO exploop1
!
!     Allow the user to close its files
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      CALL close_output_files
      IF (flg_user.NE.0) THEN
          CALL user_close (flg_user,cntr)
      ENDIF
!      
!     Deallocate allocated arrays
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      DEALLOCATE (array1,array2)
      DEALLOCATE (array5,array6)
      IF (flg_user.NE.0) THEN
          DEALLOCATE (userarray)
      END IF
      DEALLOCATE (x,z)
      IF (ALLOCATED(X_NEQ)) THEN
          DEALLOCATE(X_NEQ)
      END IF
      CALL termostop
      IF (flg_oper.NE.0) THEN
          CALL tracostop
          CALL chemcostop
      END IF
!      
!     End of the main loop on all test cases
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (flg_more.EQ.1) THEN
          cntr=cntr+1
          GOTO 1
      ENDIF
!
!     Close the input file
!     ^^^^^^^^^^^^^^^^^^^^
      CALL close_input_file
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
!     END PART 7
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
3001  FORMAT (' ===================================================',/,&
     &      ' ==> Pegase now processes input case number ',i4,' <==',/,&
     &      ' ===================================================')
3002  FORMAT (' --> Reading the input file ...')
3003  FORMAT (' --> Retrieving mixture and species data ...')
3004  FORMAT (' --> Processing calculation ',i4,' of ',i4,' ...')
3005  FORMAT (' --> Printing the output files ...')
3006  FORMAT (' --> Retrieving chemical kinetic data ...') 

!     ------------------------------------------------------------------
      WRITE (*,*)
      STOP 'Program PEGASE successfully terminated'
      END
! //////////////////////////////////////////////////////////////////////


! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //               P E G A S E   4.   M A I N   ! O D E               //
! //                                                                  //
! //  SUBROUTINE presentation                                         //
! //                                                                  //
! //  The routine prints the presentation screen                      //
! //                                                                  //
! //  Benoit Bottin, 03/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE presentation
      LOGICAL iex
      CHARACTER*72 line
!      
      INQUIRE (FILE='pegase.ttl',EXIST=iex)
!
      IF (.NOT.iex) THEN
          WRITE (*,9001)
      ELSE
          OPEN (unit=32,file='pegase.ttl')
101       CONTINUE
              READ (32,1001,end=900) line
              WRITE (*,1001) line
          GOTO 101
900       CLOSE (32)
      ENDIF 
1001  FORMAT (a72)
!     ------------------------------------------------------------------
9001  FORMAT (/,/,' P E G A S E   v4.5a',/,/,                       &
     &' PErfect GAS Equation solver for arbitrary mixtures',/,/,    &
     &' By B. Bottin, with D. Vanden Abeele and P. Barbante',/,/,   &
     &' VKI AR Dpt., January 1999',/)
!     ------------------------------------------------------------------
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

