! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  INPUT   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE explicit_input                                       //
! //                                                                  //
! //  The subroutine is used to get all the required parameters for   //
! //  the computation from the input file. It also returns whether    //
! //  a stop or a more instruction has been reached.                  //
! //                                                                  //
! //  The format of the input file is explicit, self-explanatory but  //
! //  of fixed format and rather lengthy. It is not possible to       //
! //  modify it in any way. If a short and concise input is preferred //
! //  then the subroutine implicit_input should be used instead.      //
! //                                                                  //
! //  Input: none. The input file has unit number 11.                 //    
! //                                                                  //
! //  Outs:   mixname     name of mixture file (80 char. max.)        // 
! //          casename    string to comment the actual computation    //
! //          flg_neq     non-equilibrium computation case            //
! //          flg_oper    0: single, 1: equilibrium, 2: frozen        //
! //          flg_scan    0: no, 1: single, 2: double                 //
! //          flg_mode    input mode 1: p,T 2: rho,T                  //
! //          flg_anha    0: no, 1: simple, 2: full                   //
! //          flg_out     written output file (y/n)                   //
! //          flg_dat     tabulated output file (y/n)                 //
! //          flg_da2     additional tabulated output (y/n)           //
! //          flg_cust    custom output file (y/n)                    //
! //          flg_termo   use of thermo curve fits (y/n)              //
! //          flg_traco   use of traco curve fits (y/n)               //
! //          flg_stop    stop program on library error (y/n)         //
! //          flg_more    0: stops program, 1: another case follows   //
! //          flg_user    0: do not (1:do) run through user_ routines //
! //          var1_from   value of first variable (p or rho) start    //
! //          var1_to     value of first variable (p or rho) end      //
! //          var1_step   value of first variable (p or rho) step     //
! //          var2_from   value of second variable (T) start          //
! //          var2_to     value of second variable (T) end            //
! //          var2_step   value of second variable (T) step           //
! //          epsilon     finite differences tolerance value          //
! //          t_neq(4)    array of nonequilibrium temperatures        //
! //          x_neq(sp)   array of initial compositions (mole frac)   //
! //          eps         value of epsilon for finite differences     //
! //                                                                  //
! //  The attached routine IFERROR is simply used to avoid numerous   //
! //  tests of the returned variable from SKIPLINES, so that the code //
! //  remains compact. It tests istop, prints the appropriate error   //
! //  message and stops the code.                                     //
! //                                                                  //
! //  Benoit Bottin, 03/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE explicit_input(mixname,casename,formatname,flg_neq, &
     &flg_oper,flg_scan,flg_mode,flg_anha,flg_out,flg_dat,flg_da2,   &
     &flg_cust,flg_termo,flg_traco,flg_stop,flg_more,flg_user,       &
     &t_neq,var1_from,var1_to,var1_step,var2_from,var2_to,           &
     &var2_step,epsilon,findiff,ndata,nlvl,niter,ntraco,flg_log,flg_tec)

      USE global_pegase
      IMPLICIT NONE
      CHARACTER*4 more
      CHARACTER*6 pgname
      CHARACTER*70 localinfo
      CHARACTER*80 mixname,casename,formatname,dummy
      INTEGER flg_neq,flg_oper,flg_scan,flg_mode,flg_anha,flg_out
      INTEGER flg_dat,flg_da2,flg_cust,flg_termo,flg_traco,flg_stop
      INTEGER flg_more,flg_user,istop,i,isp,ierr,ndata,nlvl,niter
      INTEGER ntraco,flg_log,flg_tec
      LOGICAL, SAVE :: firstcall = .TRUE.
      LOGICAL, SAVE :: implicit_format = .FALSE.
      REAL(kind=8) var1_from,var1_to,var1_step,var2_from,var2_to,var2_step
      REAL(kind=8) epsilon,findiff,t_neq(4)
!
      istop=0 
      pgname='pegase'
      CALL FILL_WITH_BLANKS(localinfo)   
!
!     Check for input file format on first read
!
      IF (firstcall) THEN
          READ (11,1001,err=901,end=902) dummy
          IF (dummy(1:50).NE.REPEAT('/',50)) THEN
              implicit_format = .TRUE.
          END IF
          REWIND(11)
          firstcall = .FALSE.
      END IF
!
!     Check for implicit file format on all reads
!
      IF (implicit_format) THEN
        CALL implicit_input(mixname,casename,formatname,flg_neq,      &
     &  flg_oper,flg_scan,flg_mode,flg_anha,flg_out,flg_dat,flg_da2,  &
     &  flg_cust,flg_termo,flg_traco,flg_stop,flg_more,flg_user,      &
     &  t_neq,var1_from,var1_to,var1_step,var2_from,var2_to,var2_step,&
     &  epsilon,findiff,ndata,nlvl,niter,ntraco,flg_log,flg_tec)
        RETURN
      END IF
!
!     Start reading explicit input file
!
      CALL SKIPLINES (11,5,istop)
      CALL IFERROR (istop)
      READ (11,1001,err=901,end=902) mixname
      CALL SKIPLINES (11,1,istop)
      CALL IFERROR (istop)
      READ (11,1002,err=901,end=902) casename
      CALL SKIPLINES (11,2,istop)
      CALL IFERROR (istop)
      READ (11,1003,err=901,end=902) flg_oper
      READ (11,1003,err=901,end=902) flg_mode
      READ (11,1003,err=901,end=902) flg_scan
      READ (11,1003,err=901,end=902) flg_anha
      READ (11,1003,err=901,end=902) flg_neq
      READ (11,1003,err=901,end=902) flg_termo
      READ (11,1003,err=901,end=902) flg_traco
      READ (11,1003,err=901,end=902) flg_stop
      READ (11,1003,err=901,end=902) flg_user
      CALL SKIPLINES (11,2,istop)
      CALL IFERROR (istop)
      READ (11,1004,err=901,end=902) flg_log
      READ (11,1004,err=901,end=902) flg_out
      READ (11,1004,err=901,end=902) flg_dat
      READ (11,1004,err=901,end=902) flg_da2
      READ (11,1004,err=901,end=902) flg_cust
      READ (11,1001,err=901,end=902) formatname
      CALL SKIPLINES (11,2,istop)
      CALL IFERROR (istop)
      READ (11,1005,err=901,end=902) var1_from
      READ (11,1005,err=901,end=902) var1_to
      READ (11,1005,err=901,end=902) var1_step
      READ (11,1005,err=901,end=902) var2_from
      READ (11,1005,err=901,end=902) var2_to
      READ (11,1005,err=901,end=902) var2_step
      CALL SKIPLINES (11,2,istop)
      CALL IFERROR (istop)
      READ (11,1006,err=901,end=902) epsilon
      READ (11,1006,err=901,end=902) findiff
      READ (11,1003,err=901,end=902) niter
      READ (11,1003,err=901,end=902) nlvl
      READ (11,1003,err=901,end=902) ndata
      READ (11,1003,err=901,end=902) nusr
      READ (11,1003,err=901,end=902) ntraco
      CALL SKIPLINES (11,2,istop)
      CALL IFERROR (istop)
      DO i=1,4
          READ (11,1006,err=901,end=902) t_neq(i)
      END DO
      READ (11,1007,err=901,end=902) isp
      IF (ALLOCATED(X_NEQ)) THEN
          DEALLOCATE (X_NEQ)
      END IF
      IF (isp.EQ.0) THEN
          ALLOCATE (X_NEQ(1:1),STAT=ierr)
      ELSE
          ALLOCATE (X_NEQ(1:isp),STAT=ierr)
      END IF
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of X_NEQ            '
          CALL PRINT_ERROR (37,pgname,localinfo,1)
          STOP 'PEGASE terminated due to an error in the main code'
      END IF
      DO i=1,isp
          READ (11,1008,err=901,end=902) X_NEQ(i)
      END DO
      CALL SKIPLINES (11,1,istop)
      CALL IFERROR (istop)
      READ (11,1009,err=901,end=902) more
      IF (more.EQ.'MORE') then
          flg_more=1
      ELSE
          flg_more=0
      ENDIF
      RETURN
!     ------------------------------------------------------------------
901   CALL IFERROR (3)
902   CALL IFERROR (2)
!     ------------------------------------------------------------------
1001  FORMAT (16x,a80)
1002  FORMAT (14x,a80)
1003  FORMAT (59x,i5)
1004  FORMAT (25x,i5)
1005  FORMAT (20x,f20.6)
1006  FORMAT (34x,f20.6)
1007  FORMAT (46x,i5)
1008  FORMAT (42x,f20.6)
1009  FORMAT (59x,a4)
8001  FORMAT ('           Used size:',i4,' Declared size:',i4)
!     ------------------------------------------------------------------
      END
!     ------------------------------------------------------------------
      SUBROUTINE iferror (istop)
!     1: error in the skiplines subroutine
!     2: generic EOF error
!     3: error in the data retrieval reads
!     ------------------------------------------------------------------
      INTEGER istop
      CHARACTER*6 pgname
      CHARACTER*70 localinfo
!
      pgname='pegase'
      localinfo(1:1)=' '     
      IF (istop.EQ.1) THEN
          CALL PRINT_ERROR (14,pgname,localinfo,0)
          STOP 'PEGASE terminated during input file analysis'
      ELSEIF (istop.EQ.2) THEN
          CALL PRINT_ERROR (15,pgname,localinfo,0)
          STOP 'PEGASE terminated during input file analysis'
      ELSEIF (istop.EQ.3) THEN
          CALL PRINT_ERROR (16,pgname,localinfo,0)
          STOP 'PEGASE terminated during input file analysis'
      ELSE
          RETURN
      ENDIF
      END
! //////////////////////////////////////////////////////////////////////
