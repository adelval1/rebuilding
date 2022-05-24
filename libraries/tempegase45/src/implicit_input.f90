! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  INPUT   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE implicit_input                                       //
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
! //  Benoit Bottin, 12/01/99                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE implicit_input(mixname,casename,formatname,flg_neq, &
     &flg_oper,flg_scan,flg_mode,flg_anha,flg_out,flg_dat,flg_da2,   &
     &flg_cust,flg_termo,flg_traco,flg_stop,flg_more,flg_user,       &
     &t_neq,var1_from,var1_to,var1_step,var2_from,var2_to,           &
     &var2_step,epsilon,findiff,ndata,nlvl,niter,ntraco,flg_log,flg_tec)

      USE global_pegase
      IMPLICIT NONE
      CHARACTER*6 pgname
      CHARACTER*20 com,subcom,oldcom
      CHARACTER*70 localinfo
      CHARACTER*80 mixname,casename,formatname,line,arg
      INTEGER flg_neq,flg_oper,flg_scan,flg_mode,flg_anha,flg_out
      INTEGER flg_dat,flg_da2,flg_cust,flg_termo,flg_traco,flg_stop
      INTEGER flg_more,flg_user,istop,i,isp,ierr,ndata,nlvl,niter
      INTEGER ntraco,flg_log,flg_tec
      INTEGER scan_var(1:3),cnt
      LOGICAL unknown, is_var(1:3)
      REAL(kind=8) var1_from,var1_to,var1_step,var2_from,var2_to,var2_step
      REAL(kind=8) epsilon,findiff,t_neq(4)
      REAL(kind=8) a
!
      is_var(:) = .FALSE.
      scan_var(:) = 0
      istop=0 
      cnt = 0
      pgname='pegase'
      CALL FILL_WITH_BLANKS(localinfo)   
!
!    Setting up default information
!    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
        is_var(2) = .FALSE.
        is_var(3) = .TRUE.
        mixname(1:29) = 'no_file_defined_in_input_file'
        formatname(1:29) = 'no_format_file_name_specified'
        flg_oper = 1
        flg_mode = 2
        flg_scan = 0
        flg_anha = 0
        flg_neq = 0
        flg_termo = 0
        flg_traco = 0
        flg_stop = 1
        flg_user = 0
        flg_log = 1
        flg_out = 0
        flg_dat = 1
        flg_da2 = 1
        flg_cust = 0
        flg_tec = 0
        var1_from = 1.29313d0
        var1_to = 1.3d0
        var1_step = 1.0d0
        var2_from = 2000.0d0
        var2_to = 2050.0d0
        var2_step = 100.0d0
        epsilon = 1.0d-12
        findiff = 1.0d-3
        niter = 1000
        nlvl = 20 
        ndata = 50
        nusr = 9
        ntraco = 3
        T_NEQ(:) = 0.0d0
        flg_more = 0
        isp = 0
      IF (ALLOCATED(X_NEQ)) THEN
          DEALLOCATE (X_NEQ)
      END IF
!
!     Reading the input file until a more or stop instruction
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
    DO
        READ (UNIT=11, FMT=1001, ERR=901, END=902) line
        CALL PARSE_INPUT (line, com, subcom, arg)
1001    FORMAT (a80)
1002    FORMAT (e20.12)
1003    FORMAT (a75)
1004    FORMAT (a80)
!       ---------------------
!       = Full comment line =
!       ---------------------
        IF (com(1:1).EQ.'!') THEN 
            CYCLE
        END IF
!       -------------------------------------------------------
!       = Repeated subcommand, needs to assign stored command =
!       -------------------------------------------------------
        IF (com(1:1).EQ.'"') THEN 
            com = oldcom
        END IF
        unknown = .FALSE.
        CALL FILL_WITH_BLANKS (localinfo)
!
        SELECT CASE (com)
!   ------------------------------------------------------------------
        CASE ('output              ') !------------------------ OUTPUT
            SELECT CASE (subcom)
            CASE ('log file            ')
                IF (arg(1:1).EQ.'y') THEN
                    flg_log = 1
                ELSE
                    flg_log = 0
                END IF
            CASE ('output file         ')
                IF (arg(1:1).EQ.'y') THEN
                    flg_out = 1
                ELSE
                    flg_dat = 0
                END IF
            CASE ('data file 1         ','first data file     ')
                IF (arg(1:1).EQ.'y') THEN
                    flg_dat = 1
                ELSE
                    flg_dat = 0
                END IF
            CASE ('data file 2         ','second data file    ')
                IF (arg(1:1).EQ.'y') THEN
                    flg_da2 = 1
                ELSE
                    flg_da2 = 0
                END IF
            CASE ('custom file         ')
                IF (arg(1:1).EQ.'y') THEN
                    flg_cust = 1
                ELSE
                    flg_cust = 0
                END IF
            CASE ('use custom format   ')
                READ (UNIT=arg, FMT=1003) formatname(1:75)
                formatname(76:80)='     '
                formatname=ADJUSTL(formatname)
                flg_cust = 1
            CASE ('use tecplot format  ')
                IF (arg(1:1).EQ.'y') THEN
                    flg_tec = 1
                ELSE
                    flg_tec = 0
                END IF
            CASE DEFAULT
                unknown = .TRUE.
            END SELECT
!   ------------------------------------------------------------------
        CASE ('settings            ') !---------------------- SETTINGS
            SELECT CASE (subcom)
            CASE ('stop on errors      ')
                IF (arg(1:1).EQ.'y') THEN
                    flg_stop = 1
                ELSE
                    flg_stop = 0
                END IF
            CASE ('computation type    ')
                IF (arg(1:6).EQ.'single') THEN
                    flg_oper = 0
                ELSEIF (arg(1:6).EQ.'equili') THEN
                    flg_oper = 1
                ELSEIF (arg(1:6).EQ.'frozen') THEN
                    flg_oper = 2
                END IF
            CASE ('residual tolerance  ')
                READ (UNIT=arg, FMT=1002) epsilon
            CASE ('finite difference   ')
                READ (UNIT=arg, FMT=1002) findiff
            CASE ('maximum iteration   ')
                READ (UNIT=arg, FMT=1002) a
                niter = IDINT(a)
            CASE ('electronic levels   ')
                READ (UNIT=arg, FMT=1002) a
                nlvl = IDINT(a)
            CASE ('reference file size ')
                READ (UNIT=arg, FMT=1002) a
                ndata = IDINT(a)
            CASE DEFAULT
                unknown = .TRUE.
            END SELECT
!   ------------------------------------------------------------------
       CASE ('user routines       ') !------------------ USER ROUTINES
            SELECT CASE (subcom)
            CASE ('enabled             ')
                IF (arg(1:1).EQ.'y') THEN
                    flg_user = -1
                ELSE
                    flg_user = 0
                END IF
            CASE ('disabled            ')
                IF (arg(1:1).EQ.'y') THEN
                    flg_user = 0
                ELSE
                    flg_user = -1
                END IF
            CASE ('control flag value  ')
                READ (UNIT=arg, FMT=1002) a
                flg_user = IDINT(a)
            CASE ('user array size     ')
                READ (UNIT=arg, FMT=1002) a
                nusr = IDINT(a)
            CASE DEFAULT
                unknown = .TRUE.
            END SELECT
!   ------------------------------------------------------------------
        CASE ('run with            ') !---------------------- RUN WITH
            SELECT CASE (subcom)
            CASE ('computation type    ')
                IF (arg(1:6).EQ.'single') THEN
                    flg_oper = 0
                ELSEIF (arg(1:6).EQ.'equili') THEN
                    flg_oper = 1
                ELSEIF (arg(1:6).EQ.'frozen') THEN
                    flg_oper = 2
                END IF
            CASE ('thermodynamics from ')
                IF (arg(1:6).EQ.'theory') THEN
                    flg_termo = 0
                ELSEIF (arg(1:6).EQ.'refere') THEN
                    flg_termo = 1
                END IF
            CASE ('transport from      ')
                IF (arg(1:6).EQ.'theory') THEN
                    flg_traco = 0
                ELSEIF (arg(1:6).EQ.'gupta ') THEN
                    flg_traco = 1
                END IF
            CASE ('anharmonicity       ')
                IF (arg(1:4).EQ.'none') THEN
                    flg_anha = 0
                ELSEIF (arg(1:4).EQ.'simp') THEN
                    flg_anha = 1
                ELSEIF (arg(1:4).EQ.'adva') THEN
                    flg_anha = 2
                END IF
            CASE ('thermal equilibrium ')
                IF (arg(1:1).EQ.'y') THEN
                    flg_neq = 0
                ELSE
                    flg_neq = 1
                END IF
            CASE ('sonine polynomials  ')
                READ (UNIT=arg, FMT=1002) a
                ntraco = IDINT(a)
            CASE ('mixture file        ')
                READ (UNIT=arg, FMT=1003) mixname(1:75)
                mixname(76:80)='     '
                mixname=ADJUSTL(mixname)
            CASE ('case label          ','label               ')
                READ (UNIT=arg, FMT=1003) casename(1:75)
                casename(76:80)='     '
                casename=ADJUSTL(casename)
            CASE DEFAULT
                unknown = .TRUE.
            END SELECT
!   ------------------------------------------------------------------
        CASE ('pressure            ') !---------------------- PRESSURE
            i = 1
            is_var(i) = .TRUE.
            SELECT CASE (subcom)
            CASE ('value               ')
                READ (UNIT=arg, FMT=1002) var1_from
                var1_to = var1_from + 1.0d0
                var1_step = var1_from + 2.0d0
                scan_var(i) = 0
            CASE ('from                ')
                READ (UNIT=arg, FMT=1002) var1_from
                scan_var(i) = 1
            CASE ('to                  ')
                READ (UNIT=arg, FMT=1002) var1_to
                scan_var(i) = 1
            CASE ('step                ')
                READ (UNIT=arg, FMT=1002) var1_step
                scan_var(i) = 1
            CASE ('scan type           ')
                IF (arg(1:4).EQ.'none') THEN
                    scan_var(i) = 0
                ELSEIF (arg(1:4).EQ.'line') THEN
                    scan_var(i) = 1
                ELSEIF (arg(1:4).EQ.'loga') THEN
                    scan_var(i) = 2
                END IF
            CASE DEFAULT
                unknown = .TRUE.
                is_var(i) = .FALSE.
            END SELECT
!   ------------------------------------------------------------------
        CASE ('density             ') !----------------------- DENSITY
            i = 2
            is_var(i) = .TRUE.
            SELECT CASE (subcom)
            CASE ('value               ')
                READ (UNIT=arg, FMT=1002) var1_from
                var1_to = var1_from + 1.0d0
                var1_step = var1_from + 2.0d0
                scan_var(i) = 0
            CASE ('from                ')
                READ (UNIT=arg, FMT=1002) var1_from
                scan_var(i) = 1
            CASE ('to                  ')
                READ (UNIT=arg, FMT=1002) var1_to
                scan_var(i) = 1
            CASE ('step                ')
                READ (UNIT=arg, FMT=1002) var1_step
                scan_var(i) = 1
            CASE ('scan type           ')
                IF (arg(1:4).EQ.'none') THEN
                    scan_var(i) = 0
                ELSEIF (arg(1:4).EQ.'line') THEN
                    scan_var(i) = 1
                ELSEIF (arg(1:4).EQ.'loga') THEN
                    scan_var(i) = 2
                END IF
            CASE DEFAULT
                unknown = .TRUE.
                is_var(i) = .FALSE.
            END SELECT
!   ------------------------------------------------------------------
        CASE ('temperature         ') !------------------- TEMPERATURE
            i = 3
            is_var(i) = .TRUE.
            SELECT CASE (subcom)
            CASE ('value               ')
                READ (UNIT=arg, FMT=1002) var2_from
                var2_to = var2_from + 1.0d0
                var2_step = var2_from + 2.0d0
                scan_var(i) = 0
            CASE ('from                ')
                READ (UNIT=arg, FMT=1002) var2_from
                scan_var(i) = 1
            CASE ('to                  ')
                READ (UNIT=arg, FMT=1002) var2_to
                scan_var(i) = 1
            CASE ('step                ')
                READ (UNIT=arg, FMT=1002) var2_step
                scan_var(i) = 1
            CASE ('scan type           ')
                IF (arg(1:4).EQ.'none') THEN
                    scan_var(i) = 0
                ELSEIF (arg(1:4).EQ.'line') THEN
                    scan_var(i) = 1
                ELSEIF (arg(1:4).EQ.'loga') THEN
                    scan_var(i) = 2
                END IF
            CASE DEFAULT
                unknown = .TRUE.
                is_var(i) = .FALSE.
            END SELECT
!   ------------------------------------------------------------------
        CASE ('frozen composition  ') !------------ FROZEN COMPOSITION
            SELECT CASE (subcom)
            CASE ('number of species   ')
                READ (UNIT=arg, FMT=1002) a
                isp = IDINT(a)
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
            CASE ('mole fraction       ')
                IF (cnt.EQ.isp) THEN
                    CALL PRINT_ERROR (38,pgname,localinfo,0)
                END IF
                IF (.NOT.ALLOCATED(X_NEQ)) THEN
                    CALL PRINT_ERROR (39,pgname,localinfo,1)
                    STOP 'PEGASE terminated due to an error in the main code'
                END IF
                cnt = cnt + 1
                READ (UNIT=arg, FMT=1002) X_NEQ(cnt)
            CASE DEFAULT
                unknown = .TRUE.
            END SELECT
!   ------------------------------------------------------------------
        CASE ('more                ') !--------------------- MORE/STOP
            flg_more = 1
            EXIT
        CASE ('stop                ')
            flg_more = 0
            EXIT
!   ------------------------------------------------------------------
        CASE DEFAULT
            unknown = .TRUE.
        END SELECT
!   ------------------------------------------------------------------
        IF (unknown) THEN
            localinfo(12:70)=line(1:59)
            CALL PRINT_ERROR (41,pgname,localinfo,0)
        END IF
        oldcom = com
    END DO
!
!     Performing post-processing operations (setting the scan)
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
902 i = -1
    IF (is_var(1).AND.is_var(3)) THEN   ! pressure-temperature scan
        i = 10*scan_var(1) + scan_var(3)
        flg_mode = 1
    END IF
    IF (is_var(2).AND.is_var(3)) THEN    ! density-temperature scan
        i = 10*scan_var(2) + scan_var(3)
        flg_mode = 2
    END IF
    IF (i.EQ.-1) THEN                              ! unrecognized scan
            CALL PRINT_ERROR (40,pgname,localinfo,1)
            STOP 'PEGASE terminated due to an error in the main code'
    END IF
    SELECT CASE (i)
    CASE (00)
        flg_scan = 0
    CASE (10)
        flg_scan = 1
    CASE (01)
        flg_scan = 2
    CASE (11)
        flg_scan = 3
    CASE (20,21)
        flg_scan = 4
    CASE (02,12)
        flg_scan = 5
    CASE (22)
        flg_scan = 6
    END SELECT
!
      RETURN
!     ------------------------------------------------------------------
901   CALL PRINT_ERROR (16,pgname,localinfo,0)
      STOP 'PEGASE terminated during input file analysis'
      END
! //////////////////////////////////////////////////////////////////////
