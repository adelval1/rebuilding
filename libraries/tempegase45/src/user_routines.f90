! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //                                                                  //
! //               P E G A S E   4.   U S E R   C O D E               //
! //                                                                  //
! //                   M A I N    S U B R O U T I N E                 //
! //                                                                  //
! //                                                                  //
! //  Use this routine to perform all the necessary side computations //
! //  you would like to do in PEGASE.                                 //
! //                                                                  //
! //  The call to this routine in the main code is located within the //
! //  nested DO loops, after the output files have been printed.      //
! //                                                                  //
! //  Check your user manual to know the use of each variable and     //
! //  function.                                                       //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //  Benoit Bottin, 30/10/96. Modified F90 format 7/10/97.           //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE user_compute
!
      USE global_pegase
      USE input_pegase
      USE global_thermo
      USE interf_thermo
      USE interf_traco
      IMPLICIT NONE
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
!     PART 1: DEFINITION OF VARIABLES
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
    REAL (kind=8), ALLOCATABLE :: h(:)
    REAL (kind=8) tneq(1:4), enth(1:8), T
    INTEGER i, j, nunit

!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
!     END PART 1
!
!     PART 2: USER-DEFINED PROGRAMMING
!
!     Insert here all the necessary code you want to include
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
!     This version of user routines is used to generate, species by
!     species, data with a variable number of electronic levels.
!     13-species air is considered. flg_user is used to pass the max.
!     number of electronic levels to include.
    ALLOCATE (h(1:flg_user))     
    TNEQ = 0.0d0
    T = array3(3)
    DO i = 1, nsp
        DO j = 1, flg_user
            nelvl(i) = j
            CALL species_enthalpy (i,T,TNEQ,flg_anha,flg_mode,0,1,0,enth)
            h(j) = enth(1)/8.314d0/T
        END DO
        nunit = 600 + i
        WRITE (UNIT=nunit,FMT=100) T,(h(j),j=1,flg_user)
100     FORMAT (20e15.8)
    END DO
    DEALLOCATE (h)
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
!     END PART 2
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------

END SUBROUTINE USER_COMPUTE
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //                                                                  //
! //               P E G A S E   4.   U S E R   C O D E               //
! //                                                                  //
! //                   O P E N    S U B R O U T I N E                 //
! //                                                                  //
! //                                                                  //
! //  Use this routine to perform all the necessary open instructions //
! //  you would like to do in PEGASE.                                 //
! //                                                                  //
! //  The call to this routine in the main code is located before the //
! //  computations really start, just after the case is read in the   //
! //  input file.                                                     //
! //                                                                  //
! //  Check your user manual to know the use of each unit number to   //
! //  avoid redundancy.                                               //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //  Benoit Bottin, 30/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE user_open (in_dir,out_dir,thermo_dir,traco_dir,cntr,&
     &flg_user)
      USE global_thermo
      IMPLICIT NONE
!
!     Definition of variables passed from the main code
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      CHARACTER*80 in_dir,out_dir,thermo_dir,traco_dir
      INTEGER cntr,flg_user
!
!     ------------------------------------------------------------------
!     Fill here your variable declarations
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    CHARACTER*80 filename
    CHARACTER*20 text
    CHARACTER*5 bits (1:20)
    INTEGER i, j, k, l, nunit 
!     ------------------------------------------------------------------
!     Fill here your program statements
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    IF (cntr.EQ.1) THEN
    DO i = 1, 20
        WRITE (unit=bits(i),fmt=10) i
    END DO
10  FORMAT ('"',i2,'" ')
!
    DO i = 1, nsp
        filename = REPEAT (' ',80)
        l = LEN_TRIM (out_dir)
        filename(1:l) = out_dir(1:l)
        k = LEN_TRIM (SPC_SYMBOL(i))
        filename(l+1:l+k) = SPC_SYMBOL(i)(1:k)
        filename(l+k+1:l+k+4) = '.dat'
        nunit = 600 + i
        OPEN (file=filename,unit=nunit)
        WRITE (unit=nunit,FMT=100)
        WRITE (unit=nunit,FMT=101) (bits(j),j=1,flg_user)
    END DO
100 FORMAT (' TITLE = "Effect of variation of electronic levels" ')
101 FORMAT (' VARIABLES = "T" ',20a5)
102 FORMAT (' ZONE F=POINT T=" ',a20,'"')
    END IF
!
    DO i = 1, nsp
        nunit = 600 + i
        text = SPC_SYMBOL(i)
        SELECT CASE (cntr)
        CASE (1)
            text(6:20) = '   (RRHO)      '
            WRITE (unit=nunit,FMT=102) text
        CASE (3)
            text(6:20) = '   (RVC)       '
            WRITE (unit=nunit,FMT=102) text
        CASE (5)
            text(6:20) = '   (RVEC)      '
            WRITE (unit=nunit,FMT=102) text
        END SELECT
    END DO
!     ------------------------------------------------------------------
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //                                                                  //
! //               P E G A S E   4.   U S E R   C O D E               //
! //                                                                  //
! //                  C L O S E    S U B R O U T I N E                //
! //                                                                  //
! //                                                                  //
! //  Use this routine to perform all the necessary close operations  //
! //  you would like to do in PEGASE.                                 //
! //                                                                  //
! //  The call to this routine in the main code is located after the  //
! //  nested DO loops, before the input file is read again.           //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //  Benoit Bottin, 30/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE user_close (flg_user,cntr)
      USE global_thermo
      IMPLICIT NONE 
      INTEGER flg_user,cntr
!
!     ------------------------------------------------------------------
!     Fill here your variable declarations
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    INTEGER i, nunit

!     ------------------------------------------------------------------
!     Fill here your program statements
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    IF (cntr.EQ.6) THEN
    DO i = 1, nsp
        nunit = 600 + i
        CLOSE (unit=nunit)
    END DO
    END IF
!     ------------------------------------------------------------------
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
