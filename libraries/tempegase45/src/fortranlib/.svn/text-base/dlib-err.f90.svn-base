!     //////////////////////////////////////////////////////////////////
!     ///                                                            ///
!     ///   "DUALE" FORTRAN LIBRARY                                  ///
!     ///   =======================                                  ///
!     ///                                                            ///
!     ///   SUBROUTINE print_error                                   ///
!     ///                                                            ///
!     ///  The routine prints error messages from all subroutines    ///
!     ///                                                            ///
!     ///  Input: errnb    number of the error message               ///
!     ///         pgname   the program name                          ///
!     ///         localinfo local information at error code          ///
!     ///         istop    intended action: 0 continue (warning only)///
!     ///                                   1 stop (fatal error)     ///
!     ///                                  -1 prompt for continuing  ///
!     ///                                                            ///
!     ///   B. Bottin, 13/06/97                                      ///
!     ///                                                            ///
!     //////////////////////////////////////////////////////////////////
      SUBROUTINE print_error (errnb,pgname,localinfo,istop)
!      
      IMPLICIT NONE
      LOGICAL iex
      INTEGER errnb, istop, lpn, i
      CHARACTER*(*) pgname
      CHARACTER*1 ans
      CHARACTER*12 errfile
      CHARACTER*70 msg,localinfo,line1,line2
      INTEGER LENTRIM,FILE_EXISTS,FILE_OPEN,FILE_CLOSE
!
!     Step 1. Build the errors file name
!     ----------------------------------
      lpn=LENTRIM(pgname)
      IF (lpn.GT.8) THEN
          errfile(1:8)=pgname(1:8)
          i=9
      ELSE
          errfile(1:lpn)=pgname(1:lpn)
          i=lpn+1
      ENDIF
      errfile(i:i+3)='.err'
      CALL LCASE (errfile,errfile,i)
!
!     Step 2. Check for existence and open, the errors file
!     -----------------------------------------------------
      IF (FILE_EXISTS(errfile).NE.1) THEN
          WRITE (*,1004) errfile,errnb
          STOP 'Fatal error encountered'
      ENDIF
!
      IF (FILE_OPEN(errfile,99,1).NE.1) THEN
          WRITE (*,1005) errfile,errnb
          STOP 'Fatal error encountered'
      ENDIF
!
!     Step 3. Read up to the current error number and write message
!     -------------------------------------------------------------
      DO i=1,errnb
          READ (99,1000,END=999) line1
          READ (99,1000,END=999) line2
      END DO
!
      WRITE (*,1001) line1
      WRITE (*,1001) line2
      WRITE (*,1001) localinfo
!
!     Step 4. Check if a log file (unit 66) is open, to print the error
!     -----------------------------------------------------------------
      INQUIRE (unit=66,opened=iex,err=9999)
      IF (iex) THEN
          WRITE (66,1001) line1
          WRITE (66,1001) line2
          WRITE (66,1001) localinfo
          WRITE (66,*)
      END IF
!
!     Step 5. Close the errors file and stop or exit depending on istop
!     -----------------------------------------------------------------
      IF (FILE_CLOSE(99).NE.1) THEN
          STOP 'Cannot close the error file'
      ENDIF
!
      IF (istop.EQ.-1) THEN
          WRITE (*,1002)
          READ (*,1003) ans
          IF ((ans.EQ.'y').OR.(ans.EQ.'Y')) THEN
              istop=1
          ENDIF
      ENDIF
!
      IF (istop.EQ.1) THEN
          msg(2:lpn+1)=pgname
          msg(lpn+2:lpn+27)=' terminated on fatal error'
          WRITE (*,1000) msg
          STOP '======================================================='
      ENDIF
      RETURN
!     ------------------------------------------------------------------
1000  FORMAT (a70)
1001  FORMAT (' ',a70)
1002  FORMAT (' Do you want to abort execution of this program (y/n) ?')
1003  FORMAT (a1)
1004  FORMAT (' Cannot find the errors file ',a12,' for error nb. ',i5)
1005  FORMAT (' Cannot open the errors file ',a12,' for error nb. ',i5)
1006  FORMAT (' Cannot find error nb. ',i5)
1007  FORMAT (' Error while inquiring for existence of log unit 66')
!     ------------------------------------------------------------------
999   WRITE (*,1006) errnb
      STOP 'FORTRANLIB Fatal error encountered'
9999  WRITE (*,1007)
      STOP 'FORTRANLIB Fatal error encountered'
      END
! //////////////////////////////////////////////////////////////////////
