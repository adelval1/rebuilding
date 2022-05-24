! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  INPUT   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE open_input_file                                      //
! //                                                                  //
! //  In:     inputname   full name of the input file                 // 
! //                                                                  //
! //  The routine opens the input file as unit 11.                    //
! //                                                                  //
! //  Benoit Bottin, 03/10/96. Adapted to Win95 16/06/97.             //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE open_input_file (inputname)
!
      IMPLICIT NONE
      CHARACTER*80 inputname
      CHARACTER*70 localinfo
      CHARACTER*6 pgname
      INTEGER l,FILE_OPEN
!
      pgname='pegase'
      localinfo(12:70)=inputname(1:59)
      l=FILE_OPEN(inputname,11,1)
      IF (l.NE.1) THEN
          CALL PRINT_ERROR (10,pgname,localinfo,1)
      ENDIF
      RETURN
      END
