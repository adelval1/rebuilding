! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  INPUT   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE close_input_file                                     //
! //                                                                  //
! //  The routine closes the input file (unit 11).                    //
! //                                                                  //
! //  Benoit Bottin, 03/10/96. Adapted to Win95 16/06/97.             //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE close_input_file  
!
      IMPLICIT NONE
      CHARACTER*70 localinfo
      CHARACTER*6 pgname
      INTEGER i,FILE_CLOSE
!
      pgname='pegase'
      localinfo(12:21)='Unit = 11 '
      i=FILE_CLOSE (11)
      IF (i.EQ.-1) THEN
          CALL PRINT_ERROR (13,pgname,localinfo,-1)
      ENDIF
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
