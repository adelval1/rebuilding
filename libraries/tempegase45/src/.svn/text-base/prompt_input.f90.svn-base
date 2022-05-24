! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  INPUT   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE prompt_input                                         //
! //                                                                  //
! //  Input:  in_dir      working directory for the input file        //    
! //                                                                  //
! //  Outs:   inputname   full name of the input file                 // 
! //          outname     output root without path                    //
! //                                                                  //
! //  The routine prompts for the input file to use for the           //
! //  calculation. It concatenates the input name and the input       //
! //  directory and returns the value. If the file does not exist,    //
! //  an error is generated.                                          //
! //                                                                  //
! //  The default extension if not provided is .IN                    //
! //                                                                  //
! //  The routine also returns the root name (without the extension)  //
! //  without the path, to be able to redirect output files.          //
! //                                                                  //
! //  Benoit Bottin, 03/10/96. Adapted for Win95 16/06/97             //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE prompt_input_file (in_dir,inputname,outname)
!
      IMPLICIT NONE
      CHARACTER*80 in_dir,filename,inputname,outname
      CHARACTER*70 localinfo
      CHARACTER*6 pgname
      INTEGER d,l,iex,FILE_EXISTS,LENTRIM
!     
      CALL FILL_WITH_BLANKS (outname)
      pgname='pegase'
      localinfo(1:1)=' '
!
!     Prompting on the screen
!     ^^^^^^^^^^^^^^^^^^^^^^^
      WRITE (*,3001)
      READ (*,3002) filename
      WRITE (*,*)
!
!     Getting the output root
!     ^^^^^^^^^^^^^^^^^^^^^^^
      d=INDEX(filename,'.')
      IF (d.EQ.0) THEN
          l=LENTRIM(filename)
          IF (l.EQ.0) THEN
              CALL PRINT_ERROR (7,pgname,localinfo,1)
          ELSE
              filename(l+1:l+3)='.in'
          ENDIF
      ENDIF
      d=INDEX(filename,'.')
      outname(1:d-1)=filename(1:d-1)
!
!     Form the full input name
!     ^^^^^^^^^^^^^^^^^^^^^^^^
      CALL CONCATENATE(in_dir,filename,inputname,iex)
!
!     Check the existence of the input file
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      iex=FILE_EXISTS(inputname)
      IF (iex.EQ.0) THEN
          localinfo(12:70)=inputname(1:59)
          CALL PRINT_ERROR (8,pgname,localinfo,1)
      ELSEIF (iex.EQ.-1) THEN
          localinfo(12:70)=inputname(1:59)
          CALL PRINT_ERROR (9,pgname,localinfo,1)
      ENDIF
      RETURN
!     ------------------------------------------------------------------                  
!     Formats
!     ^^^^^^^
3001  FORMAT (/,/,' Please provide your input file name:')
3002  FORMAT (a80)
!     ------------------------------------------------------------------      
      END
! //////////////////////////////////////////////////////////////////////
