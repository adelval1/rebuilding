SUBROUTINE READ_TRACOFITS (nsp,traconame)
!
!***********************************************************************
! Read Traco curve-fit coefficients from the correct data-file.
! Only for neutral-particle interactions; charged-charged interactions
! are contained in the code itself.
!***********************************************************************
      USE global_traco
      IMPLICIT NONE
!***********************************************************************
!
      INTEGER :: I, J, K ,ierr, nsp
      INTEGER FILE_EXISTS,FILE_OPEN,FILE_CLOSE
      CHARACTER * 8  :: pgname,empty
      CHARACTER * 30 :: title
      CHARACTER * 70 :: localinfo
      CHARACTER * 80 :: traconame
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 1. Open the file
!
      pgname='pegaslib'
      localinfo(1:11)='           '
      localinfo(12:70)=traconame(1:59)
      ierr=FILE_EXISTS (traconame)
      IF (ierr.EQ.0) THEN
          CALL PRINT_ERROR (8,pgname,localinfo,1)
      ELSEIF (ierr.EQ.-1) THEN
          CALL PRINT_ERROR (9,pgname,localinfo,1)
      ENDIF
      ierr=FILE_OPEN (traconame,14,1)
      IF (ierr.NE.1) THEN
          CALL PRINT_ERROR (10,pgname,localinfo,1)
      ENDIF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! 2. Read the curve-fit coefficients
!
      DO K = 1, number_fits
         READ (14,*,err=99) title
!
         DO I = 1, nneut
            READ (14,*,err=99) title
1000        FORMAT (a30)
!
            DO J = 1, nsp
!
               READ (14,*,err=99) empty, fits (K, I, J, 1), fits (K, I, J, 2), &
              & fits (K, I, J, 3), fits (K, I, J, 4)
1001           FORMAT (a8,4f15.6)
!
            END DO
         END DO
!
      END DO
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 3. Close the file
!
      ierr=FILE_CLOSE (14)
      IF (ierr.EQ.-1) THEN
          CALL PRINT_ERROR (13,pgname,localinfo,1)
      ENDIF
      RETURN
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 4. Error trapping
!
99    CALL PRINT_ERROR (46,pgname,localinfo,1)
!
!***********************************************************************
END SUBROUTINE READ_TRACOFITS

