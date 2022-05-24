SUBROUTINE Tracodef (n)

!***********************************************************************
! Set all initial things needed by the traco-routines (arrays, etc ...)
! No header: a LOT of stuff needs to be set - would become a huge
! parameter list.
!
! WARNING: tracodef calls for stuff defined by the subroutine TERMODEF.
! It is important that tracodef be called after termodef!!
!***********************************************************************
      USE global_thermo
      USE global_traco
      IMPLICIT NONE
!***********************************************************************
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      INTEGER ierr,n
!***********************************************************************
  pgname='pegaslib'
  CALL FILL_WITH_BLANKS (localinfo)
  
! Initialize the order of Sonine polynomials expansion
  norder=n
  IF ((norder.LT.2).OR.(norder.GT.4)) THEN
      CALL PRINT_ERROR (45,pgname,localinfo,0)
      norder=3
  END IF
! Initialize some arrays etc ...
  CALL Allocate_Traco
! Determine number of ions, electrons, molecules etc. Derive
! traco format and put into 'index_array' (passed via module traco).
  CALL Init_Particles
! Read traco curve-fits; reorder to traco format.
  ALLOCATE (fits(1:number_fits, 1:nneut, 1:nsp, 1:4),STAT=ierr)
  IF (ierr.NE.0) THEN
      localinfo(12:51)='Context: allocation of fits             '
      CALL PRINT_ERROR (43,pgname,localinfo,1)
      STOP 'Program ended due to PEGASE transport library error'
  END IF
  fits=0.
  CALL Read_Tracofits (nsp, traconame)

! Reorder molar masses to traco format.
  CALL Order_Data (nsp,index_array, MMOL(1:nsp), MMOL_ord(1:nsp))

! Reorder stoechiometric matrix etc to traco format. 
  CALL Order_Reactions (nr,nsp,index_array, NU, NU_ord)

! Set coefficients needed for evaluation of Devoto's q-{mp} matrix.
  CALL Init_omega_el

!***********************************************************************
END SUBROUTINE Tracodef


SUBROUTINE Allocate_Traco

!***********************************************************************
! Allocate & initialize traco-related stuff.
!***********************************************************************
      USE global_thermo
      USE global_traco
      USE traco_ord_var
      IMPLICIT NONE
      INTEGER ierr
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
!***********************************************************************
      CALL FILL_WITH_BLANKS (localinfo)
      pgname='pegaslib'
!
      C_omega_el = 0.
!
!   Allocation of ordered arrays and index arrays
!   ---------------------------------------------
      ALLOCATE (HFOR_ord(1:nsp),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of HFOR_ord         '
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE transport library error'
      END IF
      ALLOCATE (NU_ord(nr, nsp),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of NU_ord           '
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE transport library error'
      END IF
      ALLOCATE (HREACT(1:nr),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of HREACT           '
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE transport library error'
      END IF
       ALLOCATE (MMOL_ord(1:nsp),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of MMOL_ord         '
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE transport library error'
      END IF
      ALLOCATE (index_array(1:nsp),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of index_array      '
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE transport library error'
      END IF
      ALLOCATE (neutral_index_array(1:nsp),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:53)='Context: allocation of neutral_index_array'
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE transport library error'
      END IF
      neutral_index_array = 0
      index_array = 0
      MMOL_ord = 0.
      HREACT = 0.
      NU_ord = 0
      HFOR_ord = 0.
!
!   Allocation of arrays needed by the traco routines
!   -------------------------------------------------
      ALLOCATE (Omega(1:2, 1:nsp, 1:nsp),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of omega            '
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE transport library error'
      END IF
      ALLOCATE (BStar(1:nsp, 1:nsp),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of Bstar            '
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE transport library error'
      END IF
      ALLOCATE (CStar(1:nsp, 1:nsp),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of Cstar            '
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE transport library error'
      END IF
      ALLOCATE (Delta(1:2, 1:nsp, 1:nsp),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of delta            '
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE transport library error'
      END IF
      ALLOCATE (visc_factor(1:nsp),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of visc_factor      '
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE transport library error'
      END IF
      ALLOCATE (lambda_factor(1:nsp),STAT=ierr)
      IF (ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of lambda_factor    '
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE transport library error'
      END IF
      Omega = 0.
      BStar = 0.
      CStar = 0.
      Delta = 0.
      visc_factor = 0.
      lambda_factor = 0.


!   Allocation of ordered arrays for mole fractions, number densities and
!   spcific heats

      ALLOCATE(x_ord(1:nsp),STAT=ierr)
      IF(ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of x_ord            '
          CALL PRINT_ERROR(43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE transport library error'
      END IF
      ALLOCATE(n_ord(1:nsp),STAT=ierr)
      IF(ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of n_ord            '
          CALL PRINT_ERROR(43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE transport library error'
      END IF
      ALLOCATE(cpr_ord(1:nsp),STAT=ierr)
      IF(ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of cpr_ord          '
          CALL PRINT_ERROR(43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE transport library error'
      END IF
      ALLOCATE(cpv_ord(1:nsp),STAT=ierr)
      IF(ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of cpv_ord          '
          CALL PRINT_ERROR(43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE transport library error'
      END IF
      ALLOCATE(cpe_ord(1:nsp),STAT=ierr)
      IF(ierr.NE.0) THEN
          localinfo(12:51)='Context: allocation of cpe_ord          '
          CALL PRINT_ERROR(43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE transport library error'
      END IF
      x_ord=0.
      n_ord=0.
      cpr_ord=0.
      cpv_ord=0.
      cpe_ord=0.


!***********************************************************************
END SUBROUTINE Allocate_Traco

SUBROUTINE Tracostop

!***********************************************************************
! Deallocate all arrays that have been allocated in tracodef
!***********************************************************************
      USE global_thermo
      USE global_traco
      USE traco_ord_var
      IMPLICIT NONE
!
      DEALLOCATE ( Omega , BStar , CStar )
      DEALLOCATE ( Delta )
      DEALLOCATE ( visc_factor , lambda_factor )
      DEALLOCATE ( index_array , neutral_index_array )
      DEALLOCATE ( NU_ord )
      DEALLOCATE ( MMOL_ord  )
      DEALLOCATE ( HFOR_ord , HREACT )
      DEALLOCATE ( fits )
      DEALLOCATE (x_ord)
      DEALLOCATE (n_ord)
      DEALLOCATE (cpr_ord)
      DEALLOCATE (cpv_ord)
      DEALLOCATE (cpe_ord)
!
END SUBROUTINE Tracostop
!

