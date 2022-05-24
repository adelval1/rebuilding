MODULE global_traco
!               
      IMPLICIT NONE
!
!   Sizing parameters
!   -----------------
      INTEGER,SAVE :: nelectrons, nions, nneut, nmolecn, nmoleci, norder
!
!   Curve fit parameters
!   --------------------
      INTEGER,PARAMETER :: number_fits=3
      REAL(kind=8),ALLOCATABLE,SAVE :: fits (:, :, :, :)
!
!   Parameters required for the computation of transport properties
!   ---------------------------------------------------------------
      REAL(kind=8),ALLOCATABLE,SAVE :: Omega (:, :, :), BStar (:, :), CStar (:, :) 
      REAL(kind=8),ALLOCATABLE,SAVE :: Delta (:, :, :)
      REAL(kind=8),ALLOCATABLE,SAVE :: visc_factor (:), lambda_factor (:)

      REAL(kind=8),SAVE :: C_omega_el (1:2, 0:3, 0:3, 1:7)
      REAL(kind=8),SAVE :: a_average_lambda, a_average_visc
      REAL(kind=8),SAVE :: qel (0:3, 0:3)
!
!   Ordered arrays needed by the traco routines
!   -------------------------------------------
      INTEGER, ALLOCATABLE,SAVE     :: index_array (:), neutral_index_array (:)
      INTEGER, ALLOCATABLE,SAVE     :: NU_ord (:, :)
      REAL(kind=8), ALLOCATABLE,SAVE :: MMOL_ord (:) 
      REAL(kind=8), ALLOCATABLE,SAVE :: HFOR_ord (:), HREACT (:)
!
END MODULE global_traco


MODULE traco_ord_var

   IMPLICIT NONE

!  Module contains arrays necessary to compute transport properties in
!  traco order

   REAL(kind=8), ALLOCATABLE, SAVE :: x_ord(:),n_ord(:)
   REAL(kind=8), ALLOCATABLE, SAVE :: cpr_ord(:),cpv_ord(:),cpe_ord(:)


END MODULE traco_ord_var
