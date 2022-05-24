! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! // PEGASE 4. CHEMCO LIBRARY                                         //
! //                                                                  //
! // MODULE global_chemco                                             //
! //                                                                  //
! // Definition of the global variables related to the chemical       //
! // non-equilibrium calculation either in pegase or in other codes.  //
! // It must be used with the module global_termo.                    //
! //                                                                  //
! // J. P. Mellado, July 98                                           //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////

MODULE global_chemco

    IMPLICIT NONE

!   parameters for the sizing of the arrays 
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!   cne_nr     number of reactions

    INTEGER, SAVE :: cne_nr

!   array for the chemical non-equilibrium reactants stoichiometric coefficients
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!   CNE_R      reactives array        (cne_nr,nsp)

    INTEGER,ALLOCATABLE,SAVE :: CNE_R(:,:)

!   array for the chemical non-equilibrium products stoichiometric coefficients
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!   CNE_P      products array         (cne_nr,nsp)

    INTEGER,ALLOCATABLE,SAVE :: CNE_P(:,:)

!   array for the total stoichiometric coefficients
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!   STO        CNE_P-CNE_R            (cne_nr,nsp)

    INTEGER,ALLOCATABLE,SAVE :: STO(:,:) 

!   array for the forward reaction rate coefficients, Arrhenius law
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!   FRRC       coefficients array     (cne_nr,3)

    REAL(kind=8),ALLOCATABLE,SAVE :: FRRC(:,:)

END MODULE
