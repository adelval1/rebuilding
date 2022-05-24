! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //  PEGASE 4. CHEMCO LIBRARY                                        //
! //                                                                  //
! //  SUBROUTINE CHEMCODEF_DVDA                                       //
! //                                                                  //
! //  Set all initial things needed by the chemco-routines            // 
! //  (arrays, etc ...)                                               //
! //  chemcodef calls for things defined by the subroutine TERMODEF.  // 
! //  It is important that chemcodef be called after termodef!!       //
! //                                                                  //
! //  J. P. Mellado,Paolo Barbante Dec. 98                            //
! //  Modified to include thermal non-equilibrium by                  //
! //  D. Vanden Abeele, April 2000                                    //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////

SUBROUTINE chemcodef_dvda()

!   Definition of global variables
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    USE global_thermo
    USE global_chemco
    USE global_chemco_dvda

!   Definition of local variables
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    IMPLICIT NONE
    CHARACTER*70 :: localinfo
    CHARACTER*8 :: pgname
    CHARACTER*4 :: com
    CHARACTER*80 :: fullcom
    INTEGER :: ierr, flg_stop, i, j,alpha
    REAL(kind=8) :: val
    INTEGER :: transf

!   Definition of functions
!   ^^^^^^^^^^^^^^^^^^^^^^^

    INTEGER :: FILE_EXISTS, FILE_OPEN, FILE_CLOSE

!   Initilization
!   ^^^^^^^^^^^^^

    transf=0

    pgname = 'pegaslib'
    CALL FILL_WITH_BLANKS( localinfo )
    localinfo(12:70) = chemconame(1:59) 
    flg_stop = 0
    com = '    '
    cne_nr = 0

!   Open the chemical data file
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ierr = FILE_EXISTS( chemconame )
    IF( ierr.EQ.0 ) THEN
       CALL PRINT_ERROR( 8, pgname, localinfo, 1 )
    ELSE IF( ierr.EQ.-1 ) THEN
       CALL PRINT_ERROR( 9, pgname, localinfo, 1 )
    END IF

    ierr = FILE_OPEN( chemconame, 40, 1 )
    IF( ierr.NE.1 ) THEN
       CALL PRINT_ERROR( 10, pgname, localinfo, 1 )
    END IF

!   Read the data file. Allocation of the arrays
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    DO WHILE( com.NE.'stop' )
       CALL FILL_WITH_BLANKS( fullcom )
       com = '    '
       READ( 40, 1000, err=900 ) fullcom
       CALL LCASE( fullcom(1:4), com, ierr )
       IF( com.EQ.'nr  ' ) THEN
          READ( 40, 1001, err=900 ) val
          cne_nr = INT(val)
          CALL size_def_dvda 
          CALL read_reactions_dvda ()

!   Added to transform reaction rates in SI units
       ELSE IF(com.EQ.'si  ') THEN
          transf=1

       ELSE IF(com.EQ.'mks ') THEN
          transf=1

       ELSE IF(com.EQ.'cgs ') THEN
          transf=2

       END IF
    END DO 
       
!   Close the file
!   ^^^^^^^^^^^^^^

    ierr = FILE_CLOSE( 40 )
    IF( ierr.EQ.-1 ) THEN
       CALL PRINT_ERROR( 13, pgname, localinfo, 1 )
    END IF

!   Definition of chemical non-equilibrium stoichimetric coefficients
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    STO = CNE_P - CNE_R

!   Transformation of chemical reaction rates in SI units
    IF(transf.EQ.2) THEN
       DO i=1,cne_nr
          alpha=0
          DO j=1,nsp
             alpha=alpha+CNE_R(i,j)
          END DO
          alpha=alpha-1
          FRRC(i,1) = FRRC(i,1)*(1.0d6)**(-alpha)
       END DO

    ELSE IF(transf.EQ.0) THEN
       call PRINT_ERROR(47,pgname,localinfo,1)

    END IF
    

    RETURN

!   Labels
!   ^^^^^^

     900  CALL PRINT_ERROR( 50, pgname, localinfo, 1 )

    1000  FORMAT(a80)
    1001  FORMAT(f12.6)
    1002  FORMAT(3f12.6)

END SUBROUTINE chemcodef_dvda

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  SETUP   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE size_def_dvda                                        //
! //                                                                  //
! //  The subroutine is used to allocate the dynamic arrays used      //
! //  in the chemical non-equilibrium calculation                     //
! //                                                                  //
! //  J. P. Mellado, July 98                                          //
! //  Modified to include thermal non-equilibrium by                  //
! //  D. Vanden Abeele, April 2000                                    //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////

SUBROUTINE size_def_dvda()

!   Definition of global variables
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    USE global_thermo
    USE global_chemco
    USE global_chemco_dvda

!   Definition of local variables
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    IMPLICIT NONE
    INTEGER :: ierr
    CHARACTER*8 :: pgname
    CHARACTER*70 :: localinfo

!   Initialization
!   ^^^^^^^^^^^^^^

    pgname='pegaslib'
    CALL FILL_WITH_BLANKS( localinfo )

    IF( cne_nr.LE.0 ) THEN
        CALL PRINT_ERROR( 51, pgname, localinfo, 1 )
    END IF

!   Allocation
!   ^^^^^^^^^^

    ALLOCATE( CNE_R(1:cne_nr,1:nsp), STAT=ierr )
    IF( ierr.NE.0 ) THEN
       CALL FILL_WITH_BLANKS( localinfo )
       localinfo(12:39) = 'Context: allocation of CNE_R'
       CALL PRINT_ERROR( 43, pgname, localinfo, 1 )
    END IF

    ALLOCATE( CNE_P(1:cne_nr,1:nsp), STAT=ierr )
    IF( ierr.NE.0 ) THEN
       CALL FILL_WITH_BLANKS( localinfo )
       localinfo(12:39) = 'Context: allocation of CNE_P'
       CALL PRINT_ERROR( 43, pgname, localinfo, 1 )
    END IF

    ALLOCATE( STO(1:cne_nr,1:nsp), STAT=ierr )
    IF( ierr.NE.0 ) THEN
       CALL FILL_WITH_BLANKS( localinfo )
       localinfo(12:37) = 'Context: allocation of STO'
       CALL PRINT_ERROR( 43, pgname, localinfo, 1 )
    END IF

    ALLOCATE( FRRC(1:cne_nr,1:3), STAT=ierr )
    IF( ierr.NE.0 ) THEN
        CALL FILL_WITH_BLANKS( localinfo )
        localinfo(12:38) = 'Context: allocation of FRRC'
        CALL PRINT_ERROR( 43, pgname, localinfo, 1 )
    END IF

    ALLOCATE( QREAC(1:cne_nr,1:2), STAT=ierr )
    IF( ierr.NE.0 ) THEN
        CALL FILL_WITH_BLANKS( localinfo )
        localinfo(12:38) = 'Context: allocation of QREAC'
        CALL PRINT_ERROR( 43, pgname, localinfo, 1 )
    END IF

!   Initialization of the allocated arrays
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    CNE_R = 0
    CNE_P = 0
    STO = 0
    FRRC  = 0.0
    QREAC  = 0.0

END SUBROUTINE size_def_dvda

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  SETUP   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE read_reactions_dvda                                  //
! //                                                                  //
! //  The subroutine is used to read the alphanumerical reactions     //
! //  and fill the arrays CNE_R and CNE_P with the stoichiometric     //
! //  coefficients                                                    //
! //                                                                  //
! //  J. P. Mellado, July 98                                          //
! //  Modified to include thermal non-equilibrium by                  //
! //  D. Vanden Abeele, April 2000                                    //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////

SUBROUTINE read_reactions_dvda()

!   Definition of global variables
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    USE global_thermo
    USE global_chemco
    USE global_chemco_dvda

!   Definition of local variables
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    IMPLICIT NONE
    CHARACTER*70 :: localinfo
    CHARACTER*8 :: pgname           
    CHARACTER*1 :: dummy
    CHARACTER*80 :: reaction
    CHARACTER*20 :: symbol, coeff_char
    INTEGER :: i, j, k, l, isp
    INTEGER :: flg_symbol, flg_prod, coeff_value

!   Definition of functions
!   ^^^^^^^^^^^^^^^^^^^^^^^

    INTEGER LENTRIM

!   Initialization
!   ^^^^^^^^^^^^^^

    pgname = 'pegaslib'
    CALL FILL_WITH_BLANKS( localinfo )
    CALL FILL_WITH_BLANKS( reaction )
    CALL FILL_WITH_BLANKS( symbol )
    CALL FILL_WITH_BLANKS( coeff_char )
    flg_symbol = 0
    flg_prod = 0
  
!   Read the cne_nr reactions
!   ^^^^^^^^^^^^^^^^^^^^^^^^^

    DO i = 1, cne_nr
       READ( 40, 1000 ) reaction

!      Read every character of the reaction
!      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

       DO j = 1, 80 
          dummy = reaction(j:j)

!         Build the stoichiometric coefficient
!         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

          IF( ((ICHAR(dummy).GE.48).AND.(ICHAR(dummy).LE.57)).AND.&
         &    (flg_symbol.EQ.0) ) THEN
             k = LENTRIM(  coeff_char )
             coeff_char(k+1:k+1) = dummy
   
!         Build the species symbol
!         ^^^^^^^^^^^^^^^^^^^^^^^^^

          ELSE IF( ( (ICHAR(dummy).GE.65).AND.(ICHAR(dummy).LE.90) ).OR. &
         &        ( (ICHAR(dummy).GE.97).AND.(ICHAR(dummy).LE.122) ).OR.&
         &        ( ((ICHAR(dummy).GE.48).AND.(ICHAR(dummy).LE.57)).AND.&
         &          (flg_symbol.EQ.1) ) ) THEN
             flg_symbol = 1
             k = LENTRIM( symbol )
             symbol(k+1:k+1) = dummy

!         Identificate the species and get the coefficient value
!         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

          ELSE IF( (dummy.EQ.'+').OR.(dummy.EQ.'=') ) THEN

             DO isp = 1, nsp + 1
                IF( isp.EQ.(nsp+1) ) THEN
                   localinfo(12:31) = symbol 
                   CALL PRINT_ERROR( 52, pgname, localinfo, 1 )
                ELSE IF( SPC_SYMBOL(isp).EQ.symbol ) THEN
                   EXIT
                END IF
             END DO

             k = LENTRIM( coeff_char )
             IF( k.EQ.0 ) THEN
                k = 1
                coeff_char(1:1) = '1'
             END IF
             coeff_value = 0
             DO l = k, 1, -1
                coeff_value = coeff_value + &
               &               INT(10**(k-l))*(ICHAR(coeff_char(l:l))-48)
             END DO
             IF( flg_prod.EQ.0 ) THEN
                CNE_R(i,isp) = CNE_R(i,isp) + coeff_value
             ELSE
                CNE_P(i,isp) = CNE_P(i,isp) + coeff_value
             END IF

             flg_symbol = 0
             CALL FILL_WITH_BLANKS( coeff_char )
             CALL FILL_WITH_BLANKS( symbol )

!         Error: inadmissible character
!         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

          ELSE IF( dummy.NE.' ' ) THEN
             CALL PRINT_ERROR( 53, pgname, localinfo, 1 )
          END IF

!         Determine if we are in the LHS or RHS of the reaction
!         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

          IF( dummy.EQ.'=' ) THEN
             IF( flg_prod.EQ.1 ) THEN
                CALL PRINT_ERROR( 53, pgname, localinfo, 1 )
             ELSE
                flg_prod = 1
             END IF
          END IF

       END DO

!      Repeat it for the last species of the reaction
!      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

       DO isp = 1, nsp + 1
          IF( isp.EQ.(nsp+1) ) THEN
             localinfo(12:31) = symbol 
             CALL PRINT_ERROR( 52, pgname, localinfo, 1 )
          ELSE IF( SPC_SYMBOL(isp).EQ.symbol ) THEN
             EXIT
          END IF
       END DO

       k = LENTRIM( coeff_char )
       IF( k.EQ.0 ) THEN
          k = 1
          coeff_char(1:1) = '1'
       END IF
       coeff_value = 0
       DO l = k, 1, -1
          coeff_value = coeff_value + &
         &               INT(10**(k-l))*(ICHAR(coeff_char(l:l))-48)
       END DO
       CNE_P(i,isp) = CNE_P(i,isp) + coeff_value

       flg_symbol = 0
       CALL FILL_WITH_BLANKS( coeff_char )
       CALL FILL_WITH_BLANKS( symbol )

!   Reads the chemical constants for the reaction
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       READ( 40, 1002, err=900 ) ( FRRC(i,j), j=1,3 )
       READ( 40, 1003, err=900 ) ( QREAC(i,j), j=1,2 )

!   End of this reaction
!   ^^^^^^^^^^^^^^^^^^^^
       flg_prod = 0

    END DO

    RETURN

!   Labels
!   ^^^^^^

    900  CALL PRINT_ERROR( 50, pgname, localinfo, 1 )
    1000 FORMAT( a80 ) 
    1002 FORMAT(3f12.6)
    1003 FORMAT(2f12.6)

END SUBROUTINE read_reactions_dvda

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  SETUP   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE chemcostop_dvda                                      //
! //                                                                  //
! //  The subroutine is used to deallocate the dynamic arrays used by //
! //  the thermodynamic library PEGASELIB.                            //
! //                                                                  //
! //  J. P. Mellado, July 98                                          //
! //  Modified to include thermal non-equilibrium by                  //
! //  D. Vanden Abeele, April 2000                                    //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////

SUBROUTINE chemcostop_dvda

    USE global_chemco
    USE global_chemco_dvda 
 
    IMPLICIT NONE

    DEALLOCATE(CNE_R)  
    DEALLOCATE(CNE_P)  
    DEALLOCATE(STO)  
    DEALLOCATE(FRRC)  
    DEALLOCATE(QREAC)  
    
   
END SUBROUTINE chemcostop_dvda

