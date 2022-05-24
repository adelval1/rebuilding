! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE species_cp                                           //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          t           equilibrium temperature (K)                 //
! //          tneq        nonequilibrium temperatures (K)             //
! //                                                                  //
! //  Flags:  flg_anha    anharmonicity corrections flag              //
! //          mode        defines a p,T or rho,T input                //
! //          flg_neq     if 1, considers multiple temperaturess      //
! //          flg_stop    stops if errors in subroutines              //
! //          flg_termo   statistical or reftable computation         //
! //                                                                  //
! //  Output: cp(8)      array of internal energies                   //
! //              1           complete (except formation enthalpy)    //
! //              2           translation                             //
! //              3           rotation (anha=0)                       //
! //              4           vibration (anha=0)                      //
! //              5           electronic (anha=0)                     //
! //              6           coupled rotation-vibration (anha=1)     //
! //              7           coupled internal (anha=2)               //
! //              8           formation enthalpy/enthalpy             //
! //                                                                  //
! //  This subroutine provides the internal energies of the species   //
! //  isp according to the flags specified.                           //
! //                                                                  //
! //  Benoit Bottin, 11/10/96. Modified and cleaned 25/7/97.          //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE  species_cp (isp,T,TNEQ,&
      flg_anha,mode,flg_neq,flg_stop,flg_termo,cp)
!
      USE global_thermo
      IMPLICIT NONE
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      INTEGER isp,flg_anha,mode,flg_neq,flg_stop,flg_termo,istop,i,atom
      REAL(kind=8) T,TNEQ(4),cp(8),ttra,trot,tvib,tele
      REAL(kind=8) CP_TRANS_X,CP_TRANS_Y
      REAL(kind=8) CV_ROT_X,CV_ROT_Y
      REAL(kind=8) CV_VIB_X,CV_VIB_Y
      REAL(kind=8) CV_ELEC_X,CV_ELEC_Y
      REAL(kind=8) CV_RV_X,CV_RV_Y
      REAL(kind=8) CV_INT_X,CV_INT_Y
      REAL(kind=8) CP_REFTABLE_X,CP_REFTABLE_Y
!
!     Checks of correct variable values and initialize
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      pgname='pegaslib'
      istop=0
      IF (flg_neq.EQ.0) THEN
          IF (T.LE.0) THEN
              WRITE (localinfo,101) T
101           FORMAT (11x,'Routine: species_cp. Value:',d15.8)
              CALL PRINT_ERROR (24,pgname,localinfo,0)
              istop=1
          ENDIF
          ttra=T
          trot=T
          tvib=T
          tele=T      
      ELSE
          DO i=1,4
              IF (TNEQ(i).LE.0) THEN
                  WRITE (localinfo,102) i,TNEQ(i)
102   FORMAT (11x,'Routine: species_cp. Index:',i3,' Value:',d15.8)
                  CALL PRINT_ERROR (27,pgname,localinfo,0)
                  istop=1
              ENDIF
          END DO
          ttra=TNEQ(1)
          trot=TNEQ(2)
          tvib=TNEQ(3)
          tele=TNEQ(4)
      ENDIF
      DO i=1,8
          cp(i)=0.0d0
      END DO
      IF (SIG(isp).EQ.0) THEN
          atom=1
      ELSE
          atom=0
      ENDIF
!
!     Main part: the first test bears on p,T or rho,T operation.
!     The second test bears on statisvical physics or reference tables.
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (mode.EQ.1) THEN
          IF (flg_termo.EQ.1) THEN
              IF (flg_neq.NE.0) THEN
                  localinfo(1:31)='           Routine: species_cp.'
                  CALL PRINT_ERROR (28,pgname,localinfo,0)
                  t=ttra
              ENDIF
              cp(1)= CP_REFTABLE_X (isp,t,flg_stop)
          ELSE
              cp(2)=CP_TRANS_X(isp,ttra)
              IF (flg_anha.EQ.0) THEN
                  IF (atom.EQ.0) THEN
                      cp(3)=CV_ROT_X(isp,trot)
                      cp(4)=CV_VIB_X(isp,tvib)
                  ENDIF
                  cp(5)=CV_ELEC_X(isp,tele)
                  cp(1)=cp(2)+cp(3)+cp(4)+cp(5)
              ELSEIF (flg_anha.EQ.1) THEN
                  cp(5)=CV_ELEC_X(isp,tele)
                  IF (atom.EQ.0) THEN
                      cp(6)=CV_RV_X(isp,trot,tvib)
                  ENDIF
                  cp(1)=cp(2)+cp(5)+cp(6)
              ELSEIF (flg_anha.EQ.2) THEN
                  IF (atom.EQ.0) THEN
                      cp(7)=CV_INT_X(isp,trot,tvib,tele)
                  ELSE
                      cp(7)=CV_ELEC_X(isp,tele)
                  ENDIF
                  cp(1)=cp(2)+cp(7)
              ELSE
                  localinfo(1:31)='           Routine: species_cp.'
                  CALL PRINT_ERROR (29,pgname,localinfo,0)
                  istop=1
              ENDIF
          ENDIF
      ELSE
          IF (flg_termo.EQ.1) THEN
              IF (flg_neq.NE.0) THEN
                  localinfo(1:31)='           Routine: species_cp.'
                  CALL PRINT_ERROR (28,pgname,localinfo,0)
                  t=ttra
              ENDIF
              cp(1)= CP_REFTABLE_Y (isp,t,flg_stop)
          ELSE
              cp(2)=CP_TRANS_Y(isp,ttra)
              IF (flg_anha.EQ.0) THEN
                  IF (atom.EQ.0) THEN
                      cp(3)=CV_ROT_Y(isp,trot)
                      cp(4)=CV_VIB_Y(isp,tvib)
                  ENDIF
                  cp(5)=CV_ELEC_Y(isp,tele)
                  cp(1)=cp(2)+cp(3)+cp(4)+cp(5)
              ELSEIF (flg_anha.EQ.1) THEN
                  cp(5)=CV_ELEC_Y(isp,tele)
                  IF (atom.EQ.0) THEN
                      cp(6)=CV_RV_Y(isp,trot,tvib)
                  ENDIF
                  cp(1)=cp(2)+cp(5)+cp(6)
              ELSEIF (flg_anha.EQ.2) THEN
                  IF (atom.EQ.0) THEN
                      cp(7)=CV_INT_Y(isp,trot,tvib,tele)
                  ELSE
                      cp(7)=CV_ELEC_Y(isp,tele)
                  ENDIF
                  cp(1)=cp(2)+cp(7)
              ELSE
                  localinfo(1:31)='           Routine: species_cp.'
                  CALL PRINT_ERROR (29,pgname,localinfo,0)
                  istop=1
              ENDIF
          ENDIF
      ENDIF
!
!     Error handling
!     ^^^^^^^^^^^^^^
      IF ((istop.EQ.1).AND.(flg_stop.EQ.1)) THEN
          STOP 'Program ended due to PEGASE thermodynamic library error'
      ENDIF
!     ------------------------------------------------------------------
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS cp_trans_x AND cp_trans_y                             //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          tt          translation temperature (K)                 //
! //                                                                  //
! //  This function returns the translation  enthalpy per             //
! //  unit mole (p,T input) or unit mass (rho,T input).               //
! //                                                                  //
! //  Benoit Bottin, 11/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION cp_trans_x (isp,tt)
      
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp,i
      REAL(kind=8) :: tt,cp
      
      cp=2.5d0*RUNIV
      cp_trans_x=cp
!
!     isp and tt are passed but not used (keep the same look for all calls)
!     the next statement removes the compile-time error
      i=isp
      cp=tt      

      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION cp_trans_y (isp,tt)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: tt,cp
      REAL(kind=8) :: CP_TRANS_X
!      
      cp=CP_TRANS_X (isp,tt)
      cp_trans_y=cp/MMOL(isp)

      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS cp_reftable_x and cp_reftable_y                       //
! //                                                                  //
! //  Input:  T           temperature of the species (K)              //    
! //          isp         number of the concerned species             //
! //          flg_stop    stops if errors in the subroutine           //
! //                                                                  //
! //  This function returns the specific heat of the species from the //
! //  reference tables such as Gurvich or Janaf (mole or mass input)  //
! //                                                                  //
! //  Benoit Bottin, 11/12/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION cp_reftable_x (isp,T,flg_stop)
!
      USE global_thermo
      IMPLICIT NONE
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      INTEGER isp,ilow,ihigh,istop,flg_stop
      REAL(kind=8) t,x,x1,x2,y,y1,y2
!
!     Find low/high index of temperature
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      pgname='pegaslib'
      ilow=0
10    CONTINUE
          IF (REFTABLE(isp,1,ilow+1).LT.t) THEN
              ilow=ilow+1
              IF (ilow+1.GT.IDATA(isp)) THEN
                  WRITE (localinfo,102) t
102               FORMAT (11x,'Routine: species_cp. Value:',d15.8)
                  CALL PRINT_ERROR (31,pgname,localinfo,0)
                  istop=1
              ELSE
                  GOTO 10
              ENDIF
          ENDIF
      ihigh=ilow+1
      IF (ilow.EQ.0) THEN
          WRITE (localinfo,101) t
101       FORMAT (11x,'Routine: species_cp. Value:',d15.8)
          CALL PRINT_ERROR (30,pgname,localinfo,0)
          istop=1
      ENDIF
      IF ((flg_stop.EQ.1).AND.(istop.EQ.1)) THEN
          STOP 'Program terminated on PEGASE THERMO LIBRARY error'
      ENDIF
!
!     Perform the interpolation
!     ^^^^^^^^^^^^^^^^^^^^^^^^^
      x=t
      x1=REFTABLE(isp,1,ilow)
      x2=REFTABLE(isp,1,ihigh)
      y1=REFTABLE(isp,2,ilow)
      y2=REFTABLE(isp,2,ihigh)
      y=(y2-y1)/(x2-x1)*(x-x1)+y1
!
      cp_reftable_x=y
      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION cp_reftable_y (isp,T,flg_stop)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER isp,flg_stop
      REAL(kind=8) t,cp,CP_REFTABLE_X
!
      cp=CP_REFTABLE_X (isp,T,flg_stop)
      cp_reftable_y=cp/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
