! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE species_helmholtz                                    //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          var         value of pressure (Pa) or density (kg/m3)   //
! //          t           equilibrium temperature (K)                 //
! //          tneq        nonequilibrium temperatures (K)             //
! //                                                                  //
! //  Flags:  flg_anha    anharmonicity corrections flag              //
! //          mode        defines a p,T or rho,T input                //
! //          flg_neq     if 1, considers multiple temperaturess      //
! //          flg_stop    stops if errors in subroutines              //
! //          flg_termo   statistical or reftable computation         //
! //                                                                  //
! //  Output: helm(8)     array of helmholtz free energy              //
! //                                                                  //
! //  This subroutine provides the helmholtz free energy of species   //
! //  isp according to the flags specified.                           //
! //                                                                  //
! //  Benoit Bottin, 16/10/96. Modified and cleaned 25/7/97.          //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE species_helmholtz (isp,var,T,TNEQ,&
     &flg_anha,mode,flg_neq,flg_stop,flg_termo,helm)
!
      USE global_thermo
      IMPLICIT NONE
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      INTEGER isp,flg_anha,mode,flg_neq,flg_stop,flg_termo,istop,i,atom
      REAL(kind=8) var,T,TNEQ(4),helm(8),ttra,trot,tvib,tele
      REAL(kind=8) HELMHOLTZ_TRANS_X,HELMHOLTZ_TRANS_Y
      REAL(kind=8) GIBBS_ROT_X,GIBBS_ROT_Y
      REAL(kind=8) GIBBS_VIB_X,GIBBS_VIB_Y
      REAL(kind=8) GIBBS_ELEC_X,GIBBS_ELEC_Y
      REAL(kind=8) GIBBS_RV_X,GIBBS_RV_Y
      REAL(kind=8) GIBBS_INT_X,GIBBS_INT_Y
      REAL(kind=8) GIBBS_REFTABLE_X,GIBBS_REFTABLE_Y
      REAL(kind=8) MOL_TO_MASS
!
!     Checks of correct variable values and initialize
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      pgname='pegaslib'
      istop=0
      IF (flg_neq.EQ.0) THEN
          IF (T.LE.0) THEN
              WRITE (localinfo,101) T
101           FORMAT (11x,'Routine: species_helmholtz. Value:',d15.8)
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
102   FORMAT(11x,'Routine: species_helmholtz.Index:',i3,' Value:',d15.8)
                  CALL PRINT_ERROR (27,pgname,localinfo,0)
                  istop=1
              ENDIF
          END DO
          ttra=TNEQ(1)
          trot=TNEQ(2)
          tvib=TNEQ(3)
          tele=TNEQ(4)
      ENDIF
      IF (var.LE.0) THEN
          WRITE (localinfo,103) var
103       FORMAT (11x,'Routine: species_helmholtz. Value:',d15.8)
          CALL PRINT_ERROR (24+mode,pgname,localinfo,0)
          istop=1
      ENDIF
      DO i=1,8
          helm(i)=0.0d0
      END DO
      IF (SIG(isp).EQ.0) THEN
          atom=1
      ELSE
          atom=0
      ENDIF
!
!     Main part: the first test bears on p,T or rho,T operation.
!     The second test bears on statistical physics or reference tables.
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (mode.EQ.1) THEN
          IF (flg_termo.EQ.1) THEN
              IF (flg_neq.NE.0) THEN
              localinfo(1:38)='           Routine: species_helmholtz.'
                  CALL PRINT_ERROR (28,pgname,localinfo,0)
                  t=ttra
              ENDIF
              helm(1)=GIBBS_REFTABLE_X(isp,var,t,flg_stop)-RUNIV*T
              helm(8)=HFOR(isp)
          ELSE
              helm(2)=HELMHOLTZ_TRANS_X(isp,var,ttra)
              helm(8)=HFOR(isp)
              IF (flg_anha.EQ.0) THEN
                  IF (atom.EQ.0) THEN
                      helm(3)=GIBBS_ROT_X(isp,trot)
                      helm(4)=GIBBS_VIB_X(isp,tvib)
                  ENDIF
                  helm(5)=GIBBS_ELEC_X(isp,tele)
                  helm(1)=helm(2)+helm(3)+helm(4)+helm(5)+helm(8)
              ELSEIF (flg_anha.EQ.1) THEN
                  helm(5)=GIBBS_ELEC_X(isp,tele)
                  IF (atom.EQ.0) THEN
                      helm(6)=GIBBS_RV_X(isp,trot,tvib)
                  ENDIF
                  helm(1)=helm(2)+helm(5)+helm(6)+helm(8)
              ELSEIF (flg_anha.EQ.2) THEN
                  IF (atom.EQ.0) THEN
                      helm(7)=GIBBS_INT_X(isp,trot,tvib,tele)
                  ELSE
                      helm(7)=GIBBS_ELEC_X(isp,tele)
                  ENDIF
                  helm(1)=helm(2)+helm(7)+helm(8)
              ELSE
              localinfo(1:38)='           Routine: species_helmholtz.'
                  CALL PRINT_ERROR (29,pgname,localinfo,0)
                  istop=1
              ENDIF
          ENDIF
      ELSE
          IF (flg_termo.EQ.1) THEN
             IF (flg_neq.NE.0) THEN
              localinfo(1:38)='           Routine: species_helmholtz.'
                  CALL PRINT_ERROR (28,pgname,localinfo,0)
                  t=ttra
             ENDIF
             helm(1)=GIBBS_REFTABLE_Y(isp,var,t,flg_stop)-RSP(isp)*T
             helm(8)=MOL_TO_MASS(HFOR(isp),isp,0.0d0,flg_stop)
          ELSE
              helm(2)=HELMHOLTZ_TRANS_Y(isp,var,ttra)
              helm(8)=MOL_TO_MASS(HFOR(isp),isp,0.0d0,flg_stop)
              IF (flg_anha.EQ.0) THEN
                  IF (atom.EQ.0) THEN
                      helm(3)=GIBBS_ROT_Y(isp,trot)
                      helm(4)=GIBBS_VIB_Y(isp,tvib)
                  ENDIF
                  helm(5)=GIBBS_ELEC_Y(isp,tele)
                  helm(1)=helm(2)+helm(3)+helm(4)+helm(5)+helm(8)
              ELSEIF (flg_anha.EQ.1) THEN
                  helm(5)=GIBBS_ELEC_Y(isp,tele)
                  IF (atom.EQ.0) THEN
                      helm(6)=GIBBS_RV_Y(isp,trot,tvib)
                  ENDIF
                  helm(1)=helm(2)+helm(5)+helm(6)+helm(8)
              ELSEIF (flg_anha.EQ.2) THEN
                  IF (atom.EQ.0) THEN
                      helm(7)=GIBBS_INT_Y(isp,trot,tvib,tele)
                  ELSE
                      helm(7)=GIBBS_ELEC_Y(isp,tele)
                  ENDIF
                  helm(1)=helm(2)+helm(7)+helm(8)
              ELSE
              localinfo(1:38)='           Routine: species_helmholtz.'
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
! //  FUNCTIONS helmholtz_trans_x AND helmholtz_trans_y               //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          p           pressure (Pa)                               //
! //          rho         density (kg/m3)                             //
! //          tt          translation temperature (K)                 //
! //                                                                  //
! //  This function returns the translation helmholtz free energy per //
! //  unit mole (p,T input) or unit mass (rho,T input).               //
! //                                                                  //
! //  Benoit Bottin, 16/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION helmholtz_trans_x (isp,p,tt)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: p,tt,helm
!
      helm=RUNIV*tt*(-2.5d0*LOG(tt)-LOG(ST(isp))+LOG(p)-1.0d0)
      helmholtz_trans_x=helm
!
      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION helmholtz_trans_y (isp,rho,tt)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: rho,tt,helm,p
      REAL(kind=8) :: HELMHOLTZ_TRANS_X
!
      p=rho*RSP(isp)*tt
      helm=HELMHOLTZ_TRANS_X (isp,p,tt)
      helmholtz_trans_y=helm/MMOL(isp)

      RETURN
      END
! //////////////////////////////////////////////////////////////////////
