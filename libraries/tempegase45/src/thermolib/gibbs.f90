! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE species_gibbs                                        //
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
! //  Output: gibb(8)      array of gibbs free energies               //
! //                                                                  //
! //  This subroutine provides the gibbs energies of the species      //
! //  isp according to the flags specified. The routine also provides //
! //  the chemical potential of the species in the mixture.           //
! //  In order to get the temperature-dependent part only (mu-zero)   //
! //  simply use a fake pressure of 1 Pa in a mole calculation.       //
! //                                                                  //
! //  Benoit Bottin, 09/10/96. Modified to Gibbs and cleaned 25/7/97. //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE species_gibbs (isp,var,T,TNEQ,&
     &flg_anha,mode,flg_neq,flg_stop,flg_termo,gibb)
!
      USE global_thermo
      IMPLICIT NONE
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      INTEGER isp,flg_anha,mode,flg_neq,flg_stop,flg_termo,istop,i,atom
      REAL(kind=8) var,T,TNEQ(4),gibb(8),ttra,trot,tvib,tele
      REAL(kind=8) GIBBS_TRANS_X,GIBBS_TRANS_Y
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
101           FORMAT (11x,'Routine: species_gibbs. Value:',d15.8)
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
102   FORMAT (11x,'Routine: species_gibbs. Index:',i3,' Value:',d15.8)
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
103       FORMAT (11x,'Routine: species_gibbs. Value:',d15.8)
          CALL PRINT_ERROR (24+mode,pgname,localinfo,0)
          istop=1
      ENDIF
      DO i=1,8
          gibb(i)=0.0d0
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
                  localinfo(1:34)='           Routine: species_gibbs.'
                  CALL PRINT_ERROR (28,pgname,localinfo,0)
                  t=ttra
              ENDIF
              gibb(1)=GIBBS_REFTABLE_X (isp,var,t,flg_stop)
              gibb(8)=HFOR(isp)
          ELSE
              gibb(2)=GIBBS_TRANS_X(isp,var,ttra)
              gibb(8)=HFOR(isp)
              IF (flg_anha.EQ.0) THEN
                  IF (atom.EQ.0) THEN
                      gibb(3)=GIBBS_ROT_X(isp,trot)
                      gibb(4)=GIBBS_VIB_X(isp,tvib)
                  ENDIF
                  gibb(5)=GIBBS_ELEC_X(isp,tele)
                  gibb(1)=gibb(2)+gibb(3)+gibb(4)+gibb(5)
              ELSEIF (flg_anha.EQ.1) THEN
                  gibb(5)=GIBBS_ELEC_X(isp,tele)
                  IF (atom.EQ.0) THEN
                      gibb(6)=GIBBS_RV_X(isp,trot,tvib)
                  ENDIF
                  gibb(1)=gibb(2)+gibb(5)+gibb(6)
              ELSEIF (flg_anha.EQ.2) THEN
                  IF (atom.EQ.0) THEN
                      gibb(7)=GIBBS_INT_X(isp,trot,tvib,tele)
                  ELSE
                      gibb(7)=GIBBS_ELEC_X(isp,tele)
                  ENDIF
                  gibb(1)=gibb(2)+gibb(7)
              ELSE
                  localinfo(1:34)='           Routine: species_gibbs.'
                  CALL PRINT_ERROR (29,pgname,localinfo,0)
                  istop=1
              ENDIF
          ENDIF
      ELSE
          IF (flg_termo.EQ.1) THEN
              IF (flg_neq.NE.0) THEN
                  localinfo(1:34)='           Routine: species_gibbs.'
                  CALL PRINT_ERROR (28,pgname,localinfo,0)
                  t=ttra
              ENDIF
              gibb(1)=GIBBS_REFTABLE_Y (isp,var,t,flg_stop)
              gibb(8)=MOL_TO_MASS(HFOR(isp),isp,0.0d0,flg_stop)
          ELSE
              gibb(2)=GIBBS_TRANS_Y(isp,var,ttra)
              gibb(8)=MOL_TO_MASS(HFOR(isp),isp,0.0d0,flg_stop)
              IF (flg_anha.EQ.0) THEN
                  IF (atom.EQ.0) THEN
                      gibb(3)=GIBBS_ROT_Y(isp,trot)
                      gibb(4)=GIBBS_VIB_Y(isp,tvib)
                  ENDIF
                  gibb(5)=GIBBS_ELEC_Y(isp,tele)
                  gibb(1)=gibb(2)+gibb(3)+gibb(4)+gibb(5)
              ELSEIF (flg_anha.EQ.1) THEN
                  gibb(5)=GIBBS_ELEC_Y(isp,tele)
                  IF (atom.EQ.0) THEN
                      gibb(6)=GIBBS_RV_Y(isp,trot,tvib)
                  ENDIF
                  gibb(1)=gibb(2)+gibb(5)+gibb(6)
              ELSEIF (flg_anha.EQ.2) THEN
                  IF (atom.EQ.0) THEN
                      gibb(7)=GIBBS_INT_Y(isp,trot,tvib,tele)
                  ELSE
                      gibb(7)=GIBBS_ELEC_Y(isp,tele)
                  ENDIF
                  gibb(1)=gibb(2)+gibb(7)
              ELSE
                  localinfo(1:34)='           Routine: species_gibbs.'
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
! //  FUNCTIONS gibbs_trans_x AND gibbs_trans_y                       //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          p           pressure (Pa)                               //
! //          rho         density (kg/m3)                             //
! //          tt          translation temperature (K)                 //
! //                                                                  //
! //  This function returns the translation chemical gibbs per        //
! //  unit mole (p,T input) or unit mass (rho,T input).               //
! //                                                                  //
! //  Benoit Bottin, 09/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION gibbs_trans_x (isp,p,tt)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: p,tt,gibb
!
      gibb=RUNIV*tt*(-2.5d0*LOG(tt)-LOG(ST(isp))+LOG(p))
      gibbs_trans_x=gibb

      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION gibbs_trans_y (isp,rho,tt)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: rho,tt,gibb,p
      REAL(kind=8) :: GIBBS_TRANS_X
!
      p=rho*RSP(isp)*tt
      gibb=GIBBS_TRANS_X (isp,p,tt)
      gibbs_trans_y=gibb/MMOL(isp)

      RETURN
      END
! //////////////////////////////////////////////////////////////////////
      
! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS gibbs_rot_x AND gibbs_rot_y                           //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          trot        rotational temperature (K)                  //
! //                                                                  //
! //  This function returns the rotation chemical gibbs per           //
! //  unit mole (p,T input) or unit mass (rho,T input).               //
! //                                                                  //
! //  Benoit Bottin, 09/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION gibbs_rot_x (isp,trot)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: gibb,trot
!
      gibb=LOG( (trot/SIG(isp)/TR(1,isp))*(1.0d0+TR(1,isp)/3.0d0/trot) )
      gibb=-RUNIV*trot*gibb
      gibbs_rot_x=gibb
!
      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION gibbs_rot_y (isp,trot)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: trot,gibb,GIBBS_ROT_X
!
      gibb=GIBBS_ROT_X (isp,trot)
      gibbs_rot_y=gibb/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
      
! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS gibbs_vib_x AND gibbs_vib_y                           //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          tvib        vibrational temperature (K)                 //
! //                                                                  //
! //  This function returns the vibration chemical gibbs per          //
! //  unit mole (p,T input) or unit mass (rho,T input).               //
! //                                                                  //
! //  Benoit Bottin, 09/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION gibbs_vib_x (isp,tvib)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp,i
      REAL(kind=8) :: gibb,tvib,u
!
      gibb=0.0d0
      DO i=1,nvibmode(isp)
      u=TV(i,1,isp)/tvib
      IF (u.LT.307) THEN
          gibb=gibb+RUNIV*tvib*LOG(1.0d0-EXP(-u))
      ENDIF
      END DO
      gibbs_vib_x=gibb
!
      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION gibbs_vib_y (isp,tvib)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: tvib,gibb,GIBBS_VIB_X
!
      gibb=GIBBS_VIB_X (isp,tvib)
      gibbs_vib_y=gibb/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
      
! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS gibbs_elec_x AND gibbs_elec_y                          //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          tele        electronic temperature (K)                  //
! //                                                                  //
! //  This function returns the electronic chemical gibbs per         //
! //  unit mole (p,T input) or unit mass (rho,T input).               //
! //                                                                  //
! //  Benoit Bottin, 09/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION gibbs_elec_x (isp,tele)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp,i
      REAL(kind=8) :: gibb,tele,q,e
!
      gibb=0.0d0
      q=0.0d0
      DO i=1,NELVL(isp)
          e=TE(i,isp)/tele
          IF (e.LT.307) THEN
              q=q+GE(i,isp)*EXP(-e)
          ENDIF
      END DO
      gibb=-RUNIV*tele*LOG(q)
      gibbs_elec_x=gibb
!
      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION gibbs_elec_y (isp,tele)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: tele,gibb,GIBBS_ELEC_X
!
      gibb=GIBBS_ELEC_X (isp,tele)
      gibbs_elec_y=gibb/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS gibbs_rv_x AND gibbs_rv_y                             //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          trot        rotational temperature (K)                  //
! //          tvib        vibrational temperature (K)                 //
! //                                                                  //
! //  This function returns the rotation-vibration coupled chemical   //
! //  gibbs per unit mole (p,T input) or unit mass (rho,T input).     //
! //                                                                  //
! //  Benoit Bottin, 09/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION gibbs_rv_x (isp,trot,tvib)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: trot,tvib,d,g,x,u,n,mu,den
      REAL(kind=8) :: gibbs_ROT_X,gibbs_VIB_X
!
      g=G_A(1,isp)
      d=D_A(1,isp)
      x=X_A(1,isp)
      u=TV(1,1,isp)/tvib
      n=TR(1,isp)/trot
      den=EXP(u)-1.0d0
!
      mu=gibbs_ROT_X(isp,trot)+gibbs_VIB_X(isp,tvib)
      mu=mu-RUNIV*(g*trot/n+d*tvib/den+2.0d0*x*TV(1,1,isp)/(den**2))
!
      gibbs_rv_x=mu

      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION gibbs_rv_y (isp,trot,tvib)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: trot,tvib,gibb,GIBBS_RV_X
!
      gibb=GIBBS_RV_X (isp,trot,tvib)
      gibbs_rv_y=gibb/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS gibbs_rv_x AND gibbs_rv_y                             //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          trot        rotational temperature (K)                  //
! //          tvib        vibrational temperature (K)                 //
! //          tele        electronic temperature (K)                  //
! //                                                                  //
! //  This function returns the internally-coupled chemical           //
! //  gibbs per unit mole (p,T input) or unit mass (rho,T input).     //
! //                                                                  //
! //  Benoit Bottin, 09/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION gibbs_int_x (isp,trot,tvib,tele)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp,i
      REAL(kind=8) :: trot,tvib,tele,a,b,c,e,u,n,den,enern,enerd,expu
!
      enern=0.0d0
      enerd=0.0d0
      DO i=1,NELVL(isp)
          IF (TV(1,i,isp).GT.1.0d-6) THEN
             e=TE(i,isp)/tele
             u=TV(1,i,isp)/tvib
             n=TR(i,isp)/trot
             expu=EXP(u)
             den=expu-1.0d0
             a=GE(i,isp)*EXP(-e)
             b=SIG(isp)*n*(1.0d0-1.0d0/expu)
             c=1.0d0+n/3.0d0+n**2/15.0d0
             c=c+G_A(i,isp)/n+D_A(i,isp)/den+2.0d0*X_A(i,isp)*u/(den**2)
             enerd=enerd+a*c/b
          ENDIF
      END DO
!
      gibbs_int_x=-RUNIV*trot*LOG(enerd)

      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION gibbs_int_y (isp,trot,tvib,tele)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: trot,tvib,tele,gibb,GIBBS_INT_X
!
      gibb=GIBBS_INT_X (isp,trot,tvib,tele)
      gibbs_int_y=gibb/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS gibbs_reftable_x and gibbs_reftable_y                 //
! //                                                                  //
! //  Input:  T           temperature of the species (K)              //    
! //          p           pressure of the species (Pa)                //
! //          rho         density of the species (kg/m3)              //
! //          isp         number of the concerned species             //
! //          flg_stop    stops if errors in the subroutine           //
! //                                                                  //
! //  This function returns the chemical gibbs of a species from      //
! //  reference tables such as Gurvich or Janaf (mole or mass input)  //
! //                                                                  //
! //  Benoit Bottin, 11/12/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION gibbs_reftable_x (isp,p,T,flg_stop)
!
      USE global_thermo
      IMPLICIT NONE
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      INTEGER isp,ilow,ihigh,istop,flg_stop
      REAL(kind=8) p,t,x,x1,x2,h,s,y1,y2
      REAL(kind=8) ENTHALPY_REFTABLE_X,ENTROPY_REFTABLE_X
!
      s=ENTROPY_REFTABLE_X (isp,p,T,flg_stop)
      h=ENTHALPY_REFTABLE_X (isp,T,flg_stop)
      gibbs_reftable_x=h-T*s
!
      RETURN
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
102               FORMAT (11x,'Routine: species_gibbs. Value:',d15.8)
                  CALL PRINT_ERROR (31,pgname,localinfo,0)
                  istop=1
              ELSE
                  GOTO 10
              ENDIF
          ENDIF
      ihigh=ilow+1
      IF (ilow.EQ.0) THEN
          WRITE (localinfo,101) t
101       FORMAT (11x,'Routine: species_gibbs. Value:',d15.8)
          CALL PRINT_ERROR (30,pgname,localinfo,0)
          istop=1
      ENDIF
      IF ((flg_stop.EQ.1).AND.(istop.EQ.1)) THEN
          STOP 'Program ended due to PEGASE thermodynamic library error'
      ENDIF
!
!     Perform the interpolation on enthalpy
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      x=t
      x1=REFTABLE(isp,1,ilow)
      x2=REFTABLE(isp,1,ihigh)
      y1=REFTABLE(isp,4,ilow)
      y2=REFTABLE(isp,4,ihigh)
      h=(y2-y1)/(x2-x1)*(x-x1)+y1
!
!     Perform the interpolation on entropy
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      y1=REFTABLE(isp,3,ilow)
      y2=REFTABLE(isp,3,ihigh)
      s=(y2-y1)/(x2-x1)*(x-x1)+y1
!
      gibbs_reftable_x=h-T*(s-RUNIV*DLOG(p/REFPRESS(isp)))
      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION gibbs_reftable_y (isp,rho,T,flg_stop)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER isp,flg_stop
      REAL(kind=8) p,rho,t,gibb,GIBBS_REFTABLE_X
!
      gibb=0.0d0
      p=rho*RSP(isp)*T
      gibb=GIBBS_REFTABLE_X (isp,p,T,flg_stop)
      gibbs_reftable_y=gibb/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
