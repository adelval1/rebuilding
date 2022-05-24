! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE species_entropy                                      //
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
! //  Output: entro(8)      array of entropys                         //
! //              1           complete                                //
! //              2           translation                             //
! //              3           rotation (anha=0)                       //
! //              4           vibration (anha=0)                      //
! //              5           electronic (anha=0)                     //
! //              6           coupled rotation-vibration (anha=1)     //
! //              7           coupled internal (anha=2)               //
! //              8           not used                                //
! //                                                                  //
! //  This subroutine provides the entropy of the species             //
! //  isp according to the flags specified.                           //
! //                                                                  //
! //  Benoit Bottin, 11/10/96. Modified and cleaned 25/7/97.          //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE species_entropy (isp,var,T,TNEQ,&
     &flg_anha,mode,flg_neq,flg_stop,flg_termo,entro)
!
      USE global_thermo
      IMPLICIT NONE
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      INTEGER isp,flg_anha,mode,flg_neq,flg_stop,flg_termo,istop,i,atom
      REAL(kind=8) var,T,TNEQ(4),entro(8),ttra,trot,tvib,tele
      REAL(kind=8) ENTROPY_TRANS_X,ENTROPY_TRANS_Y
      REAL(kind=8) ENTROPY_ROT_X,ENTROPY_ROT_Y
      REAL(kind=8) ENTROPY_VIB_X,ENTROPY_VIB_Y
      REAL(kind=8) ENTROPY_ELEC_X,ENTROPY_ELEC_Y
      REAL(kind=8) ENTROPY_RV_X,ENTROPY_RV_Y
      REAL(kind=8) ENTROPY_INT_X,ENTROPY_INT_Y
      REAL(kind=8) ENTROPY_REFTABLE_X,ENTROPY_REFTABLE_Y
!
!     Checks of correct variable values and initialize
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      pgname='pegaslib'
      istop=0
      IF (flg_neq.EQ.0) THEN
          IF (T.LE.0) THEN
              WRITE (localinfo,101) T
101           FORMAT (11x,'Routine: species_entropy. Value:',d15.8)
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
102   FORMAT (11x,'Routine: species_entropy. Index:',i3,' Value:',d15.8)
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
103       FORMAT (11x,'Routine: species_entropy. Value:',d15.8)
          CALL PRINT_ERROR (24+mode,pgname,localinfo,0)
          istop=1
      ENDIF
      DO i=1,8
          entro(i)=0.0d0
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
                  localinfo(1:36)='           Routine: species_entropy.'
                  CALL PRINT_ERROR (28,pgname,localinfo,0)
                  t=ttra
              ENDIF
              entro(1)=ENTROPY_REFTABLE_X (isp,var,t,flg_stop)
          ELSE
              entro(2)=ENTROPY_TRANS_X(isp,var,ttra)
              IF (flg_anha.EQ.0) THEN
                  IF (atom.EQ.0) THEN
                      entro(3)=ENTROPY_ROT_X(isp,trot)
                      entro(4)=ENTROPY_VIB_X(isp,tvib)
                  ENDIF
                  entro(5)=ENTROPY_ELEC_X(isp,tele)
                  entro(1)=entro(2)+entro(3)+entro(4)+entro(5)
              ELSEIF (flg_anha.EQ.1) THEN
                  entro(5)=ENTROPY_ELEC_X(isp,tele)
                  IF (atom.EQ.0) THEN
                      entro(6)=ENTROPY_RV_X(isp,trot,tvib)
                  ENDIF
                  entro(1)=entro(2)+entro(5)+entro(6)
              ELSEIF (flg_anha.EQ.2) THEN
                  IF (atom.EQ.0) THEN
                      entro(7)=ENTROPY_INT_X(isp,trot,tvib,tele)
                  ELSE
                      entro(7)=ENTROPY_ELEC_X(isp,tele)
                  ENDIF
                  entro(1)=entro(2)+entro(7)
              ELSE
                  localinfo(1:36)='           Routine: species_entropy.'
                  CALL PRINT_ERROR (29,pgname,localinfo,0)
                  istop=1
              ENDIF
          ENDIF
      ELSE
          IF (flg_termo.EQ.1) THEN
              IF (flg_neq.NE.0) THEN
                  localinfo(1:36)='           Routine: species_entropy.'
                  CALL PRINT_ERROR (28,pgname,localinfo,0)
                  t=ttra
              ENDIF
              entro(1)=ENTROPY_REFTABLE_Y (isp,var,t,flg_stop)
          ELSE
              entro(2)=ENTROPY_TRANS_Y(isp,var,ttra)
              IF (flg_anha.EQ.0) THEN
                  IF (atom.EQ.0) THEN
                      entro(3)=ENTROPY_ROT_Y(isp,trot)
                      entro(4)=ENTROPY_VIB_Y(isp,tvib)
                  ENDIF
                  entro(5)=ENTROPY_ELEC_Y(isp,tele)
                  entro(1)=entro(2)+entro(3)+entro(4)+entro(5)
              ELSEIF (flg_anha.EQ.1) THEN
                  entro(5)=ENTROPY_ELEC_Y(isp,tele)
                  IF (atom.EQ.0) THEN
                      entro(6)=ENTROPY_RV_Y(isp,trot,tvib)
                  ENDIF
                  entro(1)=entro(2)+entro(5)+entro(6)
              ELSEIF (flg_anha.EQ.2) THEN
                  IF (atom.EQ.0) THEN
                      entro(7)=ENTROPY_INT_Y(isp,trot,tvib,tele)
                  ELSE
                      entro(7)=ENTROPY_ELEC_Y(isp,tele)
                  ENDIF
                  entro(1)=entro(2)+entro(7)
              ELSE
                  localinfo(1:36)='           Routine: species_entropy.'
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
! //  FUNCTIONS entropy_trans_x AND entropy_trans_y                   //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          p           pressure (Pa)                               //
! //          rho         density (kg/m3)                             //
! //          tt          translation temperature (K)                 //
! //                                                                  //
! //  This function returns the translation entropy per               //
! //  unit mole (p,T input) or unit mass (rho,T input).               //
! //                                                                  //
! //  Benoit Bottin, 11/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION entropy_trans_x (isp,p,tt)
!      
      USE global_thermo
      IMPLICIT NONE
      REAL(kind=8) p,tt,entro
      INTEGER isp
!      
      entro=0.0d0
      entro=2.5d0*RUNIV+RUNIV*DLOG(ST(isp))+2.5d0*RUNIV*DLOG(tt)
      entro=entro-RUNIV*DLOG(p)
      entropy_trans_x=entro
!      
      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION entropy_trans_y (isp,rho,tt)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER isp
      REAL(kind=8) rho,tt,entro,p
      REAL(kind=8) ENTROPY_TRANS_X
!      
      entro=0.0d0
      p=rho*RSP(isp)*tt
      entro=ENTROPY_TRANS_X (isp,p,tt)
      entropy_trans_y=entro/MMOL(isp)
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
      
! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS entropy_rot_x AND entropy_rot_y                       //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          trot        rotational temperature (K)                  //
! //                                                                  //
! //  This function returns the rotation entropy per                  //
! //  unit mole (p,T input) or unit mass (rho,T input).               //
! //                                                                  //
! //  Benoit Bottin, 11/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION entropy_rot_x (isp,trot)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER isp
      REAL(kind=8) entro,trot
!      
      entro=0.0d0
      entro=1.0d0+DLOG(trot/SIG(isp)/TR(1,isp))
      entro=entro+DLOG(1.0d0+TR(1,isp)/3.0d0/trot)
      entro=entro-TR(1,isp)/(TR(1,isp)+3.0d0*trot)
      entro=entro*RUNIV
      entropy_rot_x=entro
!      
      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION entropy_rot_y (isp,trot)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER isp
      REAL(kind=8) trot,entro,ENTROPY_ROT_X
!
      entro=0.0d0
      entro=ENTROPY_ROT_X (isp,trot)
      entropy_rot_y=entro/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
      
! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS entropy_vib_x AND entropy_vib_y                       //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          tvib        vibrational temperature (K)                 //
! //                                                                  //
! //  This function returns the vibration entropy per                 //
! //  unit mole (p,T input) or unit mass (rho,T input).               //
! //                                                                  //
! //  Benoit Bottin, 11/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION entropy_vib_x (isp,tvib)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER isp,i
      REAL(kind=8) entro,tvib,u
!      
      entro=0.0d0
      DO i=1,nvibmode(isp)
      u=TV(i,1,isp)/tvib
      IF (u.LT.307) THEN
          entro=entro+RUNIV*(u/(DEXP(u)-1.0d0)-DLOG(1-DEXP(-u)))
      ENDIF
      END DO 
      entropy_vib_x=entro
!      
      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION entropy_vib_y (isp,tvib)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER isp
      REAL(kind=8) tvib,entro,ENTROPY_VIB_X
!
      entro=0.0d0
      entro=ENTROPY_VIB_X (isp,tvib)
      entropy_vib_y=entro/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
      
! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS entropy_elec_x AND entropy_elec_y                     //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          tele        electronic temperature (K)                  //
! //                                                                  //
! //  This function returns the electronic entropy per                //
! //  unit mole (p,T input) or unit mass (rho,T input).               //
! //                                                                  //
! //  Benoit Bottin, 11/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION entropy_elec_x (isp,tele)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER isp,i
      REAL(kind=8) entro,tele,q,num
!      
      entro=0.0d0
      q=0.0d0
      num=0.0d0
      DO i=1,NELVL(isp)
          IF ((TE(i,isp)/tele).LT.307) THEN
              q=q+GE(i,isp)*DEXP(-TE(i,isp)/tele)
              num=num+TE(i,isp)*GE(i,isp)*DEXP(-TE(i,isp)/tele)
          ENDIF
      END DO
      entro=RUNIV*DLOG(q)+(RUNIV/tele)*(num/q)
      entropy_elec_x=entro
!      
      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION entropy_elec_y (isp,tele)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER isp
      REAL(kind=8) tele,entro,ENTROPY_ELEC_X
!
      entro=0.0d0
      entro=ENTROPY_ELEC_X (isp,tele)
      entropy_elec_y=entro/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS entropy_rv_x AND entropy_rv_y                         //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          trot        rotational temperature (K)                  //
! //          tvib        vibrational temperature (K)                 //
! //                                                                  //
! //  This function returns the rotation-vibration coupled            //
! //  entropy per unit mole (p,T input) or unit mass (rho,T input).   //
! //                                                                  //
! //  Benoit Bottin, 11/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION entropy_rv_x (isp,trot,tvib)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER isp
      REAL(kind=8) trot,tvib,d,g,x,u,n,s,den
      REAL(kind=8) ENTROPY_ROT_X,ENTROPY_VIB_X
!
      g=G_A(1,isp)
      d=D_A(1,isp)
      x=X_A(1,isp)
      u=TV(1,1,isp)/tvib
      n=TR(1,isp)/trot
      den=DEXP(u)-1.0d0
!
      s=ENTROPY_ROT_X(isp,trot)+ENTROPY_VIB_X(isp,tvib)
      s=s+2.0d0*RUNIV*g/n
      s=s+d*RUNIV/den+d*RUNIV*u*DEXP(u)/(den**2)
      s=s+(4.0d0*x*RUNIV*(u**2)*DEXP(u))/(den**3)
      entropy_rv_x=s
      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION entropy_rv_y (isp,trot,tvib)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER isp
      REAL(kind=8) trot,tvib,entro,ENTROPY_RV_X
!
      entro=0.0d0
      entro=ENTROPY_RV_X (isp,trot,tvib)
      entropy_rv_y=entro/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS entropy_int_x AND entropy_int_y                       //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          trot        rotational temperature (K)                  //
! //          tvib        vibrational temperature (K)                 //
! //          tele        electronic temperature (K)                  //
! //                                                                  //
! //  This function returns the internally-coupled                    //
! //  entropy per unit mole (p,T input) or unit mass (rho,T input).   //
! //                                                                  //
! //  Benoit Bottin, 11/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION entropy_int_x (isp,trot,tvib,tele)
      USE global_thermo
      IMPLICIT NONE
      INTEGER isp,i
      REAL(kind=8) trot,tvib,tele,a,b,c,da,db,dc,e,u,n,den,enern,enerd
!
      enern=0.0d0
      enerd=0.0d0
      DO i=1,NELVL(isp)
          IF (TV(1,i,isp).GT.1.0d-6) THEN
             e=TE(i,isp)/tele
             u=TV(1,i,isp)/tvib
             n=TR(i,isp)/trot
             den=DEXP(u)-1.0d0
             a=GE(i,isp)*DEXP(-e)
             da=GE(i,isp)*e/tele*DEXP(-e)
             b=SIG(isp)*n*(1.0d0-DEXP(-u))
             db=SIG(isp)*n/trot*(-1.0d0+DEXP(-u)*(1.0d0-u))
             c=1.0d0+n/3.0d0+n**2/15.0d0
             c=c+G_A(i,isp)/n+D_A(i,isp)/den+2.0d0*X_A(i,isp)*u/(den**2)
             dc=-n/3.0d0/trot-2.0d0*n**2/15.0d0/trot+G_A(i,isp)/n/trot
             dc=dc+(D_A(i,isp)*u*DEXP(u)-2.0d0*X_A(i,isp)*u)/trot/den**2
             dc=dc+4.0d0*X_A(i,isp)*u**2*DEXP(u)/trot/den**3
             enern=enern+((da*c+a*dc)*b-a*c*db)/b/b
             enerd=enerd+a*c/b
          ENDIF
      END DO
!
      entropy_int_x=RUNIV*DLOG(enerd)+RUNIV*trot*enern/enerd
      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION entropy_int_y (isp,trot,tvib,tele)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER isp
      REAL(kind=8) trot,tvib,tele,entro,ENTROPY_INT_X
!
      entro=0.0d0
      entro=ENTROPY_INT_X (isp,trot,tvib,tele)
      entropy_int_y=entro/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS entropy_reftable_x and entropy_reftable_y             //
! //                                                                  //
! //  Input:  T           temperature of the species (K)              //    
! //          p           pressure of the species (Pa)                //
! //          rho         density of the species (kg/m3)              //
! //          isp         number of the concerned species             //
! //          flg_stop    stops if errors in the subroutine           //
! //                                                                  //
! //  This function returns the entropy  of the species from          //
! //  reference tables such as Gurvich or Janaf (mole or mass input)  //
! //                                                                  //
! //  Benoit Bottin, 11/12/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION entropy_reftable_x (isp,p,T,flg_stop)
!
      USE global_thermo
      IMPLICIT NONE
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      INTEGER isp,ilow,ihigh,istop,flg_stop
      REAL(kind=8) p,t,x,x1,x2,y,y1,y2
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
102               FORMAT (11x,'Routine: species_entropy. Value:',d15.8)
                  CALL PRINT_ERROR (31,pgname,localinfo,0)
                  istop=1
              ELSE
                  GOTO 10
              ENDIF
          ENDIF
      ihigh=ilow+1
      IF (ilow.EQ.0) THEN
          WRITE (localinfo,101) t
101       FORMAT (11x,'Routine: species_entropy. Value:',d15.8)
          CALL PRINT_ERROR (30,pgname,localinfo,0)
          istop=1
      ENDIF
      IF ((flg_stop.EQ.1).AND.(istop.EQ.1)) THEN
          STOP 'Program ended due to PEGASE thermodynamic library error'
      ENDIF
!
!     Perform the interpolation
!     ^^^^^^^^^^^^^^^^^^^^^^^^^
      x=t
      x1=REFTABLE(isp,1,ilow)
      x2=REFTABLE(isp,1,ihigh)
      y1=REFTABLE(isp,3,ilow)
      y2=REFTABLE(isp,3,ihigh)
      y=(y2-y1)/(x2-x1)*(x-x1)+y1
!
      entropy_reftable_x=y-RUNIV*DLOG(p/REFPRESS(isp))
      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION entropy_reftable_y (isp,rho,T,flg_stop)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER isp,flg_stop
      REAL(kind=8) p,rho,t,entro,ENTROPY_REFTABLE_X
!
      entro=0.0d0
      p=rho*RSP(isp)*T
      entro=ENTROPY_REFTABLE_X (isp,p,T,flg_stop)
      entropy_reftable_y=entro/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
