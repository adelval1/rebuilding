! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE species_cv                                           //
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
! //  Output: cv(8)       array of internal energies                  //
! //              1           complete (except formation energy)      //
! //              2           translation                             //
! //              3           rotation (anha=0)                       //
! //              4           vibration (anha=0)                      //
! //              5           electronic (anha=0)                     //
! //              6           coupled rotation-vibration (anha=1)     //
! //              7           coupled internal (anha=2)               //
! //              8           formation energy/enthalpy               //
! //                                                                  //
! //  This subroutine provides the internal energies of the species   //
! //  isp according to the flags specified.                           //
! //                                                                  //
! //  Benoit Bottin, 31/11/96. Modified and cleaned 25/7/97.          //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE species_cv (isp,T,TNEQ,&
     &flg_anha,mode,flg_neq,flg_stop,flg_termo,cv)
!
      USE global_thermo
      IMPLICIT NONE
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      INTEGER isp,flg_anha,mode,flg_neq,flg_stop,flg_termo,istop,i,atom
      REAL(kind=8) T,TNEQ(4),cv(8),ttra,trot,tvib,tele
      REAL(kind=8) CV_TRANS_X,CV_TRANS_Y
      REAL(kind=8) CV_ROT_X,CV_ROT_Y
      REAL(kind=8) CV_VIB_X,CV_VIB_Y
      REAL(kind=8) CV_ELEC_X,CV_ELEC_Y
      REAL(kind=8) CV_RV_X,CV_RV_Y
      REAL(kind=8) CV_INT_X,CV_INT_Y
      REAL(kind=8) CV_REFTABLE_X,CV_REFTABLE_Y
!
!     Checks of correct variable values and initialize
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      pgname='pegaslib'
      istop=0
      IF (flg_neq.EQ.0) THEN
          IF (T.LE.0) THEN
              WRITE (localinfo,101) T
101           FORMAT (11x,'Routine: species_cv. Value:',d15.8)
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
102   FORMAT (11x,'Routine: species_cv. Index:',i3,' Value:',d15.8)
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
          cv(i)=0.0d0
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
                  localinfo(1:31)='           Routine: species_cv.'
                  CALL PRINT_ERROR (28,pgname,localinfo,0)
                  t=ttra
              ENDIF
              cv(1)= CV_REFTABLE_X (isp,t,flg_stop)
          ELSE
              cv(2)=CV_TRANS_X(isp,ttra)
              IF (flg_anha.EQ.0) THEN
                  IF (atom.EQ.0) THEN
                      cv(3)=CV_ROT_X(isp,trot)
                      cv(4)=CV_VIB_X(isp,tvib)
                  ENDIF
                  cv(5)=CV_ELEC_X(isp,tele)
                  cv(1)=cv(2)+cv(3)+cv(4)+cv(5)
              ELSEIF (flg_anha.EQ.1) THEN
                  cv(5)=CV_ELEC_X(isp,tele)
                  IF (atom.EQ.0) THEN
                      cv(6)=CV_RV_X(isp,trot,tvib)
                  ENDIF
                  cv(1)=cv(2)+cv(5)+cv(6)
              ELSEIF (flg_anha.EQ.2) THEN
                  IF (atom.EQ.0) THEN
                      cv(7)=CV_INT_X(isp,trot,tvib,tele)
                  ELSE
                      cv(7)=CV_ELEC_X(isp,tele)
                  ENDIF
                  cv(1)=cv(2)+cv(7)
              ELSE
                  localinfo(1:31)='           Routine: species_cv.'
                  CALL PRINT_ERROR (29,pgname,localinfo,0)
                  istop=1
              ENDIF
          ENDIF
      ELSE
          IF (flg_termo.EQ.1) THEN
              IF (flg_neq.NE.0) THEN
                  localinfo(1:31)='           Routine: species_cv.'
                  CALL PRINT_ERROR (28,pgname,localinfo,0)
                  t=ttra
              ENDIF
              cv(1)= CV_REFTABLE_Y (isp,t,flg_stop)
          ELSE
              cv(2)=CV_TRANS_Y(isp,ttra)
              IF (flg_anha.EQ.0) THEN
                  IF (atom.EQ.0) THEN
                      cv(3)=CV_ROT_Y(isp,trot)
                      cv(4)=CV_VIB_Y(isp,tvib)
                  ENDIF
                  cv(5)=CV_ELEC_Y(isp,tele)
                  cv(1)=cv(2)+cv(3)+cv(4)+cv(5)
              ELSEIF (flg_anha.EQ.1) THEN
                  cv(5)=CV_ELEC_Y(isp,tele)
                  IF (atom.EQ.0) THEN
                      cv(6)=CV_RV_Y(isp,trot,tvib)
                  ENDIF
                  cv(1)=cv(2)+cv(5)+cv(6)
              ELSEIF (flg_anha.EQ.2) THEN
                  IF (atom.EQ.0) THEN
                      cv(7)=CV_INT_Y(isp,trot,tvib,tele)
                  ELSE
                      cv(7)=CV_ELEC_Y(isp,tele)
                  ENDIF
                  cv(1)=cv(2)+cv(7)
              ELSE
                  localinfo(1:31)='           Routine: species_cv.'
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
! //  FUNCTIONS cv_trans_x AND cv_trans_y                             //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          tt          translation temperature (K)                 //
! //                                                                  //
! //  This function returns the translation internal energy per       //
! //  unit mole (p,T input) or unit mass (rho,T input).               //
! //                                                                  //
! //  Benoit Bottin, 31/11/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION cv_trans_x (isp,tt)
      
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp,i
      REAL(kind=8) :: tt,cv

      cv=1.5d0*RUNIV
      cv_trans_x=cv
!
!     isp and tt are passed to keep the similarity between all calls
!     following assignement removes any compiling error
      i=isp      
      cv=tt

      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION cv_trans_y (isp,tt)
      
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: tt,cv
      REAL(kind=8) :: CV_TRANS_X
      
      cv=CV_TRANS_X (isp,tt)
      cv_trans_y=cv/MMOL(isp)

      RETURN
      END
! //////////////////////////////////////////////////////////////////////
      
! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS cv_rot_x AND cv_rot_y                                 //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          trot        rotational temperature (K)                  //
! //                                                                  //
! //  This function returns the rotation internal cv per              //
! //  unit mole (p,T input) or unit mass (rho,T input).               //
! //                                                                  //
! //  Benoit Bottin, 31/11/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION cv_rot_x (isp,trot)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: cv,trot
!      
      cv=RUNIV*(1.0d0-(TR(1,isp)/(TR(1,isp)+3.0d0*trot))**2)
      cv_rot_x=cv

      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION cv_rot_y (isp,trot)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: trot,cv,CV_ROT_X
!
      cv=CV_ROT_X (isp,trot)
      cv_rot_y=cv/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
      
! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS cv_vib_x AND cv_vib_y                                 //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          tvib        vibrational temperature (K)                 //
! //                                                                  //
! //  This function returns the vibration internal cv per             //
! //  unit mole (p,T input) or unit mass (rho,T input).               //
! //                                                                  //
! //  Benoit Bottin, 31/11/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION cv_vib_x (isp,tvib)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: i,isp
      REAL(kind=8) :: cv,tvib,var,u,expu
!     
      cv=0.0d0
      do i=1,nvibmode(isp)
      u=0.5*TV(i,1,isp)/tvib
      IF (u.LT.307) THEN
          expu=EXP(u)
          var=0.5*(expu-1.0d0/expu)
          cv=cv+RUNIV*(u/var)**2
      ENDIF
      END DO 
      cv_vib_x=cv
      
!      
      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION cv_vib_y (isp,tvib)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: tvib,cv,CV_VIB_X
!
      cv=CV_VIB_X (isp,tvib)
      cv_vib_y=cv/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
      
! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS cv_elec_x AND cv_elec_y                               //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          p           pressure (Pa)                               //
! //          rho         density (kg/m3)                             //
! //          tele        electronic temperature (K)                  //
! //                                                                  //
! //  This function returns the electronic internal cv per            //
! //  unit mole (p,T input) or unit mass (rho,T input).               //
! //                                                                  //
! //  Benoit Bottin, 31/11/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION cv_elec_x (isp,tele)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp,i
      REAL(kind=8) :: cv,tele,q,num1,num2,e,ge_exp_e
!      
      cv=0.0d0
      q=0.0d0
      num1=0.0d0
      num2=0.0d0
      DO i=1,NELVL(isp)
          e=TE(i,isp)/tele
          IF (e.LT.307) THEN
              ge_exp_e=GE(i,isp)*EXP(-e)
              q=q+ge_exp_e
              num1=num1+TE(i,isp)*ge_exp_e
              num2=num2+(TE(i,isp)**2)*ge_exp_e
          ENDIF
      END DO       
      cv=RUNIV/(tele**2)*(num2/q-(num1/q)**2)
      cv_elec_x=cv
!      
      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION cv_elec_y (isp,tele)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: tele,cv,CV_ELEC_X
!
      cv=CV_ELEC_X (isp,tele)
      cv_elec_y=cv/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS cv_rv_x AND cv_rv_y                                   //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          p           pressure (Pa)                               //
! //          rho         density (kg/m3)                             //
! //          trot        rotational temperature (K)                  //
! //          tvib        vibrational temperature (K)                 //
! //                                                                  //
! //  This function returns the rotation-vibration coupled internal   //
! //  cv per unit mole (p,T input) or unit mass (rho,T input).        //
! //                                                                  //
! //  Benoit Bottin, 31/11/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION cv_rv_x (isp,trot,tvib)
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: trot,tvib,d,g,x,u,n,cv,den,cva,expu,alpha
      REAL(kind=8) :: CV_ROT_X,CV_VIB_X

      g=G_A(1,isp)
      d=D_A(1,isp)
      x=X_A(1,isp)
! **  This part is modified       
      u=TV(1,1,isp)/tvib
      n=TR(1,isp)/trot
      expu=EXP(u)
      den=expu-1.0d0

      cv=CV_ROT_X(isp,trot)+CV_VIB_X(isp,tvib)
      cva=2.0d0*g/n
      alpha=u*expu/den**2
      cva=cva+u*alpha*(-d   &
&             +(d*expu-4.0d0*x*u-8.0d0*x)/den   &
&             +12.0d0*x*alpha)

      cv_rv_x=cv+cva*RUNIV

      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION cv_rv_y (isp,trot,tvib)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: trot,tvib,cv,CV_RV_X
!
      cv=0.0d0
      cv=CV_RV_X (isp,trot,tvib)
      cv_rv_y=cv/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS cv_int_x AND cv_int_y                                 //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          p           pressure (Pa)                               //
! //          rho         density (kg/m3)                             //
! //          trot        rotational temperature (K)                  //
! //          tvib        vibrational temperature (K)                 //
! //          tele        electronic temperature (K)                  //
! //                                                                  //
! //  This function returns the internally-coupled internal           //
! //  cv per unit mole (p,T input) or unit mass (rho,T input).        //
! //                                                                  //
! //  Benoit Bottin, 31/11/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION cv_int_x (isp,trot,tvib,tele)
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp,i
      REAL(kind=8) :: trot,tvib,tele,a,b,c,da,db,dc,e,u,n,den,enern,enerd
      REAL(kind=8) :: cv1,cv2,nrot,nvib,nele,expu,exp_u,exp_e
!
!     Computation at temperatures T
      enern=0.0d0
      enerd=0.0d0
      DO i=1,NELVL(isp)
! ** This part is modified
          IF (TV(1,i,isp).GT.1.0d-6) THEN
             e=TE(i,isp)/tele
             u=TV(1,i,isp)/tvib
             n=TR(i,isp)/trot
             expu=EXP(u)
             exp_u=1.0d0/expu
             exp_e=EXP(-e)
             den=expu-1.0d0
             a=GE(i,isp)*exp_e
             da=GE(i,isp)*e/tele*exp_e
             b=SIG(isp)*n*(1.0d0-exp_u)
             db=SIG(isp)*n/trot*(-1.0d0+exp_u*(1.0d0-u))
             c=1.0d0+n/3.0d0+n**2/15.0d0
             c=c+G_A(i,isp)/n+D_A(i,isp)/den+2.0d0*X_A(i,isp)*u/(den**2)
             dc=(-n/3.0d0-2.0d0*n**2/15.0d0+G_A(i,isp)/n)/trot
             dc=dc+(D_A(i,isp)*expu-2.0d0*X_A(i,isp))*u/trot/den**2
             dc=dc+4.0d0*X_A(i,isp)*u**2*expu/trot/den**3
             enern=enern+((da*c+a*dc)*b-a*c*db)/b/b
             enerd=enerd+a*c/b
          ENDIF
      END DO
      cv1=RUNIV*trot**2*enern/enerd
!
!     Computation at temperatures T+dT
      nrot=trot*(1.0d0+1.0d-10)
      nvib=tvib*(1.0d0+1.0d-10)
      nele=tele*(1.0d0+1.0d-10)
      enern=0.0d0
      enerd=0.0d0
      DO i=1,NELVL(isp)
          IF (TV(1,i,isp).GT.1.0d-6) THEN
             e=TE(i,isp)/nele
             u=TV(1,i,isp)/nvib
             n=TR(i,isp)/nrot
             expu=EXP(u)
             exp_u=1.0d0/expu
             exp_e=EXP(-e)
             den=expu-1.0d0
             a=GE(i,isp)*exp_e
             da=GE(i,isp)*e/nele*exp_e
             b=SIG(isp)*n*(1.0d0-exp_u)
             db=SIG(isp)*n/nrot*(-1.0d0+exp_u*(1.0d0-u))
             c=1.0d0+n/3.0d0+n**2/15.0d0
             c=c+G_A(i,isp)/n+D_A(i,isp)/den+2.0d0*X_A(i,isp)*u/(den**2)
             dc=(-n/3.0d0-2.0d0*n**2/15.0d0+G_A(i,isp)/n)/nrot
             dc=dc+(D_A(i,isp)*expu-2.0d0*X_A(i,isp))*u/nrot/den**2
             dc=dc+4.0d0*X_A(i,isp)*u**2*expu/nrot/den**3
             enern=enern+((da*c+a*dc)*b-a*c*db)/b/b
             enerd=enerd+a*c/b
          ENDIF
      END DO
      cv2=RUNIV*nrot**2*enern/enerd
!
!     Computation of cv by finite differences
      cv_int_x=(cv2-cv1)/nrot/1.0d-10
!      
      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION cv_int_y (isp,trot,tvib,tele)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: trot,tvib,tele,cv,CV_INT_X
!
      cv=CV_INT_X (isp,trot,tvib,tele)
      cv_int_y=cv/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS cv_reftable_x and cv_reftable_y                       //
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
      REAL(kind=8) FUNCTION cv_reftable_x (isp,T,flg_stop)
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
102               FORMAT (11x,'Routine: species_cv. Value:',d15.8)
                  CALL PRINT_ERROR (31,pgname,localinfo,0)
                  istop=1
              ELSE
                  GOTO 10
              ENDIF
          ENDIF
      ihigh=ilow+1
      IF (ilow.EQ.0) THEN
          WRITE (localinfo,101) t
101       FORMAT (11x,'Routine: species_cv. Value:',d15.8)
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
      y1=REFTABLE(isp,2,ilow)
      y2=REFTABLE(isp,2,ihigh)
      y=(y2-y1)/(x2-x1)*(x-x1)+y1
!
      cv_reftable_x=y-RUNIV
      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION cv_reftable_y (isp,T,flg_stop)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER isp,flg_stop
      REAL(kind=8) t,cv,CV_REFTABLE_X
!
      cv=CV_REFTABLE_X (isp,T,flg_stop)
      cv_reftable_y=cv/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
