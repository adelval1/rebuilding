! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE species_energy                                       //
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
! //  Output: ener(8)      array of internal energies                 //
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
! //  Benoit Bottin, 11/10/96. Modified and cleaned 25/7/97.          //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE species_energy (isp,T,TNEQ,&
     &flg_anha,mode,flg_neq,flg_stop,flg_termo,ener)
!
      USE global_thermo
      IMPLICIT NONE
      INTERFACE
        REAL(kind=8) FUNCTION mol_to_mass (val,isp,mm,flg_stop)
            REAL(kind=8) val,mm
            INTEGER isp,flg_stop
        END FUNCTION mol_to_mass
      END INTERFACE
!
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      INTEGER isp,flg_anha,mode,flg_neq,flg_stop,flg_termo,istop,i,atom
      REAL(kind=8) T,TNEQ(4),ener(8),ttra,trot,tvib,tele
      REAL(kind=8) ENERGY_TRANS_X,ENERGY_TRANS_Y
      REAL(kind=8) ENERGY_ROT_X,ENERGY_ROT_Y
      REAL(kind=8) ENERGY_VIB_X,ENERGY_VIB_Y
      REAL(kind=8) ENERGY_ELEC_X,ENERGY_ELEC_Y
      REAL(kind=8) ENERGY_RV_X,ENERGY_RV_Y
      REAL(kind=8) ENERGY_INT_X,ENERGY_INT_Y
      REAL(kind=8) ENERGY_REFTABLE_X,ENERGY_REFTABLE_Y
!
!     Checks of correct variable values and initialize
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !WRITE(*,*)'inizio' !pietro
      pgname='pegaslib'
      istop=0
      IF (flg_neq.EQ.0) THEN
          IF (T.LE.0) THEN
              WRITE (localinfo,101) T
101           FORMAT (11x,'Routine: species_energy. Value:',d15.8)
              CALL PRINT_ERROR (24,pgname,localinfo,0)
              istop=1
          ENDIF
	  IF (T.LE.100.) THEN !pietro
	  write(*,*)'T',T
	  !T=300.!pietro
	  ENDIF !pietro
          ttra=T
          trot=T
          tvib=T
          tele=T
      ELSE
          DO i=1,4
              IF (TNEQ(i).LE.0) THEN
                  WRITE (localinfo,102) i,TNEQ(i)
102   FORMAT (11x,'Routine: species_energy. Index:',i3,' Value:',d15.8)
                  CALL PRINT_ERROR (27,pgname,localinfo,0)
                  istop=1
              ENDIF
          END DO
          ttra=TNEQ(1)
          trot=TNEQ(2)
          tvib=TNEQ(3)
          tele=TNEQ(4)
	  IF (TNEQ(1).LE.100.) THEN !pietro
	  write(*,*)'TNEQ',TNEQ!pietro
	  ENDIF !pietro
      ENDIF
      DO i=1,8
          ener(i)=0.0d0
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
                  localinfo(1:35)='           Routine: species_energy.'
                  CALL PRINT_ERROR (28,pgname,localinfo,0)
                  t=ttra
              ENDIF
              ener(1)= ENERGY_REFTABLE_X (isp,t,flg_stop)
              ener(8)=HFOR(isp)
          ELSE
              ener(2)=ENERGY_TRANS_X(isp,ttra)
              ener(8)=HFOR(isp)
              IF (flg_anha.EQ.0) THEN
                  IF (atom.EQ.0) THEN
		      !WRITE(*,*)'ciao_00' !pietro
                      ener(3)=ENERGY_ROT_X(isp,trot)
                      !WRITE(*,*)'ciao_0' !pietro
		      ener(4)=ENERGY_VIB_X(isp,tvib)
                  ENDIF
                  ener(5)=ENERGY_ELEC_X(isp,tele)
                  ener(1)=ener(2)+ener(3)+ener(4)+ener(5)
              ELSEIF (flg_anha.EQ.1) THEN
                  ener(5)=ENERGY_ELEC_X(isp,tele)
                  IF (atom.EQ.0) THEN
                      ener(6)=ENERGY_RV_X(isp,trot,tvib)
                  ENDIF
                  ener(1)=ener(2)+ener(5)+ener(6)
              ELSEIF (flg_anha.GE.2) THEN
                  IF (atom.EQ.0) THEN
                      ener(7)=ENERGY_INT_X(isp,trot,tvib,tele)
                  ELSE
                      ener(7)=ENERGY_ELEC_X(isp,tele)
                  ENDIF
                  ener(1)=ener(2)+ener(7)
              ELSE
                  localinfo(1:35)='           Routine: species_energy.'
                  CALL PRINT_ERROR (29,pgname,localinfo,0)
                  istop=1
              ENDIF
          ENDIF
      ELSE
          IF (flg_termo.EQ.1) THEN
              IF (flg_neq.NE.0) THEN
                  localinfo(1:35)='           Routine: species_energy.'
                  CALL PRINT_ERROR (28,pgname,localinfo,0)
                  t=ttra
              ENDIF
              ener(1)= ENERGY_REFTABLE_Y (isp,t,flg_stop)
              ener(8)=MOL_TO_MASS(HFOR(isp),isp,0.0d0,flg_stop)
          ELSE
              ener(2)=ENERGY_TRANS_Y(isp,ttra)
              ener(8)=MOL_TO_MASS(HFOR(isp),isp,0.0d0,flg_stop)
              IF (flg_anha.EQ.0) THEN
                  IF (atom.EQ.0) THEN
                      ener(3)=ENERGY_ROT_Y(isp,trot)
                      ener(4)=ENERGY_VIB_Y(isp,tvib)
                  ENDIF
                  ener(5)=ENERGY_ELEC_Y(isp,tele)
                  ener(1)=ener(2)+ener(3)+ener(4)+ener(5)
              ELSEIF (flg_anha.EQ.1) THEN
                  ener(5)=ENERGY_ELEC_Y(isp,tele)
                  IF (atom.EQ.0) THEN
                      ener(6)=ENERGY_RV_Y(isp,trot,tvib)
                  ENDIF
                  ener(1)=ener(2)+ener(5)+ener(6)
              ELSEIF (flg_anha.GE.2) THEN
                  IF (atom.EQ.0) THEN
                      ener(7)=ENERGY_INT_Y(isp,trot,tvib,tele)
                  ELSE
                      ener(7)=ENERGY_ELEC_Y(isp,tele)
                  ENDIF
                  ener(1)=ener(2)+ener(7)
              ELSE
                  localinfo(1:35)='           Routine: species_energy.'
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
! //  FUNCTIONS energy_trans_x AND energy_trans_y                     //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          tt          translation temperature (K)                 //
! //                                                                  //
! //  This function returns the translation internal energy per       //
! //  unit mole (p,T input) or unit mass (rho,T input).               //
! //                                                                  //
! //  Benoit Bottin, 11/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION energy_trans_x (isp,tt)

      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp,i
      REAL(kind=8) :: tt,ener

      ener=1.5d0*RUNIV*tt
      energy_trans_x=ener

!     isp is passed to keep the similarity between all calls
!     following assignement removes any compiling error
      i=isp      

      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION energy_trans_y (isp,tt)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: tt,ener
      REAL(kind=8) :: ENERGY_TRANS_X

      ener=ENERGY_TRANS_X (isp,tt)
      energy_trans_y=ener/MMOL(isp)

      RETURN
      END
! //////////////////////////////////////////////////////////////////////
      
! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS energy_rot_x AND energy_rot_y                         //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          trot        rotational temperature (K)                  //
! //                                                                  //
! //  This function returns the rotation internal energy per          //
! //  unit mole (p,T input) or unit mass (rho,T input).               //
! //                                                                  //
! //  Benoit Bottin, 11/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION energy_rot_x (isp,trot)

      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp,i
      REAL(kind=8) :: ener,trot

      ener=RUNIV*trot*(1.0d0-TR(1,isp)/(TR(1,isp)+3.0d0*trot))
      energy_rot_x=ener
!
!     isp is passed to keep the similarity between all calls
!     following assignement removes any compiling error
      i=isp

      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION energy_rot_y (isp,trot)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: trot,ener,ENERGY_ROT_X
!
      ener=ENERGY_ROT_X (isp,trot)
      energy_rot_y=ener/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
      
! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS energy_vib_x AND energy_vib_y                         //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          tvib        vibrational temperature (K)                 //
! //                                                                  //
! //  This function returns the vibration internal energy per         //
! //  unit mole (p,T input) or unit mass (rho,T input).               //
! //                                                                  //
! //  Benoit Bottin, 11/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION energy_vib_x (isp,tvib)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp,i
      REAL(kind=8) :: ener,tvib,u
!
      !WRITE(*,*)'ciao_1'!pietro
      ener=0.0d0        
      DO i=1,nvibmode(isp)
      u=TV(i,1,isp)/tvib
      !write(*,*)i,tvib,u !pietro
      !write(*,*)'tvib',tvib!pietro
      IF (u.LT.307) THEN
          ener=ener+RUNIV*TV(i,1,isp)/(EXP(u)-1.0d0)
	  !WRITE(*,*)'ciao_2'!pietro
      ENDIF
      END DO 
      energy_vib_x=ener
!
      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION energy_vib_y (isp,tvib)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: tvib,ener,ENERGY_VIB_X
      !WRITE(*,*)
      ener=ENERGY_VIB_X (isp,tvib)
      energy_vib_y=ener/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
      
! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS energy_elec_x AND energy_elec_y                       //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          p           pressure (Pa)                               //
! //          rho         density (kg/m3)                             //
! //          tele        electronic temperature (K)                  //
! //                                                                  //
! //  This function returns the electronic internal energy per        //
! //  unit mole (p,T input) or unit mass (rho,T input).               //
! //                                                                  //
! //  Benoit Bottin, 11/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION energy_elec_x (isp,tele)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp,i
      REAL(kind=8) :: ener,tele,q,num,e,ge_exp_e
!
      ener=0.0d0
      q=0.0d0
      num=0.0d0
      DO i=1,NELVL(isp)
          e=TE(i,isp)/tele
          IF (e.LT.307) THEN
              ge_exp_e=GE(i,isp)*EXP(-e)
              q=q+ge_exp_e
              num=num+TE(i,isp)*ge_exp_e
          ENDIF
      END DO
      ener=RUNIV*num/q
      energy_elec_x=ener
!
      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION energy_elec_y (isp,tele)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: tele,ener,ENERGY_ELEC_X
!
      ener=ENERGY_ELEC_X (isp,tele)
      energy_elec_y=ener/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS energy_rv_x AND energy_rv_y                           //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          p           pressure (Pa)                               //
! //          rho         density (kg/m3)                             //
! //          trot        rotational temperature (K)                  //
! //          tvib        vibrational temperature (K)                 //
! //                                                                  //
! //  This function returns the rotation-vibration coupled internal   //
! //  energy per unit mole (p,T input) or unit mass (rho,T input).    //
! //                                                                  //
! //  Benoit Bottin, 11/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION energy_rv_x (isp,trot,tvib)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: trot,tvib,d,g,x,u,n,e,den,expu
      REAL(kind=8) :: ENERGY_ROT_X,ENERGY_VIB_X
!
      g=G_A(1,isp)
      d=D_A(1,isp)
      x=X_A(1,isp)
      u=TV(1,1,isp)/tvib
      n=TR(1,isp)/trot
      expu=EXP(u)
      den=expu-1.0d0
!
      e=ENERGY_ROT_X(isp,trot)+ENERGY_VIB_X(isp,tvib)
      e=e+RUNIV*trot/n*g
      e=e+((d*expu-2.0d0*x)/den**2)*RUNIV*TV(1,1,isp)
      e=e+(4.0d0*x*expu/den**3)*RUNIV*u*TV(1,1,isp)
      !WRITE(*,*)'ciao_3'!pietro
      energy_rv_x=e

      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION energy_rv_y (isp,trot,tvib)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: trot,tvib,ener,ENERGY_RV_X
!
      ener=ENERGY_RV_X (isp,trot,tvib)
      energy_rv_y=ener/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS energy_rv_x AND energy_rv_y                           //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          p           pressure (Pa)                               //
! //          rho         density (kg/m3)                             //
! //          trot        rotational temperature (K)                  //
! //          tvib        vibrational temperature (K)                 //
! //          tele        electronic temperature (K)                  //
! //                                                                  //
! //  This function returns the internally-coupled internal           //
! //  energy per unit mole (p,T input) or unit mass (rho,T input).    //
! //                                                                  //
! //  Benoit Bottin, 11/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION energy_int_x (isp,trot,tvib,tele)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp,i
      REAL(kind=8) :: trot,tvib,tele,a,b,c,da,db,dc,e,u,n,den,enern,enerd
      REAL(kind=8) :: D,X,expu,exp_u,exp_e
!
      enern=0.0d0
      enerd=0.0d0
      DO i=1,NELVL(isp)
          IF (TV(1,i,isp).GT.1.0d-6) THEN
             !WRITE(*,*)'ciao_4'!pietro
	     e=TE(i,isp)/tele
             u=TV(1,i,isp)/tvib
             n=TR(i,isp)/trot
             expu=EXP(u)
             exp_u=1.0d0/expu
             exp_e=EXP(-e)
             den=expu-1.0d0
             D=D_A(i,isp)
             X=X_A(i,isp)
!
!             Electronic contribution
             a=GE(i,isp)*exp_e
             da=GE(i,isp)*e/tele*exp_e
!
!             Rotational contribution
             b=SIG(isp)*n*(1.0d0-exp_u)
             db=SIG(isp)*n/trot*(-1.0d0+exp_u*(1.0d0-u))
!
!             RR-HA contribution
             c=1.0d0+n/3.0d0+n**2/15.0d0
             c=c+G_A(i,isp)/n
             dc=-n/3.0d0/trot-2.0d0*n**2/15.0d0/trot
             dc=dc+G_A(i,isp)/n/trot
!
!             Mayer and Mayer Anharmonicity corrections
             c=c+D/den+2.0d0*X*u/(den**2)
             dc=dc+(D*u*expu-2.0d0*X*u)/tvib/den**2
             dc=dc+4.0d0*X*u**2*expu/tvib/den**3
!
!             Total
             enern=enern+((da*c+a*dc)*b-a*c*db)/b/b
             enerd=enerd+a*c/b
          ENDIF
      END DO
!
      energy_int_x=RUNIV*trot**2*enern/enerd

      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION energy_int_y (isp,trot,tvib,tele)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: trot,tvib,tele,ener,ENERGY_INT_X
!
      ener=ENERGY_INT_X (isp,trot,tvib,tele)
      energy_int_y=ener/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS energy_reftable_x and energy_reftable_y               //
! //                                                                  //
! //  Input:  T           temperature of the species (K)              //    
! //          isp         number of the concerned species             //
! //          flg_stop    stops if errors in the subroutine           //
! //                                                                  //
! //  This function returns the internal energy of the species from   //
! //  reference tables such as Gurvich or Janaf (mole or mass input)  //
! //                                                                  //
! //  Benoit Bottin, 11/12/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION energy_reftable_x (isp,T,flg_stop)
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
102               FORMAT (11x,'Routine: species_energy. Value:',d15.8)
                  CALL PRINT_ERROR (31,pgname,localinfo,0)
                  istop=1
              ELSE
                  GOTO 10
              ENDIF
          ENDIF
      ihigh=ilow+1
      IF (ilow.EQ.0) THEN
          WRITE (localinfo,101) t
101       FORMAT (11x,'Routine: species_energy. Value:',d15.8)
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
      y1=REFTABLE(isp,4,ilow)
      y2=REFTABLE(isp,4,ihigh)
      y=(y2-y1)/(x2-x1)*(x-x1)+y1
!
      energy_reftable_x=y*1000.0d0-RUNIV*t
      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      REAL(kind=8) FUNCTION energy_reftable_y (isp,T,flg_stop)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER isp,flg_stop
      REAL(kind=8) t,ener,ENERGY_REFTABLE_X
!
      ener=0.0d0
      ener=ENERGY_REFTABLE_X (isp,T,flg_stop)
      energy_reftable_y=ener/MMOL(isp)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
