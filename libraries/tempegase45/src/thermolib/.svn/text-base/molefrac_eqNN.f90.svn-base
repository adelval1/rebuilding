! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE molefrac_mod                                         //
! //                                                                  //
! //  Input:  x_guess     array of mass fractions (initial guess)     //
! //          p           pressure                          kg/m3     //
! //          T           equilibrium temperature           K         //
! //          flg_guess   if 1, indicates an initial guess            //
! //          flg_anha    anharmonicity corrections flag needed to    //
! //                      compute the anharmonicity corrections       //
! //          flg_stop    if 1, stops on subroutine errors            //
! //          flg_termo   needed to know the mode of computation of   //
! //                      the equilibrium constants                   //
! //                                                                  //
! //  Output: x           array of mole fractions                     //
! //                                                                  //
! //  This subroutine is called to compute the mole fraction from the //
! //  pressure and temperature. The differencies between this and the //
! //  original subroutine lies in the system of equation that is      //
! //  solved: this subroutine in fact doesn't use the condition "sum  //
! //  of the mole fractions =1", because the Xc vector is computed by //
! //  means of the nuclei fraction equations(see Bottin thesis for    //
! //  better expl).                                                   //
! //                                                                  //
! //  The array passed to the newtonsolver routine is the logarithm   //
! //  of the mass fractions, so the initial guess is converted in     //
! //  this routine (or taken as 0 = log 1).                           //
! //                                                                  //
! //  Nicola Nannipieri 11/10/02.                                     //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE MOLEFRAC_MOD(x,p,T,TNEQ,eps,flg_log,&
     & flg_anha,flg_neq,flg_stop,flg_termo,flg_guess,x_guess,niter,X_C)
!
      USE global_thermo
      USE global_equilibrium
      IMPLICIT NONE
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      INTEGER i,flg_log,flg_stop,flg_anha,flg_neq,flg_termo,flg_guess,niter
      REAL(kind=8) x(nsp),x_guess(nsp),X_C(2) 
      REAL(kind=8) p,T,TNEQ(4),eps 
!
      pgname='pegaslib'
      
      !flg_termo=0         !NANNI
      
      !write(*,*)'guarda un po'
 !     pause
     ! write(*,*)'molefrac_mod started',T,niter,X_C(1),X_C(2)
      
      IF (flg_log.EQ.1) THEN
          IF (flg_neq.EQ.0) THEN
              WRITE (66,1000) p,T
          ELSE
              WRITE (66,1001) p,(TNEQ(i),i=1,4)
          END IF
      END IF
      
!      WRITE(*,*)'I am working in MASSFRAC',xini !pietro
1000  FORMAT (' Information relative to equilibrium composition ', &
        & 'calculation:',/,' ---------------------------------------', &
        & '--------------------:',/,                                   &
        & ' Density = ',e15.9,' kg/m^3,  Temperature = ',e15.9,' K',/)
1001  FORMAT (' Information relative to equilibrium composition ', &
        & 'calculation:',/,' ---------------------------------------', &
        & '--------------------:',/,                                   &
        & ' Density = ',e15.9,' kg/m^3,  Temperatures = ',4e15.9,' K',/)

      !write(*,*)'Before system_k',TNEQ,flg_termo
!
!     Calculation of equilibrium constants
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      
      CALL SYSTEM_K_MOD (1,T,TNEQ,p,flg_anha,flg_neq,flg_stop,flg_termo)
!
!     Calculation of initial values
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (flg_guess.NE.0) THEN
          DO i=1,nsp
              IF (x_guess(i).LE.0.0d0) THEN
                  WRITE (localinfo,101) i,x_guess(i)
101   FORMAT (11x,'Routine: massfrac. Index:',i3,' Value:',d15.8)
                  CALL PRINT_ERROR (36,pgname,localinfo,0)
                  x_guess(i)=1.0d0
              ENDIF
              v(i)=LOG(x_guess(i))
          END DO
      ELSE
          v=0.0
      ENDIF
!
      !write(*,*)'Before newtonsolver'


!     Calculation of mole fractions
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      CALL NEWTONSOLVER_MOD(1,niter,eps,flg_log,T,X_C)!,T) pietro 
!
!     Retrieval of mass fractions
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^      
      !write(*,*)'retrieval of mass fraction',v
      !write(*,*)'XCONS',XCONS
      !write(*,*)'YCONS',YCONS
      
      DO i=1,nsp
          !write(*,*)'Assigning the x(i)'
	  x(i)=EXP(v(i))
      END DO
!
      !write(*,*)'Before return x(i)',x
      RETURN
      END
! //////////////////////////////////////////////////////////////////////



! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE newtonsolver_mod                                     //
! //                                                                  //
! //  Input:  v           array of unknowns to be computed            //
! //          mode        indicates a mole or mass computation        //
! //          k           array of equilibrium constants              //
! //          iniguess    initial values of the unknowns              //
! //                                                                  //
! //  Output: v           correct values of the unknowns              //
! //                                                                  //
! //  This routine executes the multidimensional Newton-Raphson       //
! //  algorithm to solve the non-linear problem of the equilibrium    //
! //  composition.                                                    //
! //                                                                  //
! //  Corrected on 8/11/96 to make the convergence check on the v(i)  //
! //  instead of on the log v(i). Write v(i) and old test declaration //
! //  left in the code can be used for further debugging.             //
! //                                                                  //
! //  Nicola Nannipieri 11/10/02.                                     //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE NEWTONSOLVER_MOD(mode,niter,eps,flg_log,T,X_C)!,T) pietro
!     
      USE global_thermo
      USE global_equilibrium
      IMPLICIT NONE
!
      INTEGER i,j,mode,cnt_up,cnt_down,niter,flg_log,iiter
      REAL(kind=8) d,resid,eps,T,X_C(2) 
!     
!     initializations
!     ^^^^^^^^^^^^^^^
      v(nsp+1)=1.0d0
      iiter=0;cnt_up=0;cnt_down=0
      truevar(1:nsp)=DEXP(v(1:nsp))
      IF (flg_log.EQ.1) THEN
          WRITE (66,1000)
1000      FORMAT (' iiter          residual')
      END IF
!
!     come-back point for the iterative procedure during convergence
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      DO
          iiter=iiter+1
!
!         computation of main matrices: Jacobian, RHS terms and then
!         reduced system matrices
!         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          
	   !write(*,*)'Before calling LOWJAC_X_MOD'     !NANNI
	  IF (mode.EQ.1) THEN       
              CALL LOWJAC_X_MOD
              !write(*,*)'Before calling F_X_MOD',X_C     !NANNI
	      CALL F_X_MOD(X_C) 
          ELSE
              CALL LOWJAC_Y 
              CALL F_Y 
          ENDIF
          !write(*,*)'After calling LOWJAC_X_MOD'     !NANNI
	  
	  CALL MAINJAC

          !write(*,*)'After calling MAINJAC'     !NANNI
!
!         Compute the residual
!         ^^^^^^^^^^^^^^^^^^^^
          resid=0.0d0
          DO i=1,nc
              resid=resid+indmn(i)**2
          END DO
          resid=DSQRT(resid)
!
!         computation of the delta variables (reduced system)
!         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          CALL LUDCMP(jacmn(1:nc,1:nc),nc,nc,luindx(1:nc),d)
          CALL LUBKSB(jacmn(1:nc,1:nc),nc,nc,luindx(1:nc),indmn(1:nc))
          DO i=1,nc
              delv(nr+i)=indmn(i)
          END DO
!
!         Exact substitution of the linear part of the system
!         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          DO i=1,nr
              delv(i)=0.0d0
              DO j=1,nc
                  delv(i)=delv(i)-jac12(i,j)*delv(nr+j)
              END DO
              delv(i)=delv(i)-indsd(i)
          END DO
!
!         computation of new unknown values and apply bounding if
!         the new value is problematic (GT.1 or LT.limit value).
!         Finally, update truevar with the real mole or mass fraction
!         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	  !IF (T.LT.400) THEN !pietro
          !write(*,'(e14.9)')T !pietro
	!	pause
	  !v=v+0.1*delv !pietro
	  !ELSE
	  v=v+delv
	  !END IF

          !write(*,*)'Step 2'     !NANNI
	  
	  DO i=1,nsp
              IF (v(i).GT.0.0d0) THEN ! 0.0d0 invece di 1.0d0
                  v(i)=0.0d0  ! 0.0d0 invece di 1.0d0
                  cnt_up=cnt_up+1
              ENDIF
              IF (v(i).LT.-227.0d0) THEN!227 invece di 527
                  v(i)=-227.0d0!227 invece di 527
                  cnt_down=cnt_down+1
              ENDIF
          END DO
          truevar(1:nsp)=DEXP(v(1:nsp))
!
!         Printout of log information
!         ^^^^^^^^^^^^^^^^^^^^^^^^^^^
          IF (flg_log.EQ.1) THEN
              WRITE (66,1001) iiter,resid
1001          FORMAT (' ',i5,e18.10)
          END IF
!
!         Perform convergence tests
!         ^^^^^^^^^^^^^^^^^^^^^^^^^
          IF ((iiter.EQ.niter).OR.(resid.LE.eps)) THEN
              IF ((iiter.EQ.niter).AND.(flg_log.EQ.1)) THEN
                  WRITE (66,*) 'Exited after ',niter,' iterations'
              END IF
              EXIT 
          ENDIF
      END DO
!
      !write(*,*)'Before last if'
      IF (flg_log.EQ.1) THEN
          WRITE (66,*)
      END IF
      RETURN
      END
! //////////////////////////////////////////////////////////////////////



! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE f_x_mod                                              //
! //                                                                  //
! //  Input:  k           array of equilibrium constants Ky           //
! //          jaclow      lower jacobian matrix for equ. computation  //
! //          v           array of unknowns                           //
! //                                                                  //
! //  Output: f1          array of upper right hand side vector       //
! //          f2          array of lower right hand side vector       //
! //                                                                  //
! //  Computation of the "function" terms (right hand side) in the    //
! //  case of mass fraction computation.                              //
! //                                                                  //
! //  Nicola Nannipieri 11/10/02.                                     //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE F_X_MOD(X_C)
!
      USE global_thermo
      USE global_equilibrium
      IMPLICIT NONE
      INTEGER i,j
      REAL (kind=8) X_C(2)
!

      DO i=1,nr
          f1(i)=0.0d0
          DO j=1,nsp
              f1(i)=f1(i)+NU(i,j)*v(j)
          END DO
          f1(i)=f1(i)-sys_keq(i)
      END DO
      DO i=1,nc
          f2(i)=0.0d0
          DO j=1,nsp
              f2(i)=f2(i)+jaclow(i,j)
          END DO
          f2(i)=f2(i)-X_C(i)
      END DO
      RETURN
      END
      
      ! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE lowjac_x                                             //
! //                                                                  //
! //  Input:  v           array of unknowns                           //
! //                                                                  //
! //  Output: jaclow      lower jacobian matrix for equ. computation  //
! //                                                                  //
! //  Computation of the lower part of the jacobian matrix in the     //
! //  case of mole fraction computation.                              //
! //                                                                  //
! //  Nicola Nannipieri 11/10/02.                                     //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE LOWJAC_X_MOD
!      
      USE global_thermo
      USE global_equilibrium
      IMPLICIT NONE
      INTEGER i,j
!      
!      DO j=1,nsp
!          truevar(j)=DEXP(v(j))!+1.0d-16
!      END DO
     
      !write(*,*)'LOWJAC_X_MOD started'
     
      DO i=1,nc
          DO j=1,nsp
              !write(*,*)'before jaclow',KHI(i,j)
	      jaclow(i,j)=KHI(i,j)*truevar(j)
          END DO
!          jaclow(i,nsp+1)=XCONS(i)
      END DO
!      DO i=1,nsp
!          jaclow(nc,i)=truevar(i)
!      END DO
!      jaclow(nc,nsp+1)=0.0d0
      !write(*,*)'Before return in lowjac_x_mod'
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE system_k                                             //
! //                                                                  //
! //  Input:  nsp         number of species of the problem            //    
! //          nr          number of reactions of the problem          //
! //          mode        type of computation p,T or rho,T            //
! //          T           temperature                       K         //
! //          varia       pressure or density               Pa,kg/m3  //
! //          flg_anha    anharmonicity corrections flag needed to    //
! //                      compute the anharmonicity corrections       //
! //          flg_stop    if 1, stops on subroutine errors            //
! //          flg_termo   needed to know the mode of computation of   //
! //                      the equilibrium constants                   //
! //                                                                  //
! //  Output: k(maxr)     array of logarithms of Kx or Ky (as mode)   //
! //                                                                  //
! //  This subroutine is called to compute the natural logarithms of  //
! //  the equilibrium constants involved in the system of equilibrium //
! //  composition solved by Pegase.                                   //
! //  Notes: nusum is the sum of the stoechiometric coefficients of   //
! //         the equation that is being processed                     //
! //         sum is used to build the sum of products NU*muo          //
! //         muo is the chemical potential divided by RUNIV T         //
! //                                                                  //
! //  Benoit Bottin, 09/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE SYSTEM_K_MOD (mode,T,TNEQ,varia,&
     &flg_anha,flg_neq,flg_stop,flg_termo)

      USE global_thermo
      USE global_equilibrium
      IMPLICIT NONE
      INTEGER mode,r,s,i,nusum,flg_anha,flg_neq,flg_stop,flg_termo
      REAL(kind=8) T,varia,summ,msum
      REAL(kind=8) TNEQ(4),po,pot(8)
      REAL(kind=8),ALLOCATABLE :: muo(:)
!      
      ALLOCATE (muo(1:nsp))
      muo=0.0
      po=1.0d0
      
!      DO j=1,4        !NANNI
!          TNEQ(j)=T
!      END DO

!      
      !write(*,*)'Inside system_k: before calling species_gibbs',nsp,&
      !&          T,TNEQ,po,pot   !NANNI
      DO s=1,nsp
          CALL SPECIES_GIBBS_MOD(s,po,T,TNEQ,flg_anha,1,flg_neq,flg_stop,&
     &    flg_termo,pot)
          muo(s)=(pot(1)+pot(8))/RUNIV/T
      END DO
      
      !write(*,*)'Inside system_k: after calling species_gibbs'  !NANNI
      
      DO r=1,nr
          summ=0
          nusum=0
          DO s=1,nsp
              nusum=nusum+NU(r,s)
              summ=summ+NU(r,s)*muo(s)
          END DO
          IF (mode.EQ.1) THEN
              !write(*,*)'Inside system_k: before sys_keq'  !NANNI
	      sys_keq(r)=-summ-nusum*LOG(varia)
          ELSE
              msum=0
              DO i=1,nsp
                  msum=msum+NU(r,i)*LOG(MMOL(i))
              END DO
              sys_keq(r)=-summ-nusum*LOG(varia*RUNIV*T)+msum
          ENDIF
      END DO
!
      !write(*,*)'Inside system_k: before deallocation'  !NANNI
      DEALLOCATE (muo)
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
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
      SUBROUTINE species_gibbs_mod (isp,var,T,TNEQ,&
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
      
      !write(*,*)'step 1'             !NANNI
      
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

      !write(*,*)'Step 1.2'          !NANNI
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
            !write(*,*)'Step 1.2.1'          !NANNI 
	      gibb(1)=GIBBS_REFTABLE_X (isp,var,t,flg_stop)
            !write(*,*)'where exactly?'          !NANNI   
	      gibb(8)=HFOR(isp)
          write(*,*)'Step 1.3'          !NANNI
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
           !write(*,*)'Step 1.4'          !NANNI   
	      ELSEIF (flg_anha.EQ.1) THEN
                  gibb(5)=GIBBS_ELEC_X(isp,tele)
                  IF (atom.EQ.0) THEN
                      gibb(6)=GIBBS_RV_X(isp,trot,tvib)
                  ENDIF
                  gibb(1)=gibb(2)+gibb(5)+gibb(6)
              !write(*,*)'Step 1.5'               !NANNI
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
      
      !write(*,*)'Step 2'           !NANNI
      
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
