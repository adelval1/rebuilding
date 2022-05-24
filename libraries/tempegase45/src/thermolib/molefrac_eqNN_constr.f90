
! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE molefrac                                             //
! //                                                                  //
! //  Input:  x_guess     array of mass fractions (initial guess)     //
! //          p           pressure                          Pa        //
! //          T           equilibrium temperature           K         //
! //          flg_guess   if 1, indicates an initial guess            //
! //          flg_anha    anharmonicity corrections flag needed to    //
! //                      compute the anharmonicity corrections       //
! //          flg_stop    if 1, stops on subroutine errors            //
! //          flg_termo   needed to know the mode of computation of   //
! //                      the equilibrium constants                   //
! //                                                                  //
! //  Output: x           array of mass fractions                     //
! //                                                                  //
! //  This subroutine is called to compute the mole fraction from the //
! //  pressure and temperature. An initial guess array can be input   //
! //  as a starting point in recurring applications (cfd code). If    //
! //  the initial array must be used, set flg_guess to 1 in the       //
! //  calling routine, if not set it to zero. If no initial guess is  //
! //  used, the initial value is 1 for every mole fraction.           //
! //                                                                  //
! //  The array passed to the newtonsolver routine is the logarithm   //
! //  of the mole fractions, so the initial guess is converted in     //
! //  this routine (or taken as 0 = log 1).                           //
! //                                                                  //
! //  Note: in mole fraction formulation, nc must be equal to the     //
! //  real nc + 1 to account for the additional equation of number of //
! //  moles per mole of cold gas. Thus nc is increased at the start   //
! //  of the routine and decreased again at the end. Note that the    //
! //  mole ratio is returned as x(nsp+1).                             //
! //                                                                  //
! //  Benoit Bottin, 11/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE MOLEFRAC_MOD_CONSTR(x,p,T,zo,TNEQ,eps,flg_log,&
     &flg_anha,flg_neq,flg_stop,flg_termo,flg_guess,x_guess,niter,XC)
!
      USE global_thermo
      USE global_equilibrium
      IMPLICIT NONE
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      INTEGER i,flg_log,flg_anha,flg_neq,flg_stop,flg_termo,flg_guess,niter
      REAL(kind=8) x(nsp),p,T,x_guess(nsp),zo,TNEQ(4),eps
      REAL(kind=8) XC(2)
!
      
      pgname='pegaslib'
      IF (flg_log.EQ.1) THEN
          IF (flg_neq.EQ.0) THEN
              WRITE (66,1000) p,T
          ELSE
              WRITE (66,1001) p,(TNEQ(i),i=1,4)
          END IF
      END IF
!      WRITE(*,*)'I am working in MOLEFRAC',xini !pietro
1000  FORMAT (' Information relative to equilibrium composition ', &
        & 'calculation:',/,' ---------------------------------------', &
        & '--------------------:',/,                                   &
        & ' Pressure = ',e15.9,' Pa,  Temperature = ',e15.9,' K',/)
1001  FORMAT (' Information relative to equilibrium composition ', &
        & 'calculation:',/,' ---------------------------------------', &
        & '--------------------:',/,                                   &
        & ' Pressure = ',e15.9,' Pa,  Temperatures = ',4e15.9,' K',/)

        !write(*,*)'Before calling system_k'

!
!     Calculation of equilibrium constants
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      
      CALL SYSTEM_K (1,T,TNEQ,p,flg_anha,flg_neq,flg_stop,flg_termo)
!
!     Calculation of initial values
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (flg_guess.EQ.0) THEN
          DO i=1,nsp
              v(i)=0.0d0
          END DO
      ELSE
          DO i=1,nsp
              IF (x_guess(i).LE.0.0d0) THEN
                  WRITE (localinfo,101) i,x_guess(i)
101   FORMAT (11x,'Routine: molefrac. Index:',i3,' Value:',d15.8)
                  CALL PRINT_ERROR (36,pgname,localinfo,0)
                  x_guess(i)=1.0d0
              ENDIF
              v(i)=LOG(x_guess(i))
          END DO
      ENDIF

      !write(*,*)'Before calling newtonsolver'
!      pause
!
!     Calculation of mole fractions
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      nc=nc+1
      CALL NEWTONSOLVER_MOD_CONSTR(1,niter,eps,flg_log,T,XC) !,T) pietro
      nc=nc-1

      !write(*,*)'After calling newtonsolver'

      
!
!     Retrieval of mass fractions
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^      
      DO i=1,nsp
          x(i)=EXP(v(i))
      END DO
      zo=1/v(nsp+1)

      !write(*,*)'concentration (5 values) and zo',x,zo

!
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
! //  Benoit Bottin, 09/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE LOWJAC_X_MOD_CONSTR(XC)
!      
      USE global_thermo
      USE global_equilibrium
      IMPLICIT NONE
      INTEGER i,j
      REAL(kind=8) XC(2)
!      
!      DO j=1,nsp
!          truevar(j)=DEXP(v(j))!+1.0d-16
!      END DO
      DO i=1,nc-1
          DO j=1,nsp
!              write(*,*)'Assigning KHI',KHI(i,j),i,j
	      jaclow(i,j)=KHI(i,j)*truevar(j)
          END DO
!         write(*,*)'Assigning -XC',XC
	  jaclow(i,nsp+1)=-XC(i)
      END DO
      DO i=1,nsp
          jaclow(nc,i)=truevar(i)
      END DO
      jaclow(nc,nsp+1)=0.0d0
      RETURN
      END
! //////////////////////////////////////////////////////////////////////


! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE f_x                                                  //
! //                                                                  //
! //  Input:  k           array of equilibrium constants Kx           //
! //          jaclow      lower jacobian matrix for equ. computation  //
! //          v           array of unknowns                           //
! //                                                                  //
! //  Output: f1          array of upper right hand side vector       //
! //          f2          array of lower right hand side vector       //
! //                                                                  //
! //  Computation of the "function" terms (right hand side) in the    //
! //  case of mole fraction computation.                              //
! //                                                                  //
! //  Benoit Bottin, 09/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE F_X_MOD_CONSTR(XC)
!
      USE global_thermo
      USE global_equilibrium
      IMPLICIT NONE
      INTEGER i,j
      REAL(kind=8) XC(2)

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
      END DO
      DO i=1,nc-1
          f2(i)=f2(i)-XC(i)*v(nsp+1)
      END DO
      f2(nc)=f2(nc)-1.0d0
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////


! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE newtonsolver                                         //
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
! //  Benoit Bottin, 09/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE NEWTONSOLVER_MOD_CONSTR(mode,niter,eps,flg_log,T,XC)!,T) pietro
!     
      USE global_thermo
      USE global_equilibrium
      IMPLICIT NONE
!
      INTEGER i,j,mode,cnt_up,cnt_down,niter,flg_log,iiter
      REAL(kind=8) d,resid,eps,T,XC(2) !,T) pietro
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
          IF (mode.EQ.1) THEN       
              !write(*,*)'Before calling LOWJAC_X_MOD_CONSTR'
	      CALL LOWJAC_X_MOD_CONSTR(XC)
              !write(*,*)'Before calling F_X_MOD_CONSTR'
	      CALL F_X_MOD_CONSTR(XC) 
          ELSE
              CALL LOWJAC_Y 
              CALL F_Y 
          ENDIF
          CALL MAINJAC

         !write(*,*)'Before computing the residuals'
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
      IF (flg_log.EQ.1) THEN
          WRITE (66,*)
      END IF
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
    
