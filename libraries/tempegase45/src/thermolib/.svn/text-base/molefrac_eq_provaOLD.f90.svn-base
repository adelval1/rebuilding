
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
      SUBROUTINE MOLEFRAC_MOD2(x,p,T,zo,TNEQ,eps,flg_log,&
     &flg_anha,flg_neq,flg_stop,flg_termo,flg_guess,x_guess,niter,X_C)
!
      USE global_thermo
      USE global_equilibrium
      IMPLICIT NONE
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      INTEGER i,flg_log,flg_anha,flg_neq,flg_stop,flg_termo,flg_guess,niter
      REAL(kind=8) x(nsp),p,T,x_guess(:),zo,TNEQ(4),eps,X_C(2)
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
!
!     Calculation of equilibrium constants
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      
      write(*,*)'Before system_k_mod2'
      CALL SYSTEM_K_MOD2 (1,T,TNEQ,p,flg_anha,flg_neq,flg_stop,flg_termo)
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
!
!     Calculation of mole fractions
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      write(*,*)'Before newtonsolver_x_mod2'
      nc=nc+1
      CALL NEWTONSOLVER_MOD2(1,niter,eps,flg_log,T,X_C) !,T) pietro
      nc=nc-1
!
!     Retrieval of mass fractions
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^      
      write(*,*)'retrieval of mass fraction',v
      write(*,*)'XCONS',XCONS
     
      DO i=1,nsp
          write(*,*)'Assigning the x(i)'
	  x(i)=EXP(v(i))
      END DO
      zo=1/v(nsp+1)
!
      write(*,*)'Before return',x
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
      SUBROUTINE SYSTEM_K_MOD2 (mode,T,TNEQ,varia,&
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
!      
      write(*,*)'Inside system_k: Before species_gibbs',flg_termo,T,TNEQ,&
     & po 
      DO s=1,nsp
          CALL SPECIES_GIBBS(s,po,T,TNEQ,flg_anha,1,flg_neq,flg_stop,&
     &    flg_termo,pot)
          muo(s)=(pot(1)+pot(8))/RUNIV/T
      END DO
      write(*,*)'After species_gibbs'
      
      DO r=1,nr
          summ=0
          nusum=0
          DO s=1,nsp
              nusum=nusum+NU(r,s)
              summ=summ+NU(r,s)*muo(s)
          END DO
        write(*,*)'Inside system_k: before IF mode'  
	  IF (mode.EQ.1) THEN
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
      DEALLOCATE (muo)
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE fixedjac                                             //
! //                                                                  //
! //  Computation of the fixed parts of the jacobian matrices in the  //
! //  equilibrium composition calculation. This is computed once and  //
! //  for all before iterating for the composition.                   //
! //                                                                  //
! //  Benoit Bottin, 09/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE FIXEDJAC_MOD2
!
      USE global_thermo
      USE global_equilibrium
      IMPLICIT NONE
      INTEGER i,j
!      
      DO i=1,nr
          jac11(i,i)=1.0d0/NU(i,i)
      END DO
!
      DO i=1,nr
          DO j=1,nc
              jac12(i,j)=jac11(i,i)*NU(i,j+nr)
          END DO
          jac12(i,nc+1)=0.0
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
! //  Benoit Bottin, 09/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE LOWJAC_X_MOD2(X_C)
!      
      USE global_thermo
      USE global_equilibrium
      IMPLICIT NONE
      INTEGER i,j
      REAL(kind=8) X_C(2) 
!      
!      DO j=1,nsp
!          truevar(j)=DEXP(v(j))!+1.0d-16
!      END DO
      DO i=1,nc-1
          DO j=1,nsp
              jaclow(i,j)=KHI(i,j)*truevar(j)
          END DO
          jaclow(i,nsp+1)=X_C(i)   !NANNI
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
      SUBROUTINE F_X_MOD2(X_C)
!
      USE global_thermo
      USE global_equilibrium
      IMPLICIT NONE
      INTEGER i,j
      REAL(kind=8) X_C(2)
!
      DO i=1,nr
          f1(i)=0.0d0
          DO j=1,nsp
              f1(i)=f1(i)+NU(i,j)*v(j)
          END DO
          f1(i)=f1(i)-sys_keq(i)
      END DO
      
      write(*,*)'Inside f_x_mod2: before do 1,nc',f1
      
      DO i=1,nc
          f2(i)=0.0d0
          DO j=1,nsp
              f2(i)=f2(i)+jaclow(i,j)
          END DO
      END DO
      DO i=1,nc-1
          f2(i)=f2(i)-X_C(i)*v(nsp+1)      !NANNI
      END DO
      f2(nc)=f2(nc)-1.0d0
      
      write(*,*)'Inside f_x_mod2: after do 1,nc',f2

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
      SUBROUTINE NEWTONSOLVER_MOD2(mode,niter,eps,flg_log,T,X_C)!,T) pietro
!     
      USE global_thermo
      USE global_equilibrium
      IMPLICIT NONE
!
      INTEGER i,j,mode,cnt_up,cnt_down,niter,flg_log,iiter
      REAL(kind=8) d,resid,eps,T,X_C(2) !,T) pietro
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
              write(*,*)'Before lowjac_x_mod2'
	      CALL LOWJAC_X_MOD2(X_C)
              write(*,*)'Before f_x_mod2'
	      CALL F_X_MOD2(X_C) 
          ELSE
              CALL LOWJAC_Y 
              CALL F_Y 
          ENDIF
          write(*,*)'Before calling mainjac'
	  
	  CALL MAINJAC
          
          write(*,*)'Before computing the residuals'
!
!         Compute the residual
!         ^^^^^^^^^^^^^^^^^^^^
          resid=0.0d0
          write(*,*)'nc and indmn',nc,indmn
	  DO i=1,nc
              write(*,*)'indmn',i
	      resid=resid+indmn(i)**2
          END DO
          resid=DSQRT(resid)
          write(*,*)'Inside newtonsolver_mod2: Before ludcmp'

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
          write(*,*)'Inside newtonsolver_mod2: After lubksb'

	  
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

          write(*,*)'Inside newtonsolver_mod2: Somewhere..',v

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

       write(*,*)'Inside newtonsolver_mod2: Before convergence tests'

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

      write(*,*)'Inside newtonsolver_mod2: Before last if'

!
      IF (flg_log.EQ.1) THEN
          WRITE (66,*)
      END IF
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
    
