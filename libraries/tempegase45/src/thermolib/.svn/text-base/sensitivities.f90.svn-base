! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE sensitivities_r                                      //
! //                                                                  //
! //  Input:  mm          molar mass of the mixture (kg/mol)          //
! //          p           pressure of the mixture (Pa)                //
! //          rho         density of the mixture (kg/m3)              //
! //          T           equilibrium temperature (K)                 //
! //          TNEQ        non-equilibrium temperature (K) (not used)  //
! //          y           array of mass fractons                      //
! //          ener(8)     array of energies per mass                  //
! //          enth(8)     array of enthalpies per mass                //
! //          epsilon     relative precision of the finite difference //
! //                                                                  //
! //  Flags:  flg_anha    anharmonicity corrections flag              //
! //          flg_neq     nonequilibrium temperatures flag            //
! //          flg_stop    stops if errors in subroutines              //
! //          flg_termo   statistical or reftable computation         //
! //                                                                  //
! //  Output: dedr (8)    de over dr at constant T   (energy)         //
! //          dhdr (8)    dh over dr at constant T   (enthalpy)       //
! //          detdr(8)    det over dr at constant T   (vol. energy)   //
! //          dpdr        dp over dr at constant T   (pressure)       //
! //          dmdr        dm over dr at constant T   (molar mass)     //
! //                                                                  //
! //  This subroutine provides partial derivatives of mass quantities //
! //  with respect to rho at constant t                               //
! //                                                                  //
! //  Benoit Bottin, 04/12/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE sensitivities_r (p,rho,T,TNEQ,ener,enth,mm,y,findiff, &
     &epsilon,flg_anha,flg_neq,flg_stop,flg_termo,flg_log,niter,       &
     &dmdr,dpdr,dedr,dhdr,detdr)                               
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER flg_anha,flg_neq,flg_stop,flg_termo,i,ierr,niter,flg_log
      REAL(kind=8) p,rho,T,TNEQ(4),ener(8),enth(8),mm,y(:),epsilon,findiff
      REAL(kind=8) newp,newrho,newe(8),newh(8),newm
      REAL(kind=8) dpdr,dhdr(8),dedr(8),dmdr,detdr(8)
      REAL(kind=8),ALLOCATABLE :: newy(:)
      CHARACTER*6 pgname
      CHARACTER*70 localinfo
      INTERFACE
          REAL(kind=8) FUNCTION mixture_pressure (rho,T,y,flg_stop)
              REAL(kind=8) rho,T,y(:)
              INTEGER flg_stop
          END FUNCTION mixture_pressure
          REAL(kind=8) FUNCTION mixture_molarmass (x,mode,flg_stop)
              REAL(kind=8) x(:)
              INTEGER flg_stop,mode
          END FUNCTION mixture_molarmass
          SUBROUTINE mixture_energy (T,TNEQ,x,&
             &flg_anha,mode,flg_neq,flg_stop,flg_termo,ener)
              INTEGER flg_anha,mode,flg_neq,flg_stop,flg_termo
              REAL(kind=8) T,TNEQ(4),ener(8),x(:)
          END SUBROUTINE mixture_energy
          SUBROUTINE mixture_enthalpy (T,TNEQ,x,&
             &flg_anha,mode,flg_neq,flg_stop,flg_termo,enth)
              INTEGER flg_anha,mode,flg_neq,flg_stop,flg_termo
              REAL(kind=8) T,TNEQ(4),enth(8),x(:)
          END SUBROUTINE mixture_enthalpy
          SUBROUTINE massfrac (y,rho,T,TNEQ,eps,flg_log,&
             &flg_anha,flg_neq,flg_stop,flg_termo,flg_guess,y_guess,niter)
              INTEGER flg_log,flg_stop,flg_anha,flg_neq,flg_termo,flg_guess,niter
              REAL(kind=8) y(:),rho,T,TNEQ(4),y_guess(:),eps
          END SUBROUTINE massfrac
      END INTERFACE
!
!     Flg_neq is not used here, but to keep the same structure it is passed
!     as well. The next statement removes the compile-time warning
      i=flg_neq
      pgname='pgaslib'
      CALL FILL_WITH_BLANKS (localinfo)
      ALLOCATE (newy(1:nsp),STAT=ierr)
      IF (ierr.NE.0) THEN
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE thermodynamic library error'
      END IF
!
!     New composition at rho+drho
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      newrho=rho*(1.0d0+findiff)
      CALL MASSFRAC (newy,newrho,t,TNEQ,epsilon,flg_log,&
     &flg_anha,flg_neq,flg_stop,flg_termo,1,y,niter)
!
!     New properties at rho+drho
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^
      CALL MIXTURE_ENERGY (t,tneq,newy,&
     &flg_anha,2,0,flg_stop,flg_termo,newe)
      CALL MIXTURE_ENTHALPY (t,tneq,newy,&
     &flg_anha,2,0,flg_stop,flg_termo,newh)
      newp=MIXTURE_PRESSURE (newrho,t,newy,flg_stop)
      newm=MIXTURE_MOLARMASS (newy,2,flg_stop)
!
!     Finite difference values
!     ^^^^^^^^^^^^^^^^^^^^^^^^
      DO i=1,8
          dedr(i)=(newe(i)-ener(i))/(findiff*rho)
          dhdr(i)=(newh(i)-enth(i))/(findiff*rho)
          detdr(i)=(newrho*newe(i)-rho*ener(i))/(findiff*rho)
      END DO
      dpdr=(newp-p)/(findiff*rho)
      dmdr=(newm-mm)/(findiff*rho)
!
      DEALLOCATE (newy)
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE sensitivities_t                                      //
! //                                                                  //
! //  Input:  mm          molar mass of the mixture (kg/mol)          //
! //          p           pressure of the mixture (Pa)                //
! //          rho         density of the mixture (kg/m3)              //
! //          T           equilibrium temperature (K)                 //
! //          TNEQ        non-equilibrium temperature (K) (not used)  //
! //          y           array of mass fractons                      //
! //          ener(8)     array of energies per mass                  //
! //          enth(8)     array of enthalpies per mass                //
! //          epsilon     relative precision of the finite difference //
! //                                                                  //
! //  Flags:  flg_anha    anharmonicity corrections flag              //
! //          flg_neq     nonequilibrium temperatures flag            //
! //          flg_stop    stops if errors in subroutines              //
! //          flg_termo   statistical or reftable computation         //
! //                                                                  //
! //  Output: dedt (8)    de over dt at constant rho (energy)         //
! //          dhdt (8)    dh over dt at constant rho (enthalpy)       //
! //          detdt(8)    det over dt at constant rho (vol. energy)   //
! //          dpdt        dp over dt at constant rho (pressure)       //
! //          dmdt        dm over dt at constant rho (molar mass)     //
! //                                                                  //
! //  This subroutine provides partial derivatives of mass quantities //
! //  with respect to t at constant rho.                              //
! //                                                                  //
! //  Benoit Bottin, 04/12/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE sensitivities_t (p,rho,T,TNEQ,ener,enth,mm,y,findiff, &
     &epsilon,flg_anha,flg_neq,flg_stop,flg_termo,flg_log,niter,       &
     &dmdt,dpdt,dedt,dhdt,detdt)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER flg_anha,flg_neq,flg_stop,flg_termo,i,ierr,niter,flg_log
      REAL(kind=8) p,rho,T,TNEQ(4),ener(8),enth(8),mm,y(:),epsilon,findiff
      REAL(kind=8) newp,newt,newe(8),newh(8),newm,newtneq(4)
      REAL(kind=8) dpdt,dhdt(8),dedt(8),dmdt,detdt(8)
      REAL(kind=8),ALLOCATABLE :: newy(:)
      CHARACTER*6 pgname
      CHARACTER*70 localinfo
      INTERFACE
          REAL(kind=8) FUNCTION mixture_pressure (rho,T,y,flg_stop)
              REAL(kind=8) rho,T,y(:)
              INTEGER flg_stop
          END FUNCTION mixture_pressure
          REAL(kind=8) FUNCTION mixture_molarmass (x,mode,flg_stop)
              REAL(kind=8) x(:)
              INTEGER flg_stop,mode
          END FUNCTION mixture_molarmass
          SUBROUTINE mixture_energy (T,TNEQ,x,&
             &flg_anha,mode,flg_neq,flg_stop,flg_termo,ener)
              INTEGER flg_anha,mode,flg_neq,flg_stop,flg_termo
              REAL(kind=8) T,TNEQ(4),ener(8),x(:)
          END SUBROUTINE mixture_energy
          SUBROUTINE mixture_enthalpy (T,TNEQ,x,&
             &flg_anha,mode,flg_neq,flg_stop,flg_termo,enth)
              INTEGER flg_anha,mode,flg_neq,flg_stop,flg_termo
              REAL(kind=8) T,TNEQ(4),enth(8),x(:)
          END SUBROUTINE mixture_enthalpy
          SUBROUTINE massfrac (y,rho,T,TNEQ,eps,flg_log,&
             &flg_anha,flg_neq,flg_stop,flg_termo,flg_guess,y_guess,niter)
              INTEGER flg_log,flg_stop,flg_anha,flg_neq,flg_termo,flg_guess,niter
              REAL(kind=8) y(:),rho,T,TNEQ(4),y_guess(:),eps
          END SUBROUTINE massfrac
      END INTERFACE
!
!     New composition at t+dt
!     ^^^^^^^^^^^^^^^^^^^^^^^
      pgname='pgaslib'
      CALL FILL_WITH_BLANKS (localinfo)
      ALLOCATE (newy(1:nsp),STAT=ierr)
      IF (ierr.NE.0) THEN
          CALL PRINT_ERROR (43,pgname,localinfo,1)
          STOP 'Program ended due to PEGASE thermodynamic library error'
      END IF
      IF (flg_neq.EQ.0) THEN
          newt=t*(1.0d0+findiff)
      ELSE
          newtneq=TNEQ*(1.0d0+findiff)
      ENDIF
      CALL MASSFRAC (newy,rho,newt,newtneq,epsilon,flg_log,&
     &flg_anha,flg_neq,flg_stop,flg_termo,1,y,niter)
!
!     New properties at t+dt
!     ^^^^^^^^^^^^^^^^^^^^^^
      CALL MIXTURE_ENERGY (newt,newtneq,newy,&
     &flg_anha,2,0,flg_stop,flg_termo,newe)
      CALL MIXTURE_ENTHALPY (newt,newtneq,newy,&
     &flg_anha,2,0,flg_stop,flg_termo,newh)
      newp=MIXTURE_PRESSURE (rho,newt,newy,flg_stop)
      newm=MIXTURE_MOLARMASS (newy,2,flg_stop)
!
!     Finite difference values
!     ^^^^^^^^^^^^^^^^^^^^^^^^
      DO i=1,8
          dedt(i)=(newe(i)-ener(i))/(findiff*t)
          dhdt(i)=(newh(i)-enth(i))/(findiff*t)
          detdt(i)=rho*dedt(i)
      END DO
      dpdt=(newp-p)/(findiff*t)
      dmdt=(newm-mm)/(findiff*t)
!
      DEALLOCATE (newy)
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
