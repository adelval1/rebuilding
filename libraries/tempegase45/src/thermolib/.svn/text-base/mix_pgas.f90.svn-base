! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTION mixture_pressure                                       //
! //                                                                  //
! //  Input:  rho         density of the mixture (kg/m3)              //
! //          T           temperature of the mixture (K)              //
! //          y           mass fractions                              //
! //          flg_stop    if 1, stops on library errors               //
! //                                                                  //
! //  This function calculates the pressure of a mixture given        //
! //  density, temperature and mass fractions.                        //
! //                                                                  //
! //  Benoit Bottin, 11/10/96. Modified and cleaned F90, 28/7/97.     //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION mixture_pressure (rho,T,y,flg_stop)
!
      USE global_thermo
      IMPLICIT NONE
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      REAL(kind=8) rho,T,y(:),out
      INTEGER i,istop,flg_stop
!      
!     Check input values and initialize variables
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      istop=0
      pgname='pegaslib'
      out=0.0d0
      IF (T.LE.0.0d0) THEN
          WRITE (localinfo,101) T
101       FORMAT (11x,'Routine: mixture_pressure. Value:',d15.8)
          CALL PRINT_ERROR (24,pgname,localinfo,0)
          istop=1
      ENDIF
      IF (rho.LE.0.0d0) THEN
          WRITE (localinfo,101) rho
          CALL PRINT_ERROR (26,pgname,localinfo,0)
          istop=1
      ENDIF
!
!     Compute pressure
!     ^^^^^^^^^^^^^^^^
      DO i=1,nsp
          out=out+(y(i)*rho)*RSP(i)*T
      END DO
      mixture_pressure=out
!
!     Error handling
!     ^^^^^^^^^^^^^^
      IF ((istop.EQ.1).AND.(flg_stop.EQ.1)) THEN
          STOP 'Program ended due to PEGASE thermodynamic library error'
      ENDIF
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTION mixture_density                                        //
! //                                                                  //
! //  Input:  p           pressure of the mixture (Pa)                //
! //          T           temperature of the mixture (K)              //
! //          x           mole fractions                              //
! //          flg_stop    if 1, stops on library errors               //
! //                                                                  //
! //  This function calculates the density of a mixture given         //
! //  pressure,temperature and mole fractions.                        //
! //                                                                  //
! //  Benoit Bottin, 11/10/96. Modified and cleaned F90, 28/7/97.     //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION mixture_density (p,T,x,flg_stop)
!
      USE global_thermo
      IMPLICIT NONE
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      REAL(kind=8) p,T,x(:),out
      INTEGER i,istop,flg_stop
!      
!     Check input values and initialize variables
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      istop=0
      pgname='pegaslib'
      out=0.0d0
      IF (T.LE.0.0d0) THEN
          WRITE (localinfo,101) T
101       FORMAT (11x,'Routine: mixture_density:',d15.8)
          CALL PRINT_ERROR (24,pgname,localinfo,0)
          istop=1
      ENDIF
      IF (p.LE.0.0d0) THEN
          WRITE (localinfo,101) p
          CALL PRINT_ERROR (25,pgname,localinfo,0)
          istop=1
      ENDIF
!
!     Compute density
!     ^^^^^^^^^^^^^^^^
      DO i=1,nsp
          out=out+(x(i)*p)/(RSP(i)*T)
      END DO
      mixture_density=out
!
!     Error handling
!     ^^^^^^^^^^^^^^
      IF ((istop.EQ.1).AND.(flg_stop.EQ.1)) THEN
          STOP 'Program ended due to PEGASE thermodynamic library error'
      ENDIF
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTION mixture_molarmass                                      //
! //                                                                  //
! //  Input:  x           mole or mass fraction                       //
! //          mode        if 1, mole fractions, if 2, mass fractions  //
! //          flg_stop    if 1, stops on library errors               //
! //                                                                  //
! //  This function calculates the molar mass of the mixture given    //
! //  the composition.                                                //
! //                                                                  //
! //  Benoit Bottin, 11/10/96. Modified and cleaned F90, 28/7/97.     //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION mixture_molarmass (x,mode,flg_stop)
!
      USE global_thermo
      IMPLICIT NONE
      INTERFACE
          REAL(kind=8) FUNCTION MIXTURE_RSPECIFIC (x,mode,flg_stop)
              INTEGER mode,flg_stop
              REAL(kind=8) x(:)
          END FUNCTION MIXTURE_RSPECIFIC
      END INTERFACE
!
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      REAL(kind=8) x(:),out,r
      INTEGER i,istop,flg_stop,mode
!      
!     Initialize variables
!     ^^^^^^^^^^^^^^^^^^^^
      pgname='pegaslib'
      out=0.0d0
      istop=0
!
!     Compute molar mass, depending on mode:
!         for a value of 1, it is mole fractions -> direct computation
!         for a value of 2, it is mass fractions -> go through spec. R
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (mode.EQ.1) THEN
          DO i=1,nsp
              out=out+x(i)*MMOL(i)
          END DO
      ELSE
          r=mixture_rspecific (x,mode,flg_stop)
          IF (r.LE.1.0d-15) THEN
              WRITE (localinfo,101) r
101           FORMAT (11x,'Routine: mixture_molarmass. Value:',d15.8)
              CALL PRINT_ERROR (32,pgname,localinfo,0)
              istop=1
          ENDIF
          out=RUNIV/r
      ENDIF          
      mixture_molarmass=out
!
!     Error handling
!     ^^^^^^^^^^^^^^
      IF ((istop.EQ.1).AND.(flg_stop.EQ.1)) THEN
          STOP 'Program ended due to PEGASE thermodynamic library error'
      ENDIF
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTION mixture_rspecific                                      //
! //                                                                  //
! //  Input:  x           mole or mass fraction                       //
! //          mode        if 1, mole fractions, if 2, mass fractions  //
! //          flg_stop    if 1, stops on library errors               //
! //                                                                  //
! //  This function calculates the specific gas constant of a mixture //
! //  given the composition.                                          //
! //                                                                  //
! //  Benoit Bottin, 11/10/96. Modified and cleaned F90, 28/7/97.     //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION mixture_rspecific (x,mode,flg_stop)
!
      USE global_thermo
      IMPLICIT NONE
      INTERFACE
          REAL(kind=8) FUNCTION MIXTURE_MOLARMASS (x,mode,flg_stop)
              INTEGER mode,flg_stop
              REAL(kind=8) x(:)
          END FUNCTION MIXTURE_MOLARMASS
      END INTERFACE
!
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      REAL(kind=8) x(:),out,m
      INTEGER i,istop,flg_stop,mode
!      
!     Initialize variables
!     ^^^^^^^^^^^^^^^^^^^^
      pgname='pegaslib'
      out=0.0d0
      istop=0
!
!     Compute specific gas constant, depending on mode:
!         for a value of 1, it is mole fractions -> go through molarmass
!         for a value of 2, it is mass fractions -> direct computation
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (mode.EQ.1) THEN
          m=MIXTURE_MOLARMASS (x,mode,flg_stop)
          IF (m.LE.1.0d-15) THEN
              WRITE (localinfo,101) m
101           FORMAT (11x,'Routine: mixture_rspecific. Value:',d15.8)
              CALL PRINT_ERROR (33,pgname,localinfo,0)
              istop=1
          ENDIF
          out=RUNIV/m
      ELSE
          DO i=1,nsp
              out=out+x(i)*RSP(i)
          END DO
      ENDIF          
      mixture_rspecific=out
!
!     Error handling
!     ^^^^^^^^^^^^^^
      IF ((istop.EQ.1).AND.(flg_stop.EQ.1)) THEN
          STOP 'Program ended due to PEGASE thermodynamic library error'
      ENDIF
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTION mixture_moleratio                                      //
! //                                                                  //
! //  Input:  x           mole or mass fraction                       //
! //          mode        if 1, mole fractions, if 2, mass fractions  //
! //          flg_stop    if 1, stops on library errors               //
! //                                                                  //
! //  This function calculates the ratio of hot/cold number of moles  //
! //  of gaz in the mixture, given composition                        //
! //                                                                  //
! //  Benoit Bottin, 11/10/96. Modified and cleaned F90, 28/7/97.     //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION mixture_moleratio (x,mode,flg_stop)
!
      USE global_thermo
      IMPLICIT NONE
      INTERFACE
           REAL(kind=8) FUNCTION MIXTURE_MOLARMASS (x,mode,flg_stop)
              INTEGER mode,flg_stop
              REAL(kind=8) x(:)
          END FUNCTION MIXTURE_MOLARMASS
         REAL(kind=8) FUNCTION MIXTURE_RSPECIFIC (x,mode,flg_stop)
              INTEGER mode,flg_stop
              REAL(kind=8) x(:)
          END FUNCTION MIXTURE_RSPECIFIC
      END INTERFACE
!
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      REAL(kind=8) x(:),out,r,ro,m,mo
      INTEGER istop,flg_stop,mode
!      
!     Initialize variables
!     ^^^^^^^^^^^^^^^^^^^^
      pgname='pegaslib'
      out=0.0d0
      istop=0
!
!     Compute specific gas constant, depending on mode:
!         for a value of 1, it is contained in the nsp+1 item of x
!         for a value of 2, it is mass fractions, use rspecific
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (mode.EQ.1) THEN
          m=MIXTURE_MOLARMASS (x,mode,flg_stop)
          IF (m.LE.1.0d-15) THEN
              WRITE (localinfo,101) m
101           FORMAT (11x,'Routine: mixture_moleratio. Value:',d15.8)
              CALL PRINT_ERROR (33,pgname,localinfo,0)
              istop=1
          ENDIF
          mo=MIXTURE_MOLARMASS (XINI,1,flg_stop)
          IF (mo.LE.1.0d-15) THEN
              WRITE (localinfo,102) mo
102           FORMAT (11x,'Routine: mixture_moleratio. Initial value:',d15.8)
              CALL PRINT_ERROR (33,pgname,localinfo,0)
              istop=1
          ENDIF
          out=mo/m
      ELSE
          r=MIXTURE_RSPECIFIC (x,mode,flg_stop)
          IF (r.LE.1.0d-15) THEN
              WRITE (localinfo,101) r
              CALL PRINT_ERROR (32,pgname,localinfo,0)
              istop=1
          ENDIF
          ro=MIXTURE_RSPECIFIC (XINI,1,flg_stop)
          IF (ro.LE.1.0d-15) THEN
              WRITE (localinfo,101) ro
              CALL PRINT_ERROR (34,pgname,localinfo,0)
              istop=1
          ENDIF
          out=r/ro
      ENDIF          
      mixture_moleratio=out
!
!     Error handling
!     ^^^^^^^^^^^^^^
      IF ((istop.EQ.1).AND.(flg_stop.EQ.1)) THEN
          STOP 'Program ended due to PEGASE thermodynamic library error'
      ENDIF
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
