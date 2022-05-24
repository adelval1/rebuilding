! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTION mol_to_mass                                            //
! //                                                                  //
! //  Input:  val         value per mole to be converted              //    
! //          isp         number of the species or zero for mixture   //
! //          mm          mixture molar mass (kg/mol)                 //
! //                                                                  //
! //  This function converts a quantity per mole to a quantity per    //
! //  unit mass. If isp is a valid species number, then mm can be     //
! //  whatever. If isp is zero, then a custom molar mass is supplied  //
! //  through mm.                                                     //
! //                                                                  //
! //  Benoit Bottin, 09/10/96. Modified and cleaned F90, 28/7/97.     //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION mol_to_mass (val,isp,mm,flg_stop)
!
      USE global_thermo
      IMPLICIT NONE
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      REAL(kind=8) val,mm,out
      INTEGER isp,istop,flg_stop
!      
      istop=0
      pgname='pegaslib'
      IF (isp.EQ.0) THEN
          IF (mm.EQ.0) THEN
              WRITE (localinfo,101) mm
101           FORMAT (11x,'Function: mol_to_mass. Value:',d15.8)
              CALL PRINT_ERROR (33,pgname,localinfo,0)
              istop=1
          ENDIF
          out=val/mm
      ELSE
          IF (MMOL(isp).EQ.0) THEN
              WRITE (localinfo,102) isp,MMOL(isp)
102   FORMAT (11x,'Function: mol_to_mass. Index:',i3,' Value:',d15.8)
              CALL PRINT_ERROR (35,pgname,localinfo,0)
              istop=1
          ENDIF
          out=val/MMOL(isp)
      ENDIF
      IF ((istop.EQ.1).AND.(flg_stop.EQ.1)) THEN      
          STOP 'Program ended due to PEGASE thermodynamic library error'
      ENDIF
!      
      mol_to_mass=out
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTION mass_to_mol                                            //
! //                                                                  //
! //  Input:  val         value per mole to be converted              //    
! //          isp         number of the species or zero for mixture   //
! //          mm          mixture molar mass (kg/mol)                 //
! //                                                                  //
! //  Flags:  flg_stop    stops on subroutine errors if set to 1      //
! //                                                                  //
! //  This function converts a quantity per mass to a quantity per    //
! //  unit mole. If isp is a valid species number, then mm can be     //
! //  whatever. If isp is zero, then a custom molar mass is supplied  //
! //  through mm. Flg_stop controls error-handling.                   //
! //                                                                  //
! //  Benoit Bottin, 09/10/96. Modified and cleaned F90, 28/7/97.     //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION mass_to_mol (val,isp,mm,flg_stop)
!
      USE global_thermo
      IMPLICIT NONE
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      REAL(kind=8) val,mm,out
      INTEGER isp,flg_stop,istop
!
      istop=0
      pgname='pegaslib'
      IF (isp.EQ.0) THEN
          IF (mm.EQ.0) THEN
              WRITE (localinfo,101) mm
101           FORMAT (11x,'Function: mass_to_mol. Value:',d15.8)
              CALL PRINT_ERROR (33,pgname,localinfo,0)
              istop=1
          ENDIF
          out=val*mm
      ELSE
          IF (MMOL(isp).EQ.0) THEN
              WRITE (localinfo,102) isp,MMOL(isp)
102   FORMAT (11x,'Function: mass_to_mol. Index:',i3,' Value:',d15.8)
              CALL PRINT_ERROR (35,pgname,localinfo,0)
              istop=1
          ENDIF
          out=val*MMOL(isp)
      ENDIF
      IF ((istop.EQ.1).AND.(flg_stop.EQ.1)) THEN      
          STOP 'Program ended due to PEGASE thermodynamic library error'
      ENDIF
!
      mass_to_mol=out
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE molefrac_to_massfrac                                 //
! //                                                                  //
! //  Input:  x           array of mole fractions to be converted     //    
! //          nsp         number of species in the mixture            //
! //          mm          mixture molar mass (kg/mol)                 //
! //                                                                  //
! //  Output: y           array of mass fractions                     //    
! //                                                                  //    
! //  This routine  converts mole fractions to mass fractions.        //
! //                                                                  //
! //  Benoit Bottin, 15/10/96. Modified and cleaned F90, 28/7/97.     //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE molefrac_to_massfrac (x,mm,y)
!      
      USE global_thermo
      IMPLICIT NONE
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      REAL(kind=8) x(:),y(:),mm
      INTEGER i
!      
      pgname='pegaslib'
      IF (mm.LE.0) THEN
          WRITE (localinfo,101) mm
101       FORMAT (11x,'Routine: molefrac_to_massfrac. Value:',d15.8)
          CALL PRINT_ERROR (33,pgname,localinfo,0)
          STOP 'Program ended due to PEGASE thermodynamic library error'
      ENDIF
!      
      DO i=1,nsp
          y(i)=x(i)*MMOL(i)/mm
      END DO
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE massfrac_to_molefrac                                 //
! //                                                                  //
! //  Input:  y           array of mass fractions to be converted     //    
! //          nsp         number of species in the mixture            //
! //          mm          mixture molar mass (kg/mol)                 //
! //                                                                  //
! //  Output: x           array of mole fractions                     //    
! //                                                                  //    
! //  This routine  converts mass fractions to mole fractions.        //
! //                                                                  //
! //  Benoit Bottin, 15/10/96. Modified and cleaned F90, 28/7/97.     //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE massfrac_to_molefrac (y,mm,x)
!      
      USE global_thermo
      IMPLICIT NONE
      REAL(kind=8) x(:),y(:),mm
      INTEGER i
!      
      DO i=1,nsp
          x(i)=y(i)*mm/MMOL(i)
      END DO
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTION number_density                                         //
! //                                                                  //
! //  Input:  p           pressure (Pa)                               //    
! //          T           temperature (K)                             //
! //                                                                  //
! //  Output:             number density in m-3                       //    
! //                                                                  //    
! //  This function computes the number density of a mixture (provide //
! //  the pressure) or of a species (provide the partial pressure).   //
! //                                                                  //
! //  Benoit Bottin, 15/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION number_density (p,T)
!      
      USE global_thermo
      IMPLICIT NONE
      REAL(kind=8) p,T,nd
!      
      nd=p/KUNIV/T
      number_density=nd
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTION concentration                                          //
! //                                                                  //
! //  Input:  p           pressure (Pa)                               //    
! //          T           temperature (K)                             //
! //                                                                  //
! //  Output:             concentration in mol/m3                     //    
! //                                                                  //    
! //  This function computes the concentration  of a mixture (provide //
! //  the pressure) or of a species (provide the partial pressure).   //
! //                                                                  //
! //  Benoit Bottin, 15/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION concentration (p,T)
!      
      USE global_thermo
      IMPLICIT NONE
      REAL(kind=8) p,T,c
!      
      c=p/RUNIV/T
      concentration=c
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTION species_density                                        //
! //                                                                  //
! //  Input:  isp         number of the species                       //    
! //          p           pressure of the species                     //
! //          T           temperature of the species                  //
! //                                                                  //
! //  This function computes the density from the equation of state.  //
! //                                                                  //
! //  Benoit Bottin, 16/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION species_density (isp,p,T)
!
      USE global_thermo
      IMPLICIT NONE
      REAL(kind=8) p,T,out
      INTEGER isp
!      
      out=p/RSP(isp)/T      
      species_density=out
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTION species_pressure                                       //
! //                                                                  //
! //  Input:  isp         number of the species                       //    
! //          rho         density  of the species                     //
! //          T           temperature of the species                  //
! //                                                                  //
! //  This function computes the pressure from the equation of state. //
! //                                                                  //
! //  Benoit Bottin, 16/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION species_pressure (isp,rho,T)
!
      USE global_thermo
      IMPLICIT NONE
      REAL(kind=8) rho,T,out
      INTEGER isp
!      
      out=rho*RSP(isp)*T      
      species_pressure=out
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
