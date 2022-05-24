! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  OUTPUT   L I B R A R Y             //
! //                                                                  //
! //  SUBROUTINE print_output                                         //
! //                                                                  //
! //  Input:  in_dir      working directory for the input file        //    
! //                                                                  //
! //  Outs:   outname     full name of the output file                // 
! //                                                                  //
! //  The routine prompts for the input file to use for the           //
! //  calculation. It concatenates the input name and the input       //
! //  directory and returns the value. If the file does not exist,    //
! //  an error is generated.                                          //
! //                                                                  //
! //                                                                  //
! //  Benoit Bottin, 22/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE print_output(cntr,cntr2,cntrmax,casename,mixname,&
     &flg_anha,flg_neq,flg_mode,flg_oper,flg_termo,flg_traco,flg_usr)
!
      USE global_pegase
      USE global_thermo
      IMPLICIT NONE
      CHARACTER*80 casename,mixname
      CHARACTER*69 titlename
      CHARACTER*20 tag
      CHARACTER*10 tag2(10),tagx,tagy
      INTEGER cntr,cntr2,cntrmax,i,j,flg_usr
      INTEGER flg_anha,flg_neq,flg_mode,flg_oper,flg_termo,flg_traco
      INTEGER lentrim
      REAL(kind=8) summ(5)
!     
!     Open output file
!     ^^^^^^^^^^^^^^^^
      CALL FILL_WITH_BLANKS (titlename)
!
!     If first calculation of the scan, print flags 
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (cntr2.EQ.1) THEN
          CALL CENTERSTRING (casename(1:69),titlename,i)
          WRITE (21,*)
          WRITE (21,2002) cntr,titlename
          IF (flg_oper.EQ.0) THEN
              tag='single species      '
          ELSEIF (flg_oper.EQ.1) THEN
              tag='equilibrium flow    '
          ELSEIF (flg_oper.EQ.2) THEN
              tag='frozen flow         '
          ENDIF
          WRITE (21,2004) tag
          IF (flg_mode.EQ.1) THEN
              tag='pressure-temperature'
          ELSEIF (flg_mode.EQ.2) THEN
              tag='density-temperature '
          ELSEIF (flg_mode.EQ.3) THEN
              tag='p-t logarithmic     '
          ELSEIF (flg_mode.EQ.4) THEN
              tag='rho-t logarithmic   '
          ENDIF
          WRITE (21,2005) tag
          IF (flg_anha.EQ.-1) THEN
              tag='RRHO, corrected temp'
          ELSEIF (flg_anha.EQ.0) THEN
              tag='not considered      '
          ELSEIF (flg_anha.EQ.1) THEN
              tag='simplified (RVAC)   '
          ELSE
              tag='advanced (RVEC)     '
          ENDIF
          WRITE (21,2006) tag
          IF (flg_neq.EQ.0) THEN
              tag='not considered      '
          ELSE
              tag='4 temperatures      '
          ENDIF
          WRITE (21,2007) tag
          IF (flg_termo.EQ.0) THEN
              tag='statistical physics '
          ELSEIF (flg_termo.EQ.1) THEN
              tag='reference tables    '
          ELSEIF (flg_termo.EQ.2) THEN
              tag='Gupta... curve fits '
          ELSE
              tag='Srinivasan... fits  '
          ENDIF
          WRITE (21,2008) tag
          IF (flg_traco.EQ.0) THEN
              tag='collision integrals '
          ELSE
              tag='Gupta... curve fits '
          ENDIF
          WRITE (21,2009) tag
          IF (flg_usr.NE.0) THEN
              WRITE (tag,100) 'run with flag =',flg_usr
100           FORMAT (a15,i5)
              WRITE (21,2010) tag
          ENDIF
!
!     Print mixture summary
!     ^^^^^^^^^^^^^^^^^^^^^
          WRITE (21,2070)
          j=MAX(LENTRIM(mixname),45)
          i=MAX(LENTRIM(mixname)-44,1)
          WRITE (21,2071) mixname(i:j)
          WRITE (21,2072) nsp
          WRITE (21,2073) nr
          WRITE (21,2074) nc
!
      ENDIF
!
!     Print calculation reference number 
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      WRITE (21,2011) cntr2,cntrmax
!
!     Printing mixture thermodynamic values
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      WRITE (21,2015) array3(1)
      WRITE (21,2016) array3(2)
      WRITE (21,2017) array3(17)
      IF (flg_neq.NE.0) THEN
          WRITE (21,2023) array3(18)
          WRITE (21,2024) array3(19)
          WRITE (21,2025) array3(20)
      ENDIF      
      WRITE (21,2018) array3(8)
      WRITE (21,2019) array3(6)
      WRITE (21,2020) array3(7)
      WRITE (21,2021) array3(4)
      WRITE (21,2022) array3(5)
!
!     Printing global energy parameter values
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      WRITE (21,2050)
      tag='Internal energy     '
      tagx=' J/mol    '
      tagy=' J/kg     '
      WRITE (21,2051) tag,array4(1,1),tagx,array4(10,1),tagy
      tag='Helmholtz energy'
      WRITE (21,2051) tag,array4(2,1),tagx,array4(11,1),tagy
      tag='Gibbs free energy   '
      WRITE (21,2051) tag,array4(3,1),tagx,array4(12,1),tagy
      tag='Enthalpy            '
      WRITE (21,2051) tag,array4(4,1),tagx,array4(13,1),tagy
      tagx=' J/(mol K)'
      tagy=' J/(kg K) '
      tag='Entropy             '
      WRITE (21,2051) tag,array4(5,1),tagx,array4(14,1),tagy
      tag='Frozen cv           '
      WRITE (21,2051) tag,array4(6,1),tagx,array4(15,1),tagy
      tag='Equilibrium cv      '
      WRITE (21,2051) tag,array4(7,1),tagx,array4(16,1),tagy
      tag='Frozen cp           '
      WRITE (21,2051) tag,array4(8,1),tagx,array4(17,1),tagy
      tag='Equilibrium cp      '
      WRITE (21,2051) tag,array4(9,1),tagx,array4(18,1),tagy
      tagx='          '
      tagy='          '
      tag='Frozen gamma        '
      WRITE (21,2051) tag,array3(9),tagx,array3(11),tagy
      tag='Equilibrium gamma   '
      WRITE (21,2051) tag,array3(10),tagx,array3(12),tagy
      tag='Isentropic exponent '
      WRITE (21,2051) tag,array3(23),tagx,array3(24),tagy
      WRITE (21,2052) array3(13)
      WRITE (21,2053) array3(14)
!
!     Printing transport coefficients
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      WRITE (21,*)
      WRITE (21,2101) array6(1,1)
      WRITE (21,2102) array6(2,2)+array6(2,3)+array6(2,4)+array6(2,5)
      WRITE (21,2103) array6(2,6)
      WRITE (21,2104) array6(2,7)
      WRITE (21,2105) array6(2,1)
      WRITE (21,2106) array6(4,1)
      WRITE (21,2107) array6(5,1)
      WRITE (21,2108) array6(3,1)
      WRITE (21,2109) array6(3,2)
!
!     Printing energy per mode values
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      tag2(1)= 'energy    '
      tag2(2)= 'Helmholtz '
      tag2(3)= 'Gibbs     '
      tag2(4)= 'enthalpy  '
      tag2(5)= 'entropy   '
      tag2(6)= 'frozen cv '
      tag2(7)= 'equil. cv '
      tag2(8)= 'frozen cp '
      tag2(9)= 'equil. cp '
      tag2(10)='k (W/(m K)'
      IF (flg_anha.LE.0) THEN      
          WRITE (21,2040)
      ELSEIF (flg_anha.EQ.1) THEN      
          WRITE (21,2043)
      ELSE
          WRITE (21,2044)
      END IF
      IF (flg_mode.EQ.1) THEN
          WRITE (21,2041)
          DO i=1,9
              IF (flg_anha.LE.0) THEN      
              WRITE (21,2045) tag2(i),array4(i,2),&
            & array4(i,3),array4(i,4),array4(i,5),array4(i,8)
              ELSEIF (flg_anha.EQ.1) THEN      
              WRITE (21,2046) tag2(i),array4(i,2),&
            & array4(i,6),array4(i,5),array4(i,8)
              ELSE
              WRITE (21,2047) tag2(i),array4(i,2),&
            & array4(i,7),array4(i,8)
              END IF
          END DO
      ELSE
          WRITE (21,2042)
          DO i=10,18
              IF (flg_anha.LE.0) THEN      
              WRITE (21,2045) tag2(i-9),array4(i,2),&
            & array4(i,3),array4(i,4),array4(i,5),array4(i,8)
              ELSEIF (flg_anha.EQ.1) THEN      
              WRITE (21,2046) tag2(i-9),array4(i,2),&
            & array4(i,6),array4(i,5),array4(i,8)
              ELSE
              WRITE (21,2047) tag2(i-9),array4(i,2),&
            & array4(i,7),array4(i,8)
              END IF
          END DO
      ENDIF
      IF (flg_anha.LE.0) THEN      
          WRITE (21,2045) tag2(10),array6(2,2),array6(2,3), &
        & array6(2,4),array6(2,5),array6(2,7)
      ELSEIF (flg_anha.EQ.1) THEN      
          WRITE (21,2046) tag2(10),array6(2,2),array6(2,3)+ &
        & array6(2,4),array6(2,5),array6(2,7)
      ELSE
          WRITE (21,2047) tag2(10),array6(2,2),array6(2,3)+ &
        & array6(2,4)+array6(2,5),array6(2,7)
      END IF
!
!     Printing species composition and energies
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      DO i=1,4
          summ(i)=0.0d0
      END DO
      WRITE (21,2030)
      WRITE (21,2031)
      DO i=1,nsp
          WRITE (21,2032) SPC_LABEL(i),array1(1,i),array1(2,i),   &
 &        array1(5,i),array6(6,i)
          summ(1)=summ(1)+array1(1,i)
          summ(2)=summ(2)+array1(2,i)
          summ(3)=summ(3)+array1(5,i)
      END DO
      tag='------> full mixture'
      WRITE (21,2032) tag,(summ(i),i=1,3)
!
      DO i=1,4
          summ(i)=0.0d0
      END DO
      WRITE (21,2035)
      WRITE (21,2036)
      DO i=1,nsp
          WRITE (21,2037) SPC_LABEL(i),array2(12,i,1)*array1(2,i),  &
        & array2(13,i,1)*array1(2,i),array6(7,i),array6(8,i)
          summ(1)=summ(1)+array2(12,i,1)*array1(2,i)
          summ(2)=summ(2)+array2(13,i,1)*array1(2,i)
      END DO
      tag='------> full mixture'
      WRITE (21,2032) tag,(summ(i),i=1,2)
!
!     Printing equilibrium constants
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      WRITE (21,2060)
      DO 105 i=1,nr
          WRITE (21,2061) i,(array5(i,j),j=1,4)
105   CONTINUE
!     ------------------------------------------------------------------
2001  FORMAT (' /////////////////////////////////////////////////////',&
     &'////////////////////////',/,' //                              ',&
     &'                                           //',/,               &
     &' //            P E G A S E  ---  formatted text output file .O',&
     &'UT            //',/,' //                        ',              &
     &'                                                 //')
2003  FORMAT (' //                        ',                           &
     &'                                                 //',/,         &
     &' /////////////////',                                            &
     &'////////////////////////////////////////////////////////////',/)
2002  FORMAT (' ======================================================',&
     &'=======================',/,' ==>              Pegase now proces',&
     &'ses input case number ',i4,':             <==',/,' ==> ',a69,    &
     &' <==',/,' ==================================================',   &
     &'===========================',/,/,' Flags setting:')
2004  FORMAT (7x,'Required operation ............... ',a20)
2005  FORMAT (7x,'Input type ....................... ',a20)
2006  FORMAT (7x,'Anharmonicity corrections ........ ',a20)
2007  FORMAT (7x,'Thermal nonequilibrium ........... ',a20)
2008  FORMAT (7x,'Thermodynamic properties data .... ',a20)
2009  FORMAT (7x,'Transport properties data ........ ',a20)        
2010  FORMAT (7x,'Access to user-defined routines .. ',a20)        
2011  FORMAT (/,' ------------------------------------------',/,       &
     &' Results from step ',i4,' of ',i4,' in the scan:',/,            &
     &' ------------------------------------------',/)
2015  FORMAT (' Mixture properties: pressure:                ',f14.6,  &
     &' Pa')
2016  FORMAT ('                     density:                   ',e12.6,&
     &' kg/m3')
2017  FORMAT ('                     translation temperature: ',f14.6,  &
     &' K')
2018  FORMAT ('                     total number density:      ',e12.6,&
     &' /m3')
2019  FORMAT ('                     mole ratio:              ',f14.6,  &
     &' mol/mol cold gas')
2020  FORMAT ('                     compressibility factor:  ',f14.6)
2021  FORMAT ('                     molar mass:              ',f14.6,  &
     &' kg/mol')
2022  FORMAT ('                     specific gas constant:   ',f14.6,  &
     &' J/(kg K)')
2023  FORMAT ('                     rotation temperature:    ',f14.6,  &
     &' K')
2024  FORMAT ('                     vibration temperature:   ',f14.6,  &
     &' K')
2025  FORMAT ('                     electron temperature:    ',f14.6,  &
     &' K')
2030  FORMAT (/,' Species                         Mole         Mass      ',&
     &' Number       m.f.p.')

2031  FORMAT ('                             fraction     fraction     ',  &
     &' Density          (m)')
2032  FORMAT (' ',a20,'     ',e11.5,'  ',e11.5,'  ',e11.5,'  ',e11.5)

2035  FORMAT (/,' Species                     Enthalpy      Entropy    Di',&
     &'ffusion        Lewis')

2036  FORMAT ('                               (J/kg)   (J/(kg K))     ',  &
     &'(m**2/s)       number')
2037  FORMAT (' ',a20,'     ',e11.5,'  ',e11.5,'  ',e11.5,'  ',e11.5)

2040  FORMAT (/,' Quantity        transl.     rotation    vibration    ', &
     &' electron     reactive')
2041  FORMAT (' (per mole)')
2042  FORMAT (' (per mass)')
2043  FORMAT (/,' Quantity        transl.      coupled  rot. + vib.   ', &
     &'  electron     reactive')
2044  FORMAT (/,' Quantity        transl.                  internal    ', &
     &'              reactive')
2045  FORMAT (' ',a10,' ',e12.5,' ',e12.5,' ',e12.5,' ',e12.5,' ',e12.5)
2046  FORMAT (' ',a10,' ',e12.5,' ',6x,' ',e12.5,' ',5x,' ',e12.5,' ',e12.5)
2047  FORMAT (' ',a10,' ',e12.5,' ',12x,' ',e12.5,' ',12x,' ',e12.5)
2050  FORMAT (/,'                            per unit mole            ',&
     &'    per unit mass')
2051  FORMAT (' ',a20,' ',f19.6,a10,f19.6,a10)
2052  FORMAT (/,' Frozen speed of sound      ',f16.6,' m/s')
2053  FORMAT (  ' Equilibrium speed of sound ',f16.6,' m/s')
2060  FORMAT (/,' Reaction rates                 ln Kp        ln Kc        ',&
     &'ln Kx        ln Ky') 
2061  FORMAT ('     Equation ',i4,'      ',4(1x,f12.6))
2070  FORMAT (/,' Mixture setting:')
2071  FORMAT (  '         Source file: .......... ',a45)
2072  FORMAT (  '         Nb. of species: ....... ',i4)
2073  FORMAT (  '         Nb. of reactions: ..... ',i4)
2074  FORMAT (  '         Nb. of nuclei types ... ',i4)
2101  FORMAT (  ' Viscosity:                     ',e12.6,' Pl')
2102  FORMAT (  ' Thermal conductivity, frozen   ',e12.6,' W/(K m)')
2103  FORMAT (  '                       electron ',e12.6,' W/(K m)')
2104  FORMAT (  '                       reactive ',e12.6,' W/(K m)')
2105  FORMAT (  '                       total    ',e12.6,' W/(K m)')
2106  FORMAT (  ' Electrical conductivity        ',e12.6,' S/m')
2107  FORMAT (  ' Debye length                   ',e12.6,' m')
2108  FORMAT (  ' Prandtl number, frozen         ',f12.6)
2109  FORMAT (  '                 equilibrium    ',f12.6)
!     ------------------------------------------------------------------
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  OUTPUT   L I B R A R Y             //
! //                                                                  //
! //  SUBROUTINE print_data                                           //
! //                                                                  //
! //  Input:  in_dir      working directory for the input file        //    
! //                                                                  //
! //  Outs:   outname     full name of the output file                // 
! //                                                                  //
! //  The routine prompts for the input file to use for the           //
! //  calculation. It concatenates the input name and the input       //
! //  directory and returns the value. If the file does not exist,    //
! //  an error is generated.                                          //
! //                                                                  //
! //                                                                  //
! //  Benoit Bottin, 23/10/96 - modified TECPLOT output 13/1/99       //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE print_data(flg_mode,cntr,cntr2,casename,flg_tec,nmax)
!
      USE global_pegase
      USE global_thermo
      IMPLICIT NONE
      CHARACTER*80 casename
      CHARACTER*69 titlename
      INTEGER flg_mode,cntr,cntr2,i,flg_tec,nmax
!
!     Open output file
!     ^^^^^^^^^^^^^^^^
      CALL FILL_WITH_BLANKS (titlename)
!
!     If first calculation of the scan, print number/casename, titles
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (cntr2.EQ.1) THEN
          IF (flg_tec.NE.1) THEN
              CALL CENTERSTRING (casename(1:69),titlename,i)
              WRITE (22,*)
              WRITE (22,2002) cntr,titlename
              IF (flg_mode.EQ.1) THEN
                  WRITE (22,2010)
              ELSE
                  WRITE (22,2011)
              ENDIF
          ELSE
              CALL FILL_WITH_BLANKS(titlename)
              i = MIN0(LEN_TRIM(casename),68)
              titlename(1:i) = casename(1:i)
              titlename(i+1:i+1) = '"'
              IF (flg_mode.EQ.1) THEN
                  WRITE (22,2020)
              ELSE
                  WRITE (22,2021)
              ENDIF
              WRITE (22,2022) nmax,titlename
          ENDIF
      ENDIF
!
!     Print output values in tab format
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (flg_mode.EQ.1) THEN
          WRITE (22,2012) array3(1),array3(2),array3(3),array4(1,1), &
     &    array4(4,1),array4(5,1),array3(7),array3(4),array3(5),     &
     &    array4(9,1),array6(1,1),array6(2,1)
      ELSE
          WRITE (22,2012) array3(1),array3(2),array3(3),array4(10,1),&
     &    array4(13,1),array4(14,1),array3(7),array3(4),array3(5),   &
     &    array4(18,1),array6(1,1),array6(2,1)
      ENDIF
!     ------------------------------------------------------------------
2001  FORMAT (' /////////////////////////////////////////////////////',&
     &'////////////////////////',/,' //                              ',&
     &'                                           //',/,               &
     &' //            P E G A S E  ---  tabulated data output file .D',&
     &'AT            //',/,' //                        ',              &
     &'                                                 //')
2003  FORMAT (' //                        ',                           &
     &'                                                 //',/,         &
     &' /////////////////',                                            &
     &'////////////////////////////////////////////////////////////',/)
2002  FORMAT (' ======================================================',&
     &'=======================',/,' ==>              Pegase now proces',&
     &'ses input case number ',i4,':             <==',/,' ==> ',a69,    &
     &' <==',/,' ==================================================',   &
     &'===========================',/,/)
2010  FORMAT ('         p (Pa)   rho (kg/m³)         T (K)',            &
     &'     E (J/mol)     H (J/mol)  S [J/(molK)]         Z (-)',       &
     &'    M (kg/mol)  R [J/(kg K)] cp [J/(molK)] µ [Pl,kg/sm)]',       &
     &'   k [W/(K m)]')
2011  FORMAT ('         p (Pa)   rho (kg/m³)         T (K)',            &
     &'      e (J/kg)      h (J/kg)  s [J/(kg K)]         Z (-)',       &
     &'    M (kg/mol)  R [J/(kg K)] cp [J/(kg K)] µ [Pl,kg/sm)]',       &
     &'   k [W/(K m)]')
2012  FORMAT (1x,12e14.6)
!
!    TECPLOT FORMATS
!
2020  FORMAT (' VARIABLES = "p [Pa]" "`r [kg/m^3]" "T [K]" "E [J/mol]"',&
     &' "H [J/mol]" "S [J/(molK)]" "Z [-]" "M [kg/mol]" "R [J/(kg K)]"',&
     &' "cp [J/(molK)]" "`m [Pl]" "k [W/(K m)]" ')
2021  FORMAT (' VARIABLES = "p [Pa]" "`r [kg/m^3]" "T [K]" "e [J/kg]"',&
     &' "h [J/kg]" "s [J/(kg K)]" "Z [-]" "M [kg/mol]" "R [J/(kg K)]"',&
     &' "cp [J/(kg K)]" "`m [Pl]" "k [W/(K m)]" ')
2022  FORMAT (' ZONE I=',i6,' F=POINT T="',a69)

!     ------------------------------------------------------------------
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  OUTPUT   L I B R A R Y             //
! //                                                                  //
! //  SUBROUTINE print_data2                                          //
! //                                                                  //
! //  Input:  in_dir      working directory for the input file        //    
! //                                                                  //
! //  Outs:   outname     full name of the output file                // 
! //                                                                  //
! //  The routine prompts for the input file to use for the           //
! //  calculation. It concatenates the input name and the input       //
! //  directory and returns the value. If the file does not exist,    //
! //  an error is generated.                                          //
! //                                                                  //
! //                                                                  //
! //  Benoit Bottin, 28/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE print_data2(flg_mode,cntr,cntr2,casename,flg_tec,nmax)
!
      USE global_pegase
      USE global_thermo
      IMPLICIT NONE
      CHARACTER*14 intname(99),item
      CHARACTER*80 casename
      CHARACTER*80 titlename
      CHARACTER*200 tectitle
      INTEGER flg_mode,cntr,cntr2,i,flg_tec,nmax,l,p
!     
!     Open output file
!     ^^^^^^^^^^^^^^^^
      CALL FILL_WITH_BLANKS (titlename)
!
!     If first calculation of the scan, print number/casename
!     Also print species names and numbers
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (cntr2.EQ.1) THEN
          IF (flg_tec.NE.1) THEN
              CALL CENTERSTRING (casename(1:69),titlename,i)
              WRITE (23,*)
              WRITE (23,2002) cntr,titlename
              i=1
              WRITE (23,2003) i, SPC_LABEL(i)
              DO 301 i=2,nsp
                  WRITE (23,2004) i, SPC_LABEL(i)
301           CONTINUE
              WRITE (23,*)
!
!     Building the titles
!     ^^^^^^^^^^^^^^^^^^^
              IF (flg_mode.EQ.1) THEN
                  item='        p (Pa)'
                  WRITE (intname(1),2006) item
                  item='         T (K)'
                  WRITE (intname(2),2006) item
                  DO 201 i=3,nsp+2
                      item='     x (sp   )'
                      WRITE (intname(i),2006) item
                      WRITE (intname(i)(12:13),2007) i-2
201               CONTINUE                  
                  WRITE (23,2010) (intname(i),i=1,nsp+2)
              ELSE
                  item='   rho (kg/m3)'
                  WRITE (intname(1),2006) item
                  item='         T (K)'
                  WRITE (intname(2),2006) item
                  DO 202 i=3,nsp+2
                      item='     y (sp   )'
                      WRITE (intname(i),2006) item
                      WRITE (intname(i)(12:13),2007) i-2
202               CONTINUE                  
                  WRITE (23,2010) (intname(i),i=1,nsp+2)
              ENDIF
          ELSE
              CALL FILL_WITH_BLANKS(tectitle)
              tectitle (1:23) = '"rho [kg/m^3]" "T [K]" '
              p = 24
              DO i = 1, nsp
                  l = LEN_TRIM(SPC_SYMBOL(i))
                  IF (flg_mode.EQ.1) THEN
                      tectitle(p:p+2) = '"x('
                  ELSE
                      tectitle(p:p+2) = '"y('
                  ENDIF
                  tectitle(p+2+1  :p+2+l  ) = SPC_SYMBOL(i)(1:l)
                  tectitle(p+2+l+1:p+2+l+3) = ')" '
                  p = p+2+l+4
              END DO
              WRITE (23,2021) tectitle
              CALL FILL_WITH_BLANKS(titlename)
              i = MIN0(LEN_TRIM(casename),68)
              titlename(1:i) = casename(1:i)
              titlename(i+1:i+1) = '"'
              WRITE (23,2022) nmax,titlename
          ENDIF
      ENDIF
!
!     Print output values in tab format
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (flg_mode.EQ.1) THEN
          WRITE (intname(1),2011) array3(1)
      ELSE
          WRITE (intname(1),2011) array3(2)
      ENDIF
      WRITE (intname(2),2011) array3(3)
      DO 103 i=1,nsp
          WRITE (intname(i+2),2011) array1(flg_mode,i)
103   CONTINUE
      WRITE (23,2010) (intname(i),i=1,nsp+2)
!
!     Closing the file
!     ^^^^^^^^^^^^^^^^
!      CLOSE (23)
!     ------------------------------------------------------------------
2000  FORMAT (' /////////////////////////////////////////////////////',&
     &'////////////////////////',/,' //                              ',&
     &'                                           //',/,               &
     &' //        P E G A S E  ---  second tabulated data output file',&
     &' .DA2         //',/,' //                        ',              &
     &'                                                 //')
2001  FORMAT (' //                        ',                            &
     &'                                                 //',/,          &
     &' /////////////////',                                             &
     &'////////////////////////////////////////////////////////////',/)
2002  FORMAT (' ======================================================',&
     &'=======================',/,' ==>              Pegase now proces',&
     &'ses input case number ',i4,':             <==',/,' ==> ',a69,    &
     &' <==',/,' ==================================================',   &
     &'===========================',/)                                  
2003  FORMAT ('Species considered: ',i2,') ',a20)
2004  FORMAT ('                    ',i2,') ',a20)
2006  FORMAT (a14)
2007  FORMAT (i2)
2010  FORMAT (1x,99a14)
2011  FORMAT (e14.6)
!
!    TECPLOT FORMATS
!
2021  FORMAT (' VARIABLES = ',a200)
2022  FORMAT (' ZONE I=',i6,' F=POINT T="',a69)
!     ------------------------------------------------------------------
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  OUTPUT   L I B R A R Y             //
! //                                                                  //
! //  SUBROUTINE open_output_file                                     //
! //                                                                  //
! //  In:     un         unit of the output file                      //
! //          out_dir    output file full path                        //
! //          outname    root of the output file name                 //
! //          extension  extension of the output file                 //
! //                                                                  //
! //  The routine opens the output file as unit un for new or append. //
! //  If fop=0, then the file has not yet been opened in this run,    //
! //  and the titles are printed.                                     //
! //  If fop=-1 then the file has already been opened in this run,    //
! //  thus it must be opened with the append status.                  //
! //  In all cases fop is set to 1.                                   //
! //                                                                  //
! //  Benoit Bottin, 18/10/96. Rewritten 18/06/97.                    //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE open_output_file (un,out_dir,outname,extension,flg_tec)
!
      USE global_pegase
      IMPLICIT NONE
      CHARACTER*120 FULL_FILE_NAME,filename
      CHARACTER*80 out_dir,outname
      CHARACTER*70 localinfo
      CHARACTER*6 pgname
      CHARACTER*4 extension
      INTEGER un,iex,fop,iop,flg_tec
      INTEGER FILE_EXISTS,FILE_OPEN
!
      pgname='pegase'
      localinfo(1:1)=' '
!
!     Building the output file name
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      filename=FULL_FILE_NAME(out_dir,outname,extension,flg_os)
      localinfo(12:70)=filename(1:59)
!
!     Checking that the file exists. If an I/O error occurs,
!     the program will be stopped by the call to print_error
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      iex=FILE_EXISTS(filename)      
      IF (iex.EQ.-1) THEN
          CALL PRINT_ERROR (9,pgname,localinfo,1)
      ENDIF
!
!     Select the file type based on extension and check for first
!     open instruction (fop=0)
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      fop=0
      IF (extension.EQ.'.log') THEN
          fop=fop_log
      ELSEIF (extension.EQ.'.out') THEN
          fop=fop_out
      ELSEIF (extension.EQ.'.dat') THEN
          fop=fop_dat
      ELSEIF (extension.EQ.'.da2') THEN
          fop=fop_da2
      ELSEIF (extension.EQ.'.res') THEN
          fop=fop_cust
      END IF
!
!     If the file has already been opened, proceed to append
!     If an error occurs, possible to continue as new (fop=0)
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (fop.EQ.-1) THEN
          IF (iex.EQ.0) THEN
              CALL PRINT_ERROR (18,pgname,localinfo,-1)
              fop=0
          ELSE
              iop=FILE_OPEN (filename,un,2)
              IF (iop.EQ.-4) THEN
                  CALL PRINT_ERROR (19,pgname,localinfo,-1)
                  fop=0
              END IF
          END IF
      END IF
!
!     Opening the file as new, overwriting if needed
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (fop.EQ.0) THEN
          IF (iex.EQ.0) THEN
              iop=FILE_OPEN (filename,un,0)
              IF (iop.EQ.0) THEN
                  CALL PRINT_ERROR (11,pgname,localinfo,1)
	        END IF
          ELSE
              iop=FILE_OPEN (filename,un,-1)
              IF (iop.EQ.-2) THEN
                  CALL PRINT_ERROR (12,pgname,localinfo,1)
	        END IF
	    END IF
      END IF
!
!     Printing the titles and setting fop to 1 (the file is now open)
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (extension.EQ.'.log') THEN
          IF (fop.EQ.0) THEN
              WRITE (un,2000)
              CALL PRINT_VERSION (un)
              WRITE (un,2005)
              WRITE (un,2006)
          END IF
          fop_log=1
      ELSEIF (extension.EQ.'.out') THEN
          IF (fop.EQ.0) THEN
              WRITE (un,2001)
              CALL PRINT_VERSION (un)
              WRITE (un,2005)
          END IF
          fop_out=1
      ELSEIF (extension.EQ.'.dat') THEN
          IF (fop.EQ.0) THEN
              IF (flg_tec.EQ.0) THEN  
                  WRITE (un,2002)
                  CALL PRINT_VERSION (un)
                  WRITE (un,2005)
              ELSE
                  WRITE (un,2011)
              END IF
          END IF
          fop_dat=1
      ELSEIF (extension.EQ.'.da2') THEN
          IF (fop.EQ.0) THEN
              IF (flg_tec.EQ.0) THEN  
                  WRITE (un,2003)
                  CALL PRINT_VERSION (un)
                  WRITE (un,2005)
              ELSE
                  WRITE (un,2011)
              END IF
          END IF
          fop_da2=1
      ELSEIF (extension.EQ.'.res') THEN
          IF (fop.EQ.0) THEN
              IF (flg_tec.EQ.0) THEN  
                  WRITE (un,2004)
                  CALL PRINT_VERSION (un)
                  WRITE (un,2005)
              ELSE
                  WRITE (un,2011)
              END IF
          END IF
          fop_cust=1
      END IF
      RETURN
!     ------------------------------------------------------------------
2000  FORMAT (' /////////////////////////////////////////////////////',&
     &'////////////////////////',/,' //                              ',&
     &'                                           //',/,               &
     &' //            P E G A S E  ---  message recording log file .L',&
     &'OG            //',/,' //                        ',              &
     &'                                                 //')
2001  FORMAT (' /////////////////////////////////////////////////////',&
     &'////////////////////////',/,' //                              ',&
     &'                                           //',/,               &
     &' //            P E G A S E  ---  formatted text output file .O',&
     &'UT            //',/,' //                        ',              &
     &'                                                 //')
2002  FORMAT (' /////////////////////////////////////////////////////',&
     &'////////////////////////',/,' //                              ',&
     &'                                           //',/,               &
     &' //            P E G A S E  ---  tabulated data output file .D',&
     &'AT            //',/,' //                        ',              &
     &'                                                 //')
2003  FORMAT (' /////////////////////////////////////////////////////',&
     &'////////////////////////',/,' //                              ',&
     &'                                           //',/,               &
     &' //        P E G A S E  ---  second tabulated data output file',&
     &' .DA2         //',/,' //                        ',              &
     &'                                                 //')
2004  FORMAT (' /////////////////////////////////////////////////////',&
     &'////////////////////////',/,' //                              ',&
     &'                                           //',/,               &
     &' //         P E G A S E  ---  user-defined custom output file', &
     &' .RES          //',/,' //                        ',             &
     &'                                                 //')
2005  FORMAT (' //                        ',                           &
     &'                                                 //',/,         &
     &' /////////////////',                                            &
     &'////////////////////////////////////////////////////////////',/)
2006  FORMAT (/,'STARTING DUMPING SCREEN OUTPUT:',/,&
                   &'------------------------------')
!     ------------------------------------------------------------------
!     TECPLOT FORMAT
!
2011  FORMAT (' TITLE = " Output file generated by PEGASE v4.5 " ')
!     ------------------------------------------------------------------
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  OUTPUT   L I B R A R Y             //
! //                                                                  //
! //  SUBROUTINE close_output_files                                   //
! //                                                                  //
! //  The routine closes all files previously open and affect -1 to   //
! //  the corresponding fop variable in light of future re-openings.  //
! //                                                                  //
! //  Benoit Bottin, 28/11/96. Rewritten 18/06/97.                    //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE close_output_files
!
      USE global_pegase
      IMPLICIT NONE
      CHARACTER*70 localinfo
      CHARACTER*6 pgname
      INTEGER FILE_CLOSE
!
      pgname='pegase'
      IF (fop_log.EQ.1) THEN
          IF (FILE_CLOSE(66).EQ.-1) THEN
              localinfo(12:27)='Log file .LOG   '
              CALL PRINT_ERROR (13,pgname,localinfo,1)        
          END IF
          fop_log=-1
      ENDIF
      IF (fop_out.EQ.1) THEN
          IF (FILE_CLOSE(21).EQ.-1) THEN
              localinfo(12:27)='Output file .OUT'
              CALL PRINT_ERROR (13,pgname,localinfo,1)        
          END IF
          fop_out=-1
      ENDIF
      IF (fop_dat.EQ.1) THEN
          IF (FILE_CLOSE(22).EQ.-1) THEN
              localinfo(12:33)='Tabulated file 1: .DAT'
              CALL PRINT_ERROR (13,pgname,localinfo,1)        
          END IF
          fop_dat=-1
      ENDIF
      IF (fop_da2.EQ.1) THEN
          IF (FILE_CLOSE(23).EQ.-1) THEN
              localinfo(12:33)='Tabulated file 2: .DA2'
              CALL PRINT_ERROR (13,pgname,localinfo,1)        
          END IF
          fop_da2=-1
      ENDIF
      IF (fop_cust.EQ.1) THEN
          IF (FILE_CLOSE(24).EQ.-1) THEN
              localinfo(12:27)='Custom file .RES'
              CALL PRINT_ERROR (13,pgname,localinfo,1)        
          END IF
          fop_cust=-1
      ENDIF
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  SETUP   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE print_version                                        //
! //                                                                  //
! //  Input:  un          unit of the file being printed to           //    
! //                                                                  //
! //  Outs:   none                                                    //
! //                                                                  //
! //  The routine is a simple routine opening, reading and closing    //
! //  the file PEGASE.VER located in the same directory as the        //
! //  executable. It retrieves the relevant version information and   //
! //  prints it in the output file headers                            //
! //                                                                  //
! //  Benoit Bottin, 04/12/96. Modified 18/06/97.                     //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE print_version (un)
      
      IMPLICIT NONE
      CHARACTER*80 vername
      CHARACTER*70 localinfo
      CHARACTER*6 pgname
      LOGICAL iex
      INTEGER size,i,un
!
      pgname='pegase'
      localinfo(1:1)=' '            
      INQUIRE (FILE='pegase.ver',EXIST=iex)
      IF (.NOT.iex) THEN
          CALL PRINT_ERROR (20,pgname,localinfo,0)
          WRITE (un,2001)
      ENDIF
      OPEN (UNIT=33,FILE='pegase.ver',ERR=801)
      READ (33,1000) size
      DO i=1,size
          READ (33,1001,err=901,end=901) vername
          WRITE (un,2002) vername
      END DO
      CLOSE (33)
      RETURN
801   CALL PRINT_ERROR (21,pgname,localinfo,0)
      WRITE (un,2001)
      RETURN
901   CALL PRINT_ERROR (22,pgname,localinfo,0)
      WRITE (un,2001)
      RETURN
!     ------------------------------------------------------------------      
1000  FORMAT (i4)
1001  FORMAT (a80)      
2001  FORMAT (' //                     (cannot find the version file.)',&
     &'                    //')
2002  FORMAT (a80)
!     ------------------------------------------------------------------      
      END
! //////////////////////////////////////////////////////////////////////
