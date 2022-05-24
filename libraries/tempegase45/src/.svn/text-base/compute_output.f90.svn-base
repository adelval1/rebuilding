! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //               P E G A S E   4.   M A I N   C O D E               //
! //                                                                  //
! //  SUBROUTINE fill_array1                                          //
! //                                                                  //
! //  Computation of species fractional parameters                    //
! //                                                                  //
! //  Inputs: nsp     number of species in the mixture                //
! //          oper    if 0, single-species, other means multi-species //
! //          mode    if 1 it is p,T input, if 2 it is rho,T input    //
! //         flg_stop allows to stop on subroutine errors             //
! //          v       pressure or density depending on mode           //
! //          T       temperature (K)                                 //
! //          TNEQ    nonequilibrium temperatures array (K)           //
! //          x       mole or mass fraction array, depending on mode  //
! //                                                                  //
! //  Output: array1  fractional species properties (10,maxsp):       //
! //          1       mole fraction                                   //
! //          2       mass fraction                                   //
! //          3       partial pressure                                //
! //          4       partial density                                 //
! //          5       number density                                  //
! //          6       concentration                                   //
! //          7       translation temperature                         //
! //          8       rotation temperature                            //
! //          9       vibration temperature                           //
! //          10      electronic temperature                          //
! //                                                                  //
! //  This subroutine computes the mentioned outputs and stores them  //
! //  in array1. If the operation is single-species, then mole/mass   //
! //  fractions must be imposed instead of converted for convenience. //
! //                                                                  //
! //  Benoit Bottin, 16/10/96. Modified and cleaned, 28/7/97.         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE fill_array1 (oper,mode,flg_neq,flg_stop,&
     &p,rho,v,T,TNEQ,x)
!      
      USE global_pegase
      USE global_thermo
      USE interf_thermo
      IMPLICIT NONE
      REAL(kind=8) p,rho,v,T,TNEQ(4),x(:),y(nsp),mm
      INTEGER oper,mode,i,flg_neq,flg_stop
!
      y=0.0
!
!     Computation of global pressure and density
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (mode.EQ.1) THEN
          IF (oper.EQ.0) THEN
              p=v
              rho=SPECIES_DENSITY(1,p,T)
          ELSE
              p=v
              rho=MIXTURE_DENSITY(p,T,x,flg_stop)
          ENDIF
      ELSE
          IF (oper.EQ.0) THEN
              rho=v
              p=SPECIES_PRESSURE(1,rho,T)
          ELSE
              rho=v
              p=MIXTURE_PRESSURE(rho,T,x,flg_stop)
          ENDIF
      ENDIF
!
!     Computing species properties
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      DO i=1,nsp
!
      IF (oper.EQ.0) THEN
          array1(1,i)=1.0d0
          array1(2,i)=1.0d0
      ELSE
          mm=MIXTURE_MOLARMASS (x,mode,flg_stop)
          IF (mode.EQ.1) THEN
              array1(1,i)=x(i)
              CALL MOLEFRAC_TO_MASSFRAC (x,mm,y)
              array1(2,i)=y(i)
          ELSE
              array1(2,i)=x(i)
              CALL MASSFRAC_TO_MOLEFRAC (x,mm,y)
              array1(1,i)=y(i)
          ENDIF
      ENDIF
      IF (mode.EQ.1) THEN
          array1(3,i)=v*array1(1,i)
          array1(4,i)=SPECIES_DENSITY(i,array1(3,i),T)
      ELSE
          array1(4,i)=v*array1(2,i)
          array1(3,i)=SPECIES_PRESSURE(i,array1(4,i),T)
      ENDIF
      array1(5,i)=NUMBER_DENSITY(array1(3,i),T)
      array1(6,i)=CONCENTRATION(array1(3,i),T)
      IF (flg_neq.EQ.0) THEN
          array1(7,i)=T
          array1(8,i)=T
          array1(9,i)=T
          array1(10,i)=T
      ELSE
          array1(7,i)=TNEQ(1)
          array1(8,i)=TNEQ(2)
          array1(9,i)=TNEQ(3)
          array1(10,i)=TNEQ(4)
      ENDIF
!
      END DO
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //               P E G A S E   4.   M A I N   ! O D E               //
! //                                                                  //
! //  SUBROUTINE fill_array2                                          //
! //                                                                  //
! //  Computation of species thermodynamic parameters                 //
! //                                                                  //
! //  Inputs: nsp     number of species in the mixture                //
! //          T       temperature (K)                                 //
! //          TNEQ    nonequilibrium temperatures array (K)           //
! //          x       mole or mass fraction array, depending on mode  //
! //          epsilon finite difference tolerance                     //
! //          array1  used for input purposes                         //
! //                                                                  //
! //  Flags:  _anha   anharmonicity correction flag                   //
! //          _neq    indicates a nonequilibrium calculation          //
! //          _stop   forces the program to stop on library errors    //
! //          _termo  defines the thermodynamic computation used      //
! //                                                                  //
! //  Output: array2  species thermodynamic properties (16,maxsp,8)   //
! //          1,9     internal energy per mole, per mass              //
! //          2,10    helmholtz free energy per mole, per mass        //
! //          3,11    gibbs free energy per mole, per mass            //
! //          4,12    enthalpy per mole, per mass                     //
! //          5,13    entropy per mole, per mass                      //
! //          6,14    chemical potential per mole, per mass           //
! //          7,15    constant volume specific heat per mole, mass    //
! //          8,16    constant pressure specific heat per mole, mass  //
! //                                                                  //
! //  This subroutine computes the mentioned outputs and stores them  //
! //  in array2.                                                      //
! //                                                                  //
! //  Benoit Bottin, 16/10/96. Modified and cleaned F90, 28/7/97.     //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE fill_array2 (flg_anha,flg_neq,flg_stop,flg_termo,&
     &T,TNEQ,p,x,mode)
!      
      USE global_pegase
      USE global_thermo
      USE interf_thermo
      IMPLICIT NONE
      REAL(kind=8) p,T,x(:),TNEQ(4)
      REAL(kind=8) ener(8),enth(8),helm(8),gibb(8),pot(8),entro(8),cv(8),cp(8)
      REAL(kind=8) ptot,mm
      REAL(kind=8) parfn(8),dparfn(8),ddparfn(8),dlnpf(8),ddlnpf(8)
      INTEGER mode,flg_anha,flg_neq,flg_stop,flg_termo,j,i,isp
!
!     computation of mixture pressure
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      ptot=p
!
!     loop on the number of species and initialization
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      DO isp=1,nsp
          DO i=1,8
              ener(i)=0.0d0
              enth(i)=0.0d0
              entro(i)=0.0d0
              helm(i)=0.0d0
              gibb(i)=0.0d0
              pot(i)=0.0d0
          END DO
          p=array1(3,isp)
!
!         computation of properties per unit mole
!         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          IF (p.GT.0.0d0) THEN
              CALL species_energy (isp,T,TNEQ,&
     &        flg_anha,1,flg_neq,flg_stop,flg_termo,ener)
              CALL species_enthalpy (isp,T,TNEQ,&
     &        flg_anha,1,flg_neq,flg_stop,flg_termo,enth)
              CALL species_entropy (isp,p,T,TNEQ,&
     &        flg_anha,1,flg_neq,flg_stop,flg_termo,entro)
              CALL species_helmholtz (isp,p,T,TNEQ,&
     &        flg_anha,1,flg_neq,flg_stop,flg_termo,helm)
              CALL species_gibbs (isp,p,T,TNEQ,&
     &        flg_anha,1,flg_neq,flg_stop,flg_termo,gibb)
              CALL species_gibbs (isp,ptot,T,TNEQ,&
     &        flg_anha,1,flg_neq,flg_stop,flg_termo,pot)
               CALL species_cv (isp,T,TNEQ,&
     &flg_anha,1,flg_neq,flg_stop,flg_termo,cv)
               CALL species_cp (isp,T,TNEQ,&
     &flg_anha,1,flg_neq,flg_stop,flg_termo,cp)
          ENDIF
          DO i=1,8
              array2(1,isp,i)=ener(i)
              array2(2,isp,i)=helm(i)
              array2(3,isp,i)=gibb(i)
              array2(4,isp,i)=enth(i)
              array2(5,isp,i)=entro(i)
              array2(6,isp,i)=pot(i)
              array2(7,isp,i)=cv(i)
              array2(8,isp,i)=cp(i)
          END DO
!
!         computation of properties per unit mass
!         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          mm=MIXTURE_MOLARMASS (x,mode,flg_stop)
          DO i=1,8
              DO j=1,8
                  array2(i+8,isp,j)=MOL_TO_MASS(array2(i,isp,j),isp,mm,flg_stop)
              END DO
          END DO
      END DO
!
!         computation of partition functions
!         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      DO i=1,nsp
          p=ptot*(x(i)+1.d-32)
          CALL partition_function (i,p,T,TNEQ,flg_anha,1,&
       &  flg_neq,flg_stop,flg_termo,parfn,dparfn,ddparfn,dlnpf,ddlnpf)
          DO j=1,8
              array2(17,i,j)=parfn(j)
              array2(18,i,j)=dparfn(j)
              array2(19,i,j)=ddparfn(j)
              array2(20,i,j)=dlnpf(j)
              array2(21,i,j)=ddlnpf(j)
          END DO
      END DO
!              
      p=ptot
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //               P E G A S E   4.   M A I N   C O D E               //
! //                                                                  //
! //  SUBROUTINE fill_array_mix                                       //
! //                                                                  //
! //  Computation of mixture properties arrays array3 and array4      //
! //                                                                  //
! //  Inputs: nsp     number of species in the mixture                //
! //          nr      number of chemical reactions                    //
! //          nc      number of nuclei conservation equations         //
! //          mode    if 1 it is p,T input, if 2 it is rho,T input    //
! //          v       pressure or density depending on mode           //
! //          T       temperature (K)                                 //
! //          TNEQ    nonequilibrium temperatures array (K)           //
! //          x       mole or mass fraction array, depending on mode  //
! //                                                                  //
! //  Output: array3  non-splittable mixture properties               //
! //          1       pressure                                        //
! //          2       density                                         //
! //          3       temperature                                     //
! //          4       molar mass                                      //
! //          5       specific gas constant                           //
! //          6       mole ratio                                      //
! //          7       compressibility factor                          //
! //          8       number density                                  //
! //          9       gamma frozen (molar)                            //
! //          10      gamma equilibrium (molar)                       //
! //          11      gamma frozen (massic)                           //
! //          12      gamma equilibrium (massic)                      //
! //          13      frozen speed of sound                           //
! //          14      equilibrium speed of sound                      //
! //          15      kappa                                           //
! //          16      khi                                             //
! //          17      translation temperature                         //
! //          18      rotation temperature                            //
! //          19      vibration temperature                           //
! //          20      electronic temperature                          //
! //          21      full dissociation enthalpy per mass             //
! //          22      dissociation enthalpy per mass of non-ions only //
! //                                                                  //
! //          array4  splittable mixture properties                   //
! //          1,10    internal energy per mole, per mass              //
! //          2,11    helmholtz free energy per mole, per mass        //
! //          3,12    gibbs free energy per mole, per mass            //
! //          4,13    enthalpy per mole, per mass                     //
! //          5,14    entropy per mole, per mass                      //
! //          6,15    cv, frozen per mole, per mass                   //
! //          7,16    cv, equilibrium per mole, per mass              //
! //          8,17    cp, frozen per mole, per mass                   //
! //          9,18    cp, equilibrium per mole, per mass              //
! //                                                                  //
! //  This subroutine computes the mentioned outputs and stores them  //
! //  in arrays array3 and array4.                                    //
! //                                                                  //
! //  Benoit Bottin, 16/10/96                                         //
! //                                                                  //
! //  Modifications to include the use of flg_oper to avoid crashing  //
! //  when using it single-specied...                                 //
! //                                                                  //
! //  Benoit Bottin, 30/10/96. Modified and cleaned F90, 28/7/97.     //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE fill_array_mix (mode,flg_anha,flg_neq,flg_log,&
     &flg_oper,flg_stop,flg_termo,p,rho,T,TNEQ,x,epsilon,findiff,niter,zo)
!      
      USE global_pegase
      USE global_thermo
      USE interf_thermo
      IMPLICIT NONE
      REAL(kind=8) p,rho,T,TNEQ(4),x(:),y(nsp),mm,epsilon,findiff
      REAL(kind=8) ener(8),enth(8),entro(8),helm(8),gibbs(8)
      REAL(kind=8) cvfrx(8),cvfry(8),cveqx(8),cveqy(8),cpfrx(8),cpfry(8)
      REAL(kind=8) cpeqx(8),cpeqy(8),a,kappa,ki,zo
      REAL(kind=8) dmdr,dpdr,dedr(8),dhdr(8),detdr(8)
      REAL(kind=8) dmdt,dpdt,dedt(8),dhdt(8),detdt(8)
      INTEGER mode,i,flg_anha,flg_neq,flg_oper
      INTEGER flg_stop,flg_termo,j,niter,flg_log
!
!     If oper=0 then artificially build up the mole/mass fraction array
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (flg_oper.EQ.0) THEN
          x(1)=1.0d0
          y(1)=x(1)
      ENDIF
!
!     Computation of basic mixture properties
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      
      IF (mode.EQ.2) THEN
          DO i=1,nsp
              y(i)=x(i)
          END DO
          array3(6)=MIXTURE_MOLERATIO (x,mode,flg_stop)
      ELSEIF (mode.EQ.1) THEN
          array3(6)=zo
      ENDIF
      array3(1)=p
      array3(2)=rho
      array3(3)=T
      array3(4)=MIXTURE_MOLARMASS (x,mode,flg_stop)
      array3(5)=MIXTURE_RSPECIFIC (x,mode,flg_stop)
      IF (flg_oper.NE.1) THEN
          array3(6)=1.0d0
      ENDIF
      array3(7)=array3(6)
      array3(8)=NUMBER_DENSITY (array3(1),T)
!
      mm=array3(4)
      IF (mode.EQ.1) THEN
          CALL MOLEFRAC_TO_MASSFRAC (x,mm,y)
!      CALL MASSFRAC (y,rho,t,TNEQ,epsilon,flg_log,&
!     &flg_anha,flg_stop,flg_termo,1,x,niter)
      ELSE
          CALL MASSFRAC_TO_MOLEFRAC (y,mm,x)
      ENDIF
!
!     Computation of mixture thermodynamic properties per mole
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      CALL mixture_energy (T,TNEQ,x,&
     &flg_anha,1,flg_neq,flg_stop,flg_termo,ener)
      CALL mixture_enthalpy (T,TNEQ,x,&
     &flg_anha,1,flg_neq,flg_stop,flg_termo,enth)
      CALL mixture_entropy (p,T,TNEQ,x,&
     &flg_anha,1,flg_neq,flg_stop,flg_termo,entro)
      CALL mixture_helmholtz (p,T,TNEQ,x,&
     &flg_anha,1,flg_neq,flg_stop,flg_termo,helm)
      CALL mixture_gibbs (p,T,TNEQ,x,&
     &flg_anha,1,flg_neq,flg_stop,flg_termo,gibbs)
!
      DO i=1,8
          array4(1,i)=ener(i)
          array4(2,i)=helm(i)
          array4(3,i)=gibbs(i)
          array4(4,i)=enth(i)
          array4(5,i)=entro(i)
      END DO
!
!     Computation of mixture thermodynamic properties per mass
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      DO i=1,8
          DO j=1,5
              array4(j+9,i)=MOL_TO_MASS(array4(j,i),0,mm,flg_stop)
          END DO
      END DO
!
!     Computation of frozen cp and cv per mole and mass
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      CALL cv_frozen (T,TNEQ,x,&
     &flg_anha,1,flg_neq,flg_stop,flg_termo,cvfrx)
      CALL cp_frozen (T,TNEQ,x,&
     &flg_anha,1,flg_neq,flg_stop,flg_termo,cpfrx)
      CALL cv_frozen (T,TNEQ,y,&
     &flg_anha,2,flg_neq,flg_stop,flg_termo,cvfry)
      CALL cp_frozen (T,TNEQ,y,&
     &flg_anha,2,flg_neq,flg_stop,flg_termo,cpfry)
      DO i=1,8
!         Filling up array values
          array4(6,i)=cvfrx(i)
          array4(8,i)=cpfrx(i)
          array4(15,i)=cvfry(i)
          array4(17,i)=cpfry(i)
!         Preparing energy per mass to pass to cv_equ and cp_equ
          ener(i)=array4(10,i)
          enth(i)=array4(13,i)
      END DO
!
!     Computation of equilibrium cp and cv per mole and mass
!     If nonequilibrium mixtures, set equilibrium as frozen
!     First, a test is made to use Gupta's equilibrium cp per mass      
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (flg_termo.EQ.2) THEN
          array4(18,1)=GUPTA_CP(p,T)
          GOTO 3000
      ENDIF
      IF (flg_oper.NE.1) THEN
          DO i=1,8
              array4(7,i)=cvfrx(i)
              array4(9,i)=cpfrx(i)
              array4(16,i)=cvfry(i)
              array4(18,i)=cpfry(i)
          END DO
      ELSE
      CALL sensitivities_r (p,rho,T,TNEQ,ener,enth,mm,y,findiff,&
     &epsilon,flg_anha,flg_neq,flg_stop,flg_termo,flg_log,niter,&
     &dmdr,dpdr,dedr,dhdr,detdr)
      CALL sensitivities_t (p,rho,T,TNEQ,ener,enth,mm,y,findiff,&
     &epsilon,flg_anha,flg_neq,flg_stop,flg_termo,flg_log,niter,&
     &dmdt,dpdt,dedt,dhdt,detdt)
          CALL cv_equilibrium (1,ener,mm,dedt,dmdt,cveqx)
          CALL cv_equilibrium (2,ener,mm,dedt,dmdt,cveqy)
          CALL cp_equilibrium (1,enth,mm,dhdt,dhdr,dmdt,dmdr,&
     &    dpdt,dpdr,cpeqx)
          CALL cp_equilibrium (2,enth,mm,dhdt,dhdr,dmdt,dmdr,&
     &    dpdt,dpdr,cpeqy)
          DO i=1,8
              array4(7,i)=cveqx(i)
              array4(9,i)=cpeqx(i)
              array4(16,i)=cveqy(i)
              array4(18,i)=cpeqy(i)
          END DO
      ENDIF
!
!     Computation of gamma
!     ^^^^^^^^^^^^^^^^^^^^      
      array3(9)=array4(8,1)/array4(6,1)
      array3(10)=array4(9,1)/array4(7,1)
      array3(11)=array4(17,1)/array4(15,1)
      array3(12)=array4(18,1)/array4(16,1)
      array3(23)=array3(10)*rho/p*dpdr
      array3(24)=array3(12)*rho/p*dpdr
!
      array3(13)=DSQRT(array3(11)*array3(5)*T)
      IF (flg_oper.NE.1) THEN
          a=array3(13)
      ELSE
          CALL sound_speed_equilibrium (enth,dpdt,dpdr,detdt,detdr,&
     &    a,kappa,ki)
      END IF
      array3(14)=a
      array3(15)=kappa
      array3(16)=ki
      IF (flg_neq.EQ.1) THEN
          DO i=1,4
              array3(16+i)=TNEQ(i)
          END DO
      ELSE
          DO i=1,4
              array3(16+i)=T
          END DO
      ENDIF
!
      IF (flg_termo.GE.2) THEN
          array3(21)=DISSOC_ENTHALPY (y,array3(6),2,3)
          array3(22)=array3(21)
      ELSE     
          array3(21)=DISSOC_ENTHALPY (y,array3(6),2,1)
          array3(22)=DISSOC_ENTHALPY (y,array3(6),2,2)
      END IF
!      
3000  RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //               P E G A S E   4.   M A I N   ! O D E               //
! //                                                                  //
! //  SUBROUTINE fill_array6                                          //
! //                                                                  //
! //  Computation of mixture transport properties                     //
! //                                                                  //
! //  Inputs: nsp     number of species in the mixture                //
! //          T       temperature (K)                                 //
! //          TNEQ    nonequilibrium temperatures array (K)           //
! //          x       mole or mass fraction array, depending on mode  //
! //          epsilon finite difference tolerance                     //
! //          array1  used for input purposes                         //
! //                                                                  //
! //  Flags:  _anha   anharmonicity correction flag                   //
! //          _neq    indicates a nonequilibrium calculation          //
! //          _stop   forces the program to stop on library errors    //
! //          _termo  defines the thermodynamic computation used      //
! //                                                                  //
! //  Output: array2  species thermodynamic properties (16,maxsp,8)   //
! //          1,9     internal energy per mole, per mass              //
! //          2,10    helmholtz free energy per mole, per mass        //
! //          3,11    gibbs free energy per mole, per mass            //
! //          4,12    enthalpy per mole, per mass                     //
! //          5,13    entropy per mole, per mass                      //
! //          6,14    chemical potential per mole, per mass           //
! //          7,15    constant volume specific heat per mole, mass    //
! //          8,16    constant pressure specific heat per mole, mass  //
! //                                                                  //
! //  This subroutine computes the mentioned outputs and stores them  //
! //  in array2.                                                      //
! //                                                                  //
! //  Benoit Bottin, 16/10/96. Modified and cleaned F90, 28/7/97.     //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE fill_array6 (p,T,T_NEQ,flg_neq, &
    & flg_termo,flg_traco,flg_anha)
!
      USE global_pegase
      USE global_thermo
      USE interf_thermo
      USE interf_traco
      USE global_traco
      USE traco_ord_var, only: n_ord
      IMPLICIT NONE
!
      REAL(kind=8) :: p,T,n_tot,mm,rho
      INTEGER flg_termo,flg_traco,i,flg_neq,flg_anha
      REAL(kind=8) :: x(1:nsp),n(1:nsp),T_NEQ(4)
      REAL(kind=8) cpr(1:nsp),cpv(1:nsp),cpe(1:nsp)
      REAL(kind=8) HINT(1:nsp), HINT_ord(1:nsp)
      REAL(kind=8) mfp(1:nsp),kt,kr,kv,ke,mu,cp,k,cp_sp(8)

!     Some initialisations
      kt=0.0
      kr=0.0
      kv=0.0
      ke=0.0
      n_tot=array3(8)
!
!     The use of GUPTA curve fits if flg_traco is set to 1
!     Missing information is obviously set to zero
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      array6 = 0.0d0
      IF (flg_traco.EQ.1) THEN
          array6(1,1)=GUPTA_VISCO(p,T)
          array6(2,1)=GUPTA_K(p,T)
          array6(3,1)=GUPTA_PRANDTL(p,T)
          RETURN
      ENDIF
!
!   Use of the true TRACO routines of Pegase
!   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!   a) initializing the non-equilibrium temperature array in case
!      of equilibrium (it should still be zero here!)
!      ----------------------------------------------
    IF (flg_neq.EQ.0) THEN
        T_NEQ=T
    END IF
!
!   b) preparing the arrays of mole fractions, number densities and
!      internal Cps to pass to Init-Traco. If the anharmonicity 
!      corrections were in use, then there is a need to compute
!      the species cps with harmonic oscillator option
!      --------------------------------------------------------
    x(1:nsp)=array1(1,1:nsp)
    n(1:nsp)=array1(5,1:nsp)
    mfp=0.0d0
    IF (flg_termo.NE.0) THEN         ! curve fits are used to get
        cpr(1:nsp)=array2(8,1:nsp,1) ! the cp so the cpr, cpv and
        cpv(1:nsp)=0.0d0             ! cpe parts are not known.
        cpe(1:nsp)=0.0d0             ! By default, equ. operation!
    ELSE  !-- down, no curve fits -- !
        IF (flg_anha.NE.0) THEN      ! harmonic data not yet computed
            DO i=1,nsp
                CALL species_cp (i,T,T_NEQ,0,1,flg_neq,1,flg_termo,cp_sp)
                cpr(i)=cp_sp(3)
                cpv(i)=cp_sp(4)
                cpe(i)=cp_sp(5)
            END DO          
        ELSE                         ! harmonic data already computed
            cpr(1:nsp)=array2(8,1:nsp,3)
            cpv(1:nsp)=array2(8,1:nsp,4)
            cpe(1:nsp)=array2(8,1:nsp,5)
        END IF
    END IF
    HINT = 0.; HINT_ord = 0.
    HINT (1:nsp) = array2(4,1:nsp,1) + array2(4,1:nsp,8)
    CALL Set_traco(x,n,cpr,cpv,cpe,HINT,T_NEQ,1)
!
!   c) simple parameters calculations
!      ------------------------------    
    array6(1,1)=viscosity()
    if (nelectrons.ne.0) then
    array6(4,1)=electrical_conductivity (n_ord(nsp),T_NEQ(4))
    array6(5,1)=Debye_Length_Devoto (n_ord(nsp),n_tot,T_NEQ(4))
    endif
!
!   d) thermal conductivity calculations
!      ---------------------------------
    CALL Yos_lambda_frozen(kt,kr,kv,ke)
    array6(2,2)=kt
    array6(2,3)=kr
    array6(2,4)=kv
    array6(2,5)=ke
    if (nelectrons.ne.0) then
    array6(2,6)=Devoto_lambda_electron (n_ord(nsp), T_NEQ(4))
    endif
    array6(2,7)=Yos_lambda_reactive (T)
    array6(2,1)=0.0d0
    DO i=2,7
        array6(2,1)=array6(2,1)+array6(2,i)
    END DO
!
!   e) other transport parameters calculations
!      ---------------------------------------
    mu=array6(1,1)
    cp=array4(17,1)
    k=array6(2,1)-array6(2,7)
    array6(3,1)=PRANDTL_NUMBER (mu,cp,k)
    k=array6(2,1)
    cp=array4(18,1)
    array6(3,2)=PRANDTL_NUMBER (mu,cp,k)
!
    CALL mfp_Yos (n_tot,mfp)
    array6(6,1:nsp)=mfp(1:nsp)
!
!   f) Binary diffusion coefficients. mfp_array is used again to 
!      retrieve the coefficients from the subroutine
!      ---------------------------------------------------------
    mm=array3(4)
    CALL diffusion_coefficients_multi(p,T_NEQ,mm,mfp)
    array6(7,1:nsp)=mfp(1:nsp)
    rho=array3(2)
    DO i=1,nsp
        array6(8,i)=LEWIS_NUMBER (k,rho,cp,mfp(i))
    END DO
!
!   g) Energy relaxation times
!      -----------------------
    array6(9,1)=DLOG10(DEVOTO_TAUEH (T, T))
    array6(9,2)=DLOG10(MILLIKAN_TAUVH (p, T))

    RETURN
    END
! //////////////////////////////////////////////////////////////////////
