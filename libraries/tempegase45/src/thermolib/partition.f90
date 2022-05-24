! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE partition_function                                   //
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
! //  Output: parfn(8)      array of partition functions              //
! //          dparfn(8)     array of first derivatives                //
! //          ddparfn(8)    array of second derivatives               //
! //              1           complete                                //
! //              2           translation                             //
! //              3           rotation (anha=0)                       //
! //              4           vibration (anha=0)                      //
! //              5           electronic (anha=0)                     //
! //              6           coupled rotation-vibration (anha=1)     //
! //              7           coupled internal (anha=2)               //
! //              8           not used                                //
! //                                                                  //
! //  This subroutine provides the entropy of the species             //
! //  isp according to the flags specified.                           //
! //                                                                  //
! //  Benoit Bottin, 11/10/96. Modified and cleaned 25/7/97.          //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE partition_function (isp,var,T,TNEQ,flg_anha,mode,&
     &flg_neq,flg_stop,flg_termo,parfn,dparfn,ddparfn,dlnpf,ddlnpf)
!
      USE global_thermo
      IMPLICIT NONE
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      INTEGER isp,flg_anha,mode,flg_neq,flg_stop,flg_termo,istop,i,atom
      REAL(kind=8) var,T,TNEQ(4),ttra,trot,tvib,tele,dlnq,ddlnq
      REAL(kind=8) parfn(8),dparfn(8),ddparfn(8),dlnpf(8),ddlnpf(8),q,dq,ddq
!
!     Checks of correct variable values and initialize
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      pgname='pegaslib'
      istop=0
      IF (flg_neq.EQ.0) THEN
          IF (T.LE.0) THEN
              WRITE (localinfo,101) T
101           FORMAT (11x,'Routine: partition_function. Value:',d15.8)
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
102   FORMAT(11x,'Routine:partition_function.Index:',i3,' Value:',d15.8)
                  CALL PRINT_ERROR (27,pgname,localinfo,0)
                  istop=1
              ENDIF
          END DO
          ttra=TNEQ(1)
          trot=TNEQ(2)
          tvib=TNEQ(3)
          tele=TNEQ(4)
      ENDIF
      IF (var.LE.0) THEN
          WRITE (localinfo,103) var
103       FORMAT (11x,'Routine: partition_function. Value:',d15.8)
          CALL PRINT_ERROR (24+mode,pgname,localinfo,0)
          istop=1
      ENDIF
      DO i=1,8
          parfn(i)=1.0d0
          dparfn(i)=0.0d0
          ddparfn(i)=0.0d0
          dlnpf(i)=0.0d0
          ddlnpf(i)=0.0d0
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
              localinfo(1:36)='           Routine: partition_function.'
                  CALL PRINT_ERROR (28,pgname,localinfo,0)
                  t=ttra
              ENDIF
!              parfn(1)=PARFN_REFTABLE_X (isp,var,t,flg_stop)
          ELSE
              CALL PARFN_TRANS_X(isp,var,ttra,q,dq,ddq,dlnq,ddlnq)
              parfn(2)=q
              dparfn(2)=dq
              ddparfn(2)=ddq
              dlnpf(2)=dlnq
              ddlnpf(2)=ddlnq
          ENDIF
      ELSE
          IF (flg_termo.EQ.1) THEN
              IF (flg_neq.NE.0) THEN
              localinfo(1:36)='           Routine: partition_function.'
                  CALL PRINT_ERROR (28,pgname,localinfo,0)
                  t=ttra
              ENDIF
!              parfn(1)=PARFN_REFTABLE_X (isp,var,t,flg_stop)
          ELSE
              CALL PARFN_TRANS_Y(isp,var,ttra,q,dq,ddq,dlnq,ddlnq)
              parfn(2)=q
              dparfn(2)=dq
              ddparfn(2)=ddq
              dlnpf(2)=dlnq
              ddlnpf(2)=ddlnq
          ENDIF
	ENDIF
      IF (flg_anha.LE.1) THEN
          IF (atom.EQ.0) THEN
              CALL PARFN_ROT(isp,trot,q,dq,ddq,dlnq,ddlnq)
              parfn(3)=q
              dparfn(3)=dq
              ddparfn(3)=ddq
              dlnpf(3)=dlnq
              ddlnpf(3)=ddlnq
              CALL PARFN_VIB(isp,tvib,q,dq,ddq,dlnq,ddlnq)
              parfn(4)=q
              dparfn(4)=dq
              ddparfn(4)=ddq
              dlnpf(4)=dlnq
              ddlnpf(4)=ddlnq
          ENDIF
          CALL PARFN_ELEC(isp,tele,q,dq,ddq,dlnq,ddlnq)
          parfn(5)=q
          dparfn(5)=dq
          ddparfn(5)=ddq
          dlnpf(5)=dlnq
          ddlnpf(5)=ddlnq
      ENDIF
      IF (flg_anha.EQ.1) THEN
          IF (atom.EQ.0) THEN
              CALL PARFN_RV(isp,trot,tvib,q,dq,ddq,dlnq,ddlnq)
              parfn(6)=q
              dparfn(6)=dq
              ddparfn(6)=ddq
              dlnpf(6)=dlnq
              ddlnpf(6)=ddlnq
          ENDIF
      ELSEIF (flg_anha.EQ.2) THEN
          IF (atom.EQ.0) THEN
              CALL PARFN_INT(isp,trot,tvib,tele,q,dq,ddq,dlnq,ddlnq)
          ELSE
              CALL PARFN_ELEC(isp,tele,q,dq,ddq,dlnq,ddlnq)
          ENDIF
          parfn(7)=q
          dparfn(7)=dq
          ddparfn(7)=ddq
          dlnpf(7)=dlnq
          ddlnpf(7)=ddlnq
      ELSEIF (flg_anha.GE.3) THEN
      localinfo(1:36)='           Routine: partition_function.'
          CALL PRINT_ERROR (29,pgname,localinfo,0)
          istop=1
      ENDIF
      DO i=2,7
          parfn(1)=parfn(1)*parfn(i)
          dparfn(1)=dparfn(1)+dparfn(i)/parfn(i)
          ddparfn(1)=ddparfn(1)+&
     &    (ddparfn(i)*parfn(i)-dparfn(i)**2)/parfn(i)**2
          dlnpf(1)=dlnpf(1)+dlnpf(i)
          ddlnpf(1)=ddlnpf(1)+ddlnpf(i)
      END DO
      ddparfn(1)=parfn(1)*(dparfn(1)**2+ddparfn(1))
      dparfn(1)=parfn(1)*dparfn(1)
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
! //  FUNCTIONS parfn_trans_x AND parfn_trans_y                       //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          p           pressure (Pa)                               //
! //          rho         density (kg/m3)                             //
! //          tt          translation temperature (K)                 //
! //                                                                  //
! //  This function returns the translation entropy per               //
! //  unit mole (p,T input) or unit mass (rho,T input).               //
! //                                                                  //
! //  Benoit Bottin, 11/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE parfn_trans_y (isp,rho,tt,q,dq,ddq,dlnq,ddlnq)
!      
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: rho,tt,q,dq,ddq,dlnq,ddlnq,sqtt
!      
      q=0.0d0
      sqtt=SQRT(tt)
      q=ST(isp)/KUNIV*tt*sqtt*MMOL(isp)/rho
      dq=1.5d0*ST(isp)/KUNIV*sqtt*MMOL(isp)/rho
      ddq=0.75d0*ST(isp)/KUNIV*MMOL(isp)/rho/sqtt
      dlnq=1.5d0/tt
      ddlnq=-1.5d0/tt/tt
!            
      RETURN
      END
!     ------------------------------------------------------------------
!     ==================================================================
!     ------------------------------------------------------------------
      SUBROUTINE parfn_trans_x (isp,p,tt,q,dq,ddq,dlnq,ddlnq)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: isp
      REAL(kind=8) :: p,tt,q,dq,ddq,dlnq,ddlnq,c,sqtt
!      
      q=0.0d0
      sqtt=SQRT(tt)
      c=ST(isp)/KUNIV*RUNIV*tt/p
      q=c*tt*sqtt
      dq=c*1.5d0*sqtt
      ddq=c*0.75d0/sqtt
      dlnq=1.5d0/tt
      ddlnq=-1.5d0/tt/tt
!            
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
      
! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS parfn_rot_x AND parfn_rot_y                           //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          trot        rotational temperature (K)                  //
! //                                                                  //
! //  This function returns the rotation entropy per                  //
! //  unit mole (p,T input) or unit mass (rho,T input).               //
! //                                                                  //
! //  Benoit Bottin, 11/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE parfn_rot (isp,trot,q,dq,ddq,dlnq,ddlnq)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER isp
      REAL(kind=8) trot,q,dq,ddq,eta,dlnq,ddlnq
!      
      eta=TR(1,isp)/trot
      q=0.0d0
      q=(1.0d0/eta+1.0d0/3.0d0+eta/15.0d0)/SIG(isp)
      dq=eta/SIG(isp)/trot*(1.0d0/eta/eta-1.0d0/15.0d0)
      ddq=2.0d0*eta/15.0d0/SIG(isp)/trot/trot
      dlnq=(1.0d0-eta*eta/15.0d0)/(1.0d0+eta/3.0d0+eta*eta/15.0d0)/trot
      ddlnq=1.0d0+2.0d0*eta/3.0d0+eta*eta/15.0d0
      ddlnq=ddlnq-4.0d0*eta**3/45.0d0-5.0d0*eta**4/225
      ddlnq=-ddlnq/(1.0d0+eta/3.0d0+eta*eta/15.0d0)**2/trot/trot
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
      
! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS parfn_vib_x AND parfn_vib_y                           //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          tvib        vibrational temperature (K)                 //
! //                                                                  //
! //  This function returns the vibration entropy per                 //
! //  unit mole (p,T input) or unit mass (rho,T input).               //
! //                                                                  //
! //  Benoit Bottin, 11/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE parfn_vib (isp,tvib,q,dq,ddq,dlnq,ddlnq)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER isp,i,j
      REAL(kind=8) tvib,u,q,dq,ddq,eu,den,dlnq,ddlnq,t_dq,t_ddq
      
      q=1.0d0
      dq=0.0d0
      ddq=0.0d0
      dlnq=0.0d0
      dlnq=0.0d0
      ddlnq=0.0d0
      DO i=1,nvibmode(isp)
         u=TV(i,1,isp)/tvib	
         eu=EXP(u)
         den=eu-1.0d0
         q=q*eu/den
      END DO    
       DO i=1,nvibmode(isp)
         u=TV(i,1,isp)/tvib	
         eu=EXP(u)
         den=eu-1.0d0
         t_dq=u*eu/tvib/den/den
        DO j=1,nvibmode (isp) 
         if (i.NE.j) THEN
           u=TV(i,1,isp)/tvib	
           eu=EXP(u)
           den=eu-1.0d0
           t_dq=t_dq*(eu/den)
         END IF 
        END DO
          dq=dq+t_dq 
        END DO   
 !             
       DO i=1,nvibmode(isp)
         u=TV(i,1,isp)/tvib	
         eu=EXP(u)
         den=eu-1.0d0
         t_ddq=2.0d0-u*(eu+1.0d0)/den
         t_ddq=-u*eu*ddq/tvib/tvib/den/den
        DO j=1,nvibmode(isp)  
         if (i.NE.j) THEN
           u=TV(i,1,isp)/tvib	
           eu=EXP(u)
           den=eu-1.0d0
           t_ddq=t_ddq*u*eu/tvib/den/den
         END IF 
        END DO
          ddq=ddq+t_ddq    
        END DO  
       DO i=1,nvibmode (isp)     
         u=TV(i,1,isp)/tvib	
         eu=EXP(u)
         den=eu-1.0d0
         dlnq=dlnq+(u/tvib/den)
         ddlnq=ddlnq+((u*eu/den-2.0d0)*u/tvib/tvib/den )
       END DO
!     
      
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
      
! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS parfn_elec_x AND parfn_elec_y                         //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          tele        electronic temperature (K)                  //
! //                                                                  //
! //  This function returns the electronic entropy per                //
! //  unit mole (p,T input) or unit mass (rho,T input).               //
! //                                                                  //
! //  Benoit Bottin, 11/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE parfn_elec (isp,tele,q,dq,ddq,dlnq,ddlnq)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER isp,i
      REAL(kind=8) entro,tele,q,dq,ddq,s1,s2,s3,ee,dlnq,ddlnq
!      
      entro=0.0d0
      s1=0.0d0
      s2=0.0d0
      s3=0.0d0
      q=0.0d0
      DO i=1,NELVL(isp)
          IF ((TE(i,isp)/tele).LT.307) THEN
                  ee=EXP(-TE(i,isp)/tele)
	    ELSE
                  ee=0.0d0
          ENDIF
          q=q+GE(i,isp)*ee
          s1=s1+TE(i,isp)*GE(i,isp)*ee/tele**2
          s2=s2+GE(i,isp)*TE(i,isp)*ee/tele**3
          s3=s3+GE(i,isp)*ee*(TE(i,isp)/tele**2)**2
      END DO
      dq=s1
      ddq=s3-2.0d0*s2
      dlnq=dq/q
      ddlnq=(ddq*q-dq*dq)/(q*q)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS parfn_rv_x AND parfn_rv_y                             //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          trot        rotational temperature (K)                  //
! //          tvib        vibrational temperature (K)                 //
! //                                                                  //
! //  This function returns the rotation-vibration coupled            //
! //  entropy per unit mole (p,T input) or unit mass (rho,T input).   //
! //                                                                  //
! //  Benoit Bottin, 11/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE parfn_rv (isp,trot,tvib,q,dq,ddq,dlnq,ddlnq)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER isp
      REAL(kind=8) trot,tvib,d,g,x,u,n,den,q,dq,ddq,lnq,dlnq,ddlnq,eu,euu
!
      G=G_A(1,isp)
      D=D_A(1,isp)
      X=X_A(1,isp)
!
      u=TV(1,1,isp)/tvib
      n=TR(1,isp)/trot
      eu=EXP(u)
      den=eu-1.0d0
      euu=eu**2
!
      lnq=G/n+D/den+2.0d0*X*u/den**2
      dlnq=G/n/trot
      dlnq=dlnq+u*(d*eu-2.0d0*X)/tvib/den**2
      dlnq=dlnq+(4.0d0*X*u*u*eu)/tvib/den**3
      ddlnq=-4*X*u*u*eu*(3.0d0*den-u*(1.0d0+2.0d0*eu))/tvib/tvib/den**4
      ddlnq=ddlnq-D*u*eu*(2.0d0-u*(eu+1.0d0)/den)/tvib/tvib/den/den
      ddlnq=ddlnq-4.0d0*X*u*(eu*(u-1.0d0)+1.0d0)/tvib/tvib/den**3
      q=EXP(lnq)
      dq=q*dlnq
      ddq=q*(dlnq*dlnq+ddlnq)
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  FUNCTIONS parfn_int_x AND parfn_int_y                           //
! //                                                                  //
! //  Input:  isp         number of the concerned species             //    
! //          trot        rotational temperature (K)                  //
! //          tvib        vibrational temperature (K)                 //
! //          tele        electronic temperature (K)                  //
! //                                                                  //
! //  This function returns the internally-coupled                    //
! //  entropy per unit mole (p,T input) or unit mass (rho,T input).   //
! //                                                                  //
! //  Benoit Bottin, 11/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE parfn_int (isp,trot,tvib,tele,q,dq,ddq,dlnq,ddlnq)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER isp,i
      REAL(kind=8) trot,tvib,tele,a,b,c,d,da,db,dc,dd,dda,ddb,ddc,ddd
      REAL(kind=8) GA,DAN,XA,u,n,e,den,eu,ee
      REAL(kind=8) q,dq,ddq,dlnq,ddlnq
!
      q=0.0d0;dq=0.0d0;ddq=0.0d0
      DO i=1,NELVL(isp)
!
          GA=G_A(i,isp)
          DAN=D_A(i,isp)    !not DA because of variable redundancy
          XA=X_A(i,isp)
          u=TV(1,i,isp)/tvib
          n=TR(i,isp)/trot
          e=TE(i,isp)/(tele*tele)
          eu=EXP(u)
          den=eu-1.0d0
          ee=1.0d0/eu
!
          a= GE(i,isp) * ee
          da= GE(i,isp) * e * ee
          dda= (1-2.0d0/tele) * (e*GE(i,isp)*ee)
!
          b= 1.0d0/(SIG(isp)*n)
          db= 1.0d0/(SIG(isp)*n*trot)
          ddb= 0.0d0
!
          c= eu / den
          dc= u*eu / (tvib*den*den)
          ddc= ( -u*eu/(tvib*tvib*den*den) ) * (2 - u*(eu+1.0d0)/den)
!
          d= 1.0d0 + n/3.0d0 + n*n/15.0d0 + GA/n + DAN/den + 2.0d0*XA*u/(den*den)
          dd= - n/(3.0d0*trot) - 2.0d0*n*n/(15.0d0*trot) + GA/(n*trot)
          dd=dd + u*(DAN*eu-2.0d0*XA)/(tvib*den*den)
          dd=dd + 4.0d0*XA*u*u*eu / (tvib*den*den*den)
          ddd= n/(3.0d0*trot*trot) + 4.0d0*n*n/(15.0d0*trot*trot)
          ddd=ddd - DAN*u*eu*(2.0d0-u*(eu+1.0d0)/den)/(tvib*den)**2
          ddd=ddd - 4.0d0*XA*u*(eu*(u-1.0d0)+1.0d0)/(tvib**2*den**3)
          ddd=ddd - (4.0d0*XA*u*u*eu/(tvib*den**2)**2) * (3.0d0*den-u*(1.0d0+2.0d0*eu))
!
          q= q + a*b*c*d
          dq=dq + da*b*c*d + a*db*c*d + a*b*dc*d + a*b*c*dd
          ddq=ddq + dda*b*c*d + a*b*ddc*d + a*b*c*ddd            &
        &         + 2.0d0 * (da*db*c*d + da*b*dc*d + da*b*c*dd   &
        &         + a*db*dc*d + a*db*c*dd + a*b*dc*dd)
!
      END DO
!
      dlnq = dq / q
      ddlnq = (ddq*q-dq*dq) / (q*q)
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

