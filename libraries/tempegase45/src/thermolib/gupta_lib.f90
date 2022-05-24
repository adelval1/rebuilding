! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  SETUP   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE load_gupta_fits                                      //
! //                                                                  //
! //  In:     traco_dir       path for the thermodynamic files        // 
! //                                                                  //
! //  The routine opens Gupta fits data files as units 14-17, loads   //
! //  fits information in COMMON blocks and closes the files.         //
! //                                                                  //
! //  Benoit Bottin, 04/12/96. Cleaned and updated DUALE, 4/8/97.     //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE load_gupta_fits (traco_dir,flg_stop)
!
      USE global_thermo
      IMPLICIT NONE
      CHARACTER*80 traco_dir,filename,openname
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      INTEGER i,j,l,CHECK_GUPTA_FILES,iv,ik,ipr,icp,flg_stop,istop
!
!     Zeroing all gupta arrays
!     ^^^^^^^^^^^^^^^^^^^^^^^^
      pgname='pegaslib'
      WRITE (*,3001)
      iv=0
      ik=0
      ipr=0
      icp=0
      istop=0
      DO 101 i=1,37
          DO 201 j=1,8
              VMAT(i,j)=0.0d0
              KMAT(i,j)=0.0d0
              PRMAT(i,j)=0.0d0
              CPMAT(i,j)=0.0d0
201       CONTINUE
101   CONTINUE
!
!     Verifying that the Gupta fit files are available
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      CALL FILL_WITH_BLANKS (openname)
      CALL FILL_WITH_BLANKS (filename)
      filename(1:11)='visco.gpt  '
      localinfo(1:11)='           '
      localinfo(12:70)=filename(1:59)
      iv=CHECK_GUPTA_FILES (filename,traco_dir)
      IF (iv.EQ.0) THEN
          CALL PRINT_ERROR (37,pgname,localinfo,0)
      ENDIF
      filename(1:11)='conduct.gpt'
      localinfo(12:70)=filename(1:59)
      ik=CHECK_GUPTA_FILES (filename,traco_dir)
      IF (ik.EQ.0) THEN
          CALL PRINT_ERROR (37,pgname,localinfo,0)
      ENDIF
      filename(1:11)='prandtl.gpt'
      localinfo(12:70)=filename(1:59)
      ipr=CHECK_GUPTA_FILES (filename,traco_dir)
      IF (ipr.EQ.0) THEN
          CALL PRINT_ERROR (37,pgname,localinfo,0)
      ENDIF
      filename(1:11)='cp.gpt     '
      localinfo(12:70)=filename(1:59)
      icp=CHECK_GUPTA_FILES (filename,traco_dir)
      IF (icp.EQ.0) THEN
          CALL PRINT_ERROR (37,pgname,localinfo,0)
      ENDIF
!
!     Opening the files
!     ^^^^^^^^^^^^^^^^^
      IF (iv.EQ.1) THEN
          filename(1:11)='visco.gpt  '
          localinfo(12:70)=filename(1:59)
          CALL CONCATENATE (traco_dir,filename,openname,l)
          l=INDEX(openname,' ')
          OPEN (file=openname(1:l-1),unit=18,err=901)
          GOTO 1
901       CALL PRINT_ERROR (38,pgname,localinfo,0)
          iv=0
          istop=1
      ENDIF
1     IF (ik.EQ.1) THEN
          filename(1:11)='conduct.gpt'
          localinfo(12:70)=filename(1:59)
          CALL CONCATENATE (traco_dir,filename,openname,l)
          l=INDEX(openname,' ')
          OPEN (file=openname(1:l-1),unit=15,err=902)
          GOTO 2
902       CALL PRINT_ERROR (38,pgname,localinfo,0)
          ik=0
          istop=1
      ENDIF
2     IF (ipr.EQ.1) THEN
          filename(1:11)='prandtl.gpt'
          localinfo(12:70)=filename(1:59)
          CALL CONCATENATE (traco_dir,filename,openname,l)
          l=INDEX(openname,' ')
          OPEN (file=openname(1:l-1),unit=16,err=903)
          GOTO 3
903       CALL PRINT_ERROR (38,pgname,localinfo,0)
          ipr=0
          istop=1
      ENDIF
3     IF (icp.EQ.1) THEN
          CALL FILL_WITH_BLANKS (openname)
          filename(1:11)='cp.gpt     '
          localinfo(12:70)=filename(1:59)
          CALL CONCATENATE (traco_dir,filename,openname,l)
          l=INDEX(openname,' ')
          OPEN (file=openname(1:l-1),unit=17,err=904)
          GOTO 4
904       CALL PRINT_ERROR (38,pgname,localinfo,0)
          icp=0
          istop=1
      ENDIF
!
!     Retrieving the values for the opened files
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
4     CONTINUE
      DO 103 i=1,37
          IF (iv.EQ.1) THEN 
              READ (18,1001) (vmat(i,j),j=1,8)
          ENDIF
          IF (ik.EQ.1) THEN
              READ (15,1001) kmat(i,1),kmat(i,2),(kmat(i,j),j=7,3,-1)&
             &,kmat(i,8)
          ENDIF
          IF (ipr.EQ.1) THEN
              READ (16,1001) (prmat(i,j),j=1,8)
          ENDIF
          IF (icp.EQ.1) THEN
              READ (17,1001) cpmat(i,1),cpmat(i,2),(cpmat(i,j),j=7,3,-1)&
             &,cpmat(i,8)
          ENDIF
103   CONTINUE
!
!     Closing the files
!     ^^^^^^^^^^^^^^^^^
      CLOSE (15)
      CLOSE (16)
      CLOSE (17)
      CLOSE (18)
!
!     Checking for error status
!     ^^^^^^^^^^^^^^^^^^^^^^^^^
      IF ((istop.EQ.1).AND.(flg_stop.EQ.1)) THEN
          STOP 'Program ended due to PEGASE transport library error'
      ENDIF
!     ------------------------------------------------------------------
1001  FORMAT (f6.1,f11.1,6f18.10)
3001  FORMAT (' --> Reading Gupta transport properties data files ...')
!     ------------------------------------------------------------------
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  SETUP   L I B R A R Y              //
! //                                                                  //
! //  INTEGER FUNCTION check_gupta_files                              //
! //                                                                  //
! //  In:     filename        simple name of the gupta data file      // 
! //          traco_dir       transport properties file directory     //
! //                                                                  //
! //  Out:    (function)      1 if file exists at root                //
! //                                                                  //
! //  The routine checks if the gupta transport directory is in the   //
! //  general transport coefficients directory. It returns 1 if they  //
! //  exist.                                                          //
! //                                                                  //
! //  Benoit Bottin, 04/12/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      INTEGER FUNCTION check_gupta_files (filename,traco_dir)
!
      IMPLICIT NONE
      LOGICAL iex
      INTEGER l
      CHARACTER*80 traco_dir,filename,openname
!      
!     Checking for the file
!     ^^^^^^^^^^^^^^^^^^^^^
      CALL FILL_WITH_BLANKS (openname)
      CALL CONCATENATE (traco_dir,filename,openname,l)
      l=INDEX(openname,' ')
      IF (l.EQ.0) THEN
          l=LEN(openname)+1
      ENDIF      
      INQUIRE (file=openname(1:l-1), EXIST=iex)
      IF (iex) THEN
          check_gupta_files=1
      ELSE
          check_gupta_files=0
      ENDIF
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  THERMO   L I B R A R Y             //
! //                                                                  //
! //  FUNCTION gupta_extract                                          //
! //                                                                  //
! //  In:     x               variable of the curve fit               // 
! //          p               pressure of the mixture (Pa)            //
! //          t               temperature of the mixture (K)          //
! //          MAT             curve fits matrix                       //
! //                                                                  //
! //  The function extracts a value in the Gurvich fit defined by     //
! //  MAT and returns it to the calling routine. Selection of the     //
! //  pressure range is made by lp then the temperature range is      //
! //  determined. The value is interpolated in terms of pressure.     //
! //                                                                  //
! //  Benoit Bottin, 05/12/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION gupta_extract(x,p,t,MAT)
!
      IMPLICIT NONE
      CHARACTER*70 localinfo
      CHARACTER*8 pgname
      REAL(kind=8) x,p,t,MAT(37,8),lp,left,right,y
      INTEGER cntr,k,lft,rgt,leftcntr,rightcntr,i
      CHARACTER*3 switch
!
      pgname='pegaslib'
      localinfo(1:1)=' '
      lp=DLOG10(p)
      IF ((lp.LT.-4).OR.(lp.GT.2)) THEN
          CALL PRINT_ERROR (39,pgname,localinfo,0)
          RETURN
      ENDIF
      k=IDINT(lp)
      IF (k.GT.0) THEN
          lft=k
          rgt=k+1
          IF (k.EQ.2) THEN
              rgt=2
          ENDIF
      ELSE
          lft=k-1
          rgt=k
          IF (k.EQ.-4) THEN
              lft=-4
          ENDIF
      ENDIF
!      
      switch='lft'
      k=lft
5     CONTINUE
          cntr=0
10        CONTINUE
              cntr=cntr+1
              IF (k.EQ.IDINT(MAT(cntr,1))) THEN
                  IF (t.LT.MAT(cntr,2)) THEN
                      GOTO 20
                  ENDIF
              ENDIF
              GOTO 10
!         END LOOP
20        CONTINUE
          IF (switch.EQ.'lft') THEN
              leftcntr=cntr
              switch='rgt'
              k=rgt
              GOTO 5
          ELSE          
              rightcntr=cntr
          ENDIF
!     END LOOP
!
      i=leftcntr
      left=MAT(i,3)+MAT(i,4)*x+MAT(i,5)*x**2+MAT(i,6)*x**3&
     &+MAT(i,7)*x**4+MAT(i,8)*x**5
      i=rightcntr
      right=MAT(i,3)+MAT(i,4)*x+MAT(i,5)*x**2+MAT(i,6)*x**3&
     &+MAT(i,7)*x**4+MAT(i,8)*x**5
      IF (lft.NE.rgt) THEN
          y=(right-left)/(rgt-lft)*(lp-lft)+left
      ELSE
          y=(right+left)/2
      ENDIF
      gupta_extract=y
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  THERMO   L I B R A R Y             //
! //                                                                  //
! //  FUNCTION gupta_visco                                            //
! //                                                                  //
! //  In:     p               pressure of the mixture (Pa)            //
! //          t               temperature of the mixture (K)          //
! //                                                                  //
! //  The function returns the value of viscosity of the 11-species   //
! //  air mixture according to Gupta fits in NASA RP-1260.            //
! //                                                                  //
! //  Benoit Bottin, 05/12/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION gupta_visco(p,T)                                          
!
      USE global_thermo
      IMPLICIT NONE
      REAL(kind=8) T,p,pp,x,gupta_extract
!
      pp=0.0
      IF (T.LT.500.0) THEN
          gupta_visco=1.4584e-6*(t**0.5)/(1+(110.33/t))
      ELSE
          x=t/1000.0
          pp=p/101325.0
          gupta_visco=GUPTA_EXTRACT(x,pp,T,VMAT)*0.1
      ENDIF
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  THERMO   L I B R A R Y             //
! //                                                                  //
! //  FUNCTION gupta_k                                                //
! //                                                                  //
! //  In:     p               pressure of the mixture (Pa)            //
! //          t               temperature of the mixture (K)          //
! //                                                                  //
! //  The function returns the value of conductivity of the           //
! //  11-species air mixture according to Gupta fits in NASA RP-1260. //
! //                                                                  //
! //  Benoit Bottin, 05/12/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION gupta_k (p,T)                                          
!
      USE global_thermo
      IMPLICIT NONE
      REAL(kind=8) T,p,pp,x,gupta_extract
!
      pp=0.0
      IF (T.LT.500.0) THEN
          gupta_k=5.9776e-6*(t**0.5)/(1+(194.4/t))*418.4
      ELSE
          x=DLOG(t/10000.0)
          pp=p/101325.0
          gupta_k=DEXP(GUPTA_EXTRACT(x,pp,T,KMAT))*418.4
      ENDIF
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  THERMO   L I B R A R Y             //
! //                                                                  //
! //  FUNCTION gupta_cp                                               //
! //                                                                  //
! //  In:     p               pressure of the mixture (Pa)            //
! //          t               temperature of the mixture (K)          //
! //                                                                  //
! //  The function returns the value of specific heat cp of the       //
! //  11-species air mixture according to Gupta fits in NASA RP-1260. //
! //                                                                  //
! //  Benoit Bottin, 05/12/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION gupta_cp(p,T)                                          
!
      USE global_thermo
      IMPLICIT NONE
      REAL(kind=8) T,p,pp,x,gupta_extract
!
      pp=0.0
      IF (T.LT.500.0) THEN
          gupta_cp=0.24*4184.0
      ELSE
          x=DLOG(t/10000.0)
          pp=p/101325.0
          gupta_cp=DEXP(GUPTA_EXTRACT(x,pp,T,CPMAT))*4184.0
      ENDIF
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  THERMO   L I B R A R Y             //
! //                                                                  //
! //  FUNCTION gupta_prandtl                                          //
! //                                                                  //
! //  In:     p               pressure of the mixture (Pa)            //
! //          t               temperature of the mixture (K)          //
! //                                                                  //
! //  The function returns the value of the Prandtl number of the     //
! //  11-species air mixture according to Gupta fits in NASA RP-1260. //
! //                                                                  //
! //  Benoit Bottin, 05/12/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      REAL(kind=8) FUNCTION gupta_prandtl(p,T)                                          
!
      USE global_thermo
      IMPLICIT NONE    
      REAL(kind=8) T,p,pp,x,gupta_extract,v,c
!
      pp=0.0
      IF (T.LT.500.0) THEN
          v=1.4584e-5*(t**0.5)/(1+(110.33/t))
          c=5.9776e-6*(t**0.5)/(1+(194.4/t))
          gupta_prandtl=0.24*v/c
      ELSE
          x=t/1000.0
          pp=p/101325.0
          gupta_prandtl=GUPTA_EXTRACT(x,pp,T,PRMAT)
      ENDIF
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
