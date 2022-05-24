! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  SETUP    L I B R A R Y             //
! //                                                                  //
! //  SUBROUTINE define_custom                                        //
! //                                                                  //
! //  Input:  in_dir      working directory for the input file        //    
! //          formatname  name of the format description file         //
! //          flg_cust    if 1, prints the custom .RES output file    //
! //                                                                  //
! //  Outs:   none                                                    // 
! //                                                                  //
! //  The routine opens and reads the format description file for     //
! //  the .RES results file.                                          //
! //                                                                  //
! //  Benoit Bottin, 28/07/97                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE define_custom (in_dir,formatname,flg_cust)
      USE global_pegase
      IMPLICIT NONE
!
      INTEGER nb,a,i,j,k,flg_cust
      INTEGER FILE_EXISTS,FILE_OPEN,FILE_CLOSE
      CHARACTER*120 filename
      CHARACTER*80 in_dir,formatname
      CHARACTER*70 localinfo
      CHARACTER*6 pgname
      CHARACTER*1 firstchar
!
!     Concatenation of the full path and name of the format file
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      pgname='pegase'
      CALL CONCATENATE (in_dir,formatname,filename,i)
      localinfo(1:70)=filename(1:70)
!
!     Checking if the file exists
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      i=FILE_EXISTS(filename)
      IF (i.EQ.-1) THEN
          CALL PRINT_ERROR (31,pgname,localinfo,1)
      ELSEIF (i.EQ.0) THEN
          CALL PRINT_ERROR (32,pgname,localinfo,-1)
          flg_cust=0
          RETURN
      ENDIF
!
!     Opening the file
!     ^^^^^^^^^^^^^^^^
      i=FILE_OPEN (filename,34,1)
      IF (i.LE.0) THEN
          CALL PRINT_ERROR (33,pgname,localinfo,1)
      ENDIF
!
!     Analyzing the format file
!     ^^^^^^^^^^^^^^^^^^^^^^^^^
      custmax=0
101   CONTINUE
      READ (34,1001,end=901) firstchar,nb,label(nb),a,i,j,k
      IF (firstchar.EQ.'@') THEN
          custmax=custmax+1
          custom(nb,1)=a
          custom(nb,2)=i
          custom(nb,3)=j
          custom(nb,4)=k
          IF (custmax.EQ.99) THEN
              CALL PRINT_ERROR (36,pgname,localinfo,1)
              GOTO 901
          END IF
      END IF
      GOTO 101
901   i=FILE_CLOSE(34)
      IF (i.EQ.-1) THEN
          CALL PRINT_ERROR (34,pgname,localinfo,1)
      ENDIF
!     ------------------------------------------------------------------
1001  FORMAT (a1,1x,i2,2x,a14,8x,i2,3x,i3,3x,i3,3x,i3)
!     ------------------------------------------------------------------
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  OUTPUT   L I B R A R Y             //
! //                                                                  //
! //  SUBROUTINE print_custom                                         //
! //                                                                  //
! //  Input:  array1      array of output data                        //    
! //          array2      array of output data                        //
! //          array3      array of output data                        //
! //          array4      array of output data                        //
! //          userarray   array of output data (user-defined)         //
! //                                                                  //
! //  Outs:   none                                                    // 
! //                                                                  //
! //  The routine prints output to the .RES file based on the         //
! //  indications contained in the custom format file.                //
! //                                                                  //
! //  Benoit Bottin, 28/10/96. Modified and cleaned, 28/7/97.         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE print_custom (cntr,cntr2,casename,flg_tec,nmax)
!
      USE global_pegase
      USE global_thermo
      IMPLICIT NONE
      CHARACTER*40 intname(99),item !pietro before 14
      CHARACTER*80 casename
      CHARACTER*80 titlename
      CHARACTER*200 tectitle
      INTEGER cntr,cntr2,i,flg_tec,nmax,p,l
      REAL(kind=8) value
!     
!     Open output file
!     ^^^^^^^^^^^^^^^^
      CALL FILL_WITH_BLANKS (titlename)
!
!     If first calculation of the scan, print number/casename
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (cntr2.EQ.1) THEN
          IF (flg_tec.NE.1) THEN
              CALL CENTERSTRING (casename(1:69),titlename,i)
              WRITE (24,*)
              WRITE (24,2002) cntr,titlename
!
!     Building the titles
!     ^^^^^^^^^^^^^^^^^^^
              DO i=1,custmax
                  item=label(i)
                  WRITE (intname(i),2006) item
              END DO
              WRITE (24,2010) (intname(i),i=1,custmax) !pietro 2010 instead of *
          ELSE
              CALL FILL_WITH_BLANKS(tectitle)
              p = 1
              DO i = 1, custmax
                  l = LEN_TRIM(label(i))
                  tectitle(p    :p    ) = '"'
                  tectitle(p+1  :p+l  ) = label(i)(1:l)
                  tectitle(p+l+1:p+l+2) = '" '
                  p = p+l+3
              END DO
              WRITE (24,2021) tectitle
              CALL FILL_WITH_BLANKS(titlename)
              i = MIN0(LEN_TRIM(casename),68)
              titlename(1:i) = casename(1:i)
              titlename(i+1:i+1) = '"'
              WRITE (24,2022) nmax,titlename
          ENDIF
      ENDIF
    !
!     Print output values in tab format
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      DO i=1,custmax
          IF (custom(i,1).EQ.1) THEN
              value=array1(int(custom(i,2)),int(custom(i,3)))
          ELSEIF (custom(i,1).EQ.2) THEN
              value=array2(int(custom(i,2)),int(custom(i,3)),int(custom(i,4)))
          ELSEIF (custom(i,1).EQ.3) THEN
              value=array3(int(custom(i,2)))
          ELSEIF (custom(i,1).EQ.4) THEN
              value=array4(int(custom(i,2)),int(custom(i,3)))
          ELSEIF (custom(i,1).EQ.5) THEN
              value=array5(int(custom(i,2)),int(custom(i,3)))
          ELSEIF (custom(i,1).EQ.6) THEN
              value=array6(int(custom(i,2)),int(custom(i,3)))
          ELSEIF (custom(i,1).EQ.0) THEN
              value=userarray(int(custom(i,2)))
          ENDIF
          WRITE (intname(i),*) value !pietro - 2011 instead of *
      END DO
      WRITE (24,2010) (intname(i),i=1,custmax) !pietro - 2010 instead of *
!     ------------------------------------------------------------------
2001  FORMAT (' /////////////////////////////////////////////////////',&
     &'////////////////////////',/,' //                              ',&
     &'                                           //',/,               &
     &' //         P E G A S E  ---  user-defined custom output file', &
     &' .RES          //',/,' //                        ',             &
     &'                                                 //')           
2003  FORMAT (' //                        ',                           &
     &'                                                 //',/,         &
     &' /////////////////',                                            &
     &'////////////////////////////////////////////////////////////',/)
2002  FORMAT (' ======================================================',&
     &'=======================',/,' ==>              Pegase now proces',&
     &'ses input case number ',i4,':             <==',/,' ==> ',a69,    &
     &' <==',/,' ==================================================',   &
     &'===========================',/)                                  
2006  FORMAT (a40)   !14 instead of 30                                                   
2010  FORMAT (1x,99a40)!99 invece di 100 - 14 instead of 30
2011  FORMAT (e30.20) !- old line !pietro
!2011  FORMAT (f24.10)
2021  FORMAT (' VARIABLES = ',a200)
2022  FORMAT (' ZONE I=',i6,' F=POINT T="',a69)
!     ------------------------------------------------------------------
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
