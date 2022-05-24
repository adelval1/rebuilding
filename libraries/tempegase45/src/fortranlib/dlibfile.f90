!     //////////////////////////////////////////////////////////////////
!     ///                                                            ///
!     ///   "DUALE" FORTRAN LIBRARY                                  ///
!     ///   =======================                                  ///
!     ///                                                            ///
!     ///   FUNCTION full_file_name                                  ///
!     ///                                                            ///
!     ///   Inputs:   path    full path of the working directory     ///
!     ///             root    name of the file root                  ///
!     ///             exten   name of the file extension             ///
!     ///             flg_sys system (0=UNIX, 1=PC)                  ///
!     ///                                                            ///
!     ///   Outputs:  f_value character*120 string of full file name ///
!     ///             if the returned string has *** as 3 first      ///
!     ///             characters, then the full name is too long     ///
!     ///                                                            ///
!     ///   B. Bottin, 10/6/97                                       ///
!     ///                                                            ///
!     //////////////////////////////////////////////////////////////////
 	CHARACTER*120 FUNCTION full_file_name (path,root,exten,flg_sys)
!
      IMPLICIT NONE
      CHARACTER*(*) path,root,exten
      INTEGER lp,lr,lex,strtex,ltot,flg_sys,i,LENTRIM
      CHARACTER*1 slash,cara
      CHARACTER*120 filename
!
      IF (flg_sys.EQ.0) THEN
          slash='/'
      ELSE
          slash='\'
      ENDIF
!
!     Step 1. Get the length of the dummy character strings
!     -----------------------------------------------------
      lp=LENTRIM(path)
      lr=LENTRIM(root)
      lex=LENTRIM(exten)
!
!     Step 2a. If the path string ends with a slash, remove it
!     --------------------------------------------------------
      cara=path(lp:lp)
      IF (cara.EQ.slash) THEN
          lp=lp-1
      ENDIF
!
!     Step 2b. If the root string ends with a slash, remove it
!     --------------------------------------------------------
      cara=root(lr:lr)
      IF (cara.EQ.slash) THEN
          lr=lr-1
      ENDIF
!
!     Step 2c. If the extension starts with a '.', remove it by
!              starting the string artificially at position 2
!     ---------------------------------------------------------
      cara=exten(1:1)
      IF (cara.EQ.'.') THEN
          strtex=2
          lex=lex-1
      ENDIF
!
!     Step 3. Check if the length does not exceed the final variable
!             If not, perform the concatenation operations
!     --------------------------------------------------------------
      ltot=lp+lr+lex+2
      IF (ltot.GT.120) THEN
          filename(1:3)='***'
          DO i=4,120
              filename(i:i)=' '
          END DO
      ELSE
          filename(1:lp)=path(1:lp)
          filename(lp+1:lp+1)=slash
          filename(lp+2:lp+1+lr)=root(1:lr)
          filename(lp+lr+2:lp+lr+2)='.'
          filename(lp+lr+3:lp+lr+2+lex)=exten(strtex:strtex+lex-1)
          DO i=lp+3+lr+lex,120
              filename(i:i)=' '
          END DO
      ENDIF
!
      full_file_name=filename
      RETURN
      END
!     //////////////////////////////////////////////////////////////////

!     //////////////////////////////////////////////////////////////////
!     ///                                                            ///
!     ///   "DUALE" FORTRAN LIBRARY                                  ///
!     ///   =======================                                  ///
!     ///                                                            ///
!     ///   FUNCTION file_exists                                     ///
!     ///                                                            ///
!     ///   Inputs:   filename full name of the file to be tested    ///
!     ///                                                            ///
!     ///   Outputs:  f_value integer function output about status:  ///
!     ///               1     the specified file exists              ///
!     ///               0     the specified file does not exist      ///
!     ///              -1     there was an error during the inquire  ///
!     ///              -2     the string is full of blanks           ///
!     ///                                                            ///
!     ///   B. Bottin, 10/6/97                                       ///
!     ///                                                            ///
!     //////////////////////////////////////////////////////////////////
	INTEGER FUNCTION file_exists (filename)
!
      IMPLICIT NONE
      CHARACTER*(*) filename
      LOGICAL iex
      INTEGER out,efflen,lentrim
!
      out=0
!
!     Test for the first blank character in the file name string
!     ----------------------------------------------------------
      efflen=LENTRIM(filename)
      IF (efflen.EQ.0) THEN
          out=-2
          file_exists=-2
          RETURN
      ENDIF
!
!     Test the existence of the file
!     ------------------------------
      INQUIRE (file=filename(1:efflen),exist=iex,err=999)
      IF (iex) THEN
          out=1
      ELSE
          out=0
      ENDIF
!
!     Exit the procedure either normal or via error label
!     ---------------------------------------------------
      file_exists=out
      RETURN
!
999   out=-1
      file_exists=out
      RETURN
      END
!     //////////////////////////////////////////////////////////////////

!     //////////////////////////////////////////////////////////////////
!     ///                                                            ///
!     ///   "DUALE" FORTRAN LIBRARY                                  ///
!     ///   =======================                                  ///
!     ///                                                            ///
!     ///   FUNCTION file_open                                       ///
!     ///                                                            ///
!     ///   Inputs:   filename full name of the file to be tested    ///
!     ///             fileunit unit to allocate to the file          ///
!     ///             filetype 0:new, 1:old, -1:overwrite, 2:append  ///
!     ///                                                            ///
!     ///   Outputs:  f_value integer function output about status:  ///
!     ///               1     the open operation succeeded           ///
!     ///               0     the new file could not be opened       ///
!     ///              -1     the old file could not be opened       ///
!     ///              -2     the unknown file could not be opened   ///
!     ///              -3     the string is full of blanks           ///
!     ///              -4     the old file could not be appended     ///
!     ///                                                            ///
!     ///   B. Bottin, 10/6/97                                       ///
!     ///                                                            ///
!     //////////////////////////////////////////////////////////////////
	INTEGER FUNCTION file_open (filename,fileunit,filetype)
!
      IMPLICIT NONE
      CHARACTER*(*) filename
      INTEGER out,efflen,fileunit,filetype,lentrim
!
      out=0
!
!     Test for the first blank character in the file name string
!     ----------------------------------------------------------
      efflen=LENTRIM(filename)
      IF (efflen.EQ.0) THEN
          out=-2
          file_open=-3
          RETURN
      ENDIF
!
!     Opens the file
!     --------------
      IF (filetype.EQ.2) THEN
          OPEN (file=filename(1:efflen),unit=fileunit,err=990, &
     &    status='old',position='append')
          out=1
      ELSEIF (filetype.EQ.1) THEN
          OPEN (file=filename(1:efflen),unit=fileunit,err=991, &
     &    status='old')
          out=1
      ELSEIF (filetype.EQ.0) THEN
          OPEN (file=filename(1:efflen),unit=fileunit,err=992, &
     &    status='new')
          out=1
      ELSE
          OPEN (file=filename(1:efflen),unit=fileunit,err=993, &
     &    status='unknown')
          out=1
      ENDIF
!
!     Exit the procedure either normal or via error label
!     ---------------------------------------------------
      file_open=out
      RETURN
!
!     Error when using status=OLD, access=APPEND
990   out=-4
      file_open=out
      RETURN
!
!     Error when using status=OLD
991   out=-1
      file_open=out
      RETURN
!
!     Error when using status=NEW
992   out=0
      file_open=out
      RETURN
!
!     Error when using status=UNKNOWN
993   out=-2
      file_open=out
      RETURN
      END
!     //////////////////////////////////////////////////////////////////

!     //////////////////////////////////////////////////////////////////
!     ///                                                            ///
!     ///   "DUALE" FORTRAN LIBRARY                                  ///
!     ///   =======================                                  ///
!     ///                                                            ///
!     ///   FUNCTION file_close                                      ///
!     ///                                                            ///
!     ///   Inputs:   fileunit unit to allocate to the file          ///
!     ///                                                            ///
!     ///   Outputs:  f_value integer function output about status:  ///
!     ///               1     the open operation succeeded           ///
!     ///              -1     the file could not be closed           ///
!     ///                                                            ///
!     ///   B. Bottin, 10/6/97                                       ///
!     ///                                                            ///
!     //////////////////////////////////////////////////////////////////
	INTEGER FUNCTION file_close (fileunit)
!
      IMPLICIT NONE
      INTEGER out,fileunit
!
      out=0
!
!     Closes the file
!     --------------
      CLOSE (unit=fileunit,err=999)
      out=1
!
!     Exit the procedure either normal or via error label
!     ---------------------------------------------------
      file_close=out
      RETURN
!
999   out=-1
      file_close=out
      RETURN
      END
!     //////////////////////////////////////////////////////////////////

!     //////////////////////////////////////////////////////////////////
!     ///                                                            ///
!     ///   "DUALE" FORTRAN LIBRARY                                  ///
!     ///   =======================                                  ///
!     ///                                                            ///
!     ///   FUNCTION get_file_root                                   ///
!     ///                                                            ///
!     ///   Inputs:   name    name of a file                         ///
!     ///   Outputs:  root    name of the file root                  ///
!     ///             exten   name of the file extension             ///
!     ///             f_value error code about how code worked:      ///
!     ///                      1  correct operation                  ///
!     ///                     -1  no dot present in string           ///
!     ///                     -2  string full of blanks              ///
!     ///                     -3  the root string is too short       ///
!     ///                     -4  the exten string is too short      ///
!     ///                                                            ///
!     ///   B. Bottin, 11/6/97                                       ///
!     ///                                                            ///
!     //////////////////////////////////////////////////////////////////
	INTEGER FUNCTION get_file_root (name,root,exten)
!
      IMPLICIT NONE
      CHARACTER*(*) name,root,exten
      INTEGER flg_err,dotpos,lr,lex,effln
      INTEGER LENTRIM
!
!     Step 1. Verify that the input string is not blank
!     -------------------------------------------------
      effln=LENTRIM(name)
      IF (effln.EQ.0) THEN
          flg_err=-2
          GOTO 9999
      END IF
!
!     Step 2. Detect the position of the dot in the string
!     ----------------------------------------------------
      dotpos=INDEX(name,'.')
      IF (dotpos.EQ.0) THEN
          flg_err=-1
          GOTO 9999
      END IF
!
!     Step 3. Copy the left and right parts in root and exten
!     -------------------------------------------------------
      lr=LEN(root)
      lex=LEN(exten)
      IF ((dotpos-1).GT.lr) THEN
          flg_err=-3
          GOTO 9999
      ELSE
          root(1:dotpos-1)=name(1:dotpos-1)
      END IF
      IF ((effln-dotpos).GT.lex) THEN
          flg_err=-4
          GOTO 9999
      ELSE
          exten(1:effln-dotpos)=name(dotpos+1:effln)
      END IF
      flg_err=1
!
9999  get_file_root=flg_err
      RETURN
      END
!     //////////////////////////////////////////////////////////////////

!     //////////////////////////////////////////////////////////////////
!     ///                                                            ///
!     ///   "DUALE" FORTRAN LIBRARY                                  ///
!     ///   =======================                                  ///
!     ///                                                            ///
!     ///   SUBROUTINE skiplines                                     ///
!     ///                                                            ///
!     ///   Inputs:   un      file unit                              ///
!     ///             nb      number of lines to be skipped in unit  ///
!     ///   Outputs:  istop   status of the operation                ///
!     ///                      1  correct operation                  ///
!     ///                     -1  an error occurred during the READ  ///
!     ///                     -2  end of file was reached            ///
!     ///                                                            ///
!     ///   B. Bottin, 11/6/97                                       ///
!     ///                                                            ///
!     //////////////////////////////////////////////////////////////////
      SUBROUTINE skiplines (un,nb,istop)
!     
      IMPLICIT NONE
      INTEGER un,nb,i,istop
      CHARACTER*1 aa
!     
      DO i=1,nb
          READ (un,1001,err=901,end=902) aa
      ENDDO   
1001  FORMAT (a1)
      RETURN
!
901   istop=1
      RETURN
!
902   istop=2
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
