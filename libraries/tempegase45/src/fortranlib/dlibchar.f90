!     //////////////////////////////////////////////////////////////////
!     ///                                                            ///
!     ///   "DUALE" FORTRAN LIBRARY                                  ///
!     ///   =======================                                  ///
!     ///                                                            ///
!     ///   SUBROUTINE lcase                                         ///
!     ///                                                            ///
!     ///   Inputs:   namein  a string to convert to lower case      ///
!     ///                                                            ///
!     ///   Outputs:  nameout the string converted to lowercase      ///
!     ///             flg_err error values corresponding to:         ///
!     ///                 1   correct operation of the subroutine    ///
!     ///                 0   nameout longer string warning          ///
!     ///                -1   nameout is shorter, operation aborted  ///
!     ///                                                            ///
!     ///   B. Bottin, 16/12/96                                      ///
!     ///                                                            ///
!     //////////////////////////////////////////////////////////////////
      SUBROUTINE lcase (namein,nameout,flg_err)
      
      IMPLICIT NONE
      CHARACTER*(*) namein,nameout  
      CHARACTER*1 letter
      INTEGER i,lin,lout,value,flg_err
      
      lin=LEN(namein)
      lout=LEN(nameout)
      IF (lout.LT.lin) THEN
          flg_err=-1
          RETURN
      ENDIF

      DO  i=1,lin
          letter=namein(i:i)
          value=ICHAR(letter)
          IF ((value.GE.65).AND.(value.LE.90)) THEN
              value=value+32
              letter=CHAR(value)
          ENDIF
          nameout(i:i)=letter
      END DO
      IF (lout.GT.lin) THEN
          flg_err=0
      ELSE
          flg_err=1
      ENDIF

      RETURN
      END
!     //////////////////////////////////////////////////////////////////

!     //////////////////////////////////////////////////////////////////
!     ///                                                            ///
!     ///   "DUALE" FORTRAN LIBRARY                                  ///
!     ///   =======================                                  ///
!     ///                                                            ///
!     ///   FUNCTION fill_with_blanks                                ///
!     ///                                                            ///
!     ///   Inputs:   var     a string to be filled with blanks      ///
!     ///   Output:   var     the same string but with all blanks    ///
!     ///                                                            ///
!     ///   B. Bottin, 04/10/96                                      ///
!     ///                                                            ///
!     //////////////////////////////////////////////////////////////////
      SUBROUTINE fill_with_blanks (var)
!      
      IMPLICIT NONE
      CHARACTER*(*) var
      INTEGER i,nb
!      
      nb=LEN(var)
      DO i=1,nb
          var(i:i)=' '
      END DO
!
      RETURN
      END
!     //////////////////////////////////////////////////////////////////

!     //////////////////////////////////////////////////////////////////
!     ///                                                            ///
!     ///   "DUALE" FORTRAN LIBRARY                                  ///
!     ///   =======================                                  ///
!     ///                                                            ///
!     ///   FUNCTION lentrim                                         ///
!     ///                                                            ///
!     ///   Inputs:   name    the string of which length is sought   ///
!     ///   Output:   f_value the length of the string without right ///
!     ///                     trailing blanks                        ///
!     ///                                                            ///
!     ///   B. Bottin, 10/06/97                                      ///
!     ///                                                            ///
!     //////////////////////////////////////////////////////////////////
      INTEGER FUNCTION lentrim (name)
!
      IMPLICIT NONE
      CHARACTER*(*) name
      INTEGER l,efflen
!
      l=LEN(name)
!
!     Get the real length of the string. Since the string can
!     have multiple spaces, this length is counted from right
!     -------------------------------------------------------     
      efflen=l
      DO WHILE (name(efflen:efflen).EQ.' ')
          efflen=efflen-1
          IF (efflen.EQ.0) THEN
              EXIT
          ENDIF
      END DO
!
      lentrim=efflen
      RETURN
      END
!     //////////////////////////////////////////////////////////////////

!     //////////////////////////////////////////////////////////////////
!     ///                                                            ///
!     ///   "DUALE" FORTRAN LIBRARY                                  ///
!     ///   =======================                                  ///
!     ///                                                            ///
!     ///   SUBROUTINE centerstring                                  ///
!     ///                                                            ///
!     ///   Inputs:   name    the string to be centered on itself    ///
!     ///   Output:   nameout the same string but with all blanks    ///
!     ///             flg_err an error code -1 string is all blank   ///
!     ///                                   -2 nameout too short     ///
!     ///                                    1 correct operation     ///
!     ///                                                            ///
!     ///   B. Bottin, 04/10/96                                      ///
!     ///                                                            ///
!     //////////////////////////////////////////////////////////////////
      SUBROUTINE centerstring (name,nameout,flg_err)
!      
      IMPLICIT NONE
      CHARACTER*(*) name,nameout
      INTEGER l,efflen,left,right,i,flg_err,LENTRIM
!
!     Step 1. Get the real length of the string. If nameout is not
!             long enough, produce an error flag and exit
!     ------------------------------------------------------------ 
      flg_err=1
      efflen=LENTRIM(name)
      IF (efflen.EQ.0) THEN
          flg_err=-1
      ENDIF
      l=LEN(nameout)
      IF (efflen.GT.l) THEN
          flg_err=-2
          RETURN
      ENDIF
!
!     Step 2. Compute left and right, the number of spaces to leave at
!             left and right of the centered string, and shift string
!     ----------------------------------------------------------------
      right=INT((l-efflen)/2)
      left=l-efflen-right
      DO i=1,left
          nameout(i:i)=' '
      END DO
      nameout(left+1:left+efflen)=name(1:efflen)
      DO i=left+efflen+1,l
          nameout(i:i)=' '
      END DO
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

!     //////////////////////////////////////////////////////////////////
!     ///                                                            ///
!     ///   "DUALE" FORTRAN LIBRARY                                  ///
!     ///   =======================                                  ///
!     ///                                                            ///
!     ///   SUBROUTINE concatenate                                   ///
!     ///                                                            ///
!     ///   Input:  str1,str2   strings to concatenate               ///
!     ///   Output: str3        concatenation of 1//2 without blanks ///
!     ///           flg_err     error code -1 str3 not long enough   ///
!     ///                                   1 correct operation      ///
!     ///                                                            ///
!     ///   B. Bottin, 02/10/96                                      ///
!     ///                                                            ///
!     //////////////////////////////////////////////////////////////////
      SUBROUTINE concatenate (str1,str2,str3,flg_err)
      
      IMPLICIT NONE
      CHARACTER*(*) str1,str2,str3
      INTEGER flg_err,l3,leff1,leff2,LENTRIM
!
!     Step 1. Compute length and effective length of the strings
!     ----------------------------------------------------------      
      l3=LEN(str3)
      leff1=LENTRIM(str1)
      leff2=LENTRIM(str2)
!
!     Step 2. Check that str3 is long enough for the concatenate
!     ----------------------------------------------------------
      IF ((leff1+leff2).GT.l3) THEN
          flg_err=-1
          RETURN
      ENDIF
!
!     Step 3. Produce the concatenation
!     ---------------------------------
      str3(1:leff1)=str1(1:leff1)
      str3(leff1+1:leff1+leff2)=str2(1:leff2)
      flg_err=1
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////

