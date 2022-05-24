!     //////////////////////////////////////////////////////////////////
!     ///                                                            ///
!     ///   "DUALE" FORTRAN LIBRARY                                  ///
!     ///   =======================                                  ///
!     ///                                                            ///
!     ///   SUBROUTINE parse_input                                   ///
!     ///                                                            ///
!     ///   Inputs:   line    full line from the input file          ///
!     ///                                                            ///
!     ///   Outputs:  command      main command string               ///
!     ///             subcommand   subaltern commant string          ///
!     ///             argument     argument string of the command    ///
!     ///                                                            ///
!     ///   Note: command returns '!' if full line is to be ignored  ///
!     ///         command returns '"' if there is no command but a   ///
!     ///                             subcommand exists.             ///
!     ///                                                            ///
!     ///   B. Bottin, 22/02/98                                      ///
!     ///                                                            ///
!     //////////////////////////////////////////////////////////////////
SUBROUTINE parse_input (line, command, subcommand, argument)
IMPLICIT NONE
!
    CHARACTER (LEN=80), INTENT (IN)  :: line
    CHARACTER (LEN=20), INTENT (OUT) :: command, subcommand
    CHARACTER (LEN=80), INTENT (OUT) :: argument
!
    CHARACTER (LEN=80) :: str1, str2, str3
    INTEGER :: i, n
!
    CALL FILL_WITH_BLANKS (command)
    CALL FILL_WITH_BLANKS (subcommand)
    CALL FILL_WITH_BLANKS (argument)
    CALL FILL_WITH_BLANKS (str1)
    CALL FILL_WITH_BLANKS (str2)
    CALL FILL_WITH_BLANKS (str3)
!
    str1 = line
!
!   =====================================
!   === PART 1: identifying a comment ===
!   =====================================
    i = INDEX (str1,'!')
    IF (i.NE.0) THEN
        n = 80 - i + 1
        str1(i:80) = REPEAT (' ', n)
    END IF
!
!   ===================================
!   === PART 2: parsing the command ===
!   ===================================
    i = INDEX (str1,'::')
    IF (i.EQ.0) THEN
        command(1:1) = '!'
        RETURN
    END IF
!
    n = 80 - (i+2) + 1
    str3 (1:n) = str1 (i+2:80)
    str1 (i:80) = REPEAT (' ',n+2)
!
!   Now, str1 contains the command and subcommand
!        str3 contains the argument
!   ======================================
!   === PART 3: parsing the subcommand ===
!   ======================================
    i = INDEX (str1,',')
    IF (i.NE.0) THEN
        n = 80 - (i+1) + 1
        str2 (1:n) = str1 (i+1:80)
        str1 (i:80) = REPEAT (' ',n+1)
    END IF
!
!   Now, str1 contains the command 
!        str2 contains the subcommand
!   ======================================
!   === PART 4: performing lowercase   ===
!   ======================================
    str1 = ADJUSTL(str1)
    str2 = ADJUSTL(str2)
    str3 = ADJUSTL(str3)
!
    CALL LCASE (str1(1:20),command,i)
    CALL LCASE (str2(1:20),subcommand,i)
    CALL LCASE (str3,argument,i)
!
    i = INDEX(command,' ')
    IF (i.EQ.1) THEN
        command(1:1) = '"'
    END IF
!
    RETURN
END SUBROUTINE parse_input