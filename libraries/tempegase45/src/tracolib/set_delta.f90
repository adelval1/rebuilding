SUBROUTINE Set_Delta (nsp,T, T_e)
!
!***********************************************************************
! Delta-{s1,s2} = Intermediate quantities to be used in Yos' mixture 
! rule.
!***********************************************************************
      USE global_traco
      IMPLICIT NONE
!***********************************************************************
!
      INTEGER :: nsp
      REAL(kind=8) :: T, T_e
!
      REAL(kind=8) :: C1, C1h, C1e, C2,PIUNIV,R
      INTEGER :: s1, s2
      INTEGER :: Delta_nr
      PARAMETER (PIUNIV=3.141592653589)
      PARAMETER (R=8.31451)
!
!***********************************************************************
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 1. Defining some often used expressions
!
      C1h = PIUNIV * R * T
      C1e = PIUNIV * R * T_e

      DO Delta_nr = 1, 2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 2. Determine the coefficient for the Delta calculation
!
         IF (Delta_nr .EQ. 1) THEN
            C2 = 8. / 3.
         ELSE
            C2 = 16. / 5.
         END IF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 3. Calculate the Delta's.
!
         C1 = C1h

         DO s1 = 1, nsp
            DO s2 = s1, nsp

               if (s2.eq.nsp) then
               C1 = C1e  
               else
               C1 = C1h
               endif
              
               Delta (Delta_nr, s1, s2) = C2 * Sqrt &
              & (2.*mmol_ord(s1)*mmol_ord(s2)/(C1*(mmol_ord(s1)+&
              & mmol_ord(s2)))) * Omega (Delta_nr, s1, s2)

            END DO
         END DO
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 4. Mirror the Delta-array.
!
         DO s1 = 2, nsp
            DO s2 = 1, s1 - 1
!
               Delta (Delta_nr, s1, s2) = Delta (Delta_nr, s2, s1)
!
            END DO
         END DO
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 5. Go on with the next delta.
!
      END DO
!
!-----------------------------------------------------------------------
!***********************************************************************
END SUBROUTINE Set_Delta