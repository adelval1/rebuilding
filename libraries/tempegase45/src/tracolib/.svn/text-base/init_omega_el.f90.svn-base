SUBROUTINE init_omega_el
!
!***********************************************************************
! Set constants to be used when deriving Devoto's q-{mp} matrix -
! which is needed to compute the electrical conductivity and the
! electron thermal conductivity.
!***********************************************************************
      USE global_traco
      IMPLICIT NONE
!***********************************************************************
!
!  C_omega_el is the array which includes the various factors occurring in the
!  expressions for the q-{mp} (see the work of Devoto).
!
!  1st : = 1 -> Q(1,s)
!        = 2 -> Q(2,s)
!  2nd : = m
!  3rd : = p
!  4th : = s
!
!***********************************************************************
!  All constants have been multiplied by 1024 !!!!!
!***********************************************************************

      C_omega_el (1, 0, 0, 1) = 1024.

      C_omega_el (1, 0, 1, 1) = 2560.
      C_omega_el (1, 0, 1, 2) = - 3072.

      C_omega_el (1, 1, 1, 1) = 6400.
      C_omega_el (1, 1, 1, 2) = - 15360.
      C_omega_el (1, 1, 1, 3) = 12288.

      C_omega_el (1, 0, 2, 1) = 4480.
      C_omega_el (1, 0, 2, 2) = - 10752.
      C_omega_el (1, 0, 2, 3) = 6144.

      C_omega_el (1, 1, 2, 1) = 11200.
      C_omega_el (1, 1, 2, 2) = - 40320.
      C_omega_el (1, 1, 2, 3) = 58368.
      C_omega_el (1, 1, 2, 4) = - 30720.

      C_omega_el (1, 2, 2, 1) = 19600.
      C_omega_el (1, 2, 2, 2) = - 94080.
      C_omega_el (1, 2, 2, 3) = 204288.
      C_omega_el (1, 2, 2, 4) = - 215040.
      C_omega_el (1, 2, 2, 5) = 92160.


      C_omega_el (1, 0, 3, 1) = 6720.
      C_omega_el (1, 0, 3, 2) = - 24192.
      C_omega_el (1, 0, 3, 3) = 27648.
      C_omega_el (1, 0, 3, 4) = - 10240.


      C_omega_el (1, 1, 3, 1) = 16800.
      C_omega_el (1, 1, 3, 2) = - 80640.
      C_omega_el (1, 1, 3, 3) = 165888.
      C_omega_el (1, 1, 3, 4) = - 163840.
      C_omega_el (1, 1, 3, 5) = 61440.
!
! Watch out for the change of sign above!
! Could this be a typing error in the article ?
! It probably is, so the signs in the 5 constants right above
! differ from the ones in the article.
!
      C_omega_el (1, 2, 3, 1) = 29400.
      C_omega_el (1, 2, 3, 2) = - 176400.
      C_omega_el (1, 2, 3, 3) = 499968.
      C_omega_el (1, 2, 3, 4) = - 770560.
      C_omega_el (1, 2, 3, 5) = 629760.
      C_omega_el (1, 2, 3, 6) = - 215040.

      C_omega_el (1, 3, 3, 1) = 44100.
      C_omega_el (1, 3, 3, 2) = - 317520.
      C_omega_el (1, 3, 3, 3) = 1124928.
      C_omega_el (1, 3, 3, 4) = - 2311680.
      C_omega_el (1, 3, 3, 5) = 2833920.
      C_omega_el (1, 3, 3, 6) = - 1935360.
      C_omega_el (1, 3, 3, 7) = 573440.

      C_omega_el (2, 1, 1, 2) = 1024.

      C_omega_el (2, 1, 2, 2) = 1792.
      C_omega_el (2, 1, 2, 3) = - 2048.

      C_omega_el (2, 2, 2, 2) = 4928.
      C_omega_el (2, 2, 2, 3) = - 7168.
      C_omega_el (2, 2, 2, 4) = 5120.

      C_omega_el (2, 1, 3, 2) = 2016.
      C_omega_el (2, 1, 3, 3) = - 4608.
      C_omega_el (2, 1, 3, 4) = 2560.

      C_omega_el (2, 2, 3, 2) = 7560.
      C_omega_el (2, 2, 3, 3) = - 16704.
      C_omega_el (2, 2, 3, 4) = 16000.
      C_omega_el (2, 2, 3, 5) = - 7680.

!   Expression below includes Qbar(4,4). It is set equal to its limit value of
!   .5*Qbar(2,2) and C_omega_el(2,3,3,2) is raised by .5*1024=512.

      C_omega_el (2, 3, 3, 2) = 14553. + 512.
      C_omega_el (2, 3, 3, 3) = - 38880.
      C_omega_el (2, 3, 3, 4) = 50080.
      C_omega_el (2, 3, 3, 5) = - 34560.
      C_omega_el (2, 3, 3, 6) = 13440.
!
!***********************************************************************
END SUBROUTINE init_omega_el
