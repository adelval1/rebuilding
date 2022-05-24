REAL(kind=8) FUNCTION Devoto_taueh (T_h, T_e)

!******************************************************************************
! Electronic translational energy relaxation. See e.g. Gupta et al.
! Originally by T. Magin (1998); reworked by D. Vanden Abeele (July 1999).
!******************************************************************************

  USE global_thermo 
  USE global_traco
  USE traco_ord_var
  IMPLICIT NONE

!******************************************************************************
 
  INTEGER :: l  
  REAL*8  :: T_h, T_e, sum, v_therm

!******************************************************************************

  sum = 0.d0
  DO l = 1, nsp - 1
    sum = sum + n_ord(l) * omega(1,nsp,l) * mmol_ord (nsp) / mmol_ord (l)  
  ENDDO

  v_therm = dsqrt( (8.d0*KUNIV*T_e) / (PIUNIV*mmol_ord(nsp)/NAUNIV) )

  Devoto_taueh = ((8.d0/3.d0) * v_therm * sum) **(-1)!pietro-before (-1)==>- 1

!******************************************************************************

END FUNCTION Devoto_taueh
