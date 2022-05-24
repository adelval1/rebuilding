REAL(kind=8) FUNCTION millikan_tauvh (p, T_h)

!******************************************************************************
! Vibrational translational energy relaxation. See e.g. original article
! by Millikan and White + Gnoffo et al.
!******************************************************************************

  USE global_thermo 
  USE global_traco
  USE traco_ord_var
  IMPLICIT NONE

!******************************************************************************

! Input variables 
  REAL*8  :: p, T_h

! Internal variables 
  INTEGER :: nhv, i, j
  REAL*8  :: n_sum, n_mol, tau (nsp), p_atm
  REAL*8  :: mu (nsp, nsp), A (nsp, nsp), thetav (nsp), num
  REAL*8  :: v_therm (nsp), sigma

!******************************************************************************

  p_atm = p / 101325.d0

  IF (nions.EQ.0) THEN
     nhv = nsp
  ELSE
     nhv = nsp - 1
  ENDIF

  n_sum = 0.d0
  DO i = 1, nhv
     n_sum = n_sum + n_ord (i)
  ENDDO

  ! Set reduced molecular weights of colliding species
  mu = 0.d0
  DO i = 1, nsp
     DO j = 1, nsp
        mu (i,j) = MMOL_ord (i) * MMOL_ord (j) / (MMOL_ord (i)+MMOL_ord (j))
     ENDDO
  ENDDO
  ! (in g / mole)
  mu (:,:) = mu (:,:) * 1000.d0

  thetav = 0.d0
  ! Set characteristic vibrational temperatures
  CALL Order_data (nsp,index_array, TV(:,1,:), thetav (:))

  ! Set A-coefficients
  A = 0.d0
  DO i = 1, nsp
     DO j = 1, nsp
        A (i,j) = 1.16d-3 * mu (i,j)**0.5d0 * thetav (i)**1.33d0
     ENDDO
  ENDDO

  tau = 0.d0
  DO i = 1, nhv
     IF (thetav(i).GT.0.d0) THEN
        DO j = 1, nhv
        tau (i) = tau (i) + n_ord (j) * dexp &
                & (  A(i,j) * (T_h**(-0.33d0) - 0.015d0 * mu(i,j)**(0.25d0)) - 18.42d0  )
!& (  A(i,j) * (T_h**-0.33d0 - 0.015d0 * mu(i,j)**0.25d0) - 18.42d0  ) pietro 
        ENDDO
        tau (i) = tau (i) / (n_sum * p_atm)
     ENDIF
  ENDDO

! Apply Park's correction here:

  v_therm = 0.d0
  DO i = 1, nhv
     v_therm (i) = dsqrt ( (8.d0*KUNIV*T_h)/(PIUNIV*MMOL_ord(i)/NAUNIV))
  ENDDO
  sigma = 1.d-21

  DO i = 1, nhv
     IF (thetav(i).GT.0.d0) THEN
        tau (i) = tau (i) + (sigma*v_therm(i)*n_sum)**(-1)
! pietro       tau (i) = tau (i) + (sigma*v_therm(i)*n_sum)**-1	
     ENDIF
  ENDDO
  
  num = 0.d0; n_mol = 0.d0
  DO i = 1, nhv
     IF (thetav(i).GT.0.d0) THEN
        num = num + n_ord (i) / tau (i)
        n_mol = n_mol + n_ord (i)
     ENDIF
  ENDDO

  IF (nmolecn.NE.0) THEN
     millikan_tauvh = n_mol / num
  ELSE
     millikan_tauvh = 0.d0
  ENDIF

!******************************************************************************

END FUNCTION millikan_tauvh
