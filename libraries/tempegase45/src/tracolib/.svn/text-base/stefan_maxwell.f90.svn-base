! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TRACO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE dvda_stefmax                                         //
! //                                                                  //
! // This subroutine calculates diffusive species mass fluxes through //
! // Sutton's iterative solution of the Stefan-Maxwell equations.     //
! // Correct 2-T formulation of Kolesnikov / Ramshaw.                 //
! // Only valid for situations of ambipolar (i.e., zero-current)      //
! // diffusion.                                                       //
! //                                                                  //
! // Input:      p                   pressure                         //
! //             x(nsp)              mole fractions                   //
! //             z(nsp)              partial pressure ratios p_i/p    //
! //             T(nsp)              translational temperatures       //
! //             d(nsp)              diffusion driving forces         //
! //             diff_ord(nsp,nsp)   binary diffusion coefficients    //
! //                                 (traco_order)                    //
! //             J_old(nsp)          initial guess for diffusive      //
! //                                 mass-fluxes                      //
! //             l_ref               Reference length                 //
! //                                 (for residual estimate)          //
! //             tresh               Residual tolerance               //
! //             n_it                Maximal number of iterations     //
! // Output:     J_new(nsp)          diffusive mass-fluxes            //
! //             E_amb               Ambipolar electric field         //
! //             resit               Stefan-Maxwell eqns. residual    //
! //             resmass             Mass conservation residual       //
! //             rescharge           Ambipolar constraint residual    //
! //             k                   Number of iterations taken       //
! //                                 to convergence                   //
! //                                                                  //
! // David Vanden Abeele, 27/4/2000                                   //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////

SUBROUTINE pDvda_stefmax (p,x,z,T,d,diff_ord,J_old,J_new,E_amb, &
 & tresh,n_it,l_ref,resit,resmass,rescharge,k)

!***************************************************************
! Calculate diffusive species mass fluxes through Sutton's
! iterative solution of the Stefan-Maxwell equations.
! Correct 2-T formulation of Kolesnikov / Ramshaw.
! Only valid for situations of ambipolar (i.e., zero-current)
! diffusion.
!***************************************************************
! By D. Vanden Abeele, 27/4/2000.
! Should be further optimised, take a look at P. Barbante's
! implementation of this routine for some useful ideas.
!***************************************************************
     USE global_thermo
     USE global_traco
     USE interf_thermo
     USE interf_traco
     IMPLICIT NONE
!***************************************************************
     ! Input / output variables
!***************************************************************
     REAL*8 :: p, x (nsp), z (nsp), T (nsp), d (nsp)
     REAL*8 :: diff_ord(nsp,nsp), J_old(nsp), J_new(nsp), E_amb
     REAL*8 :: tresh, l_ref, resit, resmass, rescharge
     INTEGER:: n_it, k
!***************************************************************
     ! Internal variables
!***************************************************************
     ! Everything needs to be reordered (in traco-order)
     REAL*8 :: x_ord (nsp), z_ord(nsp), T_ord(nsp), d_ord (nsp)
     REAL*8 :: charge_tmp(1:nsp), charge_ord(1:nsp)

     REAL*8 :: D_eff(nsp)
!    D_eff = `effective binary diffusion coefficients' (traco-order);

     REAL*8 :: N_new (nsp), N_old (nsp), N_tmp (nsp), J_tmp (nsp)
     ! N_ = number diffusive fluxes.
     ! J_tmp = temporary mass-diffusion flux

     REAL*8 :: res(nsp), res0
     ! Residual variables

     REAL*8 :: MM, sum, denom, numer
     ! MM = mixture molar mass. 

     REAL*8 :: d_ref, d_amb
!    d_ref = reference value of driving force for residual estimate
!    d_amb = Ambipolar driving force  

     INTEGER     :: i, j
!***************************************************************

! 1. Reorder into traco-format :

!    Calculate vector of charge numbers charge_ord:
     CHARGE_tmp = DBLE(CHARGE)
     CHARGE_ord = 0.d0
     CALL Order_data (nsp,index_array, CHARGE_tmp, CHARGE_ord)

     x_ord = 0.d0
     CALL Order_data (nsp,index_array,x,x_ord) 
     z_ord = 0.d0
     CALL Order_data (nsp,index_array,z,z_ord) 
     T_ord = 0.d0
     CALL Order_data (nsp,index_array,T,T_ord) 
     d_ord = 0.d0
     CALL Order_data (nsp,index_array,d,d_ord) 

! 2. Calculate effective binary diffusion coefficients (own
!    definition) in traco-order :

     DO i = 1, nsp
        D_eff (i) = 0.d0
        DO j = 1, nsp
           D_eff (i) = D_eff (i) + z_ord (j) / diff_ord (i,j)
        ENDDO 
        D_eff (i) = 1.d0 / D_eff (i)
        D_eff (i) = p * D_eff (i) / (RUNIV * T_ord (i))
     ENDDO

! 3. Calculate mixture molar mass :

     MM = MIXTURE_MOLARMASS(x, 1, 1)

! 4. Calculate initial number diffusion terms :

     N_tmp = 0.d0
     DO i = 1, nsp
        N_tmp (i) = J_old (i) / MMOL (i)
     ENDDO
     N_old = 0.d0
     CALL Order_data (nsp,index_array,N_tmp,N_old) 

     N_new = N_old

!--------------------------------------------------------------
! Sutton's iterative scheme (use only traco-order herein) :
! Estimate initial residual (starting guess: zero massfluxes, 
! estimated driving forces) :

  d_ref = 1.d0 / l_ref
  res0 = 0.d0
  DO i = 1, nsp
     res0 = res0 + D_eff (i)**2 * d_ref**2
  ENDDO
  res0 = sqrt (res0)
!--------------------------------------------------------------
  
!--------------------------------------------------------------
  DO k = 1, n_it

     resit = 0.d0; resmass = 0.d0; rescharge = 0.d0

!--------------------------------------------------------------
     ! A. Calculate residual except ambipolar driving force
     DO i = 1, nsp
        sum = 0.d0
        DO j = 1, nsp
           sum = sum + N_old (j) * RUNIV * T (j) / (p * diff_ord (i, j))   
        ENDDO
        res (i) = - N_new(i) - D_eff(i) * d_ord(i) + &
                  & z_ord(i) * D_eff(i) * sum
     ENDDO
!--------------------------------------------------------------

!--------------------------------------------------------------
     ! B. Calculate ambipolar field such as to ensure zero current
     !    for updated number fluxes
     IF (nions.GT.0) THEN
      denom = 0.d0
      numer = 0.d0
      DO i = 1, nsp
         denom=denom - charge_ord (i)*(res(i)+N_new(i))
         numer=numer + D_eff(i)*(x_ord(i)+1.d-100)*charge_ord(i)**2
      ENDDO
      d_amb = denom / numer
     ELSE 
      d_amb = 0.d0
     ENDIF
!--------------------------------------------------------------
     
!--------------------------------------------------------------
     ! C. Update number fluxes for above ambipolar field
     DO i = 1, nsp
        res (i) = res (i)+D_eff(i)*x_ord(i)*charge_ord(i)*d_amb
        resit = resit + res (i) ** 2
        N_new (i) = N_new (i) + res (i)
     ENDDO
     resit = sqrt (resit) / res0
!--------------------------------------------------------------

!--------------------------------------------------------------
     ! D. Ramshaw's correction to number fluxes
     sum = 0.d0
     DO i = 1, nsp
        sum = sum + MMOL_ord (i) * N_new (i)
     ENDDO
     ! sum now stands for non-zero net mass flux
     DO i = 1, nsp
        N_new (i) = N_new (i) - x_ord(i) / MM * sum
     ENDDO
     resmass = abs (sum) / (MM * res0) 
!--------------------------------------------------------------

!--------------------------------------------------------------
     ! E. Calculate charge residual
     sum = 0.d0
     DO i = 1, nsp
        sum = sum + charge_ord (i) * N_new (i)
     ENDDO
     rescharge = abs (sum) / res0
!--------------------------------------------------------------

  N_old = N_new
  if (resit.le.tresh.and.resmass.le.tresh.and.rescharge.le.tresh) exit

  ENDDO

! Reorder stuff :

     N_tmp = N_new
     CALL Reorder_data (nsp,index_array,N_tmp,N_new) 

! 5. Compute diffusive mass fluxes

     J_new = 0.
     DO i = 1, nsp
        J_new (i) = MMOL (i) * N_new (i)
     ENDDO

! 6. Compute ambipolar electric field
    
     E_amb = 0.d0
     DO i = 1, nsp 
        E_amb = E_amb + d_amb * x_ord (i) * KUNIV * T_ord (i) / EUNIV
     ENDDO

!***************************************************************
END SUBROUTINE pDvda_stefmax
