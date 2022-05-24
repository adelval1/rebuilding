REAL(kind=8) FUNCTION Brokaw_lambda_reactive (T)

!***********************************************************
! Reactive thermal conductivity (i.e. diffusion of
! reaction enthalpies in an equilibrium mixture) 
! computed by the rigorous formula of Butler and Brokaw.
!***********************************************************
      USE global_thermo
      USE global_traco
      USE traco_ord_var, only: x_ord
      IMPLICIT NONE
!***********************************************************

      INTEGER :: i, j, k, l
      REAL(kind=8) :: T, lambda, eps
      REAL(kind=8) :: D (nr, nr), W (nr)

!***********************************************************

eps = 1.d-16

! Compute upper part of matrix.

DO i = 1, nr
DO j = i, nr
 D (i,j) = 0.d0
 DO k=1, nsp - 1
    DO l = k+1, nsp
       D (i,j) = D (i,j) + delta(1, k, l) * (x_ord (k) + eps) * (x_ord (l) + eps) * &
     & ( nu_ord (i, k)/(x_ord (k) + eps) - nu_ord (i, l)/(x_ord (l) + eps) ) * &
     & ( nu_ord (j, k)/(x_ord (k) + eps) - nu_ord (j, l)/(x_ord (l) + eps) )   
    ENDDO 
 ENDDO 
 ! Mirror matrix accros diagonal.
 D (j,i) = D (i, j)
ENDDO
ENDDO

D = D / KUNIV

W = 0.d0
CALL invert(nr, D, HREACT/RUNIV**2/T**2, W)

lambda = 0.d0   
DO l = 1, nr
    lambda = lambda + W(l) * HREACT(l)
END DO
Brokaw_lambda_reactive = lambda

!***********************************************************
END FUNCTION Brokaw_lambda_reactive
!
SUBROUTINE Invert (n,A,B,D)

!***********************************************************
! Simple linear solver based upon Gaussian elimination.
! No tricks whatsoever.
!***********************************************************
 implicit none
 integer::i,j,n
 real(kind=8) :: A(1:n,1:n),B(1:n),C(1:n,1:n+1),D(1:n)
!***********************************************************

 C(:,1:n)=A
 C(:,n+1)=B

 do j=1,n
    C(j,:)=C(j,:)/C(j,j)
    do i=1,n
       if (i.ne.j) then
       C(i,:)=C(i,:)-C(i,j)*C(j,:)
       endif
    enddo
 enddo

 D=C(:,n+1)

!***********************************************************
end subroutine Invert
