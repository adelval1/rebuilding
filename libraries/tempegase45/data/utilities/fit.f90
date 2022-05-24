program fit

!----------------------------------------------------------
! Utility which reads collision integral data from standard
! input and fits them to the format used in PEGASE.
! Do not forget to link to the NAG library !
!----------------------------------------------------------
implicit none

integer:: i, npoints
real*8:: pi, tcheck, omcheck
real*8, allocatable:: T(:), omega(:)

! Stuff used when NAG routines are called
real*8:: ref, A(4)
!----------------------------------------------------------

!----------------------------------------------------------
! Initialize:
!----------------------------------------------------------

pi = 3.1415d0

read *,npoints
print *,'Number of data points:', npoints
print *,'                      '
allocate(T(npoints), omega(npoints))
T = 0.0d0; omega = 0.d0

!----------------------------------------------------------
! Read collision integral data from standard input:
!----------------------------------------------------------

do i = 1, npoints
print *,'Reading data point',i
read *, T(i), omega(i)
enddo

! Transform to logarithmic form:

do i = 1, npoints
T(i) = dlog(T(i))
omega(i) = dlog(omega(i))
enddo

!----------------------------------------------------------
! Calculate curve-fits
!----------------------------------------------------------

    A = 0.d0
    CALL E02ACF(T, omega, npoints, A, 4, ref)

!----------------------------------------------------------
! Check fitted integral
!----------------------------------------------------------

print *,'Check of fitted result'
do i = 1, npoints
tcheck = T(i)
omcheck = a(4)*tcheck**3+a(3)*tcheck**2+a(2)*tcheck**1+a(1)
print *,dexp(T(i)),dexp(omcheck),dexp(omega(i))
enddo
    print *,'    '
    print *,A(4),A(3),A(2),A(1)

!----------------------------------------------------------
end program fit
