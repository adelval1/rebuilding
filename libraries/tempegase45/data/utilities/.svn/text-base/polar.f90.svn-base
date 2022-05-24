program polar

!----------------------------------------------------------
! Utility which generates collision integrals for the
! polarizability-potential and which fits them into the
! PEGASE-format. Do not forget to link to the NAG library !
!----------------------------------------------------------

implicit none

! Intermediate data used in temperature loops
real*8:: omega11, omega12, omega13, omega22, bstar
real*8, allocatable:: om11(:), om22(:), ombs(:), x(:)
! Stuff used when NAG routines are called
real*8:: T, alpha, ref
real*8:: A11(4), A22(4), Abs(4), omega(3,4)

integer:: k,l, npoints

character*11:: interaction

!----------------------------------------------------------
! Open some files
!----------------------------------------------------------

    OPEN (UNIT=1, FILE='polar.dat', STATUS='unknown')
    OPEN (UNIT=2, FILE='om_11.pol', STATUS='unknown')
    OPEN (UNIT=3, FILE='om_22.pol', STATUS='unknown')
    OPEN (UNIT=4, FILE='om_bs.pol', STATUS='unknown')

!----------------------------------------------------------
! Initialize 
!----------------------------------------------------------

    npoints = 60 ! Number of data points of the collision
                 ! integral considered here
    omega = 0.d0
    allocate(x(npoints), om11(npoints), om22(npoints), ombs(npoints))
    x = 0.0d0; om11 = 0.d0; om22 = 0.d0; ombs = 0.d0

!----------------------------------------------------------
! Ask polarizability:
!----------------------------------------------------------

    print *,'Indicate interaction, e.g. N2-'
    read *,interaction
    print *,'Give value of polarizability'
    read *,alpha

!----------------------------------------------------------
! Temperature loop
!----------------------------------------------------------

    DO l = 1, npoints

    T = 300.d0 + (l-1)*200.d0

    omega11 = 425.4d0 * alpha**.5 / T**.5
    omega12 = .8333d0 * omega11
    omega13 = .7292d0 * omega11
    omega22 = .8710d0 * omega11
    bstar = (5.d0*omega12 - 4.d0*omega13) / omega11

       x(l) = dlog(T)
    om11(l) = dlog(omega11)  
    om22(l) = dlog(omega22)  
    ombs(l) = dlog(bstar) 

    ENDDO

!----------------------------------------------------------
! Calculate curve-fits
!----------------------------------------------------------

!   Omega11-fit:

    A11 = 0.d0
    CALL E02ACF(x, om11, npoints, A11, 4, ref)
    omega(1,1)=A11(4)
    omega(1,2)=A11(3)
    omega(1,3)=A11(2)
    omega(1,4)=A11(1)

!   Omega22-fit:

    A22 = 0.d0
    CALL E02ACF(x, om22, npoints, A22, 4, ref)
    omega(2,1)=A22(4)
    omega(2,2)=A22(3)
    omega(2,3)=A22(2)
    omega(2,4)=A22(1)

!   Bstar-fit:

    Abs = 0.d0
    CALL E02ACF(x, ombs, npoints, Abs, 4, ref)
    omega(3,1)=Abs(4)
    omega(3,2)=Abs(3)
    omega(3,3)=Abs(2)
    omega(3,4)=Abs(1)

!----------------------------------------------------------
! Write some files and close then
!----------------------------------------------------------

    DO k = 1, 3
       WRITE (1,1010) interaction, omega(k,1), &
       & omega(k,2), omega(k,3), omega(k,4)
    ENDDO

    DO l = 1, npoints

    T = 300.d0 + (l-1)*200.d0

    omega11 = omega(1,1)*dlog(T)**3 + &
  & omega(1,2)*dlog(T)**2 + omega(1,3)*dlog(T) + &
  & omega(1,4)

    omega22 = omega(2,1)*dlog(T)**3 + &
  & omega(2,2)*dlog(T)**2 + omega(2,3)*dlog(T) + &
  & omega(2,4)

    bstar  = omega(3,1)*dlog(T)**3 + &
  & omega(3,2)*dlog(T)**2 + omega(3,3)*dlog(T) + &
  & omega(3,4)

    write (2,*) T, dexp(omega11)
    Write (3,*) T, dexp(omega22)
    Write (4,*) T, dexp(bstar)

    ENDDO

    close(1)
    close(2)
    close(3)
    close(4)

!----------------------------------------------------------
1010 FORMAT (a11,1PE15.5,1PE15.5,1PE15.5,1PE15.5)
!----------------------------------------------------------

end program polar
