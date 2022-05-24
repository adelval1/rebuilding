program check

!----------------------------------------------------------
! Utility which checks new collision integral file air13.trc  
! Do not forget to link to the NAG library !!!
!----------------------------------------------------------

implicit none

character*80:: text
character*15, allocatable:: interac(:,:)

integer:: nsp, nneut, i, j, k, l, npoints

real*8:: T, pi

! Intermediate data used in temperature loops
real*8:: omega11_old, omega22_old, bstar_old
real*8:: omega11_new, omega12_new, omega13_new, &
       & omega22_new, bstar_new
real*8, allocatable:: om11(:), om22(:), ombs(:), x(:)

! Coefficients of old and new data in PEGASE form
real*8, allocatable:: omega_old(:,:,:,:), omega_new(:,:,:,:)

! Coefficients of new data in Capitelli's form
real*8, allocatable:: omega_cap(:,:,:,:)
real*8 :: a1,a2,a3,a4,a5,a6,a7,a8,a9

! Stuff used when NAG routines are called
real*8:: ref
real*8:: A11(4), A22(4), Abs(4)
!----------------------------------------------------------

!----------------------------------------------------------
! Open some files
!----------------------------------------------------------

OPEN (UNIT=1, FILE='air13.trc', STATUS='unknown')
OPEN (UNIT=2, FILE='air13.cap', STATUS='unknown')

!----------------------------------------------------------
! Initialize ...
!----------------------------------------------------------

pi = 3.1415d0
nsp = 13; nneut = 6
npoints = 60 ! Number of data points of the collision 
            ! integral considered here

allocate(interac(nneut,nsp)); interac = '                 '

allocate(omega_old(3,nneut,nsp,4)); omega_old = 0.d0
allocate(omega_cap(4,nneut,nsp,9)); omega_cap = 0.d0
allocate(omega_new(3,nneut,nsp,4)); omega_new = 0.d0

allocate(x(npoints), om11(npoints), om22(npoints), ombs(npoints))
x = 0.0d0; om11 = 0.d0; om22 = 0.d0; ombs = 0.d0

!----------------------------------------------------------
! Read old collision integral data file air13.trc
!----------------------------------------------------------

DO k = 1, 3
   read (1,*) text
   print *, text
   DO i = 1, nneut
      read (1,*) text
      print *, text
      DO j = 1, nsp
!     Index i in Pegase neutral order, j in Pegase-order
         read (1,*) interac(i,j), omega_old(k,i,j,1), &
       & omega_old(k,i,j,2), omega_old(k,i,j,3), omega_old(k,i,j,4)
         print *,'Reading interaction     ', interac(i,j)
      ENDDO
   ENDDO
ENDDO

omega_new = omega_old

!----------------------------------------------------------
! Read new collision integral data file air13.cap
!----------------------------------------------------------

omega_cap = 0.d0
DO k = 1, 27
   READ (2,*) text
   print *, text
   READ (2,*) i, j
!  Both i and j in Pegase-order
   IF (j.NE.13) THEN
      READ (2,*) omega_cap(1,i,j,1:6)
      READ (2,*) omega_cap(2,i,j,1:6)
      READ (2,*) omega_cap(3,i,j,1:6)
      READ (2,*) omega_cap(4,i,j,1:6)
   ELSE
      READ (2,*) omega_cap(1,i,j,1:9)
      READ (2,*) omega_cap(2,i,j,1:9)
      READ (2,*) omega_cap(3,i,j,1:9)
      READ (2,*) omega_cap(4,i,j,1:9)
   ENDIF
ENDDO

!----------------------------------------------------------
! Calculate new collision integral values between
! 300 and 12300 K:
!----------------------------------------------------------

print *,'CCCCCCCCCCC'
print *,'O2 = 1         O2=1'
print *,'N2 = 2         N2=2'
print *,'NO = 3         NO=3'
print *,'              NO+=4'
print *,'               O+=5'
print *,'               N+=6'
print *,'              Ar+=7'
print *,'              O2+=8'
print *,'              N2+=9'
print *,'O  = 4          O=10'
print *,'N  = 5          N=11'
print *,'Ar = 6         Ar=12'
print *,'                e=13'
print *,'CCCCCCCCCCC'
print *,'Give i and j values of interaction:'
print *,'Neutral order for i, full PEGASE order for j'
read *,i,j

print *,'CCCCCCCCCCC'
print *,'Chosen interaction:    ', interac(i,j)
print *,'Capitelli-fits:',omega_cap(1,i,j,:)
print *,'CCCCCCCCCCC'

DO l = 1, npoints

    T = 300.d0 + (l-1)*200.d0

    IF (j.NE.13) THEN

    omega11_new = pi * (omega_cap(1,i,j,1) + omega_cap(1,i,j,2)*T**omega_cap(1,i,j,3)) / &
                & (omega_cap(1,i,j,4) + omega_cap(1,i,j,5)*T**omega_cap(1,i,j,6))

    omega12_new = pi * (omega_cap(2,i,j,1) + omega_cap(2,i,j,2)*T**omega_cap(2,i,j,3)) / &
                & (omega_cap(2,i,j,4) + omega_cap(2,i,j,5)*T**omega_cap(2,i,j,6))

    omega13_new = pi * (omega_cap(3,i,j,1) + omega_cap(3,i,j,2)*T**omega_cap(3,i,j,3)) / &
                & (omega_cap(3,i,j,4) + omega_cap(3,i,j,5)*T**omega_cap(3,i,j,6))

    omega22_new = pi * (omega_cap(4,i,j,1) + omega_cap(4,i,j,2)*T**omega_cap(4,i,j,3)) / &
                & (omega_cap(4,i,j,4) + omega_cap(4,i,j,5)*T**omega_cap(4,i,j,6))

    bstar_new   = (5.d0*omega12_new - 4.d0*omega13_new) / &
                            & omega11_new

    ELSE
    
    a1 = omega_cap(1,i,j,1); a2 = omega_cap(1,i,j,2); a3 = omega_cap(1,i,j,3)
    a4 = omega_cap(1,i,j,4); a5 = omega_cap(1,i,j,5); a6 = omega_cap(1,i,j,6)
    a7 = omega_cap(1,i,j,7); a8 = omega_cap(1,i,j,8); a9 = omega_cap(1,i,j,9)
    omega11_new = (   a3*(dlog(T)**a6)*dexp((dlog(T)-a1)/a2)   ) &
              & / (   dexp((dlog(T)-a1)/a2) + dexp(-(dlog(T)-a1)/a2)   ) &
              & + a7*dexp(-((dlog(T)-a8)/a9)**2) + a4*(1+(dlog(T))**a5)
    omega11_new = omega11_new * pi

    a1 = omega_cap(2,i,j,1); a2 = omega_cap(2,i,j,2); a3 = omega_cap(2,i,j,3)
    a4 = omega_cap(2,i,j,4); a5 = omega_cap(2,i,j,5); a6 = omega_cap(2,i,j,6)
    a7 = omega_cap(2,i,j,7); a8 = omega_cap(2,i,j,8); a9 = omega_cap(2,i,j,9)
    omega12_new = (   a3*(dlog(T)**a6)*dexp((dlog(T)-a1)/a2)   ) &
              & / (   dexp((dlog(T)-a1)/a2) + dexp(-(dlog(T)-a1)/a2)   ) &
              & + a7*dexp(-((dlog(T)-a8)/a9)**2) + a4*(1+(dlog(T))**a5)
    omega12_new = omega12_new * pi

    a1 = omega_cap(3,i,j,1); a2 = omega_cap(3,i,j,2); a3 = omega_cap(3,i,j,3)
    a4 = omega_cap(3,i,j,4); a5 = omega_cap(3,i,j,5); a6 = omega_cap(3,i,j,6)
    a7 = omega_cap(3,i,j,7); a8 = omega_cap(3,i,j,8); a9 = omega_cap(3,i,j,9)
    omega13_new = (   a3*(dlog(T)**a6)*dexp((dlog(T)-a1)/a2)   ) &
              & / (   dexp((dlog(T)-a1)/a2) + dexp(-(dlog(T)-a1)/a2)   ) &
              & + a7*dexp(-((dlog(T)-a8)/a9)**2) + a4*(1+(dlog(T))**a5)
    omega13_new = omega13_new * pi

    a1 = omega_cap(4,i,j,1); a2 = omega_cap(4,i,j,2); a3 = omega_cap(4,i,j,3)
    a4 = omega_cap(4,i,j,4); a5 = omega_cap(4,i,j,5); a6 = omega_cap(4,i,j,6)
    a7 = omega_cap(4,i,j,7); a8 = omega_cap(4,i,j,8); a9 = omega_cap(4,i,j,9)
    omega22_new = (   a3*(dlog(T))**a6*dexp((dlog(T)-a1)/a2)   ) &
              & / (   dexp((dlog(T)-a1)/a2) + dexp(-(dlog(T)-a1)/a2)   ) &
              & + a7*dexp(-((dlog(T)-a8)/a9)**2) + a4*(1+(dlog(T))**a5)
    omega22_new = omega22_new * pi

    bstar_new   = (5.d0*omega12_new - 4.d0*omega13_new) / &
                            & omega11_new

    ENDIF 

       x(l) = dlog(T)
    om11(l) = dlog(omega11_new)  
    om22(l) = dlog(omega22_new)  
    ombs(l) = dlog(bstar_new) 

ENDDO

!----------------------------------------------------------
! Calculate curve-fits
!----------------------------------------------------------

! Omega11-fit:

A11 = 0.d0
CALL E02ACF(x, om11, npoints, A11, 4, ref)

! Omega22-fit:

A22 = 0.d0
CALL E02ACF(x, om22, npoints, A22, 4, ref)

! Bstar-fit:

Abs = 0.d0
CALL E02ACF(x, ombs, npoints, Abs, 4, ref)

! Fill in new collision integrals:

omega_new(1,i,j,1)=A11(4)
omega_new(1,i,j,2)=A11(3)
omega_new(1,i,j,3)=A11(2)
omega_new(1,i,j,4)=A11(1)

omega_new(2,i,j,1)=A22(4)
omega_new(2,i,j,2)=A22(3)
omega_new(2,i,j,3)=A22(2)
omega_new(2,i,j,4)=A22(1)

omega_new(3,i,j,1)=Abs(4)
omega_new(3,i,j,2)=Abs(3)
omega_new(3,i,j,3)=Abs(2)
omega_new(3,i,j,4)=Abs(1)

!----------------------------------------------------------
! Perform comparison with old collision integral 
! values:
!----------------------------------------------------------

print *,'writing down comparison files'

OPEN (UNIT=3, FILE='om_11.dat', STATUS='unknown')
OPEN (UNIT=4, FILE='om_22.dat', STATUS='unknown')
OPEN (UNIT=8, FILE='om_bs.dat', STATUS='unknown')

DO l = 1, 90

    T = 300.d0 + (l-1)*200.d0

    !------------------------------------------------------
    ! Calculate old and new collision integrals

    omega11_old = omega_old(1,i,j,1)*dlog(T)**3 + &
  & omega_old(1,i,j,2)*dlog(T)**2 + omega_old(1,i,j,3)*dlog(T) + &
  & omega_old(1,i,j,4)
    omega11_new = omega_new(1,i,j,1)*dlog(T)**3 + &
  & omega_new(1,i,j,2)*dlog(T)**2 + omega_new(1,i,j,3)*dlog(T) + &
  & omega_new(1,i,j,4)

    omega22_old = omega_old(2,i,j,1)*dlog(T)**3 + &
  & omega_old(2,i,j,2)*dlog(T)**2 + omega_old(2,i,j,3)*dlog(T) + &
  & omega_old(2,i,j,4)
    omega22_new = omega_new(2,i,j,1)*dlog(T)**3 + &
  & omega_new(2,i,j,2)*dlog(T)**2 + omega_new(2,i,j,3)*dlog(T) + &
  & omega_new(2,i,j,4)

    bstar_old  = omega_old(3,i,j,1)*dlog(T)**3 + &
  & omega_old(3,i,j,2)*dlog(T)**2 + omega_old(3,i,j,3)*dlog(T) + &
  & omega_old(3,i,j,4)
    bstar_new  = omega_new(3,i,j,1)*dlog(T)**3 + &
  & omega_new(3,i,j,2)*dlog(T)**2 + omega_new(3,i,j,3)*dlog(T) + &
  & omega_new(3,i,j,4)
    !------------------------------------------------------

    !------------------------------------------------------
    ! Write stuff:

    Write (3,*) T, dexp(omega11_old), dexp(omega11_new)
    Write (4,*) T, dexp(omega22_old), dexp(omega22_new)
    Write (8,*) T, dexp(bstar_old)  , dexp(bstar_new)
    
    !------------------------------------------------------

ENDDO

1000 FORMAT (f14.6,2x,f14.6,2x,f14.6)
1010 FORMAT (a15,f14.6,4x,f14.6,4x,f14.6,4x,f14.6)

CLOSE (3)
CLOSE (4)
CLOSE (8)

!----------------------------------------------------------
! Close files
!----------------------------------------------------------

CLOSE (1)
CLOSE (2)

!----------------------------------------------------------
end program check
