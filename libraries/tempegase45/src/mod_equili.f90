MODULE global_equilibrium
!
!   Old common block /jacobian/ for equilibrium calculation
!   and additional arrays once passed within the routines.
!   =======================================================
!   jac11       upper left Jacobian matrix              nr,nr
!   jac12       upper right Jacobian matrix             nr,nc+1
!   jaclow      lower Jacobian matrix                   nc+1,nsp+1
!   jacmn       Schur's complement of the Jacobian      nc+1,nc+1
!   f1          independent term of other rows          nr
!   f2          independent term of jacmn rows          nc+1
!   indsd       RHS term corresponding to upper part    nr
!   indmn       RHS term corresponding to lower part    nc+1
!   delv        delta of the variables in Newton        nsp+1
!   luindx      index array for LU matrix operations    nc+1
!   v           local array of massfrac/molefrac logs   nsp+1
!   truevar     local array of massfrac/molefrac vars   nsp+1
!   sys_keq     equilibrium constants array             nr
!   --------------------------------------------------------------
    REAL(kind=8),ALLOCATABLE,SAVE :: jac11(:,:),jac12(:,:)
    REAL(kind=8),ALLOCATABLE,SAVE :: jaclow(:,:),jacmn(:,:)
    REAL(kind=8),ALLOCATABLE,SAVE :: f1(:),f2(:),indsd(:),indmn(:),delv(:)
    REAL(kind=8),ALLOCATABLE,SAVE :: v(:),truevar(:),sys_keq(:)
    INTEGER,ALLOCATABLE,SAVE :: luindx(:)
!
END MODULE global_equilibrium
