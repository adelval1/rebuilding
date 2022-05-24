subroutine INIT_DIST

USE EXPERIMENTS_par
USE OPTIONS_par
USE ICP_par
USE CONVERGENCE_ENTH_par
USE OUTEREDGE_var
USE REBUILDING_var
USE INIT_DIST_par

implicit none

real * 8 Vs,Te,mue,rhoe,dummy
real * 8,allocatable:: concent(:)

allocate(concent(1:nspecies))

call EQUILCONDITION(he,Pstatic,Te,mue,concent(1:nspecies),rhoe)
call INITIALDISTRIB(Te)


write(*,*) 'End of the routine Initial Distribution'

end subroutine INIT_DIST
