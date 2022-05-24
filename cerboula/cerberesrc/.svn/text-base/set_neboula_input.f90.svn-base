subroutine SET_NEBOULA_INPUT

USE EXPERIMENTS_par
USE OUTEREDGE_var
USE CONVERGENCE_ENTH_par
USE REBUILDING_var
USE ICP_par
USE OPTIONS_par
USE INIT_DIST_par

implicit none

real * 8 pe,Vs,Te,mue,rhoe
real * 8,allocatable:: concent(:)

allocate(concent(1:nspecies))

call EQUILCONDITION(he,Pstatic,Te,mue,concent(1:nspecies),rhoe)

Ve_ext = Vs_ext*NDP_5_1

call NUMERICAL_PARAMETER(Te,rhoe)
call BC_WALL(Tw)
call BC_OUT(rhoe,mue,he,u1e_1,u1y_1,Ve_ext,delta_1,ve_1,concent(1:nspecies),Rm)

T_ext       = Te
rho_ext     = rhoe
mu_ext      = mue
h_ext       = he
concent_ext = concent

deallocate(concent)

end subroutine SET_NEBOULA_INPUT
