subroutine RICHARDSON_EXTRAPOLATION

USE REBUILDING_var
USE ICP_par
USE EXPERIMENTS_par
USE NAMES_par

implicit none

real*8 catalycity_1(1:N_gamma),catalycity_2(1:N_gamma),catalycity_3(1:N_gamma)
real*8 enth_1(1:N_gamma),enth_2(1:N_gamma),enth_3(1:N_gamma)
real*8 V_s_1(1:N_gamma),V_s_2(1:N_gamma),V_s_3(1:N_gamma)
real*8 beta_1(1:N_gamma),beta_2(1:N_gamma),beta_3(1:N_gamma)
real*8 dummy
real*8 extrapolated_enth(1:N_gamma)
integer i
character*80 folder,filename,filename2,filename3

folder = 'cerbereoutput/'

filename = trim(folder)//trim(probe_name1)
OPEN(20,FILE=filename,STATUS='OLD',ACCESS='SEQUENTIAL',&
& ACTION='READ')
   do i = 1,N_gamma
     read(20,*) catalycity_1(i),dummy,enth_1(i),dummy,dummy,dummy,dummy,dummy,V_s_1(i),dummy
     beta_1(i) = (u1e_1*V_s_1(i)/Rm)**1 
   enddo
close(20)

filename2 = trim(folder)//trim(probe_name2)
OPEN(21,FILE=filename2,STATUS='OLD',ACCESS='SEQUENTIAL',&
& ACTION='READ')
   do i = 1,N_gamma
     read(21,*) catalycity_2(i),dummy,enth_2(i),dummy,dummy,dummy,dummy,dummy,V_s_2(i),dummy
     beta_2(i) = (u1e_2*V_s_2(i)/Rm_2)**1 
   enddo
close(21)

filename3 = trim(folder)//trim(probe_name3)
OPEN(22,FILE=filename3,STATUS='OLD',ACCESS='SEQUENTIAL',&
& ACTION='READ')
   do i = 1,N_gamma
     read(22,*) catalycity_3(i),dummy,enth_3(i),dummy,dummy,dummy,dummy,dummy,V_s_3(i),dummy
     beta_3(i) = (u1e_3*V_s_3(i)/Rm_3)**1 
   enddo
close(22)

do i = 1,N_gamma
   extrapolated_enth(i) = beta_2(i)*beta_1(i) / ( (beta_3(i)-beta_1(i)) * (beta_3(i)-beta_2(i)) ) * enth_3(i)
   extrapolated_enth(i) = extrapolated_enth(i) - beta_3(i)*beta_1(i)/((beta_3(i)-beta_2(i))*(beta_2(i)-beta_1(i)))*enth_2(i)
   extrapolated_enth(i) = extrapolated_enth(i) + beta_3(i)*beta_2(i)/((beta_3(i)-beta_1(i))*(beta_2(i)-beta_1(i)))*enth_1(i)
enddo

OPEN(96,FILE='cerbereoutput/extrapolated_enthalpy.out',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')
   do i = 1,N_gamma
      write(96,*) catalycity_1(i),extrapolated_enth(i)
   enddo
CLOSE(96)

end subroutine RICHARDSON_EXTRAPOLATION
