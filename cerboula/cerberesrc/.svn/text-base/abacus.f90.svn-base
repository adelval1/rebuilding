subroutine ABACUS

USE ICP_par
USE OUTEREDGE_var
USE ABACUS_par

REAL * 8 gamma_min,gamma_max,gam,Qw
real * 8 Twall
INTEGER igamma,itw,iicount,III,jj,kk
character*27 outputfile
character*2 ext
character*80 plotfile,dummytxt
real*8 temp_wall,enth_wall


iicount=0
do igamma=1,ngamma
   gam=gammaarray(igamma)
   iicount=iicount+1

   outputfile(1:21) = 'cerbereoutput/out_neq'
   if(igamma.le.9) then
      write(ext,'(I1)') igamma
      ext='0'//ext
   elseif(igamma.ge.10) then  
      write(ext,'(I2)') igamma
   endif
   outputfile(22:23)=ext
   outputfile(24:27)='.dat'

   OPEN(90,FILE=outputfile,STATUS='UNKNOWN') 

   III=mod(iicount,2)

   if(III.eq.1) then
      do itw=1,ntw
         Twall=Tw_min+(Tw_max-Tw_min)/float(ntw-1)*float(itw-1)

         call PLOT(Twall,gam,Qw)
         call INITIALCOND

         OPEN(97,FILE='neboulaoutput/stagsol.dat',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',&
             &ACTION='READ')  
            read(97,*) dummytxt
            read(97,*) temp_wall, enth_wall
         CLOSE(97)  

         write(90,'(3F20.10)') Twall,enth_wall,Qw
      enddo
   elseif(III.eq.0) then
      do itw=1,ntw
         Twall=Tw_max-(Tw_max-Tw_min)/float(ntw-1)*float(itw-1)
         call PLOT(Twall,gam,Qw)
         call INITIALCOND
  
         OPEN(97,FILE='neboulaoutput/stagsol.dat',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',&
             &ACTION='READ')  
            read(97,*) dummytxt
            read(97,*) temp_wall, enth_wall
         CLOSE(97)   

         write(90,'(3F20.10)') Twall,enth_wall,Qw
      enddo
   endif

   call INITIALCOND
   close(90)
enddo

deallocate(gammaarray)
deallocate(concent_ext)

END subroutine ABACUS
