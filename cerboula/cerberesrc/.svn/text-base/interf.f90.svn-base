subroutine INTERF

USE OPTIONS_par
USE REBUILDING_var

implicit none

10 write(*,*) '                                                                        '
write(*,*) ' ______________________________________________________________________ '
write(*,*) '| ____________________________________________________________________ |'
write(*,*) '||                                                                    ||'
write(*,*) '||                                                                    ||'
write(*,*) '||    //////// ///////  ///////   ///////   /////// ///////  ///////  ||'
write(*,*) '||   //     // //       //    //  //     // //      //    // //       ||'
write(*,*) '||  //         //       //    //  //     // //      //    // //       ||'
write(*,*) '|| //          //////   ///////   ///////   //////  ///////  //////   ||'
write(*,*) '|| //          //       //  //    //    //  //      //  //   //       ||'
write(*,*) '|| //     //   //       //   //   //     // //      //   //  //       ||'
write(*,*) '|| ////////    ///////  //    //  ////////  /////// //    // ///////  ||'
write(*,*) '||                                                                    ||'      
write(*,*) '||      Catalycity and Enthalpy ReBuilding for a REference probe      ||'      
write(*,*) '||                           Version 2.2                              ||'      
write(*,*) '||                                                                    ||'      
write(*,*) '||                   A. Garcia Munoz, P. Barbante,                    ||'      
write(*,*) '||            H.W. Krassilchikoff, J. Thoemel, F. Panerai             ||'
write(*,*) '||____________________________________________________________________||'
write(*,*) '|______________________________________________________________________|'
write(*,*) '                                                                        '
write(*,*) '                                                                        '
write(*,*) 'Enthalpy rebuilding for a fixed gamma ?--------------------------type 0'
write(*,*) 'Enthalpy rebuilding S-curve (enthalpy vs gamma)?-----------------type 1' 
write(*,*) 'Sample catalycity determination ? -------------------------------type 2' 
write(*,*) 'Abacus plotting ? -----------------------------------------------type 3'
write(*,*) 'Initial distribution ? ------------------------------------------type 4'
write(*,*) 'Richardson extrapolation from 3 S-curves? -----------------------type 5'
write(*,*) 'dW/drho computation? --------------------------------------------type 6'
read(*,*)  choice

if (choice.gt.6) then
   write(*,*) ''
   write(*,*) 'Read better and try again!'
   write(*,*) ''
   pause
   goto 10
endif


if (choice.lt.3) then
   30 write(*,*) '...including Kp rebuilding in the loop? -------------------(yes=1/no=0)'
   read(*,*) Kp_choice
   if (Kp_choice.gt.1) then
      write(*,*) ''
      write(*,*) 'It is so difficult to type 0 or 1??'
      write(*,*) ''
      goto 30
   endif

   write(*,*) 'FIRST GUESS OF THE TEMPERATURE AT THE OUTEREDGE OF THE BL?'
   read(*,*) T_reb
endif

end subroutine INTERF
