!******************************************************************************
!******************************************************************************

! Utility that writes new thermo/*.mix and thermo/*.trc
! files. Starting from a parent data-set (e.g. air13)
! smaller data-sets may be generates (e.g. air5, oxygen5,
! argon3 ...).

! Sorry - this utility isn't too well programmed.
! Quite a lot of 'spaghetti-loops' etc ...
! Perhaps something better will be programmed
! the future.  
! David VDA

!******************************************************************************

 program selspec

   implicit none
   integer, parameter:: maxsp = 13
   character*80 thermo_inputname
   character*80 thermo_outputname
   character*80 traco_inputname
   character*80 traco_outputname
   character*80 thermodir
   character*80 tracodir
   character*20 species_names(1:maxsp)
   integer::    nsp
   logical      selection(1:maxsp)
   logical      Write_Mix 

   call fill_with_blanks(thermodir)
   call fill_with_blanks(tracodir)
   call fill_with_blanks(thermo_inputname)
   call fill_with_blanks(thermo_outputname)
   call fill_with_blanks(traco_inputname)
   call fill_with_blanks(traco_outputname)

   thermodir='../thermo/'
   tracodir='../traco/'

   call ask_filename(thermodir,1,thermo_inputname,Write_Mix)
   call ask_filename(thermodir,0,thermo_outputname,Write_Mix) 

   call open_files(thermodir,thermo_inputname,thermo_outputname,11,21,Write_Mix)
   call read_file(nsp,species_names,traco_inputname)

   call new_traco_name(thermo_outputname,traco_outputname)
   call open_files(tracodir,traco_inputname,traco_outputname,12,22,.true.)

   call ask_selection(nsp,species_names,selection)

   call write_new_traco_file(nsp,selection,species_names)

   if (Write_Mix) then
     call write_new_mixture_file(nsp,selection,species_names,traco_outputname)  
   end if

   call close_files

 end program selspec
 
!******************************************************************************

 subroutine ask_filename(thermodir,flg_input,filename,Write_Mix)

   implicit none
   character*80 thermodir,filename
   integer:: flg_input, flg_error
   logical   Write_Mix

   character*80 fullname
   character answer
   logical exists
   integer:: length
      
   Write_Mix=.true.
   if (flg_input.eq.1) then
     write (*,3001)
   else
     write (*,3002)
   end if
   read (*,4001) filename

   CALL CONCATENATE(thermodir,filename,fullname,flg_error)
   length=INDEX(fullname,' ')
   if (length.EQ.0) then
     length=LEN(fullname)+1
   end if

   INQUIRE (file=fullname(1:length-1), EXIST=exists)
   if ((.NOT.exists).and.(flg_input.eq.1)) then

     write(*,*)'File does not exist'
     stop 'Program terminated on SELSPEC input error'

   else if (exists.and.(flg_input.eq.0)) then

     write(*,*)'A file with this name already exists and will be overwritten'
     write(*,*)'Are you sure you want to overwrite this file?'
     write(*,*)"  If you respond with 'n' only the new curve-fit-file will"
     write(*,*)'    be generated, you have to make sure that the order in the'
     write(*,*)'    mixture-file is the same.'
     write(*,*)"  If you respond with 'q' no file will be generated."
     write(*,3003)'Do you want to overwrite the mixture file?'
     read(*,4002) answer
     if ((answer.ne.'y').and.(answer.ne.'Y')) then
       Write_Mix=.false.
     end if
     if ((answer.eq.'q').or.(answer.eq.'Q')) then
       STOP "No files generated"
     end if

   end if

3001 format(/,/,'Please give the name of the input mixture file')
3002 format(/,/,'Please give the name of the output mixture file')
3003 format($,a)
4001 format(a80)
4002 format(a1)
9001 format('ioERR 001: The specified file does not exists: ',a80)

 end subroutine ask_filename
 
!******************************************************************************

 subroutine open_files(filesdir,inputname,outputname,unitin,unitout,Write_Output)

   implicit none
   character*80 inputname,outputname
   character*80 filesdir
   integer unitin,unitout, flg_error

   logical Write_Output

   integer length
   character*80 fullname
   
   call fill_with_blanks(fullname)
   CALL CONCATENATE(filesdir,inputname,fullname,flg_error)
   length=INDEX(fullname,' ')
   if (length.EQ.0) then
     length=LEN(fullname)+1
   end if
   open(unit=unitin,file=fullname(1:length-1),access='sequential',status='old')

   call fill_with_blanks(fullname)
   CALL CONCATENATE(filesdir,outputname,fullname,flg_error)
   length=INDEX(fullname,' ')
   if (length.EQ.0) then
     length=LEN(fullname)+1
   end if
   if (Write_Output) then
     open(unit=unitout,file=fullname(1:length-1),access='sequential',status='unknown')
   end if

 end subroutine

!******************************************************************************

 subroutine read_file(nsp,species_names,traconame)
   
   implicit none
   integer, parameter:: maxsp = 13
   integer:: nsp
   character*20 species_names(1:maxsp)
   
   character*80 line
   character*80 traconame
   integer:: length
   integer:: I, flg_error

   call fill_with_blanks(traconame)

30 continue

     call fill_with_blanks(line)
     read(11,*,end=50)line
     call lcase(line(1:4), line(1:4), flg_error)

     if (line(1:4).eq.'nsp ') then
       read(11,*)nsp
       do I=1,nsp
         read(11,*)species_names(I)
       enddo 
     end if

     call lcase(line(1:4), line(1:4), flg_error)
     if (line(1:4).eq.'trac') then
       read(11,*)traconame
     end if

     goto 30

50 continue

 end subroutine read_file
 
!******************************************************************************

 subroutine new_traco_name(thermoname,traconame)

   implicit none
   character*80 thermoname,traconame
   
   integer:: length
   
   call fill_with_blanks(traconame)
   length=INDEX(thermoname,'.')
   traconame=thermoname(1:length)
   traconame(length+1:length+4)='trc '
 end subroutine

!******************************************************************************

 subroutine ask_selection(nsp,species_names,selection)
 
   implicit none
   integer, parameter:: maxsp = 13
   integer::nsp
   character*20 species_names(1:maxsp)
   logical selection(1:maxsp)
   
   character*80 input_string
   integer:: nr
   integer:: pos

60 continue

     call write_selection(nsp,species_names,selection)
     call fill_with_blanks(input_string)
     write(*,4020)
     read(*,3019) input_string
     nr=1

     do while (nr.ne.0)     

       if (index(input_string,'0').eq.1) then
         goto 70
       end if
       read(input_string,3020)nr
       if (nr.ne.0) then
         selection(nr)=.not.selection(nr)
       end if

       pos=index(input_string,',')
       if (pos.eq.0) then
         call fill_with_blanks(input_string)
         nr=0
       else
         input_string(1:80-pos)=input_string(1+pos:80)
       end if

     end do

   goto 60

70 return

3019 format(a80)
3020 format(i)
4020 format("Toggle selection by giving the numbers (if needed, separated by commas)",&
  &/,"   (0 stops and ignores the remaining input)")

 end subroutine ask_selection     

!******************************************************************************

 subroutine write_selection(nsp,species_names,selection)
!   use globalvar
 
   implicit none
   integer, parameter:: maxsp = 13
   integer:: nsp
   character*20 species_names(1:maxsp)
   logical selection(1:maxsp)
   
   integer:: I
   
   write(*,*)'The species in the input file'
   write(*,*)'-----------------------------'
   write(*,*)'    (the selected species are marked with *)'
   do 53,I=1,nsp
     if (selection(I)) then
       write(*,4003)I,species_names(I)
     else
       write(*,4004)I,species_names(I)
     end if
53 continue
4003 format ('* ',i2,'. ',a20)
4004 format ('  ',i2,'. ',a20)
 end subroutine write_selection   

!******************************************************************************

 subroutine write_new_traco_file(nsp,selection,species_names)

   implicit none
   integer, parameter:: maxsp = 13
   integer:: nsp
   logical selection(1:maxsp)
   character*20 species_names(1:maxsp)

   integer:: I,J
   integer:: old_pos,new_pos
   character*80 line
   character*5 check
   logical neutral
   integer:: Delta_nr
 
   call fill_with_blanks(line)

do Delta_nr=1,3
 read(12, 3009) line
 write(22,3009) line

   do I=1,nsp

     if (neutral(species_names(I))) then
       read(12, 3009) line
       if (selection(I)) write(22,3009) line
       do J=1,nsp
         read(12, 3009) line
         if (selection(I)) then
         if (selection(J)) write(22,3009) line
         endif
       enddo
     end if

   enddo

enddo

3009 format(a80)

 end subroutine write_new_traco_file

!******************************************************************************

 subroutine write_new_mixture_file(nsp,selection,species_names,traconame)

   implicit none
   integer, parameter:: maxsp = 13
   integer:: nsp
   logical selection(1:maxsp)
   character*20 species_names(1:maxsp)
   character*80 traconame
   
   character*80 line
   character*80 new_line
   character*4 command
   integer:: new_nsp
   real*4:: real_number
   integer integer_number
   integer:: old_pos,new_pos
   integer:: new_line_pos
   integer:: I,J
   integer:: number_array(1:maxsp)
   logical empty
   integer::flg_error

   rewind(11)
   new_nsp=0

   do I=1,nsp
     if (selection(I)) then
       new_nsp=new_nsp+1
     end if
   end do

100 continue

     read(11,3010)line
     call lcase(line(1:4),  line(1:4), flg_error)
     command=line(1:4)

     if (command.eq.'nsp ') then
       write(21,4012)'NSP'
       write(21,4015)dfloat(new_nsp)
       do I=1,nsp
         if (selection(I)) then
           write(21,4012)species_names(I)
         end if
       end do
       write(21,4011)
     else if (command.eq.'trac') then
       write(21,4012)'TRACO'
       write(21,4012)traconame
       write(21,4011)
     else if (command.eq.'xini') then
       write(21,4012)'XINI'
       do I=1,nsp
         read(11,3010)line
         if (selection(I))then
           write(21,4012)'  You will have to do this yourself'
         end if
       end do
       write(21,4011)
     else if (command.eq.'nr') then
       write(21,4012)'NR'
       write(21,4012)'  You will have to do this yourself'
       write(21,4011)
     else if (command.eq.'nc') then
       write(21,4012)'NC'
       write(21,4012)'  You will have to do this yourself'
       write(21,4011)
     else if (command.eq.'stop')then
       write(21,4012)'STOP'
       goto 110
     end if

   goto 100

110 continue

3010 format(a80)
3011 format(f12.6)
3012 format(13i6)
4011 format('--------------------')
4012 format(a)
4013 format($,' ')
4014 format($,i6)
4015 format(bn,f4.1)
4016 format($,a)

 end subroutine write_new_mixture_file

!******************************************************************************

 logical function neutral(species_name)
 
   implicit none
   character*20 species_name
   
   integer:: length
   logical ret
   
   length=INDEX(species_name,'.')
   neutral=(.not.(species_name(length+1:length+3).eq.'ion'))

 end function neutral
 
!******************************************************************************

 subroutine close_files

   close(unit=11,status='keep')
   close(unit=12,status='keep')
   close(unit=21,status='keep')
   close(unit=22,status='keep')

 end subroutine close_files

!******************************************************************************
