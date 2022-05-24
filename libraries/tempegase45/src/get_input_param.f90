! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  SETUP   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE get_input_param                                      //
! //                                                                  //
! //  Input: none.                                                    //    
! //                                                                  //
! //  Outs:   in_dir      working directory for input files           // 
! //          out_dir     working directory for output files          //
! //          thermo_dir  working directory for thermodynamic data    //
! //          traco_dir   working directory for transport data        //
! //          chemo_dir   working directory for chemical kinetics data//
! //                                                                  //
! //  The routine is a simple routine opening, reading and closing    //
! //  the file PEGASE.CTL located in the same directory as the        //
! //  executable. It retrieves the relevcant directory structure.     //
! //                                                                  //
! //  Benoit Bottin, 03/10/96                                         //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE get_input_param (in_dir,out_dir,thermo_dir,traco_dir,&
     &                       chemo_dir)
!      
      USE global_pegase
      IMPLICIT NONE
      INTEGER iex,FILE_EXISTS,FILE_OPEN,FILE_CLOSE
      CHARACTER*80 in_dir,out_dir,thermo_dir,traco_dir,chemo_dir,main,dummy
      CHARACTER*12 ctlname
      CHARACTER*6 pgname
      CHARACTER*70 localinfo
      CHARACTER*4 ostype
!            
      pgname='pegase'
      localinfo(1:1)=' '
!
      ctlname='pegase.ctl'
      iex=FILE_EXISTS(ctlname)
      IF (iex.EQ.0) THEN
          CALL PRINT_ERROR (1,pgname,localinfo,1)
      ELSEIF (iex.EQ.-1) THEN
          CALL PRINT_ERROR (2,pgname,localinfo,1)
	ENDIF
!
      iex=FILE_OPEN(ctlname,31,1)
      IF (iex.EQ.-1) THEN
          CALL PRINT_ERROR (3,pgname,localinfo,1)
      ENDIF
!
      READ (31,1000,ERR=901,END=902) ostype
      IF (ostype.EQ.'UNIX') THEN
          flg_os=0
          READ (31,1001,ERR=901,END=902) in_dir
          READ (31,1001,ERR=901,END=902) out_dir
          READ (31,1001,ERR=901,END=902) thermo_dir
          READ (31,1001,ERR=901,END=902) traco_dir
          READ (31,1001,ERR=901,END=902) chemo_dir
      ELSE
          flg_os=1
          READ (31,1001,ERR=901,END=902) main
          READ (31,1001,ERR=901,END=902) dummy
          CALL CONCATENATE(main,dummy,in_dir,iex)
          READ (31,1001,ERR=901,END=902) dummy
          CALL CONCATENATE(main,dummy,out_dir,iex)
          READ (31,1001,ERR=901,END=902) dummy
          CALL CONCATENATE(main,dummy,thermo_dir,iex)
          READ (31,1001,ERR=901,END=902) dummy
          CALL CONCATENATE(main,dummy,traco_dir,iex)
          READ (31,1001,ERR=901,END=902) dummy
          CALL CONCATENATE(main,dummy,chemo_dir,iex)
      ENDIF
      iex=FILE_CLOSE(31)
      IF (iex.EQ.-1) THEN
	    CALL PRINT_ERROR (6,pgname,localinfo,-1)
      ENDIF
      RETURN
!     ==================================================================
!     Error occurred during READ with flag ERR
!     ----------------------------------------
901   CALL PRINT_ERROR (4,pgname,localinfo,1)
902   CALL PRINT_ERROR (5,pgname,localinfo,1)
!     ==================================================================
!     Formats
!     -------
1000  FORMAT (a4)
1001  FORMAT (a80)      
!     ==================================================================
      END
! //////////////////////////////////////////////////////////////////////
