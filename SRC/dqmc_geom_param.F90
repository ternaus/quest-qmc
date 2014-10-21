module DQMC_GEOM_PARAM

implicit none
integer   ,   parameter  :: rdim=3
real*8    ,   parameter  :: toll=1.d-6, pi=acos(-1.d0)
complex*16,   parameter  :: im=(0.d0,1.d0)

integer,      parameter  :: N_Fields  = 10
integer,      parameter  :: NDIM_F    =  1
integer,      parameter  :: PRIM_F    =  2
integer,      parameter  :: SUPER_F   =  3
integer,      parameter  :: ORB_F     =  4
integer,      parameter  :: HAMILT_F  =  5
integer,      parameter  :: SYMM_F    =  6
integer,      parameter  :: PHASE_F   =  7
integer,      parameter  :: BONDS_F   =  8
integer,      parameter  :: PAIRS_F   =  9
integer,      parameter  :: DILUT_F   = 10

character*10, parameter  :: INPUT_FIELDS(N_Fields) = (/&
                              '  #NDIM   ',            &
                              '  #PRIM   ',            &
                              '  #SUPER  ',            &
                              '   #ORB   ',            &
                              ' #HAMILT  ',            &
                              '  #SYMM   ',            &
                              '  #PHASE  ',            &
                              '  #BONDS  ',            &
                              '   #PAIR  ',            &
                              '  #DILUT  '              /)

integer                  :: inpunit               !Assigned in Geom_Fill
logical                  :: Found_Field(N_fields) !Assigned in analyze_input
save

contains


  !----------------------------------------------
  ! Determines which Fields are specified
  !----------------------------------------------
  subroutine analyze_input
   integer            :: i,ios
   logical            :: stop_exe
   character(len=100) :: str
   do i=1,N_Fields
    rewind(inpunit)
    Found_Field(i)=.false.
    do
      read(inpunit,'(A)',iostat=ios)str
      if(ios/=0)exit
      if(index(str,adjustl(INPUT_FIELDS(i)))>0)Found_Field(i)=.true.
    enddo
   enddo
   stop_exe=.false.
   do i=NDIM_F,HAMILT_F
    if(.not.Found_Field(i))then
     write(*,'(A)')INPUT_FIELDS(i),' is compulsory in input'
     stop_exe=.true.
    endif
   enddo
   if(Found_Field(PAIRS_F).and..not.Found_Field(BONDS_F))then
    write(*,'(A)') '#PAIR requires #BONDS to be specified in input'
    stop_exe=.true.
   endif
   if(stop_exe)stop
  end subroutine analyze_input

end module DQMC_GEOM_PARAM
