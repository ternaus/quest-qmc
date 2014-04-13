module DQMC_Cfg

  use DQMC_util
  implicit none 


  ! 
  ! This module contains subroutines to read input parameters.
  !
  ! Parameters
  ! ==========

  ! string parameters
  character, parameter :: COMMENT = "#"
  character, parameter :: SEPARAT = "="
  character, parameter :: COMMA   = ","

  ! string length
  integer, parameter :: slen = 60
  integer, parameter :: llen = 256
  integer, parameter :: alen = 10    ! array limit

  ! status param
  integer, parameter :: STAT_EOF     = -1
  integer, parameter :: STAT_COMMENT = 0
  integer, parameter :: STAT_NORMAL  = 1

  ! type def
  character(8), parameter :: TYPE_STR(3) =(/"Real*8  ","Integer ","Char(60)"/)

  integer, parameter :: TYPE_REAL    = 1 
  integer, parameter :: TYPE_INTEGER = 2 
  integer, parameter :: TYPE_STRING  = 3 

  ! default parameters
  integer, parameter :: N_Param = 40

  ! name of parameters
  character(len=*), parameter :: PARAM_NAME(N_Param) =  &
       &(/"HSF    ", &    ! methods of how HSF are generated                     
       &  "HSFin  ", &    ! File name of HSF input, if HSF = HSF_FROM_FILE       
       &  "HSFout ", &    ! File name of HSF output                              
       &  "HSFtype", &    ! HSF type, 0 = discrete, 1 = continuous
       &  "L      ", &    ! Number of time slides                                
       &  "U      ", &    ! Parameter for potential energy                       
       &  "accept ", &    ! accept counts        
       &  "bcond  ", &    ! boundary conditions
       &  "debug  ", &    ! flag for debug information output                    
       &  "delta1 ", &    ! parameter for continous HSF
       &  "delta2 ", &    ! parameter for continous HSF
       &  "difflim", &    ! limit of tolerable difference                        
       &  "dmu    ", &    ! perturbation of mu                                   
       &  "dtau   ", &    ! discritize parameter                                 
       &  "errrate", &    ! tolerable error                                      
       &  "fixwrap", &    ! fix nwrap to input value
       &  "gamma  ", &    ! correction of Metropolis ratio                       
       &  "gfile  ", &    ! geometry definition                                  
       &  "mu_dn  ", &    ! parameter for chemical potential                     
       &  "mu_up  ", &    ! parameter for chemical potential                     
       &  "n      ", &    ! number of particles                                  
       &  "nbin   ", &    ! number of statisitical bins                          
       &  "nhist  ", &    ! print history or not
       &  "nitvl  ", &    ! number of interval, for continuous FT integration    
       &  "north  ", &    ! frequence of orthogonalization                       
       &  "npass  ", &    ! number of measurement sweeps                         
       &  "ntry   ", &    ! number of global moves per sweeps                    
       &  "nwarm  ", &    ! number of warmup sweeps                              
       &  "nwrap  ", &    ! frequence of recomputing H                           
       &  "nx     ", &    ! number of sites in x direction                       
       &  "ny     ", &    ! number of sites in y direction                       
       &  "nz     ", &    ! number of sites in y direction                       
       &  "ofile  ", &    ! prefix of output files                               
       &  "reject ", &    ! rejection counts                                     
       &  "seed   ", &    ! random seed                                          
       &  "ssxx   ", &    ! use iterative refinement during sweep
       &  "t_dn   ", &    ! parameter for kinetic energy                         
       &  "t_up   ", &    ! parameter for kinetic energy                         
       &  "tausk  ", &    ! frequence of unequal time measurement                
       &  "tdm    "/)     ! compute time dependent measurement

  ! default values
  character(len=*), parameter :: PARAM_DVAL(N_Param) =  &
       &(/"-1      ", &    ! HSF    
       &  "HSF.in  ", &    ! HSFin  
       &  "HSF.out ", &    ! HSFout 
       &  "0       ", &    ! HSFtype
       &  "12      ", &    ! L      
       &  "0.0     ", &    ! U      
       &  "0       ", &    ! accept 
       &  "0,0,0   ", &    ! bcond
       &  "0       ", &    ! debug  
       &  "1.0     ", &    ! delta1
       &  "1.0     ", &    ! delta2
       &  "0.001   ", &    ! difflim
       &  "0.0     ", &    ! dmu    
       &  "0.125   ", &    ! dtau   
       &  "0.001   ", &    ! errrate
       &  "0       ", &    ! fixwrap
       &  "0.0     ", &    ! gamma  
       &  "geom.def", &    ! gfile  
       &  "0.0     ", &    ! mu_up  
       &  "0.0     ", &    ! mu_dn  
       &  "16      ", &    ! n      
       &  "10      ", &    ! nbin   
       &  "0       ", &    ! nhist  
       &  "4       ", &    ! nitvl  
       &  "12      ", &    ! north  
       &  "5000    ", &    ! npass  
       &  "0       ", &    ! ntry   
       &  "1000    ", &    ! nwarm  
       &  "12      ", &    ! nwrap  
       &  "4       ", &    ! nx     
       &  "4       ", &    ! ny     
       &  "2       ", &    ! nz     
       &  "quest   ", &    ! ofile  
       &  "0       ", &    ! reject 
       &  "0       ", &    ! seed 
       &  "0       ", &    ! ssxx
       &  "1.0     ", &    ! t_up      
       &  "1.0     ", &    ! t_dn      
       &  "10      ", &    ! tausk  
       &  "0       "/)     ! tdm
  
  ! parameter type
  integer, parameter :: PARAM_TYPE(N_Param) = &
       &(/TYPE_INTEGER, &    ! HSF    
       &  TYPE_STRING,  &    ! HSFin  
       &  TYPE_STRING,  &    ! HSFout 
       &  TYPE_INTEGER, &    ! HSFtype
       &  TYPE_INTEGER, &    ! L      
       &  TYPE_REAL,    &    ! U      
       &  TYPE_INTEGER, &    ! accept 
       &  TYPE_REAL   , &    ! bcond
       &  TYPE_INTEGER, &    ! debug  
       &  TYPE_REAL,    &    ! delta1
       &  TYPE_REAL,    &    ! delta2
       &  TYPE_REAL,    &    ! difflim
       &  TYPE_REAL,    &    ! dmu    
       &  TYPE_REAL,    &    ! dtau   
       &  TYPE_REAL,    &    ! errrate
       &  TYPE_INTEGER, &    ! fixwrap
       &  TYPE_REAL,    &    ! gamma  
       &  TYPE_STRING,  &    ! gfile  
       &  TYPE_REAL,    &    ! mu_up     
       &  TYPE_REAL,    &    ! mu_dn     
       &  TYPE_INTEGER, &    ! n      
       &  TYPE_INTEGER, &    ! nbin   
       &  TYPE_INTEGER, &    ! nhist  
       &  TYPE_INTEGER, &    ! nitvl  
       &  TYPE_INTEGER, &    ! north  
       &  TYPE_INTEGER, &    ! npass  
       &  TYPE_INTEGER, &    ! ntry   
       &  TYPE_INTEGER, &    ! nwarm  
       &  TYPE_INTEGER, &    ! nwrap  
       &  TYPE_INTEGER, &    ! nx     
       &  TYPE_INTEGER, &    ! ny     
       &  TYPE_INTEGER, &    ! nz     
       &  TYPE_STRING,  &    ! ofile  
       &  TYPE_INTEGER, &    ! reject 
       &  TYPE_INTEGER, &    ! seed   
       &  TYPE_INTEGER, &    ! ssxx
       &  TYPE_REAL,    &    ! t_up      
       &  TYPE_REAL,    &    ! t_dn      
       &  TYPE_INTEGER, &    ! tausk  
       &  TYPE_INTEGER/)     ! tdm

  ! is array parameter
  logical, parameter :: PARAM_ARRAY(N_Param) = &
       &(/.false.,&           ! HSF    
       &  .false.,&           ! HSFin  
       &  .false.,&           ! HSFout 
       &  .false.,&           ! HSFtype
       &  .false.,&           ! L      
       &  .true. ,&           ! U      
       &  .false.,&           ! accept 
       &  .true. ,&           ! bcond
       &  .false.,&           ! debug  
       &  .false.,&           ! delta1
       &  .false.,&           ! delta2
       &  .false.,&           ! difflim
       &  .false.,&           ! dmu    
       &  .false.,&           ! dtau   
       &  .false.,&           ! errrate
       &  .false.,&           ! fixwrap
       &  .false.,&           ! gamma  
       &  .false.,&           ! gfile  
       &  .true. ,&           ! mu_up     
       &  .true. ,&           ! mu_dn     
       &  .false.,&           ! n      
       &  .false.,&           ! nbin   
       &  .false.,&           ! nhist  
       &  .false.,&           ! nitvl  
       &  .false.,&           ! north  
       &  .false.,&           ! npass  
       &  .false.,&           ! ntry   
       &  .false.,&           ! nwarm  
       &  .false.,&           ! nwrap  
       &  .false.,&           ! nx     
       &  .false.,&           ! ny     
       &  .false.,&           ! nz     
       &  .false.,&           ! ofile  
       &  .false.,&           ! reject 
       &  .false.,&           ! seed   
       &  .false.,&           ! ssxx      
       &  .true. ,&           ! t_up
       &  .true. ,&           ! t_dn
       &  .false.,&           ! tausk  
       &  .false./)           ! tdm

  !
  ! Data Type
  ! =========
  !
  type Param
     integer        :: id           ! hashcode from name
     character(len=slen)  :: pname  ! parameter name
     integer        :: ptype        ! type of parameter
     logical        :: isArray      ! is the param an array?
     logical        :: isSet        ! is the parameter been set?
     character(len=slen)  :: defaultval   ! default value
     type(Param), pointer :: next 

     ! values
     integer        :: ival
     real(wp)       :: rval
     integer, pointer     :: iptr(:) 
     real(wp), pointer    :: rptr(:) 
  end type Param

  type config
     type(Param), pointer :: record(:)    ! head of the linked list
     integer        :: nParam
     logical        :: hasDef = .false.
  end type config

  ! ================
  ! Access functions
  ! ================
  interface CFG_Set
     module procedure DQMC_Config_SetI, DQMC_Config_SetR
     module procedure DQMC_Config_SetS, DQMC_Config_SetPR
     module procedure DQMC_Config_SetPI
  end interface

  interface CFG_Get
     module procedure DQMC_Config_GetI, DQMC_Config_GetR
     module procedure DQMC_Config_GetS, DQMC_Config_GetPR
     module procedure DQMC_Config_GetPI
  end interface


contains

  !---------------------------------------------------------------------!
  subroutine DQMC_ReadLn(str, IPT, status)
    !
    ! Purpose
    ! =======
    !    This subrotine reads in a line from file.
    !    It will get rid of #.
    !
    ! Arguments
    ! =========
    !
    character(len=llen)  :: str
    integer              :: status
    intent(inout)        :: str, status
    integer, intent(in)  :: IPT

    ! ... Local Variables ...
    integer                :: ios, pos

    ! ... Executable ...

    read (unit=IPT, FMT="(a)", iostat=ios)  str
    status = STAT_COMMENT

    ! end of file
    if (ios .ne. 0) then
       status = STAT_EOF
       return
    end if

    ! find comment # and get rid of the tailing part
    pos = scan(str, COMMENT, .false.)
    if (pos .ne. 0) then
       ! find the comment sign
       if (pos .ge. 2) then
          str = str(1:pos-1)
       else
          str = ""
       end if
    end if

    if (len_trim(str) .gt. 0) then
       status = STAT_NORMAL
    end if

  end subroutine DQMC_ReadLn

  !---------------------------------------------------------------------!

  subroutine DQMC_Default_Def(cfg)
    !
    ! Purpose
    ! =======
    !    This subrotine initializes default configuration def.
    !    when the config.def is missing.
    !
    ! Arguments
    ! =========
    !
    type(config), intent(inout)  :: cfg          ! configuration
    
    ! ... Local Variable ...
    type(Param),pointer    :: curr 
    integer :: i

    ! ... Executable ...
    cfg%nParam = N_Param
    allocate(cfg%record(cfg%nParam))

    do i = 1, N_Param
       curr => cfg%record(i)
       curr%id         = i
       curr%pname      = PARAM_NAME(i)
       curr%ptype      = PARAM_TYPE(i)
       curr%isArray    = PARAM_ARRAY(i)
       curr%defaultval = PARAM_DVAL(i)

       if (curr%ptype .eq. TYPE_REAL) then
          read(curr%defaultval,*) curr%rval
       elseif (curr%ptype .eq. TYPE_INTEGER) then
          read(curr%defaultval,*) curr%ival
       end if
       
       curr%isSet      = .false.
       nullify(curr%iptr)
       nullify(curr%rptr)
       nullify(curr%next)
    end do

  end subroutine DQMC_Default_Def

  !---------------------------------------------------------------------!

  subroutine DQMC_Config_Free(cfg)
    !
    ! Purpose
    ! =======
    !    This subrotine initializes default configuration def.
    !    when the config.def is missing.
    !
    ! Arguments
    ! =========
    !
    type(config), intent(inout)  :: cfg          ! configuration
    
    ! ... Local ...
    integer :: i
    
    ! ... Executable ...

    do i = 1, cfg%nParam
       if (associated(cfg%record(i)%rptr)) then
          deallocate(cfg%record(i)%rptr)
       end if
       if (associated(cfg%record(i)%iptr)) then
          deallocate(cfg%record(i)%iptr)
       end if
    end do

    deallocate(cfg%record)

  end subroutine DQMC_Config_Free

  !---------------------------------------------------------------------!

  subroutine DQMC_Read_Def(cfg, IPT)
    !
    ! Purpose
    ! =======
    !    This subrotine reads in parameters from a config file.
    !
    ! Arguments
    ! =========
    !
    type(config), intent(inout)  :: cfg          ! configuration
    integer, intent(in)          :: IPT          ! Input file handle

    ! ... Local Variable ...
    integer                :: stat, i, cnt
    character(len=llen)    :: str
    type(Param),pointer    :: head 
    type(Param),pointer    :: curr 
    type(Param),pointer    :: tmp 
    logical                :: found
    ! ... Executable ...

    ! satinize
    nullify(curr)
    nullify(head)
    cnt  = 0
    stat = STAT_COMMENT

    ! read config def
    ! for fast access, sort records by name
    ! using insertion sort
    do while (stat .ne. STAT_EOF)
       call DQMC_ReadLn(str, IPT, stat)

       ! read in a parameter definition 
       if (stat .eq. STAT_NORMAL) then

          ! allocate space
          if (cnt .eq. 0) then
             allocate(head)
             curr => head
          else
             allocate(curr)
          end if
          nullify(curr%next)

          cnt = cnt + 1
          ! read in [name][type][is array][is critical][default value]
          read(str, *)  curr%pname, curr%ptype,  curr%isArray, curr%defaultval

          if (cnt .gt. 1) then
             ! insertion sort
             ! if curr < head, put it as the first one
             if (LGT(head%pname, curr%pname)) then
                curr%next => head
                head      => curr
             else  ! curr >= head
                ! find a record tmp, which is > curr, but its next is < curr.
                tmp => head
                found = .false.
                do while (.not. found .and. associated(tmp%next))
                   ! tmp > curr
                   if (LGT(tmp%next%pname, curr%pname)) then
                      curr%next => tmp%next
                      tmp%next  => curr
                      found = .true.
                   else
                      tmp => tmp%next
                   end if
                end do
                ! curr is the largest
                if (.not. found) then
                   tmp%next => curr
                end if
             end if
          end if
       end if
    end do

    ! allocate space for records
    cfg%nParam = cnt
    allocate(cfg%record(cnt))

    tmp  => head
    do i = 1, cnt
       curr => cfg%record(i)
       curr%id         = i
       curr%pname      = tmp%pname
       curr%ptype      = tmp%ptype
       curr%isArray    = tmp%isArray
       curr%defaultval = tmp%defaultval
       if (curr%ptype .eq. TYPE_REAL) then
          read(curr%defaultval,*) curr%rval
       elseif (curr%ptype .eq. TYPE_INTEGER) then
          read(curr%defaultval,*) curr%ival
       end if

       curr%isSet      = .false.
       nullify(curr%iptr)
       nullify(curr%rptr)
       nullify(curr%next)

       ! free allocated space
       tmp => head%next
       deallocate(head)
       head => tmp
    end do

    cfg%hasDef = .true.

  end subroutine DQMC_Read_Def

  !---------------------------------------------------------------------!

  subroutine DQMC_Print_Def(cfg, OPT)
    !
    ! Purpose
    ! =======
    !    This subrotine prints condifuration definitions.
    !
    ! Arguments
    ! =========
    !
    type(config), intent(inout)  :: cfg          ! configuration
    integer, intent(in)          :: OPT          ! Input file handle

    ! ... Local Variable ...
    integer                :: i
    type(Param), pointer   :: curr 

    ! ... Executable ...

    write(OPT, 200) "ID", " Name","Default", "Type" 
    write(OPT, "(53('='))")
    do i = 1, cfg%nParam
       curr => cfg%record(i)
       if (curr%isArray) then
          write(OPT, 100) curr%id, curr%pname, curr%defaultval, &
               TYPE_STR(curr%ptype), "Array  "
       else
          write(OPT, 100) curr%id, curr%pname, curr%defaultval, &
               TYPE_STR(curr%ptype), "Scalar "
       end if
    end do

100 format(i5,2X,a10,2X,a10,2X,a10,2X,a10)
200 format(a5,1X,a5,5X,a10,1X,a10)
  end subroutine DQMC_Print_Def

  !---------------------------------------------------------------------!

  function DQMC_Find_Param(cfg, pname)result(id)
    !
    ! Purpose
    ! =======
    !    This subrotine returns the id of given parameter name.
    !    It will return -1 if no match is found.
    !
    ! Arguments
    ! =========
    !
    type(config), intent(in)     :: cfg       ! configuration
    character(*), intent(in)  :: pname        ! Parameter name
    integer :: id
    character(len(pname)) :: uppname, uprecord

    ! ... Local Variable ...
    integer :: high, low
    logical :: found
    ! ... Executable ...

    ! binary search
    low   = 1
    high  = cfg%nParam
    uppname = uppercase(pname)

    found = .false.
    do id = 1, high
       uprecord = uppercase(cfg%record(id)%pname)
       if (uppname .eq. uprecord) then
          found = .true.
          exit
       endif
    enddo

    !do while (.not. found .and. (low .le. high))
    !   id = (low+high)/2
    !   ! pname > param(id)
    !   if (LGT(pname, cfg%record(id)%pname)) then
    !      low  = id + 1
    !      ! pname < param(id)	 
    !   elseif (LLT(pname, cfg%record(id)%pname)) then
    !      high = id - 1
    !      !   
    !   else
    !      found = .true.
    !   end if
    !end do

    if (.not. found) then
       id = -1
    end if

  contains
   
    function uppercase(string) result(newstring)

       character(len=*), intent(in) :: string
       character(len=len(string)) :: newstring
       integer :: j

       do j = 1,len(string)
          if(string(j:j) >= "a" .and. string(j:j) <= "z") then
              newstring(j:j) = achar(iachar(string(j:j)) - 32)
          else
              newstring(j:j) = string(j:j)
          end if
       end do
       newstring = adjustl(newstring)

    end function uppercase

  end function DQMC_Find_Param

  !---------------------------------------------------------------------!

  subroutine DQMC_Read_Config(cfg)
    !
    ! Purpose
    ! =======
    !    This subrotine reads in parameters from a config file.
    !
    ! Arguments
    ! =========
    !
    type(config), intent(inout)  :: cfg          ! configuration

    ! ... Local Variable ...
    integer                :: ios, pos, line, j, id
    character(len=llen)    :: str, attr, val
    logical                :: found
    real(wp)               :: tmp(alen)          ! for reading t
    type(Param), pointer   :: curr 
    integer, parameter     :: funit = 10
    character(len=60)      :: iname
    integer                :: IPT, status

    ! ... Executable ...

    ! Fetch input file name from command line
    call get_command_argument(1, iname, STATUS=status)
    if (status > 0) then
       call DQMC_Error("failed to retrieve input file argument", 0)
    elseif (status == -1) then
       call DQMC_Error("String 'iname' is too small to hold input file name, recompile me with a larger ifile!", 0)
    end if

    ! Open input file
    call DQMC_open_file(iname, 'old', IPT)

    ! read def first
    if (.not. cfg%hasDef) then
       ! read def file
       inquire(file="config.def", exist=found)
       if (found) then
          open(unit = funit, file="config.def")
          call DQMC_Read_Def(cfg, funit)
          close(funit)
       else
          ! use default def
          call DQMC_Default_Def(cfg) 
       end if
    end if

    ! read real config file
    line = 0
    do
       line = line + 1
       read (unit=IPT, FMT="(a)", iostat=ios)  str

       ! end of file
       if (ios .ne. 0) then
          exit
       end if

       ! find comment  # and get rid of the tailing part
       pos = scan(str, COMMENT, .false.)
       if (pos .ne. 0) then
          ! find the comment sign
          if (pos .ge. 2) then
             str = str(1:pos-1)
          else
             str = ""
          end if
       end if

       ! trim the read in string
       if (len_trim(str) .gt. 0) then
          ! find separator = 
          pos = scan(str, SEPARAT, .false.)

          if (pos .ne. 0) then
             ! read name and data 
             attr = adjustl(str(1:pos-1))
             val  = adjustl(str(pos+1:llen))

             ! search parameter definition
             id = DQMC_Find_Param(cfg, attr)

             ! found it
             if (id .gt. 0) then
                curr => cfg%record(id)
                if (curr%isArray) then
                   ! array case
                   if  (curr%ptype .eq. TYPE_REAL .or. &
                        curr%ptype .eq. TYPE_INTEGER) then
                      j = 1
                      pos = scan(val, COMMA, .false.)
                      
                      ! For more than one t
                      do while(pos .gt. 0)
                         read(val(1:pos-1), *)  tmp(j)
                         val = val(pos+1:llen)
                         j = j + 1
                         pos = scan(val, COMMA, .false.)
                      end do
                      
                      ! the last one
                      read(val,*) tmp(j)
                      
                      ! copy to new allocated PR
                      if (curr%ptype .eq. TYPE_REAL) then
                         allocate(curr%rptr(j))
                         curr%ival = j
                         curr%rptr = tmp(1:j)
                      else
                         allocate(curr%iptr(j))
                         curr%ival = j
                         curr%iptr = int(tmp(1:j))
                      end if
                   else
                      call DQMC_Warning("Array only for real and integer",1)
                   end if
                else
                   ! scalar case
                   select case(curr%ptype)
                   case (TYPE_REAL)
                      read(val,*) curr%rval
                   case (TYPE_INTEGER)
                      read(val,*) curr%ival
                   case (TYPE_STRING)
                      curr%defaultval=val(1:slen)
                   end select
                end if

                ! mark the flag
                curr%isSet = .true.

             else
                call DQMC_Warning("Warning: unknown input:"//trim(str),1) 
             end if
          else
             call DQMC_Warning("cannot recog input line :", line)
          end if
       end if
    end do

  end subroutine DQMC_Read_Config

  !---------------------------------------------------------------------!
  ! Access functions 
  !---------------------------------------------------------------------!

  subroutine DQMC_Config_SetI(cfg, name, value)
    !
    ! Purpose
    ! =======
    !    This subrotine set configurations.
    !
    ! Arguments
    ! =========
    !
    type(config), intent(inout)  :: cfg          ! configuration
    character(len=*), intent(in) :: name
    integer, intent(in)          :: value  ! 

    ! ... Local variables...
    integer :: id

    ! ... Executable ...
    
    
    id = DQMC_Find_Param(cfg, name)
    if (id .gt. 0) then
       cfg%record(id)%ival = value
       cfg%record(id)%isSet = .true.
    else
       call DQMC_Error("cannot find parameter "//name, 0)
    end if

  end subroutine DQMC_Config_SetI

  !---------------------------------------------------------------------!

  subroutine DQMC_Config_SetR(cfg, name, value)
    !
    ! Purpose
    ! =======
    !    This subrotine set configurations.
    !
    ! Arguments
    ! =========
    !
    type(config), intent(inout)  :: cfg          ! configuration
    character(len=*), intent(in) :: name
    real(wp), intent(in)         :: value

    ! ... Local variables...
    integer :: id

    ! ... Executable ...
        
    id = DQMC_Find_Param(cfg, name)
    if (id .gt. 0) then
       cfg%record(id)%rval = value
       cfg%record(id)%isSet = .true.
    else
       call DQMC_Error("cannot find parameter "//name, 0)
    end if

  end subroutine DQMC_Config_SetR

  !---------------------------------------------------------------------!

  subroutine DQMC_Config_SetPR(cfg, name, n, value)
    !
    ! Purpose
    ! =======
    !    This subrotine set configurations.
    !
    ! Arguments
    ! =========
    !
    type(config), intent(inout)  :: cfg          ! configuration
    character(len=*), intent(in) :: name
    real(wp), intent(in)         :: value(n)
    integer, intent(in)          :: n

    ! ... Local variables...
    integer :: id

    ! ... Executable ...
    id = DQMC_Find_Param(cfg, name)
    if (id .gt. 0) then
       cfg%record(id)%ival = n
       if(associated(cfg%record(id)%rptr)) then
          deallocate(cfg%record(id)%rptr)
       end if
       allocate(cfg%record(id)%rptr(n))
       cfg%record(id)%rptr = value
       cfg%record(id)%isSet = .true.
    else
       call DQMC_Error("cannot find parameter "//name, 0)
    end if
    
  end subroutine DQMC_Config_SetPR

  !---------------------------------------------------------------------!

  subroutine DQMC_Config_SetPI(cfg, name, n, value)
    !
    ! Purpose
    ! =======
    !    This subrotine set configurations.
    !
    ! Arguments
    ! =========
    !
    type(config), intent(inout)  :: cfg          ! configuration
    character(len=*), intent(in) :: name
    integer, intent(in)          :: value(n)
    integer, intent(in)          :: n

    ! ... Local variables...
    integer :: id

    ! ... Executable ...
    id = DQMC_Find_Param(cfg, name)
    if (id .gt. 0) then
       cfg%record(id)%ival = n
       if(associated(cfg%record(id)%iptr)) then
          deallocate(cfg%record(id)%iptr)
       end if
       allocate(cfg%record(id)%iptr(n))
       cfg%record(id)%iptr = value
       cfg%record(id)%isSet = .true.
    else
       call DQMC_Error("cannot find parameter "//name, 0)
    end if
    
  end subroutine DQMC_Config_SetPI

  !---------------------------------------------------------------------!

  subroutine DQMC_Config_SetS(cfg, name, value)
    !
    ! Purpose
    ! =======
    !    This subrotine set configurations.
    !
    ! Arguments
    ! =========
    !
    type(config), intent(inout)  :: cfg          ! configuration
    character(*), intent(in)     :: name
    character(*), intent(in)     :: value

    ! ... Local variables...
    integer :: id

    ! ... Executable ...
    id = DQMC_Find_Param(cfg, name)
    if (id .gt. 0) then
       cfg%record(id)%defaultval = value
       cfg%record(id)%isSet = .true.
    else
       call DQMC_Error("cannot find parameter "//name, 0)
    end if

  end subroutine DQMC_Config_SetS

  !---------------------------------------------------------------------!

  function DQMC_Config_isSet(cfg, name) result(isSet)
    !
    ! Purpose
    ! =======
    !    This subrotine set configurations.
    !
    ! Arguments
    ! =========
    !
    type(config), intent(in)  :: cfg          ! configuration
    character(*), intent(in)  :: name
    logical                   :: isSet        ! 

    ! ... local variables
    integer :: id

    ! ... Executable ...

    id = DQMC_Find_Param(cfg, name)
    isSet = .false.
    if (id .gt. 0) then
       isSet = cfg%record(id)%isSet
    else
       call DQMC_Error("cannot find parameter "//name, 0)
    end if

  end function DQMC_Config_isSet

  !---------------------------------------------------------------------!

  subroutine DQMC_Config_GetI(cfg, name, value)
    !
    ! Purpose
    ! =======
    !    This subrotine set configurations.
    !
    ! Arguments
    ! =========
    !
    type(config), intent(in)  :: cfg          ! configuration
    character(*), intent(in)  :: name
    integer, intent(out)      :: value        ! 

    ! ... Local variables...
    integer :: id

    ! ... Executable ...
    id = DQMC_Find_Param(cfg, name)
    if (id .gt. 0) then
       if (.not.cfg%record(id)%isSet) then
          call DQMC_Warning(name//" wasn't initialized,&
               & used default setting.",1)
          read(cfg%record(id)%defaultval,*) value
       else
          value = cfg%record(id)%ival   
       end if
    else
       call DQMC_Error("cannot find parameter "//name, 0)
    end if
    

  end subroutine DQMC_Config_GetI

  !---------------------------------------------------------------------!

  subroutine DQMC_Config_GetR(cfg, name, value)
    !
    ! Purpose
    ! =======
    !    This subrotine set configurations.
    !
    ! Arguments
    ! =========
    !
    type(config), intent(in)  :: cfg          ! configuration
    character(*), intent(in)  :: name
    real(wp), intent(out)     :: value        ! 

    ! ... local variables
    integer :: id

    ! ... Executable ...

    id = DQMC_Find_Param(cfg, name)
    if (id .gt. 0) then
       if (.not.cfg%record(id)%isSet) then
          call DQMC_Warning(name//" wasn't initialized, &
               & used default setting.", 1)
          read(cfg%record(id)%defaultval,*) value
       else
          value = cfg%record(id)%rval   
       end if
    else
       call DQMC_Error("cannot find parameter "//name, 0)
    end if

  end subroutine DQMC_Config_GetR

  !---------------------------------------------------------------------!

  subroutine DQMC_Config_GetPR(cfg, name, n, value)
    !
    ! Purpose
    ! =======
    !    This subrotine set configurations.
    !
    ! Arguments
    ! =========
    !
    type(config), intent(in)         :: cfg          ! configuration
    character(*), intent(in)         :: name
    real(wp), pointer, intent(inout) :: value(:) 
    integer, intent(out)             :: n

    ! ... local variables
    integer :: id

    ! ... Executable ...
    id = DQMC_Find_Param(cfg, name)
!    write(*,*) "in DQMC_Config_GetPR, id=",id
!    write(*,*) "in DQMC_Config_GetPR, name=",name
!    write(*,*) "in DQMC_Config_GetPR, value associated?",associated(value)
    if (id .gt. 0) then
       if (.not.cfg%record(id)%isSet) then
          call DQMC_Warning(name//" wasn't initialized,&
               & used default setting.",1)
          n = 1
          if (associated(value)) then
             deallocate(value)
          end if
          allocate(value(n))
          read(cfg%record(id)%defaultval,*) value(1)
       else
          n = cfg%record(id)%ival
          if (associated(value)) then
             deallocate(value)
          end if
          allocate(value(n))
          value(1:n) = cfg%record(id)%rptr(1:n)
       end if
    else
       call DQMC_Error("cannot find parameter "//name, 0)
    end if
    
  end subroutine DQMC_Config_GetPR

  !---------------------------------------------------------------------!

  subroutine DQMC_Config_GetPI(cfg, name, n, value)
    !
    ! Purpose
    ! =======
    !    This subrotine set configurations.
    !
    ! Arguments
    ! =========
    !
    type(config), intent(in)        :: cfg          ! configuration
    character(*), intent(in)        :: name
    integer, pointer, intent(inout) :: value(:) 
    integer, intent(out)            :: n

    ! ... local variables
    integer :: id

    ! ... Executable ...
    id = DQMC_Find_Param(cfg, name)
!    write(*,*) "in DQMC_Config_GetPI, id=",id
!    write(*,*) "in DQMC_Config_GetPI, name=",name
!    write(*,*) "in DQMC_Config_GetPI, value associated?",associated(value)
    if (id .gt. 0) then
       if (.not.cfg%record(id)%isSet) then
          call DQMC_Warning(name//" wasn't initialized, &
               & used default setting.", 1)
          n = 1
          read(cfg%record(id)%defaultval,*) value(1)
       else
          n = cfg%record(id)%ival
          if (associated(value)) then
             deallocate(value)
          end if
          allocate(value(n))
          value(1:n) = cfg%record(id)%iptr(1:n)
       end if
    else
       call DQMC_Error("cannot find parameter "//name, 0)
    end if
    
  end subroutine DQMC_Config_GetPI

  !---------------------------------------------------------------------!

  subroutine DQMC_Config_GetS(cfg, name, value)
    !
    ! Purpose
    ! =======
    !    This subrotine set configurations.
    !
    ! Arguments
    ! =========
    !
    type(config), intent(in)  :: cfg          ! configuration
    character(*), intent(in)  :: name
    character(len=slen)       :: value

    ! ... local variables
    integer :: id

    ! ... Executable ...
    id = DQMC_Find_Param(cfg, name)
    if (id .gt. 0) then
       if (.not.cfg%record(id)%isSet) then
          call DQMC_Warning(name//" wasn't initialized, &
               & used default setting.", 1)
       end if
       value = cfg%record(id)%defaultval   
    else
       call DQMC_Error("cannot find parameter "//name, 0)
    end if


  end subroutine DQMC_Config_GetS

  !---------------------------------------------------------------------!
  
end module DQMC_Cfg
