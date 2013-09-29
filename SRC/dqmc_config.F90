module DQMC_Config

  use DQMC_UTIL
  implicit none 

  ! 
  ! This module contains subroutines to read input parameters.
  !
  ! Data Type
  ! =========
  !

  ! int data
  integer, parameter :: PARAM_NX     = 1
  integer, parameter :: PARAM_NY     = 2
  integer, parameter :: PARAM_NZ     = 3
  integer, parameter :: PARAM_N      = 4
  integer, parameter :: PARAM_L      = 5
  integer, parameter :: PARAM_NWARM  = 6
  integer, parameter :: PARAM_NPASS  = 7
  integer, parameter :: PARAM_NMEAS  = 8
  integer, parameter :: PARAM_NBIN   = 9
  integer, parameter :: PARAM_NHIST  = 10
  integer, parameter :: PARAM_IDUM   = 11
  integer, parameter :: PARAM_NORTH  = 12
  integer, parameter :: PARAM_NWRAP  = 13
  integer, parameter :: PARAM_TAUSK  = 14
  integer, parameter :: PARAM_NITVL  = 15
  integer, parameter :: PARAM_HSF    = 16
  integer, parameter :: PARAM_NTRY   = 17
  integer, parameter :: PARAM_ACCEPT = 18
  integer, parameter :: PARAM_REJECT = 19
  integer, parameter :: PRN_T        = 20
  integer, parameter :: PRN_U        = 21
  integer, parameter :: PRN_MU       = 22

  ! real data
  real, parameter :: PARAM_DTAU  = 1.0
  real, parameter :: PARAM_ERRAT = 2.0
  real, parameter :: PARAM_DIFF  = 3.0
  real, parameter :: PARAM_GAMMA = 4.0

  ! string data
  character(1), parameter :: PARAM_FNAME = char(1)
  character(1), parameter :: PARAM_HSFIN = char(2)
  character(1), parameter :: PARAM_HSFOUT= char(3)
  character(1), parameter :: PARAM_GFILE = char(4)

  ! pointer data
  integer, parameter :: PARAM_T     = 1
  integer, parameter :: PARAM_U     = 2
  integer, parameter :: PARAM_MU    = 3

  ! parameter counts
  integer, parameter :: TOTAL_INT  = 24
  integer, parameter :: TOTAL_INTPARAM  = 17
  integer, parameter :: TOTAL_REAL = 4
  integer, parameter :: TOTAL_STR  = 4
  integer, parameter :: TOTAL_PTR  = 5

  ! string parameters
  character, parameter :: COMMENT = "#"
  character, parameter :: SEPARAT = "="
  character, parameter :: COMMA   = ","

  ! HSF parameter
  integer, parameter:: HSF_OUTPUT_UNIT = 28
  integer, parameter:: HSF_INPUT_UNIT  = 27
  integer, parameter:: HSF_FROM_FILE   =  1
  integer, parameter:: HSF_FROM_MEMORY =  0
  integer, parameter:: HSF_RANDOM_GEN  =  -1

  character(len=*), parameter :: INTSTR(TOTAL_INT) = &
       & (/"nx    ", "ny    ", "nz    ", "n     ", "L     ", &
       &   "nwarm ", "npass ", "nmeas ", "nbin  ", "nhist ", &
       &   "seed  ", "north ", "nwrap ", "tausk ", "nitvl ", &
       &   "HSF   ", "numtry", "accept", "reject", "nt_up ", &
       &   "nt_dn ", "n_U   ", "nmu_up", "nmu_dn" /)

  character(len=*), parameter :: REALSTR(TOTAL_REAL) = &
        & (/"dtau   ", "errrate", "difflim", "gamma  "/)

  character(len=*), parameter :: STRSTR(TOTAL_STR) = &
       (/"fname ", "HSFin ", "HSFout", "gfile "/)

  character(len=*), parameter :: PTRSTR(TOTAL_PTR) = (/"t_up ", "t_dn ","  U  ", "mu_up", "mu_dn"/)

  ! Input file type
  integer, parameter :: CONFIG_LEGACY    = 1
  integer, parameter :: CONFIG_CONFIG    = 2
  integer, parameter :: CONFIG_XML       = 3
  integer, parameter :: CONFIG_GUI       = 4

  ! ===========
  ! type define
  ! ===========

  integer, parameter :: namelen = 50

  type Config
     ! integer parameters
     integer   :: I(TOTAL_INT)
     integer   :: IFlag(TOTAL_INT)

     ! real number parameters
     real(wp)  :: R(TOTAL_REAL)
     integer   :: RFlag(TOTAL_REAL)

     ! string parameters
     character(namelen) :: S(TOTAL_STR)
     integer   :: SFlag(TOTAL_STR)

     ! pointer parameters
     real(wp), pointer :: t(:)
     real(wp), pointer :: U(:)
     real(wp), pointer :: mu(:)
     integer   :: PRFlag(TOTAL_PTR)

  end type Config

  ! ================
  ! Access functions
  ! ================

  interface CFG_Set
     module procedure DQMC_Config_SetI, DQMC_Config_SetR
     module procedure DQMC_Config_SetS, DQMC_Config_SetPR
  end interface

  interface CFG_Get
     module procedure DQMC_Config_GetI, DQMC_Config_GetR
     module procedure DQMC_Config_GetS, DQMC_Config_GetPR
  end interface

contains

  !---------------------------------------------------------------------!

  subroutine DQMC_Config_Free(cfg)
    ! 
    ! Releases the PTR in cfg
    !
    type(config), intent(inout)  :: cfg          ! configuration

    call DQMC_Free(cfg%t)
    call DQMC_Free(cfg%U)
    call DQMC_Free(cfg%mu)

  end subroutine DQMC_Config_Free
  
  !---------------------------------------------------------------------!

  subroutine DQMC_Config_default(cfg)
    !
    ! Purpose
    ! =======
    !    This subrotine sets default values for parameters.
    !
    ! Arguments
    ! =========
    !
    type(config), intent(inout)  :: cfg          ! configuration

    ! ... Executable ...


    ! initeger parameters
    cfg%I(PARAM_NX    ) = 0
    cfg%I(PARAM_NY    ) = 0
    cfg%I(PARAM_NZ    ) = 0
    cfg%I(PARAM_N     ) = 0
    cfg%I(PARAM_L     ) = 0
    cfg%I(PARAM_NWARM ) = 500
    cfg%I(PARAM_NPASS ) = 1000 
    cfg%I(PARAM_NMEAS ) = 12
    cfg%I(PARAM_NBIN  ) = 10
    cfg%I(PARAM_NHIST ) = 50
    cfg%I(PARAM_IDUM  ) = 0
    cfg%I(PARAM_NORTH ) = 12
    cfg%I(PARAM_NWRAP ) = 12
    cfg%I(PARAM_TAUSK ) = 10
    cfg%I(PARAM_NITVL ) = 4
    cfg%I(PARAM_HSF   ) = -1
    cfg%I(PARAM_NTRY  ) = 0
    cfg%I(PARAM_ACCEPT) = 0
    cfg%I(PARAM_REJECT) = 0

    ! real parameters
    cfg%R(int(PARAM_DTAU)) = 0
    cfg%R(int(PARAM_ERRAT))= 0.001D0
    cfg%R(int(PARAM_DIFF)) = 0.001D0
    cfg%R(int(PARAM_GAMMA))= 0.5D0


    ! pointer parameters
    nullify(cfg%t)
    nullify(cfg%U)
    nullify(cfg%mu)
    
    ! string parameter
    cfg%S(ichar(PARAM_FNAME ))  = "quest"
    cfg%S(ichar(PARAM_HSFIN ))  = ""
    cfg%S(ichar(PARAM_HSFOUT))  = "HSFout"
    cfg%S(ichar(PARAM_GFILE))   = ""

    ! Flag of set
    cfg%IFlag  = 0
    cfg%RFlag  = 0
    cfg%SFlag  = 0
    cfg%PRFlag  = -1

    ! The following vars must be initialized
    cfg%IFlag(PARAM_NX) = -1
    cfg%IFlag(PARAM_NY) = -1
    cfg%IFlag(PARAM_NZ) = -1
    cfg%IFlag(PARAM_N)  = -1
    cfg%IFlag(PARAM_L)  = -1
    cfg%RFlag(int(PARAM_DTAU)) = -1
    cfg%SFlag(ichar(PARAM_HSFIN )) = -1

  end subroutine DQMC_Config_default

  !---------------------------------------------------------------------!

  subroutine DQMC_Config_Read(cfg, IPT, fmt)
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
    integer, intent(in)          :: fmt          ! Format of config file

    ! ... Executable ...

    call DQMC_Config_default(cfg)

    select case(fmt)
    case (CONFIG_CONFIG)
       call dqmc_config_config(cfg, IPT)
    case (CONFIG_XML)
    case (CONFIG_GUI)
    case default
       call DQMC_Error("Config file format does not support", fmt)       
    end select

  end subroutine DQMC_Config_Read

 !---------------------------------------------------------------------!
  
  subroutine DQMC_Config_Config(cfg, IPT)
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
    integer, parameter     :: strlen = 256
    integer                :: ios, pos, line, i, j
    character(len=strlen)  :: str, attr, val
    logical                :: found
    real(wp)               :: tmp(strlen)        ! for reading t

    ! ... Executable ...

    line = 0
    do 
       line = line + 1
       read (unit=IPT, FMT="(a)", iostat=ios)  str

       ! end of file
       if (ios .ne. 0) then
          exit
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

       ! trim the read in string
       if (len_trim(str) .gt. 0) then
          ! find separator = 
          pos = scan(str, SEPARAT, .false.)

          if (pos .ne. 0) then
             ! read name and data 
             attr = adjustl(str(1:pos-1))
             val  = adjustl(str(pos+1:strlen))
             
             ! =============================================
             ! search the integer parameters
             found = .false.
             do i = 1, TOTAL_INTPARAM 
                if (trim(attr) .eq. trim(INTSTR(i))) then
                   read(val, *) cfg%I(i)
                   cfg%IFlag(i) = 1
                   found = .true.
                   exit
                end if
             end do
             
             ! search the real parameters
             if (.not. found) then
                do i = 1, TOTAL_REAL
                   if (trim(attr) .eq. trim(REALSTR(i))) then
                      read(val, *) cfg%R(i)
                      cfg%RFlag(i) = 1
                      found = .true.
                      exit
                   end if
                end do
             end if

             ! search the string parameters
             if (.not. found) then
                do i = 1, TOTAL_STR
                   if (trim(attr) .eq. trim(STRSTR(i))) then
                      read(val,*) cfg%S(i)
                      cfg%SFlag(i) = 1
                      found = .true.
                      exit
                   end if
                end do
             end if
             
             ! search for the pointer parameters (t)
             if (.not. found) then
                do i = 1, TOTAL_PTR
                   if (trim(attr) .eq. trim(PTRSTR(i))) then
                      ! initialize
                      j = 1
                      pos = scan(val, COMMA, .false.)
                      
                      ! For more than one t
                      do while(pos .gt. 0)
                         read(val(1:pos-1), *)  tmp(j)
                         val = val(pos+1:strlen)
                         j = j + 1
                         pos = scan(val, COMMA, .false.)
                      end do
                      
                      ! the last one
                      read(val,*) tmp(j)
                      
                      ! copy to new allocated PR
                      select case(i)
                      case (PARAM_t)
                         call DQMC_Reshape(j, cfg%t)
                         cfg%t  = tmp(1:j)
                         call CFG_Set(cfg, PRN_T, j)
                      case (PARAM_U)
                         call CFG_Set(cfg, PRN_U, j)
                         call DQMC_Reshape(j, cfg%U)
                         cfg%U  = tmp(1:j)
                      case (PARAM_MU)
                         call CFG_Set(cfg, PRN_MU, j)
                         call DQMC_Reshape(j, cfg%mu)
                         cfg%mu = tmp(1:j)
                      end select
                      
                      ! set flags
                      cfg%PRFlag(i) = 1
                      found = .true.
                   end if
                end do
             end if


             ! =============================================

             ! give warning message
             if (.not. found) then
                write(STDERR,*) "Warning: unknown input:", trim(str) 
             end if
          else
             call DQMC_Error("Cannot recognize input line", line)
          end if
       end if
    end do

  end subroutine DQMC_Config_Config

 !---------------------------------------------------------------------!

  subroutine DQMC_Config_Print(cfg, OPT)
    !
    ! Purpose
    ! =======
    !    This subrotine reads in parameters.
    !
    ! Arguments
    ! =========
    !
    type(config), intent(in)  :: cfg          ! configuration
    integer, intent(in)       :: OPT          ! Input file handle

    ! ... Local Variable ...
    integer :: i
    character(*), parameter :: FMT_STRINT  = "(a10, i10)"
    character(*), parameter :: FMT_STRDBL  = "(a10, f12.6)"

    ! ... Executable ...

    do i = 1, TOTAL_INTPARAM 
       write(OPT, FMT_STRINT) INTSTR(i), cfg%I(i)
    end do

    do i = 1, TOTAL_REAL
       write(OPT, FMT_STRDBL) REALSTR(i), cfg%R(i)
    end do

    do i = 1, TOTAL_STR
       write(OPT, *) STRSTR(i), cfg%S(i)
    end do

    write(OPT, *) PTRSTR(PARAM_T), cfg%t
    write(OPT, *) PTRSTR(PARAM_U), cfg%U
    write(OPT, *) PTRSTR(PARAM_MU), cfg%mu

  end subroutine DQMC_Config_Print

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
    integer, intent(in)          :: name, value  ! 

    ! ... Executable ...

    cfg%I(name) = value
    cfg%IFlag(name) = 1

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
    real, intent(in)             :: name
    real(wp), intent(in)         :: value

    ! ... Executable ...

    cfg%R(int(name)) = value
    cfg%RFlag(int(name)) = 1

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
    integer, intent(in)          :: name, n
    real(wp), intent(in)         :: value(n)

    ! ... Executable ...

    select case(name)
    case (PARAM_t)
       call DQMC_Reshape(n, cfg%t)
       cfg%t = value
       cfg%I(PRN_T)  = n
    case (PARAM_U)
       call DQMC_Reshape(n, cfg%U)
       cfg%U = value
       cfg%I(PRN_U)  = n
    case (PARAM_mu)
       call DQMC_Reshape(n, cfg%mu)
       cfg%mu = value
       cfg%I(PRN_MU)  = n
    end select

    cfg%PRFlag(name) = 1
    
  end subroutine DQMC_Config_SetPR

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
    character, intent(in)        :: name
    character(*), intent(in)     :: value

    ! ... Executable ...

    cfg%S(ichar(name)) = value
    cfg%SFlag(ichar(name)) = 1

    print *, "File name = ",  cfg%S(ichar(name)) 

  end subroutine DQMC_Config_SetS

  !---------------------------------------------------------------------!

  function DQMC_Config_GetI(cfg, name) result(value)
    !
    ! Purpose
    ! =======
    !    This subrotine set configurations.
    !
    ! Arguments
    ! =========
    !
    type(config), intent(in)  :: cfg          ! configuration
    integer, intent(in)       :: name
    integer                   :: value        ! 

    ! ... local variables
    character(10) :: strval
    
    ! ... Executable ...

    value = cfg%I(name)
    if (cfg%IFlag(name) .le. 0) then
       if (cfg%IFlag(name) .eq. -1) then
          call DQMC_Error(INTSTR(name)//" must be initialized", 1)
       else
          write(strval, "(i10)")  cfg%I(name)
          call DQMC_Warning(INTSTR(name)//" is not initialized. &
               & Use default "//strval, 0)
       end if
    end if

  end function DQMC_Config_GetI

  !---------------------------------------------------------------------!

  function DQMC_Config_GetR(cfg, name) result(value)
    !
    ! Purpose
    ! =======
    !    This subrotine set configurations.
    !
    ! Arguments
    ! =========
    !
    type(config), intent(in)  :: cfg          ! configuration
    real, intent(in)          :: name
    real(wp)                  :: value

    ! ... local variables
    integer :: idx
    character(20) :: strval

    ! ... Executable ...

    idx = int(name)
    value = cfg%R(idx)
    if (cfg%RFlag(idx).le. 0) then
       if (cfg%RFlag(idx) .eq. -1) then
          call DQMC_Error(REALSTR(idx)//" must be initialized", 1)
       else
          write(strval, "(f20.10)")  cfg%R(idx)
          call DQMC_Warning(REALSTR(idx)//" is not initialized. &
               & Use default"//strval, 0)
       end if
    end if

  end function DQMC_Config_GetR

  !---------------------------------------------------------------------!

  function DQMC_Config_GetPR(cfg, name, n) result(value)
    !
    ! Purpose
    ! =======
    !    This subrotine set configurations.
    !
    ! Arguments
    ! =========
    !
    type(config), intent(in)  :: cfg          ! configuration
    integer, intent(in)       :: name
    real(wp), pointer         :: value(:)
    integer, intent(out)      :: n

    ! ... Executable ...

    if (cfg%PRFlag(name) .ne. 1) then
       call DQMC_Error(PTRSTR(name)//" must be initialized", 1)
    end if

    select case(name)
    case(PARAM_t)
       n = cfg%I(PRN_t)
       value => cfg%t
    case(PARAM_U)
       n = cfg%I(PRN_U)
       value => cfg%U
    case(PARAM_mu)
       n = cfg%I(PRN_mu)
       value => cfg%mu
    end select

  end function DQMC_Config_GetPR

  !---------------------------------------------------------------------!

  function DQMC_Config_GetS(cfg, name, str) result(valid)
    !
    ! Purpose
    ! =======
    !    This subrotine set configurations.
    !
    ! Arguments
    ! =========
    !
    type(config), intent(in)  :: cfg          ! configuration
    character, intent(in)     :: name
    character(len=namelen)    :: str
    logical                   :: valid

    ! ... local variables
    integer :: idx

    ! ... Executable ...

    idx = ichar(name)
    str = cfg%S(idx)
    valid = (cfg%SFlag(idx) .gt. 0)
    if (cfg%SFlag(idx).le. 0) then
       if (cfg%SFlag(idx) .eq. -1) then
          call DQMC_Error(STRSTR(idx)//" must be initialized", 1)
!       else
!          call DQMC_Warning(STRSTR(idx)//" is not initialized. &
!               & Use default "//cfg%S(idx), 0)
       end if
    end if

  end function DQMC_Config_GetS

  !---------------------------------------------------------------------!

end module DQMC_Config
