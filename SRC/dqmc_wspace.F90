module DQMC_WSpace

  use DQMC_UTIL

  implicit none 
  
  ! 
  ! This module is used for memory mangement.
  !
  ! Data Type
  ! =========
  !
  type WSpace
     integer :: lw(5)
     real(wp), pointer :: R1(:,:)
     real(wp), pointer :: R2(:,:)
     real(wp), pointer :: R3(:,:)
     real(wp), pointer :: R4(:,:)
     real(wp), pointer :: R5(:)  
     real(wp), pointer :: R6(:)  
     real(wp), pointer :: R7(:)  
     real(wp), pointer :: R8(:,:)
     integer,  pointer :: I1(:)  
     integer,  pointer :: I2(:)  
  end type WSpace

  !
  ! Parameters
  ! ==========
  !
  integer, parameter  :: LA_SYEV   = 1     ! Index of LAPACK subroutines.
  integer, parameter  :: LA_GEQRF  = 2     ! This is used to specify 
  integer, parameter  :: LA_ORGQR  = 3     ! the size of working space.
  integer, parameter  :: LA_ORMQR  = 4     !  
  integer, parameter  :: LA_GETRI  = 5     !  

contains

  !---------------------------------------------------------------------!

  subroutine DQMC_WSpace_Allocate(n, mx_nbr, WS)
    !
    ! Purpose
    ! =======
    !    This subrotine allocates memory for MPool
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)           :: n, mx_nbr  ! 
    type(WSpace), intent(inout)   :: WS         ! working space

    ! ... Local scalar ...
    integer  :: nb, info

    ! ... Executable ....

    info = 0
    ! Initialize working space
    allocate(WS%R1(n,n))
    allocate(WS%R2(n,n))
    allocate(WS%R3(n,n))

    nb = max(n, mx_nbr)
    !allocate(WS%R4(nb,nb))
    allocate(WS%R4(n,n))
    allocate(WS%R5(n))
    allocate(WS%R6(n))
    allocate(WS%I1(n))
    allocate(WS%I2(n))

    ! Get optimal block size for working space
    ! Part of codes are copied from LAPACK
    call lapack_dsyev('V', 'L', n, WS%R1, n, WS%R2, WS%R3, -1, info)
    WS%lw(LA_SYEV) = int(WS%R3(1,1)) 
    call lapack_dgeqp3(n, n, WS%R1, n, WS%I1, WS%R5, WS%R3, -1, info)
    WS%lw(LA_GEQRF) = int(WS%R3(1,1))
    call lapack_dorgqr(n, n, n, WS%R1, n, WS%R5, WS%R3, -1, info)
    WS%lw(LA_ORGQR) = int(WS%R3(1,1))
    call lapack_dormqr('R', 'N', n, n, n, WS%R1, n, WS%R5, WS%R2, n, WS%R3, -1, info)
    WS%lw(LA_ORMQR) = int(WS%R3(1,1))
    call lapack_dgetri(n, WS%R1, n, WS%I1, WS%R3, -1, info)
    WS%lw(LA_GETRI) = int(WS%R3(1,1))

    nb = maxval(WS%lw)
    allocate(WS%R7(nb))
    allocate(WS%R8(n,n))

  end subroutine DQMC_WSpace_Allocate

  !---------------------------------------------------------------------!
  
  subroutine DQMC_WSpace_Free(WS)
    !
    ! Purpose
    ! =======
    !    This subrotine frees working space
    !
    ! Arguments
    ! =========
    !
    type(WSpace), intent(inout)   :: WS         ! working space

    ! ... Executable ....

    deallocate(WS%R1,WS%R2,WS%R3,WS%R4,WS%R5,WS%R6,WS%R7,WS%R8)
    deallocate(WS%I1, WS%I2)

  end subroutine DQMC_WSpace_Free

end module DQMC_Wspace
