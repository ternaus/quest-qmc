module DQMC_CheckerBoard
#include "dqmc_include.h"

  use DQMC_UTIL
  use DQMC_WSPACE
!  use BLAS_MOD
!  use LAPACK_MOD

  implicit none 

  ! 
  ! This module defines the type and subroutines for propagator B.
  ! B is defined by the checkerboard method. See [1]
  !
  ! For a 2-dimensional model
  !
  ! B = exp(K) = exp(K_{4})exp(K_{3})exp(K_{2})exp(K_{1})
  !            = B_4B_3B_2B_1
  !            
  !    [1] W. Hanke and Yu.V. Kopaev, "Stable Numerical Simulation 
  !        of Models of Interating Electrons in Condensed-Matter Physics."
  !        Chapter 4. Elsevier Science Pub. 
  !
  !  Data Type
  !  =========
  !
  type matB
     integer  :: n                      ! dim of B
     integer  :: m                      ! number of neighbors of lattice
     integer, pointer  :: A(:,:)        ! Adjacency info
     real(wp), pointer :: sinht(:)      ! sinh(t), cosh(t)
     real(wp), pointer :: cosht(:)      ! sinh(t), cosh(t)
     real(wp), pointer :: exptaumu(:)   ! parameters for checkerboard method
     real(wp), pointer :: work(:)       ! parameters for checkerboard method
     character(12) :: name

     logical :: exactb  !unused in ckb
  end type matB

contains


  !-------------------------------------------------------------------------!

  subroutine DQMC_B_Init(n, B, WS, Adj, ckb, t, mu, dtau)
    !
    ! Purpose
    ! =======
    !    This subroutine initiliazes the data type of Green function.
    !
    ! Pre-assumption
    ! ==============
    !    This module only used for one band Hubbard's model.
    !
    ! Arguments
    ! =========
    !

    use DQMC_STRUCT

    integer, intent(in)       :: n                     ! Number of sites
    type(MatB), intent(out)   :: B                     ! MatB
    type(CCS), intent(in)     :: adj                   ! adjacent info
    type(CCS), intent(in)     :: ckb
    real(wp), intent(in)      :: t(*)              ! model parameter
    real(wp), intent(in)      :: mu(n), dtau  
    type(WSpace), intent(in), target  :: WS                    ! shared working space

    ! ... local variables    ...
    integer :: h, k, i, j, nt, nckb
    integer :: ckbmat(n,n)
    integer :: hopmat(n,n)
    real(wp), pointer :: dum(:,:) 
    
    ! ... Executable ...

    ! So that the compiler does not warn for WS being unused
    dum   => WS%R1

    if(sum((ckb%row-adj%row)**2)+sum((ckb%cstart-adj%cstart)**2)>0)then
       write(*,*)'ckb and adj do not conform. Stop.'
       stop
    endif

    B%n    = n
    B%m    = ckb%nnz / 2
    B%name = "Checkerboard"
    nckb   = maxval(ckb%A)
    nt     = maxval(Adj%A)

    allocate(B%cosht(nt))
    allocate(B%sinht(nt))
    allocate(B%A(3,B%m))
    allocate(B%exptaumu(n))

    if(nt .gt. 0)then
       B%sinht(1:nt) = sinh(dtau*0.5*t(1:nt))
       B%cosht(1:nt) = cosh(dtau*0.5*t(1:nt))
    endif

    call dqmc_ccs_fill(n, ckbmat, ckb)
    call dqmc_ccs_fill(n, hopmat, Adj)
    h = 0
    do k = 1, nckb
       do j = 2, n 
          do i = 1, j-1
             if (ckbmat(i,j) .eq. k) then
                h = h + 1
                B%A(1, h) = i
                B%A(2, h) = j
                B%A(3, h) = hopmat(i,j)
             endif  
          enddo
       enddo
    enddo

    do i=1,n
       B%exptaumu(i) = exp(dtau*mu(i))
    enddo

    allocate(B%work(n))

  end subroutine DQMC_B_Init

  !-------------------------------------------------------------------------!

  subroutine DQMC_B_Free(B)
    !
    ! Purpose
    ! =======
    !    This subroutine frees memory of B.
    !
    ! Arguments
    ! =========
    !
    type(MatB), intent(inout)  :: B  ! MatB

    ! ... Executable ...

    deallocate(B%A, B%work)

  end subroutine DQMC_B_Free


  !----------------------------------------------------------------------!
  
  subroutine DQMC_MultB_Left(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine uses checkerboard method to multiply B from
    !    leftside.
    !
    !      B_i = V_i*B
    !
    ! Pre-assumption
    ! ==============
    !    1. Matrix M, B and C are square with the same order n.
    !    2. Vector V_i is of length n. 
    !    3. Flag lr is either 'l' or 'r'.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)              :: n              ! The order of matrix M  
    type(MatB), intent(inout)        :: B              !
    real(wp), intent(in)             :: V_i(n)         !
    real(wp), intent(inout), target  :: M(n,n)         !  
    real(wp), intent(inout), target  :: C(n,n)         ! 

    ! ... Local variables ...
    integer  :: i, j, k, h
    real(wp) :: a

    ! ... executable ...

    call DQMC_Trans(n, C, M)

    do k = B%m, 1, -1
       i = B%A(1,k)
       j = B%A(2,k)
       h = B%A(3,k)
       a = B%sinht(h)/B%cosht(h)
       call daxpy(n, a, C(:,j), 1, C(:,i), 1)
       a = B%cosht(h)*B%sinht(h)
       call daxpy(n, a, C(:,i), 1, C(:,j), 1)
       a = B%cosht(h)
       call dscal(n, a, C(:,i), 1)
       a = 1.d0/B%cosht(h)
       call dscal(n, a, C(:,j), 1)
    enddo

    do k = 1, B%m
       i = B%A(1,k)
       j = B%A(2,k)
       h = B%A(3,k)
       a = B%sinht(h)/B%cosht(h)
       call daxpy(n, a, C(:,j), 1, C(:,i), 1)
       a = B%cosht(h)*B%sinht(h)
       call daxpy(n, a, C(:,i), 1, C(:,j), 1)
       a = B%cosht(h)
       call dscal(n, a, C(:,i), 1)
       a = 1.d0/B%cosht(h)
       call dscal(n, a, C(:,j), 1)
    enddo


    ! Combin V_i with the coefficient of B
    B%work(1:n) = B%exptaumu(1:n) * V_i(1:n)
    call DQMC_ScaleCol(n, C, B%work)

    call DQMC_Trans(n, M, C)

  end subroutine DQMC_MultB_Left

  !----------------------------------------------------------------------!
  
  subroutine DQMC_MultB_Right(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine uses checkerboard method to multiply B from
    !    righthand side.
    !
    !      B_i = V_i*B
    !
    ! Pre-assumption
    ! ==============
    !    1. Matrix M, B and C are square with the same order n.
    !    2. Vector V_i is of length n. 
    !    3. Flag lr is either 'l' or 'r'.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)              :: n              ! The order of matrix M  
    type(MatB), intent(inout)        :: B              !
    real(wp), intent(in)             :: V_i(n)         !
    real(wp), intent(inout), target  :: M(n,n)         !  
    real(wp), intent(inout), target  :: C(n,n)         !    

    ! ... Local variables ...
    integer  :: i, j, k, h
    real(wp) :: a

    ! ... executable ...

    ! Combin V_i with the coefficient of B
    B%work(1:n) = B%exptaumu(1:n) * V_i(1:n)
    call DQMC_ScaleCol(n, M, B%work)

    do k = 1, B%m
       i = B%A(1,k)
       j = B%A(2,k)
       h = B%A(3,k)
       a = B%sinht(h)/B%cosht(h)
       call daxpy(n, a, M(:,j), 1, M(:,i), 1)
       a = B%cosht(h)*B%sinht(h)
       call daxpy(n, a, M(:,i), 1, M(:,j), 1)
       a = B%cosht(h)
       call dscal(n, a, M(:,i), 1)
       a = 1.d0/B%cosht(h)
       call dscal(n, a, M(:,j), 1)
    enddo
    
    do k = B%m, 1, -1
       i = B%A(1,k)
       j = B%A(2,k)
       h = B%A(3,k)
       a = B%sinht(h)/B%cosht(h)
       call daxpy(n, a, M(:,j), 1, M(:,i), 1)
       a = B%cosht(h)*B%sinht(h)
       call daxpy(n, a, M(:,i), 1, M(:,j), 1)
       a = B%cosht(h)
       call dscal(n, a, M(:,i), 1)
       a = 1.d0/B%cosht(h)
       call dscal(n, a, M(:,j), 1)
    enddo

  end subroutine DQMC_MultB_Right

  !--------------------------------------------------------------------------!
  
  subroutine DQMC_MultBi_Left(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine uses checkerboard method to multiply inv(B) from
    !    leftside.
    !
    ! Pre-assumption
    ! ==============
    !    1. Matrix M, B and C are square with the same order n.
    !    2. Vector V_i is of length n. 
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)              :: n              ! The order of matrix M  
    type(MatB), intent(inout)        :: B              !
    real(wp), intent(in)             :: V_i(n)         !
    real(wp), intent(inout), target  :: M(n,n)         !  
    real(wp), intent(inout), target  :: C(n,n)         !
    
    ! ... Local variables ...
    integer  :: i, j, k, h
    real(wp) :: a

    ! ... executable ...

    call DQMC_Trans(n, C, M)

    B%work(1:n) = 1.d0 / (B%exptaumu(1:n) * V_i(1:n))
    call DQMC_ScaleCol(n, C, B%work)

    do k = 1, B%m
       i = B%A(1,k)
       j = B%A(2,k)
       h = B%A(3,k)
       a = -B%sinht(h)/B%cosht(h)
       call daxpy(n, a, C(:,j), 1, C(:,i), 1)
       a = -B%cosht(h)*B%sinht(h)
       call daxpy(n, a, C(:,i), 1, C(:,j), 1)
       a = B%cosht(h)
       call dscal(n, a, C(:,i), 1)
       a = 1.d0/B%cosht(h)
       call dscal(n, a, C(:,j), 1)
    enddo

    do k = B%m, 1, -1
       i = B%A(1,k)
       j = B%A(2,k)
       h = B%A(3,k)
       a = -B%sinht(h)/B%cosht(h)
       call daxpy(n, a, C(:,j), 1, C(:,i), 1)
       a = -B%cosht(h)*B%sinht(h)
       call daxpy(n, a, C(:,i), 1, C(:,j), 1)
       a = B%cosht(h)
       call dscal(n, a, C(:,i), 1)
       a = 1.d0/B%cosht(h)
       call dscal(n, a, C(:,j), 1)
    enddo

    ! Combin V_i with the coefficient of B

    call DQMC_Trans(n, M, C)

  end subroutine DQMC_MultBi_Left

  !--------------------------------------------------------------------------!
  
  subroutine DQMC_MultBi_Right(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine uses checkerboard method to multiply inv(B) from
    !    right hand side.
    !
    ! Pre-assumption
    ! ==============
    !    1. Matrix M, B and C are square with the same order n.
    !    2. Vector V_i is of length n. 
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)              :: n              ! The order of matrix M  
    type(MatB), intent(inout)        :: B              !
    real(wp), intent(in)             :: V_i(n)         !
    real(wp), intent(inout), target  :: M(n,n)         !  
    real(wp), intent(inout), target  :: C(n,n)         !
    
    ! ... Local variables ...
    integer  :: i, j, k, h
    real(wp) :: a

    ! ... executable ...

    do k = B%m, 1, -1
       i = B%A(1,k)
       j = B%A(2,k)
       h = B%A(3,k)
       a = -B%sinht(h)/B%cosht(h)
       call daxpy(n, a, M(:,j), 1, M(:,i), 1)
       a = -B%cosht(h)*B%sinht(h)
       call daxpy(n, a, M(:,i), 1, M(:,j), 1)
       a = B%cosht(h)
       call dscal(n, a, M(:,i), 1)
       a = 1.d0/B%cosht(h)
       call dscal(n, a, M(:,j), 1)
    enddo

    do k = 1, B%m
       i = B%A(1,k)
       j = B%A(2,k)
       h = B%A(3,k)
       a = -B%sinht(h)/B%cosht(h)
       call daxpy(n, a, M(:,j), 1, M(:,i), 1)
       a = -B%cosht(h)*B%sinht(h)
       call daxpy(n, a, M(:,i), 1, M(:,j), 1)
       a = B%cosht(h)
       call dscal(n, a, M(:,i), 1)
       a = 1.d0/B%cosht(h)
       call dscal(n, a, M(:,j), 1)
    enddo

    B%work(1:n) = 1.d0 / (B%exptaumu(1:n) * V_i(1:n))
    call DQMC_ScaleCol(n, M, B%work)

  end subroutine DQMC_MultBi_Right

  !----------------------------------------------------------------------!

  subroutine DQMC_MultrtB0_Left(n, M, B, C)
    !
    ! Purpose
    ! =======
    !    This subroutine mulitplies matrix M by the matrix sqrt(B0) from
    !    leftside.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)      :: n              ! The order of matrix M  
    real(wp), intent(inout)  :: M(n,n)         !  
    type(MatB), intent(in)   :: B              !
    real(wp), intent(inout)  :: C(n,n)         ! working space   

    integer :: i

    ! Dummy executable to avoid warning
    C = M
    i = B%n

  end subroutine DQMC_MultrtB0_Left

  !----------------------------------------------------------------------!

  subroutine DQMC_MultrtB0i_Right(n, M, B, C)
    !
    ! Purpose
    ! =======
    !    This subroutine mulitplies matrix M by the matrix sqrt(B0_i) from
    !    rightside.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)      :: n              ! The order of matrix M  
    real(wp), intent(inout)  :: M(n,n)         !  
    type(MatB), intent(in)   :: B              !
    real(wp), intent(inout)  :: C(n,n)         ! working space  

    integer :: i

    ! Dummy executable to avoid warning
    C = M
    i = B%n

  end subroutine DQMC_MultrtB0i_Right

  !-------------------------------------------------------------------------!

  subroutine DQMC_GetB(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine returns M = V_iB
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)      :: n              ! The order of matrix M  
    real(wp), intent(inout)  :: M(n,n)         !  
    type(MatB), intent(inout):: B              !
    real(wp), intent(in)     :: V_i(n)         !
    real(wp), intent(inout)  :: C(n,n)         !
    ! ... Executable ...

    call DQMC_Eye(n, M)
    call DQMC_MultB_Left(n, M, B, V_i, C)

  end subroutine DQMC_GetB

  !-----------------------------------------------------------------------!

  subroutine DQMC_GetBi(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine returns M = inv(B)inv(V_i)
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)      :: n              ! The order of matrix M  
    real(wp), intent(inout)  :: M(n,n)         !  
    type(MatB), intent(inout):: B              ! MatB 
    real(wp), intent(in)     :: V_i(n)         !
    real(wp), intent(inout)  :: C(n,n)         !
    ! ... Executable ...

    call DQMC_Eye(n, M)
    call DQMC_MultBi_Left(n, M, B, V_i, C)

  end subroutine DQMC_GetBi

  !-----------------------------------------------------------------------!
  
end module DQMC_CheckerBoard
