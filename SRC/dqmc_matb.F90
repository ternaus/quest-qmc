module DQMC_MATB
#include "dqmc_include.h"

  use DQMC_UTIL
  use DQMC_Struct
  use DQMC_WSPACE

  implicit none 

  ! 
  ! This module defines the type and subroutines for propagator B.
  !
  !  Data Type
  !  =========
  !
  type matB
     integer  :: n                   ! dim of B
     real(wp), pointer :: K(:,:)     ! The Hamiltonian
     real(wp), pointer :: B(:,:)     ! Matrix B 
     real(wp), pointer :: Bi(:,:)    ! Inverse of matrix B
     real(wp), pointer :: rtB(:,:)   ! Matrix B 
     real(wp), pointer :: rtBi(:,:)  ! Inverse of matrix B
     character(12)     :: name
     logical           :: exactb
  end type matB
  !----------------------------------------------------------------------!

#ifdef DQMC_PROFILE

#define PROF_BEGIN if (matb_profile) call get_time(matb_t1)
#define PROF_END(c, t) if (matb_profile) call matb_end(c, t)

    real(wp) :: MultB_Left_time = 0, MultBi_Left_time = 0 
    real(wp) :: MultB_Right_time =0, MultBi_Right_time = 0
    double precision :: matb_t1, matb_t2
    logical  :: matb_profile = .false.
    integer  :: MultB_Left_count = 0, MultBi_Left_count = 0 
    integer  :: MultB_Right_count = 0, MultBi_Right_count = 0

contains

    subroutine matb_end(c, t)
      integer  :: c
      real(wp) :: t
      call get_time(matb_t2)
      c = c + 1
      t = t + matb_t2 - matb_t1
    end subroutine matb_end

    subroutine matb_print()
      write(*,*) "MultB_Left    ", MultB_Left_count, MultB_Left_time
      write(*,*) "MultBi_Left   ", MultBi_Left_count, MultBi_Left_time
      write(*,*) "MultB_Right   ", MultB_Right_count, MultB_Right_time
      write(*,*) "MultBi_Right  ", MultBi_Right_count, MultBi_Right_time
    end subroutine matb_print

#else

#define PROF_BEGIN
#define PROF_END(c, t)

contains

#endif

  subroutine DQMC_B_Init(n, B, WS, Adj, ckb, t, mu, dtau)
    !
    ! Purpose
    ! =======
    !    This subroutine getnerates exponentional matrices of [-]dtau*T.
    !
    !         Ek  = exp(dtau*T)
    !         Eki = exp(-dtau*T)
    !
    !    The steps are
    !    1. Compute the Spectral Decomposition of T = USU^T
    !    2. Exponent the diagonal matrix S and -S
    !    3. Assemble B = U*exp(S)*U^T, Bi=U*exp(-S)*U^T
    !
    !    The parameters for checkboard method are also initialized
    !    in this subroutine. For the checkboard method, see [1] for
    !    more detail.
    !
    !    [1] W. Hanke and Yu.V. Kopaev, "Stable Numerical Simulation 
    !        of Models of Interating Electrons in Condensed-Matter Physics."
    !        Chapter 4. Elsevier Science Pub. 
    !
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)       :: n            ! Number of sites
    type(MatB), intent(inout) :: B            ! MatB
    type(WSpace), intent(inout), target  :: WS   ! shared working space
    type(CCS), intent(in)     :: adj          ! adjacent info
    type(CCS), intent(out)    :: ckb
    real(wp), intent(in)      :: t(:)         ! model parameter
    real(wp), intent(in)      :: mu(n), dtau  

    ! ... local scalars    ...
    integer  :: i, j                          ! iterator
    real(wp), pointer :: H(:,:) 
    integer, pointer  :: start(:) 
    integer, pointer  :: r(:) 
    integer, pointer  :: A(:) 

    ! ... Executable ...

    B%n      = n
    ckb%n    = 0
    ckb%nnz  = 0
    B%exactb = .false.

    ! test if allocated
    allocate(B%K(n,n))
    allocate(B%B(n,n))
    allocate(B%Bi(n,n))
    allocate(B%rtB(n,n))
    allocate(B%rtBi(n,n))

    B%name = "Dense MatB"

    ! Compute the B matrix
    H   => WS%R1

    H = ZERO
    start => Adj%cstart
    r     => Adj%row
    A     => Adj%A
    
    do i = 1, n
       do j = start(i), start(i+1)-1
          H(r(j),i) = t(A(j))*dtau
       end do
       H(i,i) = mu(i)*dtau
    end do

    B%K = H

    call DQMC_B_ExpInit(B, H, WS)

  end subroutine DQMC_B_Init

  !----------------------------------------------------------------------!

  subroutine DQMC_B_ExpInit(B, K, WS)
    !
    ! Purpose
    ! =======
    !    This subroutine getnerates exponentional matrices of [-]dtau*T.
    !
    !         Ek  = exp(dtau*T)
    !         Eki = exp(-dtau*T)
    !
    !    The steps are
    !    1. Compute the Spectral Decomposition of T = USU^T
    !    2. Exponent the diagonal matrix S and -S
    !    3. Assemble B = U*exp(S)*U^T, Bi=U*exp(-S)*U^T
    !
    !    The parameters for checkboard method are also initialized
    !    in this subroutine. For the checkboard method, see [1] for
    !    more detail.
    !
    !    [1] W. Hanke and Yu.V. Kopaev, "Stable Numerical Simulation 
    !        of Models of Interating Electrons in Condensed-Matter Physics."
    !        Chapter 4. Elsevier Science Pub. 
    !
    !
    ! Arguments
    ! =========
    !
    type(MatB), intent(inout) :: B              ! MatB
    real(wp), intent(inout)   :: K(:,:)         ! model parameter
    type(WSpace), intent(inout), target  :: WS     ! shared working space

    ! ... local scalars    ...
    integer  :: i, info, n       ! iterator
    real(wp), pointer :: W5(:) 
    real(wp), pointer :: W4(:) 

    ! ... Executable ...

    W5  => WS%R5
    W4  => WS%R7
    n   = B%n

    ! Compute the Spectral Decomposition of W1
    call lapack_dsyev('V', 'L', n, K, n, W5, W4, WS%lw(LA_SYEV), info)

    !! Check error message
    if (info .ne. 0) then
       call DQMC_Error("Error: from dsyev, dqmc_getek", info)
    end if
  
    ! Exponent the diagonal elements and scale the right eigenmatrix
    WS%R2 = K
    !! Scale the right eigenmatrix
    do i = 1, n
       call blas_dscal(n, dexp(W5(i)), WS%R2(:,i), 1)
    end do
    ! Assemble B
    call blas_dgemm('N', 'T', n, n, n, ONE, K, n, WS%R2, n, ZERO, B%B, n)


    ! For the inverse part, do the same thing
    WS%R2 = K
    do i = 1, n
       call blas_dscal(n, dexp(-W5(i)), WS%R2(:,i), 1)
    end do
    ! Assemble Bi
    call blas_dgemm('N', 'T', n, n, n, ONE, K, n, WS%R2, n, ZERO, B%Bi, n)


    ! square root of B
    WS%R2 = K
    !! Scale the right eigenmatrix
    do i = 1, n
       call blas_dscal(n, sqrt(dexp(W5(i))), WS%R2(:,i), 1)
    end do
    ! Assemble B
    call blas_dgemm('N', 'T', n, n, n, ONE, K, n, WS%R2, n, ZERO, B%rtB, n)


    ! square root of inv(B)
    WS%R2 = K
    !! Scale the right eigenmatrix
    do i = 1, n
       call blas_dscal(n, sqrt(dexp(-W5(i))), WS%R2(:,i), 1)
    end do
    ! Assemble B
    call blas_dgemm('N', 'T', n, n, n, ONE, K, n, WS%R2, n, ZERO, B%rtBi, n)

  end subroutine DQMC_B_ExpInit

  !----------------------------------------------------------------------!


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
    
    deallocate(B%B, B%Bi, B%rtB, B%rtBi)

  end subroutine DQMC_B_Free

  !----------------------------------------------------------------------!
  
  subroutine DQMC_MultB_Left(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine mulitplies matrix M by the matrix B_i from
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
    integer, intent(in)      :: n              ! The order of matrix M  
    real(wp), intent(inout)  :: M(n,n)         !  
    type(MatB), intent(in)   :: B              !
    real(wp), intent(in)     :: V_i(n)         !
    real(wp), intent(inout)  :: C(n,n)         ! working space   

    integer :: i
    integer :: lwork
    real(wp), pointer :: W1(:)   
    real(wp), pointer :: W2(:)   
    real(wp), pointer :: W3(:,:) 
    ! ... executable ...

    PROF_BEGIN
    if(B%exactb) then

       lwork = 3*n
       allocate(W1(lwork), W2(n), W3(n,n))

       ! Hamiltonian at time slice
       C = B%K
       do i = 1, n
         C(i,i) = C(i,i) + log(V_i(i))
       enddo

       ! Compute the Spectral Decomposition 
       call lapack_dsyev('V', 'L', n, C, n, W2, W1, lwork, i)

       !! Check error message
       if (i .ne. 0) then
          call DQMC_Error("Error: from dsyev, dqmc_getek", i)
       end if

       ! Compute M^T U
       call blas_dgemm('T', 'N', n, n, n, ONE, M, n, C, n, ZERO, W3, n)
  
       ! Compute M^T U D
       do i = 1, n
          call blas_dscal(n, dexp(W2(i)), W3(:,i), 1)
       end do

       ! Compute U (M^T U D)^T = U D U^T M = B M
       call blas_dgemm('N', 'T', n, n, n, ONE, C, n, W3, n, ZERO, M, n)

       deallocate(W1, W2, W3)

    else

       ! Copy M to C
       C = M

       ! M = B*C
       call blas_dgemm('N', 'N', n, n, n, ONE, B%B, n, C, n, ZERO, M, n)

       ! M = Vi*M 
       call DQMC_ScaleRow(n, M, V_i)

    endif
    PROF_END(MultB_Left_count, MultB_Left_time)    
    
  end subroutine DQMC_MultB_Left

  !-----------------------------------------------------------------------!

  subroutine DQMC_MultBi_Left(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine mulitplies matrix M by the matrix B_i from
    !    leftside.
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
    integer, intent(in)      :: n              ! The order of matrix M  
    real(wp), intent(inout)  :: M(n,n)         !  
    type(MatB), intent(in)   :: B              ! MatB 
    real(wp), intent(in)     :: V_i(n)         !
    real(wp), intent(inout)  :: C(n,n)         ! working space  

    integer :: i
    integer :: lwork
    real(wp), pointer :: W1(:)   
    real(wp), pointer :: W2(:)   
    real(wp), pointer :: W3(:,:) 
    ! ... executable ...

    PROF_BEGIN
    if(B%exactb) then

       lwork = 3*n
       allocate(W1(lwork), W2(n), W3(n,n))

       ! Hamiltonian at time slice
       C = B%K
       do i = 1, n
         C(i,i) = C(i,i) + log(V_i(i))
       enddo

       ! Compute the Spectral Decomposition 
       call lapack_dsyev('V', 'L', n, C, n, W2, W1, lwork, i)

       !! Check error message
       if (i .ne. 0) then
          call DQMC_Error("Error: from dsyev, dqmc_getek", i)
       end if

       ! Compute M^T U
       call blas_dgemm('T', 'N', n, n, n, ONE, M, n, C, n, ZERO, W3, n)
  
       ! Compute M^T U D
       do i = 1, n
          call blas_dscal(n, dexp(-W2(i)), W3(:,i), 1)
       end do

       ! Compute U (M^T U D)^T = U D U^T M = B M
       call blas_dgemm('N', 'T', n, n, n, ONE, C, n, W3, n, ZERO, M, n)

       deallocate(W1, W2, W3)

    else

       ! Multiply B from left-hand-side

       ! M = inv(V_i)*M
       call DQMC_ScaleRowInv(n, M, V_i)

       ! Copy M to C
       C = M

       ! M = B*C
       call blas_dgemm('N', 'N', n, n, n, ONE, B%Bi, n, C, n, ZERO, M, n)

    endif
    PROF_END(MultBi_Left_count, MultBi_Left_time)    

  end subroutine DQMC_MultBi_Left

  !----------------------------------------------------------------------!
  
  subroutine DQMC_MultB_Right(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine mulitplies matrix M by the matrix B_i from
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
    integer, intent(in)      :: n              ! The order of matrix M  
    real(wp), intent(inout)  :: M(n,n)         !  
    type(MatB), intent(in)   :: B              !
    real(wp), intent(in)     :: V_i(n)         !
    real(wp), intent(inout)  :: C(n,n)         ! working space   

    integer :: i
    integer :: lwork
    real(wp), pointer :: W1(:)   
    real(wp), pointer :: W2(:)   
    real(wp), pointer :: W3(:,:) 
    ! ... executable ...

    PROF_BEGIN
    if(B%exactb) then

       lwork = 3*n
       allocate(W1(lwork), W2(n), W3(n,n))

       ! Generate H at time slice
       C = B%K
       do i = 1, n
         C(i,i) = C(i,i) + log(V_i(i))
       enddo

       ! Compute the Spectral Decomposition of W1
       call lapack_dsyev('V', 'L', n, C, n, W2, W1, lwork, i)

       !! Check error message
       if (i .ne. 0) then
          call DQMC_Error("Error: from dsyev, dqmc_getek", i)
       end if

       ! Compute M U
       call blas_dgemm('N', 'N', n, n, n, ONE, M, n, C, n, ZERO, W3, n)
  
       ! Compute M U D
       do i = 1, n
          call blas_dscal(n, dexp(W2(i)), W3(:,i), 1)
       end do
       ! Compute M U D U^T = M B
       call blas_dgemm('N', 'T', n, n, n, ONE, W3, n, C, n, ZERO, M, n)
 
       deallocate(W1, W2, W3)
    
    else

       ! M = M*V_i
       call DQMC_ScaleCol(n, M, V_i)
       
       ! Copy M to C
       C = M
       
       ! M = C*B
       call blas_dgemm('N', 'N', n, n, n, ONE, C, n, B%B, n, ZERO, M, n)

    endif
    PROF_END(MultB_Right_count, MultB_Right_time)    
    
  end subroutine DQMC_MultB_Right

  !-----------------------------------------------------------------------!

  subroutine DQMC_MultBi_Right(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine mulitplies matrix M by the matrix B_i from
    !    rightside.
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
    integer, intent(in)      :: n              ! The order of matrix M  
    real(wp), intent(inout)  :: M(n,n)         !  
    type(MatB), intent(in)   :: B              ! MatB 
    real(wp), intent(in)     :: V_i(n)         !
    real(wp), intent(inout)  :: C(n,n)         ! working space  

    integer :: i
    integer :: lwork
    real(wp), pointer :: W1(:)   
    real(wp), pointer :: W2(:)   
    real(wp), pointer :: W3(:,:) 
    ! ... executable ...

    PROF_BEGIN
    if(B%exactb) then

       lwork = 3*n
       allocate(W1(lwork), W2(n), W3(n,n))

       ! Generate H at time slice
       C = B%K
       do i = 1, n
         C(i,i) = C(i,i) + log(V_i(i))
       enddo

       ! Compute the Spectral Decomposition of W1
       call lapack_dsyev('V', 'L', n, C, n, W2, W1, lwork, i)

       !! Check error message
       if (i .ne. 0) then
          call DQMC_Error("Error: from dsyev, dqmc_getek", i)
       end if

       ! Compute M U
       call blas_dgemm('N', 'N', n, n, n, ONE, M, n, C, n, ZERO, W3, n)
  
       ! Compute M U D
       do i = 1, n
          call blas_dscal(n, dexp(-W2(i)), W3(:,i), 1)
       end do

       ! Compute M U D U^T = M B
       call blas_dgemm('N', 'T', n, n, n, ONE, W3, n, C, n, ZERO, M, n)
 
       deallocate(W1, W2, W3)
    
    else

       ! Multiply B from right-hand-side
       ! Copy M to C
       C = M
       
       ! M = C*B
       call blas_dgemm('N', 'N', n, n, n, ONE, C, n, B%Bi, n, ZERO, M, n)
       
       ! M = M*V_i
       call DQMC_ScaleColInv(n, M, V_i)
    PROF_END(MultBi_Right_count, MultBi_Right_time)    
    
    endif
    
  end subroutine DQMC_MultBi_Right

  !-----------------------------------------------------------------------!

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
    ! ... executable ...

    ! Multiply B from left-hand-side

    ! Copy M to C
    C = M

    ! M = B*C
    call blas_dgemm('N', 'N', n, n, n, ONE, B%rtB, n, C, n, ZERO, M, n)

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
    ! ... executable ...

    ! Copy M to C
    C = M

    ! M = C*B
    call blas_dgemm('N', 'N', n, n, n, ONE, C, n, B%rtBi, n, ZERO, M, n)

  end subroutine DQMC_MultrtB0i_Right

  !-----------------------------------------------------------------------!

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
    type(MatB), intent(in)   :: B              !
    real(wp), intent(in)     :: V_i(n)         !
    real(wp), intent(inout)  :: C(n,n)         ! working space  

    integer :: i
    ! ... executable ...

    if(B%exactb) then

       M = ZERO
       do i = 1, n
         M(i,i) = ONE
       enddo
       call DQMC_MultB_Right(n, M, B, V_i, C)

    else

       C(n,n) = 0.d0
       M = B%B
       ! M = Vi*M 
       call DQMC_ScaleRow(n, M, V_i)

    endif

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
    type(MatB), intent(in)   :: B              ! MatB 
    real(wp), intent(in)     :: V_i(n)         !
    real(wp), intent(inout)  :: C(n,n)         ! working space  

    integer :: i
    ! ... executable ...

    if(B%exactb) then

       M = ZERO
       do i = 1, n
         M(i,i) = ONE
       enddo
       call DQMC_MultBi_Right(n, M, B, V_i, C)

    else

       ! ... executable ...

       C(n,n) = 0.d0
       M = B%Bi
       ! M = Vi*M 
       call DQMC_ScaleColInv(n, M, V_i)

    endif

  end subroutine DQMC_GetBi

  !----------------------------------------------------------------------!

end module DQMC_MATB
