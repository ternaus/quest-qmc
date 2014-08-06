module DQMC_SEQB
#include "dqmc_include.h"

#ifdef DQMC_CKB
  use DQMC_CheckerBoard
#else 
  use DQMC_MATB
#endif

  use DQMC_UTIL
  use DQMC_WSPACE

  implicit none 

  ! 
  ! This module implements multiplication of a sequent B. 
  !
  ! Data Types
  ! ==========
  !
  !  UDT decomposition
  !

  type SeqB
     integer  :: n
     integer  :: L
     integer  :: nOrth
     
     ! B matrix
     type(MatB), pointer :: B      

     ! For UDT decomposition
     real(wp),   pointer :: U(:,:) 
     real(wp),   pointer :: D(:)   
     real(wp),   pointer :: T(:,:) 

     ! Working space
     real(wp),   pointer :: W1(:,:)
     real(wp),   pointer :: W2(:,:)
     real(wp),   pointer :: W3(:,:)
     real(wp),   pointer :: rw(:)  
     real(wp),   pointer :: tau(:) 
     integer,    pointer :: piv(:) 
     integer,    pointer :: lw(:)  

  end type SeqB
  
contains
  
  !----------------------------------------------------------------------!

  subroutine DQMC_SeqB_Init(n, L, nOrth, B, SB, WS)
    !
    ! Purpose
    ! =======
    ! This subroutine initializes the intermediate results
    !
    ! Arguments
    ! =========
    !
    integer, intent(in) :: n                 ! order of matrix
    integer, intent(in) :: L                 ! time slice
    integer, intent(in) :: nOrth             ! number of safe multiplication
    type(MatB), intent(in), target   :: B    ! Data structure of B matrix
    type(SeqB), intent(inout) :: SB          ! intermediate results
    type(WSpace), intent(in), target :: WS   ! intermediate results
    
    ! ... Executable ...
    SB%n = n
    SB%L = L
    SB%nOrth = nOrth
    
    ! B matrix
    SB%B => B
    
    ! working spaces
    SB%U  => WS%R1
    SB%T  => WS%R2
    SB%W1 => WS%R3
    SB%W2 => WS%R4
    SB%D  => WS%R5
    SB%tau=> WS%R6
    SB%rw => WS%R7
    SB%W3 => WS%R8
    SB%piv=> WS%I1
    SB%lw => WS%lw

  end subroutine DQMC_SeqB_Init

  !----------------------------------------------------------------------!

  subroutine DQMC_SeqB_Init2(n, L, nOrth, B, SB, U, D, T, WS)
    !
    ! Purpose
    ! =======
    ! This subroutine initializes the intermediate results
    !
    ! Arguments
    ! =========
    !
    integer, intent(in) :: n                 ! order of matrix
    integer, intent(in) :: L                 ! time slice
    integer, intent(in) :: nOrth             ! number of safe multiplication
    real(wp), intent(in), target   :: U(:,:) ! 
    real(wp), intent(in), target   :: D(:)   ! 
    real(wp), intent(in), target   :: T(:,:) ! 
    type(MatB), intent(in), target   :: B    ! Data structure of B matrix
    type(SeqB), intent(inout) :: SB          ! intermediate results
    type(WSpace), intent(in), target :: WS   ! intermediate results
    
    ! ... Executable ...
    SB%n = n
    SB%L = L
    SB%nOrth = nOrth
    
    ! B matrix
    SB%B => B
    
    ! working spaces
    SB%U  => U
    SB%D  => D
    SB%T  => T
    SB%W1 => WS%R3
    SB%W2 => WS%R4
    SB%rw => WS%R7
    SB%tau=> WS%R6
    SB%piv=> WS%I1
    SB%lw => WS%lw

  end subroutine DQMC_SeqB_Init2

  !----------------------------------------------------------------------!

  subroutine DQMC_SeqB_Free(SB)
    !  
    ! Purpose
    ! =======
    !    This subroutine frees dynamically allocated memory of SB.
    !
    ! Arguments
    ! =========
    !
    type(SeqB), intent(inout) :: SB  ! intermediate results
    SB%n = 0

  end subroutine DQMC_SeqB_Free

  !----------------------------------------------------------------------!

  subroutine DQMC_SeqB_Update(SB)
    !
    ! This subroutine initializes the intermediate results
    !
    type(SeqB), intent(inout) :: SB         ! Data structure of B matrix
    SB%n = 0

  end subroutine DQMC_SeqB_Update

  !----------------------------------------------------------------------!

  subroutine DQMC_UDTD(n, U, D, T, W1, W2, rw, tau, piv, lwork)
    ! 
    ! Purpose
    ! =======
    !    This subroutine computes (updates) UDT decomposition,
    !             
    !       A = UDT
    !
    !    where U is orthonormal, D is diagonal, and T is normalized
    !    is some way.
    !   
    !    In input, U is not orthonormal. Therefore,
    !
    !    1. QR-decomposition with pivoting on U
    !
    !          [U, R, P] = QRP(U*D)
    !
    !    2. Normalize R by its diagonal elements and set them to D.
    !
    !          R = D*R
    !
    !    3. Apply P to T,  W = P*T.
    !
    !    4. Multiply R to W to get the new T, T = R*W = R*P*T.
    !    
    !
    ! Arguments
    ! =========
    ! 
    integer,  intent(in)     :: n
    real(wp), intent(inout)  :: U(n,n)
    real(wp), intent(inout)  :: D(n)
    real(wp), intent(inout)  :: T(n,n)
    real(wp), intent(inout)  :: W1(n,n)       ! R-factor in QR factor
    real(wp), intent(inout)  :: W2(n,n)       ! working array in pivoting
    real(wp), intent(inout)  :: rw(:)         ! working array in QR factor
    real(wp), intent(inout)  :: tau(n)        ! working array in QR factor
    integer,  intent(inout)  :: piv(:)        ! pivoting array in QRD 
    integer,  intent(in)     :: lwork(:)      ! working array in QR
    
     ! ... Local variables ...
    integer  :: info, i

    ! ... Executable ...
    
    !! Compute U = U*D
    call DQMC_ScaleCol(n, U, D)

    !! Initial parameters for dgeqp3
    piv = 0
    info = 0

    !! QR factorization with column pivoting
    call lapack_dgeqp3(n, n, U, n, piv, tau, rw, lwork(LA_GEQRF), info)

    if (info .ne. 0) then
       call DQMC_Error("Error: dgeqp3 in dqmc_UDTD.", info)
    end if

    !! dgegp3 returns R-factor on the upper triangle of G.
    !! The lower triangle of G stores "reflectors",
    !! which is used to reconstruct the Q-factor.

    !! Move R-factor to R, and normalize it by diagonal elements.
    ! W1 = U
    call blas_dcopy(n*n, U, 1, W1, 1)

    do i = 1, n
       !! make T upper triangular.
       W1(i,1:i-1) = ZERO
       !! D = diag(T).
       D(i) = W1(i,i)
       if(D(i) .eq. ZERO) then
          call DQMC_Error("Error: R-factor is singular: dqmc_UDTD.", i)
       else
          !! Normalize R's row by its diagonal. R = inv(D)*R 
          call blas_dscal(n-i+1, ONE/D(i), W1(i,i), n)
       endif
    end do

    !! Compute V = P*V. (W is used as an temporary variable.)
    do i = 1, n
       W2(i,1:n) = T(piv(i), 1:n)
    end do

    !! Compute V = R*W = R*P*V
    call blas_dtrmm('L', 'U', 'N', 'U', n, n, ONE, W1, n, W2, n)

    ! T = W2
    call blas_dcopy(n*n, W2, 1, T, 1)

    ! Generate Q-factor
    call lapack_dorgqr(n, n, n, U, n, tau, rw, lwork(LA_ORGQR), info)

  end subroutine DQMC_UDTD

  !----------------------------------------------------------------------!

  subroutine DQMC_SeqMultB(il, ir, SB, V)
    !
    ! Purpose
    ! =======
    !    This subroutine computes A = B_{il}B_{il-1}...B_{ir}
    !    and returns A's UDT decomposition.
    !
    !    ir is the index of right most B, adn il is the index for the
    !    left most B. 
    !
    ! Arguments
    ! =========
    !
    integer,  intent(in)      :: il, ir           ! start/end slice 
    type(SeqB), intent(inout) :: SB               ! SeqB
    real(wp), intent(in)      :: V(SB%n,SB%L)


    ! ... local scalars    ...
    integer :: i, si, interval          ! iterator
    integer :: n, L                     ! alias 
    real(wp), pointer :: U(:,:) 
    real(wp), pointer :: T(:,:) 
    real(wp), pointer :: D(:)   
    real(wp), pointer :: W1(:,:)
    real(wp), pointer :: W2(:,:)
    real(wp), pointer :: tau(:) 
    real(wp), pointer :: rw(:)  
    integer,  pointer :: piv(:) 
    integer,  pointer :: lw(:)  

    ! ... Executable ...

    U   => SB%U
    D   => SB%D
    T   => SB%T
    W1  => SB%W1
    W2  => SB%W2
    tau => SB%tau

    n   = SB%n
    L   = SB%L

    rw  => SB%rw
    lw  => SB%lw
    piv => SB%piv

    !! Initially, Q = B_{i} = V_i*B

    ! computing the interval between i1 and i2
    if (il .ge. ir) then
       interval = il - ir + 1
    else
       interval = il + L - ir + 1
    end if

    si = ir
    if (si .gt. L) si = si - L
    
    ! Let U be B_{i1}
    call DQMC_GetB(n, U, SB%B, V(:,si), SB%W1)

    !! T = I, T will be the R-factor of the QDR factorization
    call DQMC_Eye(n, T)
    D = ONE

    ! Loop over the rest B_i
    do i = 1, interval - 1
       !! Compute the index of B_{i}
       si = si + 1
       if (si .gt. L) si = 1
  
       !! The UDT decomposition is performed at every nOrth step, and
       !! at the last step. In other steps, we just multiply B_i to 
       !! the Q-factor     
       if (mod(i, SB%nOrth) .eq. 0) then
          ! UDV decomposition
          call DQMC_UDTD(n, U, D, T, W1, W2, rw, tau, piv, lw)
       end if

       !! multiply B_i to the Q-factor
       call DQMC_MultB_Left(n, U, SB%B, V(:,si), SB%W1)
    end do

    ! before leave, make the decomposition form
    call DQMC_UDTD(n, U, D, T, W1, W2, rw, tau, piv, lw) 

  end subroutine DQMC_SeqMultB

  !----------------------------------------------------------------------!

  subroutine DQMC_SeqMultBi(il, ir, SB, V)
    !
    ! Purpose
    ! =======
    !    This subroutine computes A = inv(B_{il}B_{il-1}...B_{ir})
    !    and returns A's UDT decomposition.
    !
    ! Pre-assumption
    ! ==============
    !   ir is the index of right most B, adn il is the index for the
    !   left most B. Both i1, i2 are in [1,..L].
    !
    ! Arguments
    ! =========
    !
    integer,  intent(in)    :: il, ir          ! start/end slice 
    type(SeqB), intent(inout)  :: SB           ! MatB
    real(wp), intent(in)      :: V(SB%n,SB%L) 

    ! ... local scalars    ...
    integer :: i, si, interval          ! iterator
    integer :: n, L                     ! alias 
    real(wp), pointer :: U(:,:) 
    real(wp), pointer :: T(:,:) 
    real(wp), pointer :: D(:)   
    real(wp), pointer :: W1(:,:)
    real(wp), pointer :: W2(:,:)
    real(wp), pointer :: tau(:) 
    real(wp), pointer :: rw(:)  
    integer,  pointer :: piv(:) 
    integer,  pointer :: lw(:)  

    ! ... Executable ...

    U   => SB%U
    D   => SB%D
    T   => SB%T
    W1  => SB%W1
    W2  => SB%W2
    tau => SB%tau

    n  = SB%n
    L  = SB%L

    rw  => SB%rw
    lw  => SB%lw
    piv => SB%piv

    ! computing the interval between i1 and i2
    if (il .ge. ir) then
       interval = il - ir + 1
    else
       interval = il + L - ir + 1
    end if

    si = il
    if (si .gt. L) si = si - L

    ! Let U be B_{i2}^{-1}
    call DQMC_GetBi(n, U, SB%B, V(:,si), SB%W1)

    !! R = I, R will be the R-factor of the QDR factorization
    call DQMC_Eye(n, T)
    D = ONE
    
    ! Loop over the rest B_i
    do i = 1, interval - 1
       si = si - 1
       if (si .le. 0) si = L
       
       !! The UDT decomposition is performed at every nOrth step, and
       !! at the last step. In other steps, we just multiply B_i to 
       !! the Q-factor     
       if (mod(i, SB%nOrth) .eq. 0) then
          ! UDV decomposition
          call DQMC_UDTD(n, U, D, T, W1, W2, rw, tau, piv, lw)
       end if

       !! multiply B_i to the Q-factor
       call DQMC_MultBi_Left(n, U, SB%B, V(:,si), SB%W1)
    end do

    ! before leave, make the decomposition form
    call DQMC_UDTD(n, U, D, T, W1, W2, rw, tau, piv, lw) 

  end subroutine DQMC_SeqMultBi
  
  !----------------------------------------------------------------------!

end module DQMC_SEQB
