module DQMC_SEQB
#include "dqmc_include.h"

  use _DQMC_MATB
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
     real(wp),   pointer :: V(:,:) 

     ! For UDT decomposition
     real(wp),   pointer :: U(:,:) 
     real(wp),   pointer :: D(:)   
     real(wp),   pointer :: T(:,:) 

     ! Working space
     real(wp),   pointer :: W1(:,:)
     real(wp),   pointer :: W2(:,:)
     real(wp),   pointer :: rw(:)  
     real(wp),   pointer :: tau(:) 
     integer,    pointer :: piv(:) 
     real(wp),   pointer :: lw(:)  
  end type SeqB
  
contains
  
  !----------------------------------------------------------------------!

  subroutine DQMC_SeqB_Init(n, L, nOrth, B, V, SB, WS)
    !
    ! Purpose
    ! =======
    ! This subroutine initializes the intermediate results
    !
    ! Arguments
    ! =========
    !
    integer, intent(in) :: n            ! order of matrix
    integer, intent(in) :: L            ! time slice
    integer, intent(in) :: nOrth        ! number of safe multiplication
    real(wp), intent(in):: V(:,:)       ! Hubbard-Stratonivich Field 
    type(MatB), intent(in)    :: B      ! Data structure of B matrix
    type(SeqB), intent(inout) :: SB     ! intermediate results
    
    ! ... Executable ...
    SB%n = n
    SB%L = L
    SB%nOrth = nOrth
    
    ! B matrix
    SB%B => B
    SB%V => V
    
    ! working spaces
    SB%U  => WS%R1
    SB%D  => WS%R5
    SB%T  => WS%R2
    SB%W1 => WS%R3
    SB%W2 => WS%R4
    SB%rw => WS%R7
    SB%tau=> WS%R6
    SB%piv=> WS%I1
    SB%lw => WS%lw

  end subroutine DQMC_SeqB_Init

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
    
  end subroutine DQMC_SeqB_Free

  !----------------------------------------------------------------------!

  subroutine DQMC_SeqB_Update(SB)
    !
    ! This subroutine initializes the intermediate results
    !
    type(SqeB), intent(in) :: SB         ! Data structure of B matrix

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
    real(wp), intent(inout)  :: U(:,:)
    real(wp), intent(inout)  :: D(:)
    real(wp), intent(inout)  :: T(:,:)
    real(wp), intent(inout)  :: W1(:,:)       ! R-factor in QR factor
    real(wp), intent(inout)  :: W2(:,:)       ! working array in pivoting
    real(wp), intent(inout)  :: rw(:)         ! working array in QR factor
    real(wp), intent(inout)  :: tau(:)        ! working array in QR factor
    integer,  intent(inout)  :: piv(:)        ! pivoting array in QRD 
    integer,  intent(in)     :: lwork(:)      ! working array in QR
    
     ! ... Local variables ...
    integer  :: info, i

    ! ... Executable ...
    
    !! Compute U = U*D
    call DQMC_ScaleCol(n, U, D, .false.)

    !! Initial parameters for dgeqp3
    piv = 0

    !! QR factorization with column pivoting
    call lapack_dgeqp3(n, n, U, n, piv, tau, rw, lwork(LA_GEQRF), info)

    if (info .ne. 0) then
       call DQMC_Error("Error: dgeqp3 in dqmc_UDTD.", info)
    end if

    !! dgegp3 returns R-factor on the upper triangle of G.
    !! The lower triangle of G stores "reflectors",
    !! which is used to reconstruct the Q-factor.

    !! Move R-factor to R, and normalize it by diagonal elements.
    W1 = U
    do i = 1, n
       !! make T upper triangular.
       W1(i,1:i-1) = ZERO
       !! D = diag(T).
       D(i) = W1(i,i)
       if(D(i) .eq. ZERO) then
          call DQMC_Error("Error: R-factor is singular: dqmc_UDTD.", &
               i)
       else
          !! Normalize R's columns by its diagonal. R = inv(D)*R 
          call blas_dscal(n-i+1, ONE/D(i), W1(i,i), n)
       endif
    end do

    !! Compute V = P*V. (W is used as an temporary variable.)
    do i = 1, n
       W2(i,1:n) = T(piv(i), 1:n)
    end do

    !! Compute V = R*W = R*P*V
    call blas_dtrmm('L', 'U', 'N', 'U', n, n, ONE, W1, n, W2, n)

    T = W2

    ! Generate Q-factor
    call lapack_dorgqr(n, n, n, U, n, tau, rw, lwork(LA_ORGQR), info)
    
  end subroutine DQMC_UDTD

  !----------------------------------------------------------------------!

  subroutine DQMC_SeqMultB(il, ir, SB)
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

    ! ... local scalars    ...
    integer :: info                     ! parameters for lapack's sub
    integer :: i, si, interval          ! iterator
    
    ! ... Executable ...

    !! Initially, Q = B_{i} = V_i*B

    ! computing the interval between i1 and i2
    if (i2 .ge. i1) then
       interval = i2 - i1 + 1
    else
       interval = i2 + L - i1 + 1
    end if

    si = i1
    if (si .gt. L) si = 1
    if (si .le. 0) si = L
    
    ! Let U be B_{i1}
    U = B%B
    call DQMC_ScaleRow(n, U, V(:,si), .false.)
    !! R = I, R will be the R-factor of the QDR factorization
    call DQMC_Eye(n, T)
    D = ONE
    
    ! Loop over the rest B_i
    do i = 1, interval 
       !! Compute the index of B_{i}
       si = si + 1
       if (si .gt. L) si = 1
       
       !! The UDT decomposition is performed at every nOrth step, and
       !! at the last step. In other steps, we just multiply B_i to 
       !! the Q-factor     
       if (mod(i, B%nOrth) .eq. 0 .or. i .eq. interval) then
          ! UDV decomposition
          call DQMC_UDTD(n, U, D, T, W1, W2, rw, tau, piv, lwork)
          
          ! Apply UDV
          if (i .eq. interval) then
             !! The last step, just generate U-factor
             call lapack_dorgqr(n, n, n, U, n, tau, rw, lwork(LA_ORGQR), info)
          else
             !! Compute U = B_i*U = (V_i*B)*U.
             !! (W1 is used as an temporary variable.)
             W1 = B%B
             call lapack_dormqr('R', 'N', n, n, n, U, n, tau, W1, n, rw, &
                  lwork(LA_ORMQR), info)
             U  = W1
             call DQMC_ScaleRow(n, U, V(:,si), .false.)
          end if
          
       else
          !! The other cases: Just multiply B_i to the Q-factor
          call DQMC_MultB(n, U, B, V(:,si), 'left', W1)
       end if
       
    end do

  end subroutine DQMC_SeqMultB

  !----------------------------------------------------------------------!

  subroutine DQMC_SeqMultBi(n, L, i1, i2, &
                           B, V, U, D, T, W1, W2, rw, lwork, tau, piv)
    !
    ! Purpose
    ! =======
    !    This subroutine computes A = inv(B_{i2}B_{i2-1}...B_{i1})
    !    and returns A's UDT decomposition.
    !
    ! Pre-assumption
    ! ==============
    !   i1 is the index of right most B, adn i2 is the index for the
    !   left most B. Both i1, i2 are in [1,..L].
    !
    ! Arguments
    ! =========
    !
    integer,  intent(in)    :: n, L            ! parameter of Hub
    integer,  intent(in)    :: i1, i2          ! start/end slice 
    type(MatB), intent(in)  :: B               ! MatB 
    real(wp), intent(in)    :: V(n,L)          ! 
    real(wp), intent(inout) :: U(n,n)          ! U-factor
    real(wp), intent(inout) :: D(n)            ! D-factor
    real(wp), intent(inout) :: T(n,n)          ! T-factor
    real(wp), intent(inout) :: W1(n,n)         ! working space
    real(wp), intent(inout) :: W2(n,n)         ! working space
    real(wp), intent(inout) :: rw(:)           ! working space
    real(wp), intent(inout) :: tau(n)          ! working space
    integer,  intent(in)    :: lwork(:)        !
    integer,  intent(inout) :: piv(n)

    ! ... local scalars    ...
    integer :: info                     ! parameters for lapack's sub
    integer :: i, si, interval          ! iterator
    
    ! ... Executable ...

    !! Initially, Q = B_{i} = V_i*B

    ! computing the interval between i1 and i2
    if (i2 .ge. i1) then
       interval = i2 - i1 + 1
    else
       interval = i2 + L - i1 + 1
    end if

    si = i2
    if (si .gt. L) si = 1
    if (si .le. 0) si = L
    
    ! Let U be B_{i2}^{-1}
    U = B%Bi
    call DQMC_ScaleCol(n, U, V(:,si), .true.)
    !! R = I, R will be the R-factor of the QDR factorization
    call DQMC_Eye(n, T)
    D = ONE
    
    ! Loop over the rest B_i
    do i = 1, interval 
       si = si - 1
       if (si .le. 0) si = L
       
       !! The UDT decomposition is performed at every nOrth step, and
       !! at the last step. In other steps, we just multiply B_i to 
       !! the Q-factor     
       if (mod(i, B%nOrth) .eq. 0 .or. i .eq. interval) then
          ! UDV decomposition
          call DQMC_UDTD(n, U, D, T, W1, W2, rw, tau, piv, lwork)
          
          ! Apply UDV
          if (i .eq. l) then
             !! The last step, just generate U-factor
             call lapack_dorgqr(n, n, n, U, n, tau, rw, lwork(LA_ORGQR), info)
          else
             !! Compute U = B_i*U = (V_i*B)*U.
             !! (R is used as an temporary variable.)
             W1 = B%Bi
             call DQMC_ScaleCol(n, W1, V(:,si), .true.)
             call lapack_dormqr('R', 'N', n, n, n, U, n, tau, W1, n, rw, &
                  lwork(LA_ORMQR), info)
             U  = W1
          end if
          
       else
          !! The other cases: Just multiply B_i to the Q-factor
          call DQMC_MultBi(n, U, B, V(:,si), 'left', W1)
       end if
       
    end do
       
  end subroutine DQMC_SeqMultBi
  
  !----------------------------------------------------------------------!

end module DQMC_SEQB
