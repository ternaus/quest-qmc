module DQMC_GTAU
#include "dqmc_include.h"

  use DQMC_UTIL
  use DQMC_WSPACE
#ifdef DQMC_CKB
  use DQMC_CheckerBoard
#else
  use DQMC_MATB
#endif
  use DQMC_SEQB
  use DQMC_GFUN

  implicit none 
  
  !
  ! This module is designed for the computation of time dependent 
  ! measurement (TDM), which requires the unequal time Green's 
  ! function Gup_tau (Gdn_tau).
  ! Mathematically, Gup_tau(Gdn_tau) is the inverse of the following
  ! matrix.
  !
  !                  [  I                   B_L] 
  !                  [-B_1   I                 ]
  !              M = [     -B_2  I             ]
  !                  [          ...   ...      ]
  !                  [             -B_{L-1}  I ]
  !                   
  !
  ! The current implementation only considers to return the ith row
  ! of Gup_tau (Gdn_tau). Let Gup_ij be the (i,j) block matrix of 
  ! Gup_tau.
  !
  !     Gup_ii = inv(I+B_iB_{i-1}...B_1B_L...B_{i+1})   
  !     Gup_ij = -Gup_ii(B_iB_{i-1}...B_1B_L...B_{j+1})  for j = i...L  
  !     Gup_ij =  Gup_ii(B_iB_{i-1}...B_1B_L...B_{j+1})  for j = 1...i-1  
  !
  ! The corresponding Gdn_ij has the same structure.
  !
  ! [1] Z. Bai, W.Chen, R. Scalettar, I. Yamazaki, "Lecture Notes 
  !     on Advances of Numerical Methods for Hubbard Quantum Monte
  !     Carlo Simulation." 
  !
  type Gtau
     integer  :: n       ! number of sites
     integer  :: L       ! number of imag-times
     integer  :: nb      ! number of imag-times in A_up, A_dn
     integer  :: nnb     ! size of A_up, A_dn
     integer  :: it_up   ! time indices for upt0, up0t, dnt0 and dn0t
     integer  :: i0_up   ! time indices for upt0, up0t, dnt0 and dn0t
     integer  :: it_dn   ! time indices for upt0, up0t, dnt0 and dn0t
     integer  :: i0_dn   ! time indices for upt0, up0t, dnt0 and dn0t
     integer  :: north   ! Number of imag-times
     integer  :: which   ! part of gtau to compute
     integer  :: sfc     ! safe count

     ! Green's functions
     real(wp), pointer :: upt0(:,:)  ! Gup(t,0)
     real(wp), pointer :: up0t(:,:)  ! Gup(0,t)
     real(wp), pointer :: dnt0(:,:)  ! Gdn(t,0)
     real(wp), pointer :: dn0t(:,:)  ! Gdn(0,t)
     real(wp), pointer :: up00(:,:)  ! Gup(0,0)
     real(wp), pointer :: uptt(:,:)  ! Gup(t,t)
     real(wp), pointer :: dn00(:,:)  ! Gdn(0,0)
     real(wp), pointer :: dntt(:,:)  ! Gdn(t,t)

     ! Full Green's functions (nb imag-times. Not L!)
     real(wp), pointer :: A_up(:,:) 
     real(wp), pointer :: A_dn(:,:) 
     ! Indices of imag-times stored in the A's
     integer, pointer  :: itau_up(:) 
     integer, pointer  :: itau_dn(:) 
     
     ! Pointers to B matrices (no fields)
     type(MatB), pointer :: B_up 
     type(MatB), pointer :: B_dn 
     ! Pointers to B matrices with fields
     real(wp), pointer   :: V_up(:,:) 
     real(wp), pointer   :: V_dn(:,:) 
 
     ! Pointers to phase
     real(wp), pointer   :: P(:) 

     ! Control variables for dn computation
     logical :: comp_dn
     logical :: neg_u
    
     ! Signs
     real(wp), pointer   :: sgnup 
     real(wp), pointer   :: sgndn 

     ! Workspace variables
     integer             :: lw
     integer, pointer    :: IW(:)   
     real(wp), pointer   :: W1(:)   
     real(wp), pointer   :: W2(:,:) 
     real(wp), pointer   :: W3(:,:) 

     ! Workspace for g0 computation
     logical           :: g0_stored = .false.
     real(wp), pointer :: e0up(:)   
     real(wp), pointer :: e0dn(:)   
     real(wp), pointer :: U0up(:,:) 
     real(wp), pointer :: U0dn(:,:) 

  end type Gtau
  
  integer, parameter :: TAU_T0   =  0   ! Column
  integer, parameter :: TAU_BOTH =  1   ! column and row
  integer, parameter :: TAU_0T   =  2   ! ROW

  integer, parameter :: TAU_UP   =  1   ! Spin up
  integer, parameter :: TAU_DN   = -1   ! Spin down

  integer, parameter :: TPLUS    =  1  ! Spin up
  integer, parameter :: TMINUS   =  2  ! Spin down
  integer, parameter :: ZPLUS    =  3  ! Spin up
  integer, parameter :: ZMINUS   =  4  ! Spin down
contains

  ! Subroutines
  ! ==================================================================

  subroutine DQMC_Gtau_Init(Hub, tau)
    use DQMC_HUBBARD
    !
    ! Purpose
    ! =======
    !    This subroutine initializes tau.
    !
    ! Arguments
    ! =========
    !
    type(Gtau), intent(inout)        :: tau     
    type(Hubbard), target, intent(in)    :: Hub

    integer  :: n, nnb, info
    real(wp) :: query(1)

    ! ... Executable ...
    tau%it_up = 0
    tau%i0_up = 0
    tau%it_dn = 0
    tau%i0_dn = 0

    ! Copy variable from Hub
    tau%n       = Hub%n
    tau%L       = Hub%L
    tau%north   = Hub%SB_up%north
    tau%nb      = Hub%L / tau%north
    tau%comp_dn = Hub%comp_dn
    tau%neg_u   = Hub%neg_u

    if (mod(tau%L,tau%north) /= 0) &
       call dqmc_error("L must be an exact multiple of north. Stop.",1)

    tau%which = TAU_BOTH

    ! Allocate space for measurement
    tau%nnb = tau%n * tau%nb
    n   = tau%n
    nnb = tau%nnb
    allocate(tau%upt0(n, n))
    allocate(tau%up0t(n, n))
    allocate(tau%up00(n, n))
    allocate(tau%uptt(n, n))
    allocate(tau%A_up(nnb, nnb))
    if (tau%comp_dn .or. .not.tau%neg_u) then
       allocate(tau%dnt0(n, n))
       allocate(tau%dn0t(n, n))
       allocate(tau%dn00(n, n))
       allocate(tau%dntt(n, n))
       allocate(tau%A_dn(nnb, nnb))
    else
       tau%dnt0    => tau%upt0
       tau%dn0t    => tau%up0t
       tau%dn00    => tau%up00
       tau%dntt    => tau%uptt
       tau%A_dn    => tau%A_up
    endif

    ! sgn and itaus are identical when .not.comp_dn
    allocate(tau%sgnup)
    allocate(tau%itau_up(tau%nb))
    if (tau%comp_dn) then
       allocate(tau%sgndn)
       allocate(tau%itau_dn(tau%nb))
    else
       tau%sgndn   => tau%sgnup
       tau%itau_dn => tau%itau_up
    endif

    tau%B_up => Hub%B_up
    tau%B_dn => Hub%B_dn
    tau%V_up => Hub%V_up
    tau%V_dn => Hub%V_dn

    tau%P    => Hub%S%P

    ! Create working space 
    allocate(tau%IW(nnb))
    call lapack_dgetri(nnb, tau%A_up, nnb, tau%IW, query, -1, info)
    tau%lw = nint(query(1))
    allocate(tau%W1(tau%lw))
    allocate(tau%W2(n,n))
    allocate(tau%W3(n,n))

    tau%g0_stored = .false.

  end subroutine DQMC_Gtau_Init

  !--------------------------------------------------------------------!

  subroutine DQMC_Gtau_Free(tau)
    !
    ! Purpose
    ! =======
    !    This subroutine frees TDM.
    !
    ! Arguments
    ! =========
    !
    type(Gtau), intent(inout) :: tau      ! TDM to be initialized

    ! ... Executable ...

    deallocate(tau%A_up)
    deallocate(tau%upt0, tau%up0t)
    deallocate(tau%up00, tau%uptt)
    deallocate(tau%itau_up)
    deallocate(tau%sgnup)

    if (associated(tau%A_dn)) then
       deallocate(tau%A_dn)
       deallocate(tau%dnt0, tau%dn0t)
       deallocate(tau%dn00, tau%dntt)
    endif
    if (associated(tau%sgndn)) then
       deallocate(tau%itau_dn)
       deallocate(tau%sgndn)
    endif

    deallocate(tau%W1, tau%IW)

    if (tau%g0_stored) then
       deallocate(tau%e0up)
       deallocate(tau%e0dn)
       deallocate(tau%U0up)
       deallocate(tau%U0dn)
    endif
 
  end subroutine DQMC_Gtau_Free

  !--------------------------------------------------------------------!

  subroutine DQMC_Gtau_LoadA(tau, spin, slice, sgn)
    !
    ! This solves the gtau explicitly by Lapack.
    !

    type(gtau), target, intent(inout) :: tau
    integer,    intent(in)            :: spin
    integer,    intent(in)            :: slice
    real(wp),   intent(out)           :: sgn
    
    ! ... Local var ...
    real(wp),   pointer :: A(:,:)     
    real(wp),   pointer :: V(:,:)     
    type(MatB), pointer :: B          
    real(wp),   pointer :: s          
    integer,    pointer :: t(:)       
    real(wp),   pointer :: gtau1(:,:) 
    real(wp),   pointer :: gtau2(:,:) 

    integer  :: i, j, k, h, isl, jsl
    integer  :: n, nb, nnb, nor, L
    integer  :: info 

    ! ... Executable ...
    
    n     = tau%n
    nb    = tau%nb
    nnb   = tau%nnb
    nor   = tau%north
    L     = tau%L
    isl   = mod(slice-1,nor)+1

    if (spin==TAU_UP .or. tau%comp_dn) then

       call dqmc_gtau_setAlias(spin, tau, A=A, V=V, B=B, sgn=s, itau=t)
       ! making A
       A = ZERO
       do i = 1, nnb
          A(i,i) = ONE
       end do
       ! Compute B products:
       !    i=0  B_{isl}BB...B_{isl-north+1}
       !    i=1  B_{isl+north}BB...B_{isl+1}
       !    i=2  B_{isl+2*north}BB...B_{isl+north+1}
       !    ....
       do i = 0, nb - 1
          ! Pointers to sub-diagonal blocks (except when i=0)
          jsl  =  mod(isl+i*nor-1, L) + 1
          t(i+1) = jsl
          call DQMC_GetB(n, tau%W3, B, V(:,jsl), tau%W2)
          do j = 2, nor
             k = mod(jsl-j+L, L) + 1
             call DQMC_MultB_Right(n, tau%W3, B, V(:,k), tau%W2)
          enddo
          ! Subdiagonal blocks need negative sign
          if (i > 0) tau%W3 = -tau%W3
          j = mod(i+nb-1, nb) 
          A(i*n+1:(i+1)*n, j*n+1:(j+1)*n) = tau%W3
       end do
       ! Inversion
       call lapack_dgetrf(nnb, nnb, A, nnb, tau%IW, info)
       s = ONE
       do i = 1, nnb
         if (tau%IW(i)/= i)  s = -s
         if (A(i,i) < ZERO) s = -s
       enddo
       call lapack_dgetri(nnb, A, nnb, tau%IW, tau%W1, tau%lw, info)

    elseif (.not.tau%neg_u) then

       ! Use p-h symmetry to fill A_dn
       s => tau%sgndn
       do i = 0, nb - 1
          do j = 0, nb - 1
             gtau1  => tau%A_up(i * n + 1: (i + 1) * n, j * n + 1:(j + 1) * n)
             gtau2 => tau%A_dn(j * n + 1: (j + 1) * n, i * n + 1:(i + 1) * n)
             do h = 1, n
                do k = 1, n
                   gtau2(h,k) = -tau%P(k) * tau%P(h) * gtau1(k,h)
                enddo
             enddo
          enddo
       enddo
       ! Equal time, equal site entries need correction
       do i = 1, nnb
          tau%A_dn(i,i) = tau%A_dn(i,i) + ONE
       enddo

    else
       ! Negative U case with identical non-interacting parts
       s => tau%sgndn
    endif

    sgn = s

  end subroutine DQMC_Gtau_LoadA

  !--------------------------------------------------------------------!

  subroutine DQMC_Gtau_CopyUp(tau)
    !
    ! Load content of A in upt0, up0t, dnt0 and dn0t.
    ! Equal time upt0=G_up(t+,t), up0t=G_up(t,t+) 
    ! 
    type(Gtau), intent(inout) :: tau

    integer :: h, k
    real(wp) :: PP

    do h = 1, tau%n
       do k = 1, tau%n
          PP = -tau%P(k)*tau%P(h)
          tau%dnt0(h,k) = PP*tau%up0t(k,h)
          tau%dn0t(h,k) = PP*tau%upt0(k,h)
          tau%dn00(h,k) = PP*tau%up00(k,h)
          tau%dntt(h,k) = PP*tau%uptt(k,h)
       enddo
       tau%dn00(h,h) = ONE + tau%dn00(h,h)
       tau%dntt(h,h) = ONE + tau%dntt(h,h)
    enddo
    tau%it_dn = tau%it_up
    tau%i0_dn = tau%i0_up

  end subroutine DQMC_Gtau_CopyUp
 
  !--------------------------------------------------------------------!

  subroutine DQMC_Gtau_DumpA(tau, spin, it, i0)
    !
    ! Load content of A in upt0, up0t, dnt0 and dn0t.
    ! Equal time upt0=G_up(t+,t), up0t=G_up(t,t+) 
    ! 
    type(Gtau), intent(inout) :: tau
    integer, intent(in)       :: it, i0, spin

    integer  :: n, i, j0, jt
    integer,  pointer :: itptr    
    integer,  pointer :: i0ptr    
    integer,  pointer :: t(:)     
    real(wp), pointer :: gt0(:,:) 
    real(wp), pointer :: g0t(:,:) 
    real(wp), pointer :: g00(:,:) 
    real(wp), pointer :: gtt(:,:) 
    real(wp), pointer :: A(:,:)   

    n = tau%n

    call dqmc_gtau_setAlias(spin, tau, A=A, gt0=gt0, g0t=g0t, g00=g00, gtt=gtt, it=itptr, i0=i0ptr, itau=t)
    
    ! Save current slice in tau
    itptr = t(it)
    i0ptr = t(i0)

    ! Fill gt0 and g0t
    jt = (it-1) * n
    j0 = (i0-1) * n
    gt0 = A(jt+1:jt+n, j0+1:j0+n)
    g0t = A(j0+1:j0+n, jt+1:jt+n)
    if (i0 == it) then
       do i = 1, n
          g0t(i,i) = g0t(i,i) - ONE 
       enddo
    endif
    g00 = A(j0+1:j0+n, j0+1:j0+n)
    gtt = A(jt+1:jt+n, jt+1:jt+n)

   end subroutine DQMC_Gtau_DumpA

  !--------------------------------------------------------------------!

  subroutine DQMC_MakeGtau(tau, it, i0)
    !
    ! Purpose
    ! =======
    !    This subroutine generates Gtau.
    !
    ! Arguments
    ! =========
    !
    type(Gtau), intent(inout)    :: tau
    integer, intent(in)          :: it, i0
    
    ! ... local scalar
    integer  :: n, idx, dt, d0, L

    ! ... executable ...
    !
    !  meaning of indices
    !     ii: the ii-th block row or column
    !     ib: the block offset
    !
    !     id: 

    ! initialization
    n = tau%n
    L = tau%L
    
    ! Find increment
    d0 = i0 - tau%i0_up
    if (abs(d0) > L/2) d0 = d0 - sign(L, d0)
    dt = it - tau%it_up
    if (abs(dt) > L/2) dt = dt - sign(L, dt)

    if (abs(dt) + abs(d0) ==  1) then
       ! reduce safe count
       tau%sfc = tau%sfc - 1
       ! Map the increment in a direction of motion
       if (dt == 1) then
          idx = TPLUS
       elseif (dt == -1) then
          idx = TMINUS
       elseif (d0 ==  1) then
          idx = ZPLUS
       elseif (d0 == -1) then
          idx = ZMINUS
       endif
    else
       ! recompute if cannot use update
       tau%sfc = 0
    end if

    ! compute Gtau
    if (tau%sfc /= 0) then 
       ! Update Gtau 
       call DQMC_change_gtau_time(tau, idx, TAU_UP)
       call DQMC_change_gtau_time(tau, idx, TAU_DN)
    else
       ! Recompute Gtau from scratch
       call DQMC_GetGtau(it, i0, TAU_UP, tau%upt0, tau%up0t, tau%upt0, tau%up0t, tau)
       call DQMC_GetGtau(it, i0, TAU_DN, tau%dnt0, tau%dn0t, tau%dnt0, tau%dn0t, tau)
       !tau%sfc = tau%nWrap
    end if

  end subroutine DQMC_MakeGtau

  !-----------------------------------------------------!

  subroutine DQMC_GetGtau(it, i0, spin, gt0, g0t, g00, gtt, tau)
    !
    ! Derived from DQMC_GetGtau with a few bug fixes from Simone
    !
    ! Purpose
    ! =======
    !    
    !    This subroutine computes the (i,j) submatrix of Gtau if which
    !    equals to'R'ow or 'B'oth. and computes the (j,i) submatrix of 
    !    Gtau if which equals to 'C'olumn or 'B'oth 
    !
    ! Mathematically, Gtau(Gdn_tau) is the inverse of 
    !
    !                  [  I                   B_1] 
    !                  [-B_2   I                 ]
    !              M = [     -B_3  I             ]
    !                  [          ...   ...      ]
    !                  [             -B_{L}    I ]
    !                   
    !
    ! The (i,j) submatrix of Gtau is given as
    !
    !     G_ii =    inv(I+B_iB_{i-1}...B_1B_L...B_{i+1})   
    !     gt0 = -Gtau_ii(B_iB_{i-1}...B_1B_L...B_{j+1})  for j = i+1...L  
    !     gt0 =  Gtau_ii(B_iB_{i-1}...B_{j+1})           for j = 1...i-1 
    !
    ! In general, we can write gt0 as
    !
    !         gt0 = (+/-) inv(I+A_1A_2)A_1
    !              = (+/-) inv(inv(A_1)+A_2)
    ! where
    !
    !          A_1 = B_{i}...B_{j+1} and
    !          A_2 = B_{j}...B_{i+1}
    !
    ! The following procedure compute gt0 in a stable way
    ! 
    ! 1. Perform UDT decomposition on inv(A_1) and A_2
    !    
    !       inv(A_1) = U_1D_1T_1
    !           A_2  = U_2D_2T_2
    !
    !    See the DQMC_UDTD in DQMC_B.F90 for detail of UDT decomposition.
    !
    ! 2. Decompose D_1 = barD_1*hatD_1
    !              D_2 = barD_2*hatD_2
    !    where
    !           barD_1(i,i) = max(1, D_1(i,i)) and
    !           hatD_1(i,i) = min(1, D_1(i,i))
    !
    ! 3. Compute
    !
    !    C = hatD_2*T_2*inv(T_1)*inv(barD_1)+inv(barD_2)*inv(U_2)*U_1*hatD_1
    !    
    ! 4. Assemble G as 
    !
    !    G = inv(T_1)*inv(barD_1)*inv(C)*inv(barD_2)*inv(U_2)
    !      = inv(T_1)*inv(D_2T_2inv(T_1)+inv(U_2)U_1D_1)*inv(U_2)
    !      = inv(U_1D_1T_1+U_2D_2T_2)
    !      = inv(inv(A_1)+A_2)
    !      = inv(I+A_1A_2)A_1
    !      = inv(I+B_{i}...B_1*B_l...B_{i-1})B_i...B_{j+1}
    !
    ! Matrix g0t has very similar structure with gt0.
    !
    !     G_jj =    inv(I+B_jB_{j-1}...B_1B_L...B_{j+1})   
    !     g0t = -Gtau_jj(B_jB_{j-1}...B_1B_L...B_{i+1})  for i = j+1...L  
    !     g0t =  Gtau_jj(B_jB_{j-1}...B_{i+1})           for i = 1...j-1 
    !
    ! For a fixed i and j,
    !
    !         g0t = (+/-) inv(I+A_2A_1)A_2
    !              = (+/-) inv(inv(A_2)+A_1)
    ! 
    ! where A_1 and A_2 are as defined before.
    ! Therefore, 
    !
    !         g0t = inv(inv(U_1D_1T_1)+inv(U_2D_2T_2))
    !              = inv(inv(T_1)inv(D_1)inv(U_1)+inv(T_2)inv(D_2)inv(U_2))
    !              = U_2*inv(inv(D_1)inv(U_1)U_2+T_1*inv(T_2)*inv(D_2))*T_1
    !
    ! The same trick of bar and hat is also applied to inv(D_1) and inv(D_2).
    !        
    !         g0t = U_2*inv(barD_2)*inv(...)*inv(barD_1)*T_1
    !
    ! where (...) = hatD_1*inv(U_1)U_2*inv(barD_2)+
    !               inv(barD_1)T_1*inv(T_2)hatD_2
    !
    ! NOTE: the hatD_1, barD_1, hatD_2, and barD_2 here are different from
    !       the previous ones.
    !
    ! See working notes for more detail.
    !
    ! Arguments
    ! =========
    !
    integer,  intent(in)       :: it, i0        ! block indices
    integer,  intent(in)       :: spin
    real(wp), intent(inout)    :: gt0(:,:)     ! submatrix of Gtau
    real(wp), intent(inout)    :: g0t(:,:)     ! submatrix of Gtau
    real(wp), intent(inout)    :: g00(:,:)     ! submatrix of Gtau
    real(wp), intent(inout)    :: gtt(:,:)     ! submatrix of Gtau
    type(Gtau), intent(inout), target  :: tau

    ! ... local scalars    ...
    integer :: info           ! parameters for lapack's sub
    integer :: i              ! iterator
    integer :: n
    integer :: which

    real(wp), pointer :: U1(:,:)        ! 
    real(wp), pointer :: D1(:)          ! 
    real(wp), pointer :: T1(:,:)        ! 
    real(wp), pointer :: U2(:,:)        ! 
    real(wp), pointer :: D2(:)          ! 
    real(wp), pointer :: T2(:,:)        ! 
    real(wp), pointer :: W1(:,:)        ! working space
    real(wp), pointer :: W2(:,:)        !
    real(wp), pointer :: rw(:)          ! working space
    real(wp), pointer :: V(:,:)         ! HSF
    integer,  pointer :: lw(:)          !
    integer,  pointer :: pvt1(:)        !
    integer,  pointer :: pvt2(:)        !

    real(wp), pointer :: bar1i(:)       ! 
    real(wp), pointer :: bar2i(:)       ! 
    real(wp), pointer :: hat1(:)        ! 
    real(wp), pointer :: hat2(:)        ! 
    
    type(SeqB), pointer :: SB1 
    type(SeqB), pointer :: SB2 

    ! ... Executable ...

    write(*,*) "Need to implement g00 and gtt in GetGtau before using it"
    g00 = 0
    gtt = 0
    stop

    ! STEP 0. Initialization
    n = tau%n
    !bar1i => tau%v1
    !bar2i => tau%v2
    !hat1  => tau%v3
    !hat2  => tau%v4

    !if(spin == TAU_UP) then
    !   SB1 => tau%SB1_up
    !   SB2 => tau%SB2_up
    !   V   => tau%V_up
    !else
    !   SB1 => tau%SB1_dn
    !   SB2 => tau%SB2_dn
    !   V   => tau%V_dn
    !endif
    
    !U1 => SB1%U
    !D1 => SB1%D
    !T1 => SB1%T

    !U2 => SB2%U
    !D2 => SB2%D
    !T2 => SB2%T

    !W1 => SB1%W1
    !W2 => SB1%W2
    !rw => SB1%rw
    !lw => SB1%lw
    !pvt1 => SB1%piv
    !pvt2 => SB2%piv

    info = 0

    if (it<1 .or. it>tau%L .or. i0<1 .or. i0>tau%L) then
      write(*,'(A)')"GetGtau can only work with indices in [1,L]. Stop."
      stop
    endif
    

    ! STEP 1. Cmpute UDT decomposition of 
    !         inv(A_1) = inv(B_{i}...B_{j+1})
    !         and A_2  = B_j...B_{i+1}.
    ! ==========================================
    ! W1, W2, rw, lwork, tau, pvt1 can be reused.

    call DQMC_SeqMultB (i0, it+1, SB1, V)
    if ( it /= i0 ) then
       call DQMC_SeqMultBi(it, i0+1, SB2, V)
       which = tau%which
    else
       D2 = ONE
       call DQMC_Eye(n, U2)
       call DQMC_Eye(n, T2)
       which = TAU_T0
    endif
    
    if (which == TAU_T0 .or. which == TAU_BOTH) then
       !
       ! STEP 2.  D_1 = inv(barD_1)*hatD_1
       !          D_2 = inv(barD_2)*hatD_2
       ! ==================================
       
       do i = 1, n
          bar1i(i) = ONE / max(ONE, abs(D1(i)))
          hat1(i)  = D1(i) * bar1i(i)
          bar2i(i) = ONE / max(ONE, abs(D2(i)))
          hat2(i)  = D2(i) * bar2i(i)
       end do

       !   
       ! STEP 3. Compute C = hatD_2*T_2*inv(T_1)*inv(barD_1)+
       !                     inv(barD_2)*inv(U_2)*U_1*hatD_1
       ! =======================================================   
       
       !! Compute  T_2*inv(T_1)
       ! copy T_1 to W_2, because we may need T_1 later
       call blas_dcopy(n*n,T1(:,1),1,W2(:,1),1)
       
       ! W_1 = T_2'
       call DQMC_trans(n, W1, T2)

       ! W_1 = inv(W_2')*W_1 = inv(T_1')*T_2'
       call lapack_dgetrf(n, n, W2, n, pvt1, info)
       call lapack_dgetrs('T', n, n, W2, n, pvt1, W1, n, info)
       if (info /= 0) then
          call DQMC_Error("Error: dgetrs(1) in dqmc_getgtau.", info)
       end if

       ! T_2 = transpose(W_1) = transpose(inv(T_1')*T_2') = T_2*inv(T_1)
       call DQMC_trans(n, T2, W1)
       
       if (tau%which==TAU_T0) then
          ! if only Row is computed, then T1 is not reference, reuse it 
          call blas_dcopy(n*n,gt0(:,1),1,T1(:,1),1)
       end if
       
       ! U_1 = gt0 = U_2'*U_1
       ! ** gt0 here is used as a temp variable
       call blas_dgemm('T', 'N', n, n, n, ONE, U2, n, U1, n, ZERO, gt0, n)
       call blas_dcopy(n*n,gt0,1,U1,1)
       
       !! *** We need to keep T2 and U1 for later use.
       
       ! compute U_1 = barD_2*U_2'*U_1*hatD_1
       call DQMC_ScaleRow(n, gt0, bar2i)
       call DQMC_ScaleCol(n, gt0, hat1)
       
       ! compute W_1 = hatD_2*T_2*inv(T_1)*barD_1
       call blas_dcopy(n*n,T2(:,1),1,W1(:,1),1)
       call DQMC_ScaleRow(n, W1, hat2)
       call DQMC_ScaleCol(n, W1, bar1i)
       
       ! W_1 = W_1 + gt0 (This is called "C" where STEP 3 is defined)
       call blas_daxpy(n*n, ONE, gt0, 1, W1, 1)
       
       !   
       ! STEP 4. Compute inv(T_1)*inv(barD_1)*inv(C)*inv(barD_2)*inv(U_2)
       ! =================================================================   
       
       ! Let gt0 = inv(barD_2) * inv(U2)
       call DQMC_trans(n, gt0, U2)
       call DQMC_ScaleRow(n, gt0, bar2i)

       ! Straight inversion of "C" using LU (dgesv calls dgetrf).
       ! To be modified using the safer alternative of a further
       ! UDT decomposition followed by inversion?
       ! gt0 = inv(W_1)*inv(barD_2)*inv(U_2)
       call lapack_dgesv(n, n, W1, n, pvt2, gt0, n, info)

       if (info /= 0) then
          call DQMC_Error("Error: dgesv(2) in dqmc_getgtau.", info)
       end if
       
       ! gt0 = inv(barD_1)*gt0 = inv(barD_1)*inv(W_1)*inv(barD_2)*inv(U_2)
       call DQMC_ScaleRow(n, gt0, bar1i)
       
       ! gt0 = inv(T_1)*gt0
       !      = inv(T_1)*inv(barD_1)*inv(C)*inv(barD_2)*inv(U_2)
       call lapack_dgetrs('N', n, n, W2, n, pvt1, gt0, n, info)
       if (info /= 0) then
          call DQMC_Error("Error: dgetrs(1) in dqmc_getgtau.", info)
       end if

    end if

    !
    ! Compute g0t, repeat step 2, 3, 4 for Gji
    ! ==========================================

    if (which==TAU_0T .or. which == TAU_BOTH) then       
       !
       ! STEP 5.  inv(D_1) = barD_1*hatD_1
       !          inv(D_2) = barD_2*hatD_2
       ! ======================================
       !
       do i = 1, n
          if (D1(i) == ZERO) then
             call DQMC_Error("Error: in dqmc_getgtau, D1(i)=0.0, i=", i)
          end if
          D1(i) = ONE / D1(i)
          bar1i(i) = ONE / max(ONE, D1(i))
          hat1(i) = D1(i) * bar1i(i)

          if (D2(i) == ZERO) then
             call DQMC_Error("Error: in dqmc_getgtau, D2(i)=0.0, i=", i)
          end if
          D2(i) = ONE / D2(i)
          bar2i(i) = ONE / max(ONE, D2(i))
          hat2(i) = D2(i) * bar2i(i)
       end do
       
       !   
       ! STEP 6. Compute g0t = hatD_1*inv(U_1)U_2*inv(barD_2)+
       !                        inv(barD_1)T_1*inv(T_2)hatD_2
       ! =======================================================   
       if (tau%which == TAU_BOTH) then
          ! Previously, T_2 = T_2*inv(T_1)
          !             U_1 = inv(U_2)*U_1
          ! Therefore, we only need to invert them.
          
          ! first, compute inv(barD_1)T_1*inv(T_2)hatD_2          
          call lapack_dgetrf(n, n, T2, n, pvt1, info)
          if (info /= 0) then
             call DQMC_Error("Error: dgetrf(1) in dqmc_getgtau.", info)
          end if
          call lapack_dgetri(n, T2, n, pvt1, rw, lw(LA_GETRI), info)
          
          ! W1 = U1' = inv(inv(U_2)*U_1) = inv(U_1)*U_2
          call DQMC_Trans(n, W1, U1)
           
       else
          ! No previous computed results. Compute them from scratch.
          !
          ! (1) Compute T_1*inv(T_2) 
          !     Let W_1 = T_1'
          call DQMC_trans(n, W1, T1)

          !     W_1 = inv(T_2')*W_1 = inv(T_2')*T_1'
          call lapack_dgetrf(n, n, T2, n, pvt1, info)
          call lapack_dgetrs('T', n, n, T2, n, pvt1, W1, n, info)
          if (info /= 0) then
             call DQMC_Error("Error: dgetrs(1) in dqmc_getgtau.", info)
          end if
          !     T_2 = W_1' = (inv(T_2')*T_1')' = T_1*inv(T_2)
          call DQMC_trans(n, T2, W1)
          
          ! (2) Compute W_1 = U_2'*U_1
          call blas_dgemm('T', 'N', n, n, n, ONE, U1, n, U2, n, ZERO, W1, n)
       end if


       ! Compute inv(barD_1)T_1*inv(T_2)hatD_2
       call DQMC_ScaleRow(n, T2, bar1i)
       call DQMC_ScaleCol(n, T2, hat2)
          
       ! Compute hatD_1*inv(U_1)U_2*inv(barD_2)
       call DQMC_ScaleRow(n, W1, hat1)
       call DQMC_ScaleCol(n, W1, bar2i)
       

       ! W1 = W1 + T2
       call blas_daxpy(n*n, ONE, T2, 1, W1, 1)

       !
       ! STEP 7. Compute U_2*inv(barD_2)*inv(...)*inv(barD_1)*T_1
       ! =========================================================
       !
       call DQMC_ScaleRow(n, T1, bar1i)
       call lapack_dgesv(n, n, W1, n, pvt1, T1, n, info)
       if (info /= 0) then
          call DQMC_Error("Error: dgesv(3) in dqmc_getgtau.", info)
       end if
       
       ! inv(barD_2)*inv(...)*inv(barD_1)*T_1
       call DQMC_ScaleRow(n, T1, bar2i)

       ! copy the previous result
       call blas_dcopy(n*n,gt0(:,1),1,W2(:,1),1)
       ! multiply -U2, the sign is negative
       call blas_dgemm('N', 'N', n, n, n, -ONE, U2, n, T1, n, ZERO, g0t, n)

    end if

    if( it == i0 .and. tau%which > TAU_T0 )then
      g0t = gt0
    endif

    if ( i0 > it ) then
      if ( tau%which < TAU_0T ) gt0 = -gt0
      if ( tau%which > TAU_T0 ) g0t = -g0t
    endif

    if (spin == TAU_UP) then
       tau%it_up = it
       tau%i0_up = i0
    else
       tau%it_dn = it
       tau%i0_dn = i0
    endif
    
  end subroutine DQMC_GetGtau

  !--------------------------------------------------------------------!

  subroutine DQMC_change_gtau_time(tau, idir, spin)
    !
    ! Purpose
    ! =======
    ! This subroutine computes a new Gtau which is adjacent 
    ! to the one stored in tau using the recursion relations. 
    ! idir specifies which one of the adjacent four G has to be 
    ! computed. 
    !
    ! tau contains (or may contain) two blocks : 
    ! G(i,j) and G(j,i) where i and j are time indices. tau%ii
    ! and tau%ib are assumed to contain the indices i and j
    ! (Note that somewhere else tau%ib contains the displacement
    ! from tau%ii instead). The variable tau%which says whether
    ! G(i,j) and/or G(j,i) are stored. 
    !
    ! This routine applies the following transformation:
    !   if (idir == TPLUS [1]) G(i,j)=>G(i+1,j) and/or G(j,i)=>G(j,i+1)
    !   if (idir == TMINUS[2]) G(i,j)=>G(i-1,j) and/or G(j,i)=>G(j,i-1)
    !   if (idir == ZPLUS [3]) G(i,j)=>G(i,j+1) and/or G(j,i)=>G(j+1,i)
    !   if (idir == ZMINUS[4]) G(i,j)=>G(i,j-1) and/or G(j,i)=>G(j-1,i)
    ! keeping correctly track of the case where i==j (either 
    ! initially or after the transformation). 
    ! i and j are always kept between 1 and L.
    !
    ! Arguments
    ! =========
    !
    type(Gtau),  intent(inout) :: tau
    integer,     intent(in)    :: idir
    integer,     intent(in)    :: spin
  
    ! ... local ...
    integer :: i, j, id, n, L

    ! ... aliases ...
    type(matB), pointer :: B        
    real(wp),   pointer :: gt0(:,:) 
    real(wp),   pointer :: g0t(:,:) 
    real(wp),   pointer :: g00(:,:) 
    real(wp),   pointer :: gtt(:,:) 
    real(wp),   pointer :: W(:,:)   
    real(wp),   pointer :: V(:,:)   
    integer,    pointer :: it       
    integer,    pointer :: i0       

    n = tau%n
    L = tau%L
    W => tau%W2

    call dqmc_gtau_setAlias(spin, tau, B=B, V=V, gt0=gt0, g0t=g0t, g00=g00, gtt=gtt, it=it, i0=i0)

    if(tau%which <= TAU_BOTH) then

       select case (idir)

       case (TPLUS) ! G(i,j)=> G(i+1,j) 
          i = it + 1
          if(i > L) i = 1
          !Multiply by B_{ii+1} 
          call DQMC_MultB_Left  (n, gt0, B, V(:,i), W)
          !Time wrapped through beta. Need to change sign.
          if(i == 1)then
            gt0 = -gt0
          endif
          !Final G is equal time. Handle G(i,j) properly.
          if (i0 == i) then
             do id = 1, n
                gt0(id,id) = 1.d0 + gt0(id,id)
             enddo
          endif

       case (TMINUS) ! G(i,j)=> G(i-1,j) 
          i = it
          !Initial G is equal time. Handle G(i,j) properly.
          if (i0 == it) then
             do id = 1, n
                gt0(id,id) = -1.d0 + gt0(id,id)
             enddo
          endif
          call DQMC_MultBi_Left(n, gt0, B, V(:,i), W)
          !Time wrapped through zero. Need to change sign.
          if(i == 1)then
            gt0 = -gt0
            i = L
          else
            i = i -1
          endif

       case (ZPLUS) !G(i,j)=> G(i,j+1) 
          j = i0 + 1
          if (j > L) j = 1
          ! Initial G is equal time.
          if (i0 == it) then
             do id = 1, n
                gt0(id,id) = -1.d0 + gt0(id,id)
             enddo
          endif
          call DQMC_MultBi_Right(n, gt0, B, V(:,j), W)
          !Time wrapped through beta. Need to change sign.
          if(j == 1)then
            gt0 = -gt0
          endif

       case(ZMINUS) !G(i,j)=> G(i,j-1) 
          j = i0
          call DQMC_MultB_Right(n, gt0, B, V(:,j), W)
          !Time wrapped through zero. Need to change sign.
          if(j == 1)then
            gt0 = -gt0
            j = L
          else
            j = j - 1
          endif
          !Final G is equal time. Treat G(i,j) properly.
          if(it == j)then
             do id = 1, n
                gt0(id,id) = 1.d0 + gt0(id,id)
             enddo
          endif

       end select

    endif

    if(tau%which >= TAU_BOTH) then

       select case (idir)

       case (TPLUS) ! G(j,i)=>G(j,i+1)
          i = it+1
          if(i > L) i = 1
          !Multiply by B_{i+1} and its inverse
          call DQMC_MultBi_Right(n, g0t, B, V(:,i), W)
          !Time wrapped through beta. Need to change sign.
          if(i == 1)then
            g0t = -g0t
          endif
          ! Update equal time G at t 
          call DQMC_MultB_Left  (n, gtt, B, V(:,i), W)
          call DQMC_MultBi_Right(n, gtt, B, V(:,i), W)
          !Final G is equal time. Handle G(j,i) properly.
          if (i0 == i) then
             do id = 1, n
                g0t(id,id) = -1.d0 + g0t(id, id)
             enddo
          endif

       case(TMINUS) ! G(j,i)=>G(j,i-1)
          i = it
          !Initial G is equal time. Handle G(j,i) properly.
          if(i0 == i)then
             do id = 1, n
                g0t(id,id) = 1.d0 + g0t(id, id)
             enddo
          endif
          call DQMC_MultB_Right(n, g0t, B, V(:, i), W)
          ! Update equal time G at t 
          call DQMC_MultBi_Left  (n, gtt, B, V(:, i), W)
          call DQMC_MultB_Right(n, gtt, B, V(:, i), W)
          !Time wrapped through zero. Need to change sign.
          if(i == 1)then
            g0t = -g0t
            i = L
          else
            i = i -1
          endif
          ! Update gtt

       case(ZPLUS) !G(j,i)=>G(j+1,i)
          j = i0 + 1
          if(j > L) j = 1
          !Initial G is equal time. Handle G(j,i) properly.
          if(it == i0)then
             do id = 1, n
                g0t(id,id) = 1.d0 + g0t(id,id)
             enddo
          endif
          call DQMC_MultB_Left  (n, g0t, B, V(:,j), W)
          ! Update equal time G at 0
          call DQMC_MultB_Left  (n, g00, B, V(:,j), W)
          call DQMC_MultBi_Right(n, g00, B, V(:,j), W)
          !Time wrapped through beta. Need to change sign.
          if(j == 1)then
            g0t = -g0t
          endif

       case(ZMINUS) ! G(j,i)=>G(j-1,i)
          j = i0
          call DQMC_MultBi_Left(n, g0t, B, V(:,j), W)
          ! Update equal time G at 0
          call DQMC_MultBi_Left(n, g00, B, V(:,j), W)
          call DQMC_MultB_Right(n, g00, B, V(:,j), W)
          !Time wrapped through zero. Need to change sign.
          if(j == 1)then
             g0t = -g0t
             j = L
          else
             j = j - 1
          endif
          if (it == j) then
             do id = 1, n
                g0t(id,id) = -1.d0 + g0t(id,id)
             enddo
          endif

       end select

    endif

    !Update block index
    select case (idir)
      case(1, 2)
         it = i
      case(3, 4)
         i0 = j
    end select 

  end subroutine DQMC_change_gtau_time

  !---------------------------------------------------------------------!

  subroutine DQMC_Gtau_SetAlias(spin, tau, A, B, V, gt0, g0t, g00, gtt, it, i0, itau, sgn)
     
     integer,    intent(in) :: spin
     type(gtau), target, intent(in) :: tau

     type(MatB), pointer, optional :: B        
     integer,    pointer, optional :: itau(:)  
     real(wp),   pointer, optional :: A(:,:)   
     real(wp),   pointer, optional :: V(:,:)   
     real(wp),   pointer, optional :: gt0(:,:) 
     real(wp),   pointer, optional :: g0t(:,:) 
     real(wp),   pointer, optional :: g00(:,:) 
     real(wp),   pointer, optional :: gtt(:,:) 
     integer,    pointer, optional :: it       
     integer,    pointer, optional :: i0       
     real(wp),   pointer, optional :: sgn      

     select case (spin)
       case (TAU_UP)
          if(present(A))    A    => tau%A_up
          if(present(B))    B    => tau%B_up
          if(present(V))    V    => tau%V_up
          if(present(gt0))  gt0  => tau%upt0
          if(present(g0t))  g0t  => tau%up0t
          if(present(gtt))  gtt  => tau%uptt
          if(present(g00))  g00  => tau%up00
          if(present(it))   it   => tau%it_up
          if(present(i0))   i0   => tau%i0_up
          if(present(sgn))  sgn  => tau%sgnup
          if(present(itau)) itau => tau%itau_up
       case (TAU_DN)
          if(present(A))    A    => tau%A_dn
          if(present(B))    B    => tau%B_dn
          if(present(V))    V    => tau%V_dn
          if(present(gt0))  gt0  => tau%dnt0
          if(present(g0t))  g0t  => tau%dn0t
          if(present(gtt))  gtt  => tau%dntt
          if(present(g00))  g00  => tau%dn00
          if(present(it))   it   => tau%it_dn
          if(present(i0))   i0   => tau%i0_dn
          if(present(sgn))  sgn  => tau%sgndn
          if(present(itau)) itau => tau%itau_dn
      end select
  end subroutine DQMC_Gtau_SetAlias

  !--------------------------------------------------------------------!

  subroutine DQMC_Gtau_GetG0(n, tau, spin, slice, g0tau)
    
     integer, intent(in)       :: n
     type(gtau), intent(inout) :: tau
     integer, intent(in)       :: spin
     integer, intent(in)       :: slice
     real(wp), intent(out)     :: g0tau(n, n)

     integer  :: L, lw, info
     real(wp) :: alpha, beta
     real(wp) :: g(n)
     real(wp), allocatable :: work(:)

    !
    ! This solves the gtau explicitly by Lapack.
    !
    ! ... Executable ...

    alpha = 1.0_wp
    beta  = 1.0_wp
    g0tau = 0.0_wp
    L     = tau%L

    if (.not. tau%g0_stored) then

       call dsyev('V', 'L', n, tau%W2, n, tau%W1, tau%W1, -1, info)
       lw = nint(tau%W1(1))
       allocate(work(lw))

       allocate(tau%e0up(n))
       allocate(tau%e0dn(n))
       allocate(tau%U0up(n, n))
       allocate(tau%U0dn(n, n))

       tau%W1 = 1.0_wp

       call DQMC_GetB(n, tau%U0up, tau%B_up, tau%W1, tau%W2)
       call DQMC_GetB(n, tau%U0dn, tau%B_dn, tau%W1, tau%W2)
       call dsyev('V', 'L', n, tau%U0up, n, tau%e0up, work, lw, info)
       call dsyev('V', 'L', n, tau%U0dn, n, tau%e0dn, work, lw, info)

       tau%g0_stored = .true.

    endif

    if (spin == TAU_UP .or. spin == 0) then
       g = tau%e0up**slice / (1.0_wp + tau%e0up**L)
       call dcopy(n*n, tau%U0up(:,1), 1, tau%W2(:,1), 1)
       call dqmc_scaleCol(n, tau%W2, g)
       call dgemm('N','C', n, n, n, alpha, tau%W2, n, tau%U0up, n, beta, g0tau, n)
    endif

    if (spin == TAU_DN .or. spin == 0) then
       g = tau%e0dn**slice / (1.0_wp + tau%e0dn**L)
       call dcopy(n * n, tau%U0dn(:, 1), 1, tau%W2(:, 1), 1)
       call dqmc_scaleCol(n, tau%W2, g)
       call dgemm('N','C', n, n, n, alpha, tau%W2, n, tau%U0dn, n, beta, g0tau, n)
    endif

    if (spin == 0) g0tau = g0tau / 2

  end subroutine DQMC_Gtau_GetG0

  !--------------------------------------------------------------------!

end module DQMC_GTAU
