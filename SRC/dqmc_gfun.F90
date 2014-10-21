module DQMC_GFun
#include "dqmc_include.h"

  use DQMC_UTIL
  use DQMC_WSPACE
#ifdef DQMC_CKB
  use DQMC_CheckerBoard
#else
  use DQMC_MATB
#endif
  use DQMC_SEQB

  implicit none 

  ! 
  ! This module defines the data type and subroutines for 
  ! constructiong and manipulating the Green's function.
  ! The Green's function in Hubbard's model is defined as
  ! 
  !     G = inv((I+B_{i-1}...B_1*B_l...B_{i})),   
  !
  ! where B_{i} = V_{i}*B.  Matrix V_{i} is diagonal and B is
  ! defined as exp(-t*K), where K is an adjacent matrix for
  ! underline lattice. 
  !
  !     G = inv((I+inv(B_{i-1}...B_1*B_l...B_{i}))),   
  !
  ! For more details, see the working note.
  !
  ! [1] Z. Bai, W.Chen, R. Scalettar, I. Yamazaki, "Lecture Notes 
  !     on Advances of Numerical Methods for Hubbard Quantum Monte
  !     Carlo Simulation." 
  !
  !  Subroutine List
  !  ===============
  !    DQMC_GFun_Init(n, l, GF) : initialize the data type.
  !    DQMC_GetV(G_up, G_dn, explook, hub) : construct matrix V
  !    DQMC_GetG(sl, G, R, Q, T, rw, D, 
  !              tau, piv1, piv2, lwork, nOrth) : construct matrix G
  !
  !  Data Type
  !  =========
  !
  type G_fun
     integer  :: n                 ! Number of sites in Hubbard's model
     integer  :: L                 ! Number of slice
     real(wp) :: sgn               ! Sign of det(G)
     real(wp), pointer :: G(:,:)   ! Matrix of Green's function
     logical  :: owns_G            ! does this G_fun own its G's storage?
     real(wp), pointer :: V(:,:)   ! Matrix info of V(1)...V(l)
     integer  :: ilb               ! index of left most B
     real(wp) :: det

     ! working space, in addition to SB's
     integer, pointer  :: pvt(:) 
     real(wp),pointer  :: tmp(:,:) 
     
     ! data pointer for the C++ module
#if defined(DQMC_ASQRD)
     integer*8   :: cpp_data
#endif

     ! For numerical stab
     integer  :: nWrap, wps, lastwr, maxwrap, fixwrap
     real(wp) :: difflim, errrate
     integer  :: redo, noredo

     ! Block size of delayed update
     integer  :: nBlk   
     integer  :: blkSz
     real(wp),pointer  :: U(:,:) 
     real(wp),pointer  :: W(:,:) 
     integer  :: nModify

     real(wp), pointer :: GS(:,:) 
     real(wp), pointer :: WS(:,:) 

     logical :: sxx

  end type G_fun

  logical, parameter :: GMAT_UP = .true.
  logical, parameter :: GMAT_DN = .false.

#ifdef DQMC_PROFILE

#define PROF_BEGIN if (gfun_profile) call get_time(gfun_t1)
#define PROF_END(c, t) if (gfun_profile) call gfun_end(c, t)

    real(wp) :: GetG_time = 0, ComputeG_time = 0, Getjj_time = 0 
    real(wp) :: UpdateG_time = 0, ApplyUpdate_time = 0
    double precision :: gfun_t1, gfun_t2
    logical  :: gfun_profile = .false.
    integer  :: GetG_count = 0, ComputeG_count = 0, Getjj_count = 0 
    integer  :: UpdateG_count = 0, ApplyUpdate_count = 0 

contains

    subroutine gfun_end(c, t)
      integer  :: c
      real(wp) :: t
      call get_time(gfun_t2)
      c = c + 1
      t = t + gfun_t2 - gfun_t1
    end subroutine gfun_end

    subroutine gfun_print()
      write(*,*) "GetG          ", GetG_count, GetG_time
      write(*,*) "ComputeG      ", ComputeG_count, ComputeG_time
      write(*,*) "Getjj         ", Getjj_count, Getjj_time
      write(*,*) "UpdateG       ", UpdateG_count, UpdateG_time
      write(*,*) "ApplyUpdate   ", ApplyUpdate_count, ApplyUpdate_time
    end subroutine gfun_print

#else

#define PROF_BEGIN
#define PROF_END(c, t)

contains

#endif

  !-------------------------------------------------------------------------!

  subroutine DQMC_GFun_Init(n, L, G, V, WS, nWrap, difflim, errrate, up, ssxx, fixw)
    !
    ! Purpose
    ! =======
    !    This subroutine initiliazes the data type of Green function.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)         :: n              ! Number of sites
    integer, intent(in)         :: L              ! Number of slice
    type(G_fun), intent(inout)  :: G              ! Green's function
    real(wp), intent(in), target:: V(n,L)         ! V matrix
    type(Wspace), intent(in), target :: WS        ! working space
    real(wp), intent(in)        :: difflim        ! limit of mat diff
    real(wp), intent(in)        :: errrate        ! toleratble error
    integer, intent(in)         :: nwrap          ! safe wrap number
    logical, intent(in)         :: up             ! for G_up or G_dn
    integer, intent(in)         :: ssxx
    integer, intent(in)         :: fixw

    ! ... LAPACK subroutine
    integer  :: ILAENV

    ! ... Local variables ...
    integer  :: nb

    ! ... Executable ...

    G%n   = n
    G%L   = L
    G%sgn = ONE
    G%pvt => WS%I2
    G%V   => V

    G%ilb = -L-1
    G%difflim = difflim
    G%nWrap   = nWrap
    G%fixwrap = fixw
    G%maxWrap   = 3*nWrap
    G%wps     = 1
    G%errrate = errrate

    G%redo   = 0
    G%noRedo = 1
    G%lastwr = 0
    G%blkSz  = 0
    G%nModify = 0
    G%det    = 0
    
    nb = ILAENV( 1, "DGETRI", " ", n, -1, -1, -1 )
    G%nblk = min(nb, n)

    allocate(G%G(n,n), G%tmp(n,n))
    G%owns_G = .true.
    G%sgn = ONE

    ! U and W are used when updating G during the sweep over sites
    ! while Gs and WS are used before measuring. Hence, there should
    ! be no conflict in having them pointing to the same location.
    if (up) then
       G%U   => WS%R1
       G%W   => WS%R2
       G%GS  => WS%R1
       G%WS  => WS%R2
    else
       G%U   => WS%R3
       G%W   => WS%R4
       G%GS  => WS%R3
       G%WS  => WS%R4
    end if

    G%sxx = ssxx.gt.0
    
    ! C++ module initialization
#if defined(DQMC_ASQRD)
    call cpp_gfun_init(G%cpp_data, n, L, nWrap, fixw)
#endif
    
  end subroutine DQMC_GFun_Init

  !-------------------------------------------------------------------------!

  subroutine DQMC_GFun_Init_original(n, L, G, V, WS, nWrap, difflim, errrate, up, ssxx, fixw)
    !
    ! Purpose
    ! =======
    !    This subroutine initiliazes the data type of Green function.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)         :: n              ! Number of sites
    integer, intent(in)         :: L              ! Number of slice
    type(G_fun), intent(inout)  :: G              ! Green's function
    real(wp), intent(in), target:: V(n,L)         ! V matrix
    type(Wspace), intent(in), target :: WS        ! working space
    real(wp), intent(in)        :: difflim        ! limit of mat diff
    real(wp), intent(in)        :: errrate        ! toleratble error
    integer, intent(in)         :: nwrap          ! safe wrap number
    logical, intent(in)         :: up             ! for G_up or G_dn
    integer, intent(in)         :: ssxx
    integer, intent(in)         :: fixw

    ! ... LAPACK subroutine
    integer  :: ILAENV

    ! ... Local variables ...
    integer  :: nb

    ! ... Executable ...

    G%n   = n
    G%L   = L
    G%sgn = ONE
    G%pvt => WS%I2
    G%V   => V

    G%ilb = -L-1
    G%difflim = difflim
    G%nWrap   = nWrap
    G%fixwrap = fixw
    G%wps     = nWrap
    G%errrate = errrate

    G%redo   = 0
    G%noRedo = 1
    G%lastwr = 0
    G%blkSz  = 0
    G%nModify = 0
    G%det    = 0
    
    nb = ILAENV( 1, "DGETRI", " ", n, -1, -1, -1 )
    G%nblk = min(nb, n)

    allocate(G%G(n,n), G%tmp(n,n))
    G%sgn = ONE

    ! U and W are used when updating G during the sweep over sites
    ! while Gs and WS are used before measuring. Hence, there should
    ! be no conflict in having them pointing to the same location.
    if (up) then
       G%U   => WS%R1
       G%W   => WS%R2
       G%GS  => WS%R1
       G%WS  => WS%R2
    else
       G%U   => WS%R3
       G%W   => WS%R4
       G%GS  => WS%R3
       G%WS  => WS%R4
    end if

    G%sxx = ssxx.gt.0
    
  end subroutine DQMC_GFun_Init_original

  !-------------------------------------------------------------------------!

  subroutine DQMC_GFun_Clone(G1, G2)
    !
    ! Purpose
    ! =======
    !    G1 = G2
    !
    ! Arguments
    ! =========
    !
    type(G_fun), intent(inout)  :: G1              ! Green's function
    type(G_fun), intent(in)     :: G2              ! Green's function 

    ! ... Executable ...

    G1%n   = G2%n
    G1%L   = G2%L
    G1%det = G2%det

    G1%sgn = G2%sgn
    G1%G   => G2%G
    G1%owns_G = .false.
    G1%V   => G2%V
    G1%GS  => G2%GS

    G1%sxx = G2%sxx

#if defined(DQMC_ASQRD)
    G1%cpp_data = G2%cpp_data
#endif

  end subroutine DQMC_GFun_Clone

  !-------------------------------------------------------------------------!

  subroutine DQMC_GFun_Duplicate(G1, G2)
    !
    ! Purpose
    ! =======
    !    G1 = G2
    !
    ! Arguments
    ! =========
    !
    type(G_fun), intent(inout)  :: G1              ! Green's function
    type(G_fun), intent(in)     :: G2              ! Green's function 

    ! ... Executable ...

    G1%n    = G2%n
    G1%L    = G2%L
    G1%sgn  = G2%sgn
    G1%pvt => G2%pvt
    G1%V   => G2%V

    G1%ilb     = G2%ilb
    G1%difflim = G2%difflim
    G1%nWrap   = G2%nWrap
    G1%wps     = G2%wps
    G1%errrate = G2%errrate

    G1%redo    = G2%redo   
    G1%noRedo  = G2%noRedo 
    G1%lastwr  = G2%lastwr 
    G1%blkSz   = G2%blkSz  
    G1%nModify = G2%nModify
    G1%det     = G2%det    
    
    G1%nblk = G2%nblk

    allocate(G1%G(G2%n,G2%n), G1%tmp(G2%n,G2%n))
    G1%owns_G = .true.

    G1%G    = G2%G
    G1%tmp  = G2%tmp
    G1%U   => G2%U 
    G1%W   => G2%W  
    G1%GS  => G2%GS 
    G1%WS  => G2%WS 

    G1%sxx = G2%sxx

#if defined(DQMC_ASQRD)
    call cpp_gfun_init(G1%cpp_data, G2%n, G2%L, G2%nWrap, G2%fixwrap)
#endif

  end subroutine DQMC_GFun_Duplicate

  !-------------------------------------------------------------------------!

  subroutine DQMC_GFun_CopyUp(G1, G2, P)
    !
    ! Purpose
    ! =======
    !    G1 = G2
    !
    ! Arguments
    ! =========
    !
    type(G_fun), intent(inout)  :: G1              ! Green's function
    type(G_fun), intent(in)     :: G2              ! Green's function 
    real(wp),    intent(in)     :: P(G2%n)

    integer :: i, j

    ! ... Executable ...

    G1%ilb = G2%ilb
    G1%sgn = G2%sgn

    do i = 1, G2%n
       do j = 1, G2%n
          G1%G(i,j) = -P(i)*P(j)*G2%G(j,i)
       end do
       G1%G(i,i) = G1%G(i,i) + ONE 
    end do

    G1%det = G2%det
    do i = 1, G2%L
       do j = 1, G2%n
          G1%det = G1%det + log(G2%V(j,i))
       end do
    end do

  end subroutine DQMC_GFun_CopyUp

  !-------------------------------------------------------------------------!

  subroutine DQMC_Gfun_Free(G)
    !
    ! Purpose
    ! =======
    !    This subroutine frees memory of G.
    !
    ! Arguments
    ! =========
    !
    type(G_fun), intent(inout)  :: G  ! Green's function

    ! ... Executable ...

    if (G%owns_G) then
       deallocate(G%G, G%tmp)
#if defined(DQMC_ASQRD)     
       call cpp_gfun_free(G%cpp_data)
#endif
    end if

  end subroutine DQMC_Gfun_Free
  
  !--------------------------------------------------------------------------!

  subroutine DQMC_GetG(il, G, SB)
    !
    ! Purpose
    ! =======
    !    This subroutine returns 
    !    
    !         G = inv((I+B_{il}...B_1*B_L...B_{il+1})),
    !
    !    where B_{i} = V_{i}*B. Matrix V_{i} is diagonal,
    !    whose elements are stored in a vector V.
    !    It also returns the sign of det(G).
    !
    ! Arguments
    ! =========
    !
    integer,  intent(in)       :: il            ! starting slice
    type(G_fun), intent(inout) :: G
    type(SeqB), intent(inout)  :: SB

    ! ... local scalars    ...
    logical  :: match
    real(wp) :: diff, sgn

    ! ... Executable ...

    match = ((G%ilb .eq. G%L .and. il .eq. 1) .or. (G%ilb .eq. il-1)) 

    if (match) then
       ! Swap the slice of G_up
       PROF_BEGIN
#if defined(DQMC_ASQRD)
       call cpp_gfun_swapg(G%cpp_data, SB%B%B, SB%B%Bi, G%V(:,il), G%G)
#else
       call DQMC_MultB_Left  (G%n, G%G, SB%B, G%V(:,il), SB%W1)
       call DQMC_MultBi_Right(G%n, G%G, SB%B, G%V(:,il), SB%W1)
#endif
       PROF_END(GetG_count,GetG_time)

       ! check if the swapping is safe
       G%wps = G%wps - 1

       ! need to recompute
       if (G%wps .eq. 0) then
          ! store current result
          ! G%tmp = G%G
          call blas_dcopy(G%n*G%n, G%G, 1, G%tmp, 1)
          ! recompute G
          PROF_BEGIN
#if defined(DQMC_ASQRD)
          call cpp_gfun_computeg(G%cpp_data, il, sgn, G%G, G%V, SB%B%B, SB%nOrth, G%det)  
#else
          call DQMC_ComputeG(il, G%n, sgn, G%G, G%V, SB, G%pvt, .true., G%det, G%sxx)
#endif
          PROF_END(ComputeG_count,ComputeG_time)

          ! evaluate the difference
          diff = DQMC_MatDiff(G%n, G%tmp, G%G)
          
          ! increase the statistics
          if(diff .gt. G%difflim .or. abs(sgn-G%sgn) .gt. 1.0d-3) then
             G%redo = G%redo + 1
             write(*,'(3x,A, f12.6)')'Error in updating the Green''s function :', diff
          else
             G%noredo = G%noredo + 1
          end if          

          G%wps = G%nWrap
          G%sgn = sgn
       end if
    else
       ! cannot use swap, recompute
       PROF_BEGIN
#if defined(DQMC_ASQRD)
       call cpp_gfun_computeg(G%cpp_data, il, G%sgn, G%G, G%V, SB%B%B, SB%nOrth, G%det)  
#else
       call DQMC_ComputeG(il, G%n, G%sgn, G%G, G%V, SB, G%pvt, .true., G%det, G%sxx)
#endif
       PROF_END(ComputeG_count,ComputeG_time)
       G%wps = G%nWrap
    end if

    G%ilb = il
       
  end subroutine DQMC_GetG

 !--------------------------------------------------------------------------!

  subroutine DQMC_GetG_2nd_order(G, B)
     !
     ! Purpose
     ! =======
     !    This subroutine computes the GF correct to dtau^2.
     ! 
     ! Argument
     ! ========

     type(G_fun), intent(inout) :: G
     type(matB),  intent(in) :: B

     integer :: n

     n = B%n

     G%GS = G%G

     ! 
     ! Temporary fix for the inconsistency b/t equal-time
     ! and time-dependent tau=0 measurements. The inconsistency
     ! is caused by using 2nd and 1st order Trotter breakup
     ! when generating B matrices.
     !
     !call DQMC_MultrtB0i_Right(n, G%GS, B, G%WS)
     !call DQMC_MultrtB0_Left(n, G%GS, B, G%WS)

  end subroutine

  !--------------------------------------------------------------------------!

  subroutine DQMC_ComputeG(il, n, sgn, G, V, SB, pvt2, compDet, det, sxx)

    integer,  intent(in)       :: il            ! starting slice
    integer, intent(in)        :: n             ! dim of G
    real(wp), intent(inout)    :: sgn           ! sign(det(G)) 
    real(wp), intent(in)       :: V(:,:)        ! HSF
    real(wp), intent(inout)    :: G(:,:)        ! Green's function
    type(SeqB), intent(inout)  :: SB            
    integer, intent(inout)     :: pvt2(:)       ! working space 
    logical, intent(in)        :: compDet
    real(wp), intent(inout)    :: det           ! det(G)
    logical, intent(in)        :: sxx

    ! ... local scalars    ...
    integer :: info           ! parameters for lapack's sub
    integer :: i, j           ! iterator

#   ifdef _SXX
       ! Error estimates
       real(wp) :: rcond
       real(wp) :: berr(n)
       ! Additional working space
       real(wp) :: work(4*n)
       integer  :: iwork(n)
       !More parameters for dgesvxx
       integer, parameter :: nerr = 2, maxerr = 3
       real(wp) :: errn(n, maxerr), errc(n, maxerr)
       integer, parameter ::  nparam = 0, maxparam = 3
       real(wp) :: param(maxparam)
#   endif

    !Alias to SB's working space
    real(wp), pointer :: U(:,:)  
    real(wp), pointer :: T(:,:)  
    real(wp), pointer :: D(:)    
    real(wp), pointer :: D2(:)   
    real(wp), pointer :: R(:)    
    real(wp), pointer :: C(:)    
    real(wp), pointer :: W1(:,:) 
    real(wp), pointer :: W2(:,:) 
    real(wp), pointer :: W3(:,:) 
    integer,  pointer :: pvt1(:) 

    !real(wp) :: W3(n,n)

    ! ... Executable ...
#ifdef DQMC_PROFILE
    call profile_enable()
#endif
    U    => SB%U
    D    => SB%D
    T    => SB%T
    W1   => SB%W1
    W2   => SB%W2
    pvt1 => SB%piv
    D2   => SB%tau
    W3   => SB%W3
    C    => SB%tau
    R    => SB%rw


    !
    ! STEP 0. Initialization
    ! ======================
    !! Setup index for B_{i}
    sgn  = ONE
    det  = ZERO
    info = 0
    !
    ! STEP 1. Cmpute A = B_{il}...B_1*B_l...B_{il+1}.
    !         and its UDT decomposition, A = UDT
    ! =============================================
    call DQMC_SeqMultB(il, il+1, SB, V)

    ! STEP 2. Compute C = inv(U)*inv(T) + D.
    ! ======================================

    !! Compute the LU decomposition of T first.
    call lapack_dgetrf(n, n, T, n, pvt1, info)
    if (info .ne. 0) then
       T(info, info) = epsilon(ONE)
       call DQMC_Warning("T is singular, dgetrf in dqmc_ComputeG.", info)
    end if

    ! determine the sign of det(T)
    call DQMC_DetSgn(n, T, pvt1, sgn)
    if (compDet) then
       call DQMC_DetLog(n, T, det)
    end if

    !! Solve T'W = U. W = transpose(inv(U)inv(T))
    ! W1 = U
    call blas_dcopy(n*n, U, 1, W1, 1)
    call lapack_dgetrs('T', n, n, T, n, pvt1, W1, n, info)
    if (info .ne. 0) then
       call DQMC_Error("Error: dgetrs in dqmc_getgp.", info)
    end if

    ! STEP 2.5 Here we need to decompose D = D2/D1=D1\D2.
    !       if   abs(D(i,i))>1,
    !            D1(i,i) = 1/abs(D(i,i))
    !            D2(i,i) = sign(D(i,i))
    !       else
    !            D1(i,i) = 1
    !            D2(i,i) = D(i,i) 
    ! ======================================

    !! Compute W*D1+D2  
    !!   = transpose(inv(U)inv(T))*D1+D2
    !!   = transpose(D1*inv(U)inv(T) + D2)
  
    do i = 1,n
       if (abs(D(i)) .gt. ONE) then
          if (D(i) .gt. ZERO) then
             D2(i) = ONE
          else
             D2(i) = -ONE
          end if
          D(i)  = ONE/abs(D(i))
       else
          D2(i) = D(i)
          D(i)  = ONE
       end if

      ! W' = W'*D1 
       do j = 1, n
          W1(j,i) = W1(j,i)*D(i)
       end do
       ! W' = W' + D2 = W'*D1+D2
       W1(i,i) = W1(i,i) + D2(i)
    end do
    
    !   
    ! STEP 3. Compute G = inv(T)*D1*inv(C)*inv(U)
    ! ===========================================   

    if (sxx) call blas_dcopy(n*n, W1, 1, W3, 1)
    
    !! Compute the LU decomposition of W1.
    call lapack_dgetrf(n, n, W1, n, pvt2, info)

    if (info .ne. 0) then
       call DQMC_Warning("matrix D+inv(U)inv(T) is singular. &
            & dgetrf(2) in dqmc_computeG.", info)
       W1(info, info) = epsilon(ONE)
    end if

    ! determine the sign of det(W)
    call DQMC_DetSgn(n, W1, pvt2, sgn)
    if (compDet) then
       call DQMC_DetLog(n, W1, det)
    end if

    !! G = inv(U), which is transpose of U.
    call DQMC_Trans(n, G, U)

    !! G = D1*G = D1*inv(U)
    call DQMC_ScaleRow(n, G, D)
    call blas_dcopy(n*n, G, 1, U, 1)  

    !! Solve W'G = inv(U), G = inv(W')*G
    call lapack_dgetrs('T', n, n, W1, n, pvt2, G, n, info)
    if (info .ne. 0) then
       call DQMC_Error("Error: dgetrs(2) in dqmc_getgp.", info)
    end if
    !Refine the solution
    if (sxx) then
#      ifdef _SXX
          call lapack_dgerfsx('T', 'N', n, n, W3, n, W1, N, pvt2, R, C,   &
            & U, n, G, n, rcond, berr, nerr, errn, errc, &
            & nparam, param, work, iwork, info)
#      endif
    endif

    !! Solve TG = inv(W')*inv(U), G = inv(T)inv(W')*D1*inv(U)
    call lapack_dgetrs('N', n, n, T, n, pvt1, G, n, info)
    if (info .ne. 0) then
       call DQMC_Error("Error: dgetrs(3) in dqmc_getgp.", info)
    end if

    !
    ! STEP 4. Compute the sign of det(G)
    ! ==================================
    
    !! Compute the LU decomposition of U first.
    call lapack_dgetrf(n, n, U, n, pvt2, info)
    if (info .ne. 0) then
       call DQMC_Error("Error: dgetrf(4) in dqmc_getgp.", info)
    end if
    call DQMC_DetSgn(n, U, pvt2, sgn)

    ! det = det(D1)/(det(W)*det(T))
    if (compDet) then
       det = -det
       do i = 1, n
          det = det + log(abs(D(i)))
       end do
    end if

#ifdef DQMC_PROFILE
    call profile_disable()
#endif

  contains

    !--------------------------------------------------------!
    
    subroutine DQMC_DetSgn(n, A, pvt, sgn)
      !
      ! Purpose
      ! =======
      !    This subroutine computes sign(DET(A))
      !    where A is the content after dgetrf;
      !    pvt is the vipoting vector.
      ! 
      ! Pre-assumption
      ! ==============
      !    A contains the content of dgetrf
      !
      ! Argument
      ! ========
      integer, intent(in)     :: n
      real(wp), intent(in)    :: A(n,n)
      integer, intent(in)     :: pvt(n)
      real(wp), intent(inout) :: sgn
      
      ! ... local scalar 
      integer :: i
      
      !! decide the sgn of det(Q)
      do i = 1, n
         if (pvt(i).ne. i) then
            sgn = -sgn
         end if

         !! ?? stable ??
         if (A(i,i) .lt. ZERO) then
            sgn = -sgn
         end if
      end do
      
    end subroutine DQMC_DetSgn

    !--------------------------------------------------------!
    
    subroutine DQMC_DetLog(n, A, det)
      !
      ! Purpose
      ! =======
      !    This subroutine computes log(abs(DET(A)))
      !    where A is the content after dgetrf;
      !    pvt is the vipoting vector.
      ! 
      ! Pre-assumption
      ! ==============
      !    A contains the content of dgetrf
      !
      ! Argument
      ! ========
      integer, intent(in)     :: n
      real(wp), intent(in)    :: A(n,n)
      real(wp), intent(inout) :: det
      
      ! ... local scalar 
      integer :: i
      
      !! decide the sgn of det(Q)
      do i = 1, n
         det = det + log(abs(A(i,i)))
      end do
      
    end subroutine DQMC_DetLog

  end subroutine DQMC_ComputeG

  !-----------------------------------------------------------------------!

  subroutine DQMC_UpdateWraps(G)
    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !    This subroutine updates G%nWraps depends on previous 
    !    execution results
    !
    ! Arguments
    ! =========
    !
    type(G_fun), intent(inout) :: G
        
    ! ... paremeters ...
    integer, parameter  :: DQMC_REDO_F      = 10
    integer, parameter  :: DQMC_NOREDO_F    = 100
    real(wp), parameter :: DQMC_REDO_RATE   = 0.20_wp
  

    ! ... local scalar ...
    real(wp) :: redorat
    integer  :: redo, noredo
#   ifdef _QMC_MPI
       integer  :: send_cnt(2), recv_cnt(2), err
#   endif

    if (G%fixwrap .le. 0) then

#      ifdef _QMC_MPI
          send_cnt(1) = G%redo
          send_cnt(2) = G%noredo
          call mpi_allreduce(send_cnt, recv_cnt, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD, err)           
          redo   = recv_cnt(1)
          noredo = recv_cnt(2)
#      else  
          redo   = G%redo   
          noredo = G%noredo 
#      endif

       redorat = dble(redo)/dble(redo+noredo)
       !write(*,'(3x,A,i4,A,i6)') 'Total number of errors : ', redo,' out of ', noredo+redo

       ! ... Executable ...
       if(redo .gt. DQMC_REDO_F) then
          if(redorat.gt. G%errrate)then
             !write(*,'(3x,A)')'Lowering nwrap and resetting counter'
             G%nwrap = G%nwrap - 1
             G%maxwrap = G%nwrap
             G%redo  = 0
             G%noredo = 1
             !write(*,'(3x,A,i4)')'New nwrap :', G%nwrap
          endif
       endif
       
       if(NoRedo .gt. DQMC_NOREDO_F) then
          if(redorat .lt. DQMC_REDO_RATE*G%errrate) then
             if (G%nwrap .ge. G%lastwr) then
                if (G%nwrap < G%maxwrap) then
                   !write(*,'(3x,A)')'Raising nwrap and resetting counter'
                   G%nwrap = G%nwrap + 1
                   G%redo = 0
                   G%noredo = 1
                   G%lastwr = G%nwrap
                   !write(*,'(3x,A,i4)')'New nwrap :', G%nwrap
                endif
             end if
          endif
       endif

    endif

  end subroutine DQMC_UpdateWraps

  !-----------------------------------------------------------------------!

  subroutine DQMC_SyncWraps(G1, G2)
    !
    ! Purpose
    ! =======
    !    This subroutine syncronize the value of nwrap
    !    after this was updated for one of the G's
    !
    ! Arguments
    ! =========
    !
    type(G_fun), intent(inout) :: G1, G2

    if (G1%redo==0 .and. G1%noredo==1) then
       G2%redo   = 0
       G2%noredo = 1
       G2%nwrap   = G1%nwrap
       G2%maxwrap = G1%maxwrap
       G2%lastwr  = G1%lastwr
    elseif (G2%redo==0 .and. G2%noredo==1) then
       G1%redo   = 0
       G1%noredo = 1
       G1%nwrap   = G2%nwrap
       G1%maxwrap = G2%maxwrap
       G1%lastwr  = G2%lastwr
    end if

  end subroutine DQMC_SyncWraps


  !-----------------------------------------------------------------------!
  ! The following three subroutines are used for delayed update.          !
  !-----------------------------------------------------------------------!
  
  function DQMC_Gfun_Getjj(n, j, blksz, G, U, W) result(gjj)
    !
    ! Purpose
    ! =======
    !    This function returns the (j,j)th element of G to gjj,
    !    in which 
    !
    !       G = G_1 + UV'
    !
    !    where U,V are n*BlkSz.
    !    
    !    Therefore, 
    !
    !       gjj = G_1(j,j)+U(j,:)*V(j,:)'
    !
    ! Arguments
    ! =========
    !
    real(wp)                :: gjj
    integer , intent(in)    :: n, j, blksz
    real(wp), intent(in)    :: G(n,n)
    real(wp), intent(in)    :: U(n,n)
    real(wp), intent(in)    :: W(n,n)

    ! ... BLAS function ...
    real(wp), external :: blas_ddot

    ! ... Executable ...

    PROF_BEGIN
    gjj = G(j,j)
    if (blkSz .gt. 0) then
       gjj = gjj + blas_ddot(blkSz, U(j,1), n, W(j,1), n)
    end if
    PROF_END(Getjj_count,Getjj_time)

  end function DQMC_Gfun_Getjj
  
 !--------------------------------------------------------------------------!
  
  subroutine DQMC_UpdateG(j, gamma, G)
    !
    ! Purpose
    ! =======
    !    This subroutine updates U, V, and D, which accumulate
    !    rank-1 updates of G.
    !
    !    Matrix G has the form
    !
    !       G_k = G_1 + UV'
    !
    !    where U,V are n*BlkSz.
    !    
    !    The new matrix is 
    !
    !       G_{k+1} = G_{k} + gamma*xy'
    !    
    !    where 
    !    
    !       x = G_1(:,j) + UV'(:,j)
    !       y' = G_1(j,:) - e_j' + U(j,:)V'
    !
    !    The update of U, V, and D is as follows.
    !
    !       U = [U x], V = [V gamma*y]
    !
    !    If blkSz==nBlk, then apply the update
    ! Arguments
    ! =========
    !

    integer , intent(in)       :: j
    real(wp), intent(in)       :: gamma
    type(g_fun), intent(inout), target :: G

    ! ... local variables ...
    real(wp), pointer :: x(:) 
    real(wp), pointer :: y(:) 
    real(wp)  :: xx(G%n)
    integer   :: n, blksz

    ! ... Executable ...

    PROF_BEGIN
    n     = G%n
    blksz = G%blkSz
    x => G%U(1:n, blkSz+1)
    y => G%W(1:n, blkSz+1)

    x = G%G(1:n,j)
    y = G%G(j,1:n)
    y(j) = y(j) - ONE

    if (blkSz .gt. 0) then
       ! if U, V are not empty, add their effects
       xx = G%W(j,1:n)
       call blas_dgemv('N', n, blkSz, ONE, G%U, n, xx, 1, ONE, x, 1)
       xx = G%U(j,1:n)
       call blas_dgemv('N', n, blkSz, ONE, G%W, n, xx, 1, ONE, y, 1)
    end if
    call blas_dscal(n, gamma, y, 1)
    
    G%blkSz = G%blkSz + 1
 
    ! apply the update when necessary
    PROF_END(UpdateG_count,UpdateG_time)
    call DQMC_ApplyUpdate(G, .false.)
    !call DQMC_ApplyUpdate(G, .true.)
        
  end subroutine DQMC_UpdateG

  !--------------------------------------------------------------------------!
  
  subroutine DQMC_ApplyUpdate(G, forced)
    !
    ! Purpose
    ! =======
    !    This subroutine updates G with U, V, and D.
    !    Matrix G has the form
    !
    !       G_k = G_1 + UV'
    !
    !    where U,V are n*BlkSz.
    !
    ! Arguments
    ! =========
    !
    type(g_fun), intent(inout) :: G
    logical, intent(in)        :: forced

    ! ... local variables ...
    integer :: n

    ! ... Executable ...

    if (forced .or. G%blkSz .eq. G%nBlk) then
       PROF_BEGIN
       n = G%n
       ! apply the update when necessary
       call blas_dgemm('N','T', n, n, G%blkSz, ONE, G%U, n, G%W, n, ONE, G%G, n)
       
       ! reset the block size
       G%blkSz = 0
       PROF_END(ApplyUpdate_count,ApplyUpdate_time)
    end if

  end subroutine DQMC_ApplyUpdate

  !--------------------------------------------------------------------------!

end module DQMC_GFun
