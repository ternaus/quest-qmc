module DQMC_PREFIX
#include "dqmc_include.h"

  use DQMC_UTIL
  use _DQMC_MATB
  use DQMC_SEQB

  implicit none 

  ! 
  ! This module contains data structure for statble computation of 
  !     B_{i1}...B{i2}
  ! where i1>i2. 
  ! The major data struture is the UDT decomposition, where
  ! U is an orthogonal matrix, D is diagonal, and T is a 'good-conditioned'
  ! matrix.
  !
  ! The tree of B-segment is arranged as follows.
  ! Let the number of leaf node be some power of two.
  ! The number of parent node is the 'prefix' of the numbers of its 
  ! children nodes. Use an eight leaf tree as an example.
  ! Each leaf node contains 12 matrix product. 
  ! 
  !                          1 (1-96) or (49-96-1-48)
  !                    +-----------------------------+     
  !                    |                             |
  !                10(1-48)                      11(49-96)
  !           +---------------+               +--------------+
  !           |               |               |              |
  !       100(1-24)      101(25-48)      110(49-72)     111(73-96) 
  !       +-------+      +--------+      +--------+     +--------+
  !       |       |      |        |      |        |     |        |
  !     1000    1001    1010    1011    1100    1101   1110    1111    
  !    (1-12) (13-24) (25-36) (37-48) (49-60) (61-72) (73-84) (85-96)
  !
  !  Types
  ! =========
  !
  !  UDT decomposition
  !
  type UDT
     integer  :: n                      ! dim of B
     integer  :: il, ir                 ! indices of B_i
     real(wp),pointer :: U(:,:)    ! 
     real(wp),pointer :: D(:)      ! 
     real(wp),pointer :: T(:,:)    ! 
  end type UDT

  integer, parameter :: POS_LEFT  = 0
  integer, parameter :: POS_RIGHT = 1
  integer, parameter :: POS_ROOT  = 2

  !
  ! Intermediate results
  !
  type PrefixB
     integer   :: n                     ! dim of Bs
     integer   :: L                     ! number of Bs
     integer   :: nOrth                 ! number of Bs per node
     integer   :: nSeg                  ! number of leaf segments
     type(UDT), pointer :: Bseg(:)      ! nodes (UDT)
     integer, pointer   :: queue(:)     ! queue of product
     integer   :: leaflevel             !
     integer   :: maxq                  !
     real(wp), pointer  :: W3(:,:)
     logical   :: inv
  end type PrefixB
  
contains

  !----------------------------------------------------------------------!

  subroutine DQMC_PrefixB_Init(n, L, nOrth, V, SB, PB, inv)
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
    real(wp), intent(in):: V(n,L)       ! Hubbard-Stratonivich Field 
    type(SeqB), intent(inout) :: SB     ! Seq B
    type(PrefixB), intent(inout) :: PB  ! Prefix B
    logical, intent(in) :: inv          ! true if inv(B_i) is requred

    ! ... Local scalar ....
    integer  :: i, nSeg, il, ir
    logical  :: re_alloc

    ! ... Executable ...

    PB%n     = n
    PB%L     = L
    PB%nOrth = nOrth
    PB%inv   = inv

    ! nSeg is the total number of leaf nodes
    nSeg = L/nOrth
    if (nSeg .eq. 1) return      ! quick return for a pecial case 

    ! test if nLeaf is power of 2
    if (iand(nSeg, nSeg-1) .ne. 0) then
       call DQMC_Error("DQMC_PrefixB_Init: nLeaf is not power of two,", nSeg)
    end if
    ! the number of starting leafnode
    PB%leaflevel = nSeg
        
    ! nSeg is the total number of leaf nodes.
    ! PB%nSeg is the total number of nodes, internal and leaf
    PB%nSeg  = 2*nSeg-1

    ! allocate space for BSeg
    re_alloc = .true.
    if (associated(PB%BSeg)) then
       if (size(PB%BSeg) .eq. PB%nSeg) then
          re_alloc = .false.
       else
          deallocate(PB%BSeg)
       end if
    end if
    if (re_alloc) then
       allocate(PB%BSeg(PB%nSeg))
    end if
    
    ! We do not need to allocate space for root node
    do i = 2, PB%nSeg
       call DQMC_Reshape(n, n, PB%BSeg(i)%U)
       call DQMC_Reshape(n,    PB%BSeg(i)%D)
       call DQMC_Reshape(n, n, PB%BSeg(i)%T)
    end do

    ! Construct the leaf node and propagage (merge) them up  
    il = nOrth
    ir = 1
    if (PB%nSeg .gt. 1) then
       do i = PB%leaflevel, PB%nSeg
          if (inv) then
             call DQMC_SeqMultBi(il, ir, SB, V)
          else
             call DQMC_SeqMultB(il, ir, SB, V)
          end if
          
          ! copy the result
          PB%BSeg(i)%U = SB%U
          PB%BSeg(i)%D = SB%D
          PB%BSeg(i)%T = SB%T
          PB%BSeg(i)%il = il
          PB%BSeg(i)%ir = ir

          ir = il + 1          
          il = il + nOrth
       end do
    end if

    ! from the left most non-leaf node
    call DQMC_Reshape(n, n, PB%W3)
    do i = PB%leaflevel - 1, 2, -1
       ! NOTE the index !!!
       il = DQMC_Left_Child(i)
       ir = DQMC_Right_Child(i)
       ! Copy i1 (right child) to i (parent)
       PB%Bseg(i)%U = PB%Bseg(il)%U
       PB%Bseg(i)%D = PB%Bseg(il)%D
       PB%Bseg(i)%T = PB%Bseg(il)%T
       PB%Bseg(i)%il = PB%Bseg(ir)%il
       PB%Bseg(i)%ir = PB%Bseg(il)%ir

       call DQMC_UDT_Merge(n, PB%Bseg(i)%U, PB%Bseg(i)%D, PB%Bseg(i)%T,&
            PB%Bseg(ir)%U, PB%Bseg(ir)%D, PB%Bseg(ir)%T, & 
            SB%W1, SB%W2, PB%W3, SB%rw, SB%lw, SB%tau, SB%piv)
    end do

    ! .................................................
    ! allocate queue, whose size is the number of level
    ! maxq = log_2(nSeg)
    PB%maxq  = 1
    do while (nSeg .ne. 0)
       nSeg = nSeg/2
       PB%maxq  = PB%maxq + 1
    end do

    call DQMC_Reshape(PB%maxq, PB%queue)

  end subroutine DQMC_PrefixB_Init

  !----------------------------------------------------------------------!

  subroutine DQMC_PrefixB_Free(PB)
    !  
    ! Purpose
    ! =======
    !    This subroutine frees dynamically allocated memory of MP.
    !
    ! Arguments
    ! =========
    !
    type(PrefixB), intent(inout) :: PB  ! intermediate results
    
    ! ... local scalar ....
    integer  :: i

    ! ... Executable ...
    do i = 2, PB%nSeg
       call DQMC_Free(PB%BSeg(i)%U)
       call DQMC_Free(PB%BSeg(i)%D)
       call DQMC_Free(PB%BSeg(i)%T)
    end do
    
    deallocate(PB%Bseg)
    deallocate(PB%queue)
    call DQMC_Free(PB%W3)

  end subroutine DQMC_PrefixB_Free

  !----------------------------------------------------------------------!

  subroutine DQMC_UDT_Merge(n, U1, D1, T1, U2, D2, T2, W1, W2, W3, &
       rw, lwork, tau, piv)
    !
    ! U_1D_1T_1 = U_2(D_2T_2*U_1)D_1T_1
    !
    ! steps:
    ! ======
    !    1. U_1 = D_2T_2U_1
    !    2. update U_1D_1T_1 
    !    3. combine U_2U_1
    ! 
    ! Arguments
    ! =========
    integer,  intent(in)    :: n               ! parameter of Hub
    real(wp), intent(inout) :: U1(n,n)         ! U-factor
    real(wp), intent(inout) :: D1(n)           ! D-factor
    real(wp), intent(inout) :: T1(n,n)         ! T-factor
    real(wp), intent(in)    :: U2(n,n)         ! U-factor
    real(wp), intent(in)    :: D2(n)           ! D-factor
    real(wp), intent(in)    :: T2(n,n)         ! T-factor
    real(wp), intent(inout) :: W1(n,n)         ! working space
    real(wp), intent(inout) :: W2(n,n)         ! working space
    real(wp), intent(inout) :: W3(n,n)         ! working space
    real(wp), intent(inout) :: rw(:)           ! working space
    real(wp), intent(inout) :: tau(n)          ! working space
    integer,  intent(in)    :: lwork(:)        !
    integer,  intent(inout) :: piv(n)
    
    ! ... Executable ...
    
    ! W_3 = D_2*T_2*U_1
    call blas_dgemm('N', 'N', n, n, n, ONE, T2, n, U1, n, ZERO, W3, n)
    call DQMC_ScaleRow(n, W3, D2)

    ! update UDT decomposition for W_3
    call DQMC_UDTD(n, W3, D1, T1, W1, W2, rw, tau, piv, lwork)

    ! U_1 = U_2*W_3
    call blas_dgemm('N', 'N', n, n, n, ONE, U2, n, W3, n, ZERO, U1, n)

  end subroutine DQMC_UDT_Merge

  !----------------------------------------------------------------------!

  subroutine DQMC_PrefixB_Update(n, L, idx, PB, SB, V)
    !
    ! This subroutine recomputes the nodes that affected by tbe cbange of V.
    !
    ! Arguments 
    ! =========
    integer, intent(in) :: n, L         ! order of matrix
    integer, intent(in) :: idx          ! 
    type(SeqB), intent(inout) :: SB     ! Data structure of B matrix
    real(wp), intent(in) :: V(n,L)      ! Hubbard-Stratonivich Field 
    type(PrefixB), intent(inout) :: PB  ! intermediate results

    ! ... local scalars    ...
    integer :: g1, g2, g3, root         ! iterator
    integer :: il, ir, nDiv
    
    ! ... Executable ...

    ! quick returne for a special case
    if (PB%nSeg .le. 0) return

    ! the normal case
    nDiv = PB%nOrth
    g1   = PB%leaflevel + (idx-1)/nDiv

    il = PB%BSeg(g1)%il
    ir = PB%BSeg(g1)%ir
    if (PB%inv) then
       call DQMC_SeqMultBi(il, ir, SB, V)
    else
       call DQMC_SeqMultB(il, ir, SB, V)
    end if
    PB%BSeg(g1)%U = SB%U
    PB%BSeg(g1)%D = SB%D
    PB%BSeg(g1)%T = SB%T

    root = DQMC_Root()
    g2 = DQMC_Parent(g1)

    ! Update the ancestors.
    do while(g2 .ne. root)
       if (DQMC_Is_Left(g1)) then
          g3 = DQMC_Right_Child(g2)

          ! Copy g1 (left child) to g2 (parent)
          PB%Bseg(g2)%U = PB%Bseg(g1)%U
          PB%Bseg(g2)%D = PB%Bseg(g1)%D
          PB%Bseg(g2)%T = PB%Bseg(g1)%T

          ! Merge g1 g3 to g2
          call DQMC_UDT_Merge(n, PB%Bseg(g2)%U, PB%Bseg(g2)%D, &
               PB%Bseg(g2)%T, PB%Bseg(g3)%U, PB%Bseg(g3)%D, PB%Bseg(g3)%T, & 
               SB%W1, SB%W2, PB%W3, SB%rw, SB%lw, SB%tau, SB%piv)

       else ! g1 is on right
          g3 = DQMC_Left_Child(g2)

          ! Copy g3 (left child) to g2 (parent)
          PB%Bseg(g2)%U = PB%Bseg(g3)%U
          PB%Bseg(g2)%D = PB%Bseg(g3)%D
          PB%Bseg(g2)%T = PB%Bseg(g3)%T

          ! Merge g1 g3 to g2
          call DQMC_UDT_Merge(n, PB%Bseg(g2)%U, PB%Bseg(g2)%D, &
               PB%Bseg(g2)%T, PB%Bseg(g1)%U, PB%Bseg(g1)%D, PB%Bseg(g1)%T, & 
               SB%W1, SB%W2, PB%W3, SB%rw, SB%lw, SB%tau, SB%piv)
       end if

       g1 = g2
       g2 = DQMC_Parent(g1)
    end do

  end subroutine DQMC_PrefixB_Update

  !----------------------------------------------------------------------!

  subroutine DQMC_PrefixMultB(il, ir, PB, SB, V)
    !
    ! Purpose
    ! =======
    !    This subroutine computes A = B_{i2}B_{i2-1}...B_{i1}
    !    and returns A's UDT decomposition.
    !
    ! Pre-assumption
    ! ==============
    !   i1 is the index of right most B, and i2 is the index for the
    !   left most B. Both i1, i2 are in [1,..L].
    !
    ! Arguments
    ! =========
    !

    integer,  intent(in)    :: il, ir          ! start/end slice 
    type(PrefixB), intent(inout)  :: PB        ! Mid products
    type(SeqB), intent(inout)     :: SB        ! 
    real(wp), intent(in)    :: V(:,:)          ! 

    ! ... local scalars    ...
    integer :: nq                      ! Number of Bs in queue
    integer :: i, idx
    integer :: fil, fir

    ! ... Executable ...

    nq = 0
    ! special case, when there is only one block
    if (PB%nSeg .eq. 0) then
       if (PB%inv) then
          call DQMC_SeqMultBi(il, ir, SB, V)
       else
          call DQMC_SeqMultB(il, ir, SB, V)
       end if
    else
       ! fix indeces. Note it is in reverse order.
       fil = ir
       fir = il
       if (fil .gt. PB%L) then
          fil = fil - PB%L
       end if
       if (fir .le. 0) then
          fir = 1
       end if

       ! Find minimum nodes to be merged
       call DQMC_UDT_Enqueue(fil, fir, PB%nOrth, PB%queue, &
            PB%leaflevel, PB%nSeg, PB%maxq, nq)

       idx = PB%queue(1)
       SB%U = PB%Bseg(idx)%U
       SB%D = PB%Bseg(idx)%D
       SB%T = PB%Bseg(idx)%T
       
       do i = 2, nq
          idx = PB%queue(i)
          call DQMC_UDT_Merge(PB%n, SB%U, SB%D, SB%T, PB%Bseg(idx)%U, &
               PB%Bseg(idx)%D, PB%Bseg(idx)%T, SB%W1, SB%W2, PB%W3, & 
               SB%rw, SB%lw, SB%tau, SB%piv)
       end do
    end if

  end subroutine DQMC_PrefixMultB

  !----------------------------------------------------------------------!

  subroutine DQMC_UDT_Enqueue(il, ir, nDiv, queue, lbnd, rbnd, maxq, nUDT)
    !
    ! Purpose
    ! =======
    !    This subroutine decides the minimum UDT to be merged 
    !    for given i1 and i2. The indices of UDTs are put into
    !    queue and the function returns the number of UDTs in queue.
    !
    !    B_{i2}....B_{i1}
    !   
    ! Assumption
    ! ==========
    !    i2 is assumed larger to i1. If not, that means the indices
    !    are round. 
    !    
    !    B_{i2}...B_{1}B_{L}...B_{i1}
    !
    ! Arguments
    ! =========
    integer, intent(in)    :: il, ir        ! start and end index
                                            ! requested interval
    integer, intent(in)    :: nDiv          ! Group size
    integer, intent(in)    :: lbnd, rbnd    ! left/right boundary of leaf nodes
    integer, intent(in)    :: maxq          ! maximum nodes in queue
    integer, intent(inout) :: queue(maxq)   ! result queue
    integer                :: nUDT

    ! ... Local scalar ...
    integer :: gl, gr, lb, rb
    integer :: i, root
    integer :: temp_queue(maxq), temp_nUDT

    ! We need to put UDTs from g2 into a temp queue, and 
    ! merge it later.

    ! Initialization
    nUDT = 0
    temp_nUDT = 0
    lb = lbnd
    rb = rbnd
    ! g1 and g2 are leaf node id of given index i1, i2
    gl   = lbnd + (il-1)/nDiv
    gr   = lbnd + (ir-1)/nDiv

    do while (.true.)
       ! If g1 and g2 are in the same group, push it into the queue.
       ! and return.

       if (gl .eq. gr) then
          ! special treatment for root node
          if (DQMC_Is_Root(gl)) then
             nUDT = 2
             root = DQMC_Root()
             if (ir .gt. il) then
                queue(1) = DQMC_Left_Child(root)
                queue(2) = DQMC_Right_Child(root)
             else
                queue(1) = DQMC_Right_Child(root)
                queue(2) = DQMC_Left_Child(root)
             end if
             return
          end if
          
          call DQMC_Push_Int(queue, gl, nUDT)
          ! Before return, we need to check the temp queue
          ! If it is not empty, merge it to the queue
          do i = temp_nUDT, 1, -1
             call DQMC_Push_Int(queue, temp_queue(i), nUDT)
          end do
          return
       end if

       ! Two are siblings. Push both and return
       if (DQMC_Right_Sibling(gl,lb, rb) .eq. gr) then
          if (DQMC_Parent(gl).ne.DQMC_Parent(gr)) then
             call DQMC_Push_Int(queue, gl, nUDT)
             call DQMC_Push_Int(queue, gr, nUDT)
             ! Before return, we need to check the temp queue
             ! If it is not empty, merge it to the queue
             do i = temp_nUDT, 1, -1
                call DQMC_Push_Int(queue, temp_queue(i), nUDT)
             end do
             return
          end if
       end if
        
       ! If gl is the left node, then push it into the queue and 
       ! use its right sibling to continue the work
       if (DQMC_Is_Right(gl)) then
          ! push gl into queue
          call DQMC_Push_Int(queue, gl, nUDT)
          ! Get gl's right sibling's parent
          gl = DQMC_Parent(DQMC_Right_Sibling(gl, lb, rb))
       else
          gl = DQMC_Parent(gl)
       end if
       
       ! If gr is the right node, then push it into the temp queue and 
       ! use its left sibling to continue the work
       if (DQMC_Is_Left(gr)) then
          ! push gr into queue
          call DQMC_Push_Int(temp_queue, gr, temp_nUDT)
          ! Get gr's right sibling's parent
          gr = DQMC_Parent(DQMC_Left_Sibling(gr, lb, rb))
       else
          ! Otherwise, just go to gr's parent
          gr = DQMC_Parent(gr)
       end if

       ! update lBound and rBound
       rb = lb - 1
       lb = lb/2

    end do

  end subroutine DQMC_UDT_Enqueue

  !----------------------------------------------------------------------!
  ! The following functions should be inlined by compiler.
  
  function DQMC_Is_Left(idx) result(isleft)
    integer, intent(in) :: idx
    logical             :: isleft

    isleft = .not. btest(idx, 0)
  end function DQMC_Is_Left

  function DQMC_Is_Right(idx) result(isright)
    integer, intent(in) :: idx
    logical             :: isright

    isright = btest(idx, 0)
  end function DQMC_Is_Right

  function DQMC_Is_Root(idx) result(isroot)
    integer, intent(in) :: idx
    logical             :: isroot

    isroot = (idx .eq. 1)
  end function DQMC_Is_Root

  function DQMC_Parent(idx) result(parent)
    integer, intent(in) :: idx
    integer             :: parent

    parent = ishft(idx, -1) 
  end function DQMC_Parent

  function DQMC_Right_Sibling(idx, lbnd, rbnd) result(rsib)
    integer, intent(in) :: idx, lbnd, rbnd
    integer             :: rsib

    rsib = idx + 1 
    if (rsib .gt. rbnd) rsib = lbnd
  end function DQMC_Right_Sibling
  

  function DQMC_Left_Sibling(idx, lbnd, rbnd) result(lsib)
    integer, intent(in) :: idx, lbnd, rbnd
    integer             :: lsib

    lsib = idx - 1 
    if (lsib .lt. lbnd) lsib = rbnd
       
  end function DQMC_Left_Sibling

  function DQMC_Root() result(root)
    integer             :: root

    root = 1
  end function DQMC_Root
  
  function DQMC_Right_Child(idx) result(rchild)
    integer, intent(in) :: idx
    integer             :: rchild

    rChild = idx*2 + 1
  end function DQMC_Right_Child
  

  function DQMC_Left_Child(idx) result(lChild)
    integer, intent(in) :: idx
    integer             :: lChild

    lChild = idx*2
  end function DQMC_Left_Child

  !----------------------------------------------------------------------!

  subroutine DQMC_Push_Int(queue, data, pt)
    integer, intent(inout) :: queue(:), pt 
    integer, intent(in)    :: data

    pt = pt + 1
    queue(pt) = data
  end subroutine DQMC_Push_Int

end module DQMC_PREFIX
