module DQMC_Phy2
#include "dqmc_include.h"

  use DQMC_UTIL
  use DQMC_STRUCT
  use DQMC_WSPACE
  use DQMC_BONDS
  use DQMC_SYMM

  implicit none 
  
  !
  ! This module contains data structure and subroutines for 
  ! pair measurement.
  !
  ! List of subroutines
  ! ===================
  !    1. DQMC_Phy2_Init(P2, nBin, nb)
  !    2. DQMC_Phy2_Avg(P2, W, T, ldt)
  !    3. DQMC_Phy2_Print(P2, OPT)
  !    4. DQMC_Phy2_GetErr(P2)
  !    5. DQMC_Phy2_Meas(n, P2, G_up, G_dn, sgn, S)
  !
  ! Data Structure and Parameters
  ! =============================
  !    The data type Phy2 performs pair measurement. For each paired sites, 
  !    it computes the convolution of their correlation and their neighbors' 
  !    correlation. (CHECK) The results are accumulated into bins and then
  !    transformed by the 'wave' function. 
  !
  type Phy2

     ! Measurement of pair
     integer  :: nb                          ! Number of bonds
     integer  :: nWave                       ! Number of wave functions
     integer  :: nBin                        ! Number of bins
     integer  :: idx                         ! current bin index
     integer  :: cnt                         ! Number of measurement for 
     integer  :: ncell
                                             ! current bin 
     integer  :: avg, err                    ! Index for average and error
                                             ! bins.
     real(wp), pointer :: sgn(:)             ! sign  
     real(wp), pointer :: M1(:,:)            ! pair measurement
     real(wp), pointer :: M2(:,:)            ! accumulated pair measurement
     real(wp), pointer :: M3(:,:)            ! averaged pair measurement
     real(wp), pointer :: M4(:,:)            ! averaged pair measurement when waves are not specified
     real(wp), pointer :: M5(:,:)            ! averaged pair measurement connected

     ! working space
     real(wp), pointer :: T(:,:) 
     integer  :: ldt
     integer  :: nData
     logical  :: compute

     logical  :: connected   
     logical  :: diagonalize=.false.

  end type Phy2

  interface DQMC_Phy2_avg 
     module procedure DQMC_Phy2_Avg_Wave, DQMC_Phy2_Avg_Symm
  end interface

  interface DQMC_Phy2_GetIrrep
     module procedure DQMC_Phy2_GetIrrep_Connected, DQMC_Phy2_GetIrrep_Full
  end interface
contains

  ! Subroutines
  ! ==================================================================
  
  subroutine DQMC_Phy2_Init(P2, nBin, S, WS, meas)
    !
    ! Purpose
    ! =======
    !    This subroutine initializes Phy2.
    !
    ! Arguments
    ! =========
    !
    type(Phy2), intent(inout) :: P2          ! Phy2 to be initialized
    integer, intent(in)       :: nBin        ! No of bin
    type(Struct), intent(in)  :: S
    type(Wspace), intent(in), target :: WS   ! working space
    logical, intent(out)      :: meas

    ! ... Local Vars ...
    integer :: nb, nWave
    real(wp), pointer :: dum(:,:) 

    ! ... Executable ...

    P2%compute = S%checklist(STRUCT_BOND)
    meas = P2%compute

    if (P2%compute) then
       nWave    = S%nWave
       nb     = S%n_b

       P2%nBin  = nBin
       P2%nb    = nb
       P2%nwave = nWave
       P2%ncell = S%ncell

       ! Allocate storages
       !! Allocate two additional bins for storing 
       !! average value and error
       allocate(P2%M1(nb,nb))
       allocate(P2%M2(nb,nb))
       allocate(P2%sgn(nBin+2))
       if(S%checklist(STRUCT_WAVE))then
         allocate(P2%M3(nWave,nBin+2))
         nullify(P2%M4, P2%M5)
       else
         allocate(P2%M3(S%nClass_b,nBin+2))
         allocate(P2%M4(nWave,nbin+2))
         allocate(P2%M5(nWave,nbin+2))
         P2%diagonalize=.true.
       endif
       
       ! Initialize 
       P2%M2      = ZERO
       P2%sgn     = ZERO
       P2%avg     = nBin+1
       P2%err     = nBin+2
       P2%cnt     = 0
       P2%idx     = 1
       
       dum  =>  WS%R4
       !P2%ldt     = size(WS%R4, 1)
       allocate(P2%T(nb, nWave))
       P2%ldt     = nb
       P2%nData   = nWave + 1
    else
       P2%nData = 0
    end if

  end subroutine DQMC_Phy2_Init

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy2_Free(P2)
    !
    ! Purpose
    ! =======
    !    This subroutine frees Phy2.
    !
    ! Arguments
    ! =========
    !
    type(Phy2), intent(inout) :: P2      ! Phy0 to be initialized

    ! ... Executable ...

    if (P2%compute) then
       deallocate(P2%sgn, P2%M1, P2%M2, P2%M3, P2%T)
       if(associated(P2%M4)) deallocate(P2%M4)
       if(associated(P2%M5)) deallocate(P2%M5)
    end if
 
  end subroutine DQMC_Phy2_Free

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy2_Meas(n, M1, M2, P2, Bond, G_up, G_dn, sgn)
    !
    ! Purpose
    ! =======
    !    This subroutine performs physics measurement on
    !    the neighbors of pairs of sites 
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)          :: n             ! Number of sites
    real(wp), intent(inout)      :: M1(:, :)      ! temp variables
    real(wp), intent(inout)      :: M2(:, :)      !
    type(Phy2), intent(inout)    :: P2            ! Phy2
    type(CCS), intent(in)        :: Bond
    real(wp), intent(in)         :: G_up(n,n)     ! Green's function
    real(wp), intent(in)         :: G_dn(n,n)     ! for spin up and down
    real(wp), intent(in)         :: sgn           ! Sgnup*sgndn
 
    ! ... local scalar ...
    integer  :: nb
    real(wp) :: factor

    if (P2%compute) then

       ! Compute the pair measurement       
       call DQMC_Phy2_Pair(n, M1, G_up, G_dn, Bond)

       ! Averaging and accumulation
       factor = sgn / P2%ncell
       nb = P2%nb
       call blas_daxpy(nb*nb, factor, M1, 1, M2, 1)
       P2%sgn(P2%idx) = P2%sgn(P2%idx) + sgn
       !write(*,*) 'Phy2_Meas',P2%ncell   
       
       ! increase counter
       p2%cnt = P2%cnt + 1
    end if

  end subroutine DQMC_Phy2_Meas

  !--------------------------------------------------------------------!
  
  subroutine DQMC_Phy2_Pair(n, M, up, dn, Bond)
    !
    ! Purpose
    ! =======
    !    This subroutine computes pair measurements.
    !    For each pair (i,j), compute the effect of 
    !    all their neighbor pairs.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)          :: n             ! Number of sites
    real(wp), intent(inout)      :: M(:, :)       ! pair measurement
                                                  ! (P_dd')
    real(wp), intent(in)         :: up(:,:)       ! Green's function
    real(wp), intent(in)         :: dn(:,:)       ! for spin up and down
    type(CCS), intent(in)        :: Bond          ! bond info
 
    ! ... local variables ...
    integer  :: i, j, inbr, jnbr, ni, nj, nci, ncj
    integer, pointer  :: start(:) 
    integer, pointer  :: r(:) 
    integer, pointer  :: A(:) 

    ! ... Executable ...

    ! alias
    start => Bond%cstart
    r     => Bond%row
    A     => Bond%A
    
    ! clean up
    M = ZERO

    do i = 1, n
       do inbr = start(i), start(i+1)-1
          ni  = r(inbr)              ! bond of site i
          nci = A(inbr)              ! bond class of (i,ri)
          do j = 1, n
             do jnbr = start(j), start(j+1)-1
                nj  = r(jnbr)        ! bond of site j
                ncj = A(jnbr)        ! bond class of (j,rj)
                M(ncj, nci) = M(ncj, nci) + up(nj, ni)*dn(j, i)
             end do
          end do
       end do
    end do

  end subroutine DQMC_Phy2_Pair

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy2_Avg_Wave(P2, W)
    !
    ! Purpose
    ! =======
    !    This subroutine averges the pair measurements.
    !    which is stored in P2%M2.
    !
    ! Arguments
    ! =========
    !
    type(Phy2), intent(inout) :: P2                 ! phy2
    real(wp), intent(in)      :: W(P2%nb,P2%nwave)! wave

    ! ... local scalar ...
    real(wp) :: factor

    ! ... Executable ...

    if (P2%compute) then
       factor = ONE/P2%cnt       
       ! CC
       !write(*,*) 'Avg_Wave', P2%M2(:,:)
       call DQMC_Wave_Avg(P2%nb, P2%nWave, W, P2%T, P2%M2, &
           P2%M3(1:P2%nWave,P2%idx), factor, P2%ldt)
       P2%sgn(P2%idx) = P2%sgn(P2%idx) * factor
    
       ! Reset counter and change bins
       P2%idx = P2%idx + 1
       P2%cnt = 0
       P2%M2  = ZERO
    end if

  end subroutine DQMC_Phy2_Avg_Wave

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy2_Avg_Symm(P2, S)
    !
    ! Purpose
    ! =======
    !    This subroutine averges the pair measurements.
    !    which is stored in P2%M2.
    !
    ! Arguments
    ! =========
    !
    type(Phy2), intent(inout) :: P2                 ! phy2
    type(struct), intent(in)  :: S

    ! ... local scalar ...
    real(wp) :: factor

    ! ... Executable ...

    if (P2%compute) then
       factor = ONE/P2%cnt       
       call DQMC_Pair_Symm(P2%M2, P2%M3(1:S%nclass_b,P2%idx), &
            S%class_b, S%size_b, S%n_b, S%nclass_b, factor)
       P2%sgn(P2%idx) = P2%sgn(P2%idx) * factor
    
       ! Reset counter and change bins
       P2%idx = P2%idx + 1
       P2%cnt = 0
       P2%M2  = ZERO
    end if

  end subroutine DQMC_Phy2_Avg_Symm

  !--------------------------------------------------------------------!

  subroutine DQMC_Wave_Avg(nb, nWave, W, T, I, O, factor, ldt)
    !
    ! Purpose
    ! =======
    !    This subroutine computes the waves from input data.
    !
    !    The averaging process runs as follows. 
    !    1. Let T = I*W
    !    2. Compute O(i) = W(:,i)'*T(:,i)
    !    3. Averaging O  = O*factor
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)     :: nB              ! number of bonds
    integer, intent(in)     :: nWave           ! number of waves
    real(wp), intent(in)    :: W(nb, nWave)    ! wave matrix
    real(wp), intent(in)    :: T(ldt, nWave)   ! temp matrix
    real(wp), intent(in)    :: I(nb, nb)       ! input matrix
    real(wp), intent(inout) :: O(:)            ! output matrix
    real(wp), intent(in)    :: factor
    integer, intent(in)     :: ldt

    ! ... Local Variables ...
    integer :: j

    ! ... BLAS function ...
    real(wp), external :: blas_ddot

    ! ... Executable ....

    call blas_dgemm('N','N', nb, nWave, nb, factor, I, nb, W, nb, &
         ZERO, T, ldt)
    do j = 1, nWave
       O(j) = blas_ddot(nb, W(1:nb,j), 1, T(1:nb,j), 1)
    end do

  end subroutine DQMC_Wave_Avg

  !--------------------------------------------------------------------!

  subroutine DQMC_Pair_Symm(M2, M3, class_b, size_b, n, nclass, factor)
    !
    ! Purpose
    ! =======
    !    This subroutine symmetrizes the pairing matrix.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)         :: n, nclass
    integer, intent(in)         :: size_b(nclass)
    integer, intent(in)         :: class_b(n,n)
    real(wp), intent(in)        :: M2(n,n)
    real(wp), intent(in)        :: factor
    real(wp), intent(out)       :: M3(nclass)

    ! ... Local Variables ...
    integer :: i, j, ij

    ! ... Executable ...

    !Symmetrize
    M3=0.d0
    do i=1,n
      do j=1,n
        ij=class_b(i,j)
        M3(ij)=M3(ij)+M2(i,j)
      enddo
    enddo

    !Normalize
    M3=factor*M3/size_b

  end subroutine

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy2_GetErr(P2)
    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !    This subroutine computes the average and errors
    !    of measurements. After this routine is called
    !    bins are modified to contain data from the 
    !    jack-knived distribution.
    !
    ! Argument
    ! ========
    type(Phy2), intent(inout) :: P2

    ! ... Local Scalar ...
    integer  :: i
    integer  :: n, avg, err, m
    integer  :: nproc

    ! ... Local Array
    real(wp) :: sgn(P2%nBin) , sum_sgn, y(P2%nBin), data(P2%nBin)
    
    ! ... Executable ...

    if (P2%compute) then
       n   = P2%nBin
       avg = P2%avg
       err = P2%err
       m   = size(P2%M3,1)
       nproc = qmc_sim%size
       
       if (nproc == 1) then
          ! Averge the other terms
          data = P2%sgn(1:n)
          call DQMC_JackKnife(n, P2%sgn(avg), P2%sgn(err), data, y, sgn, sum_sgn)
          do i = 1, m
             data = P2%M3(i, 1:n)
             call DQMC_SignJackKnife(n, P2%M3(i, avg), P2%M3(i, err), &
                  data, y, sgn, sum_sgn)
          end do

          !Store Jackknife in bins
          P2%sgn(1:n) = sgn(1:n) / dble(n-1)
          do i = 1, m
             P2%M3(i,1:n) = (sum_sgn*P2%M3(i,avg) - P2%M3(i,1:n)) / dble(n-1)
          enddo

       else

#         ifdef _QMC_MPI
             !Sum sign
             call mpi_allreduce(P2%sgn(1), P2%sgn(avg), 1, mpi_double, &
                mpi_sum, mpi_comm_world, n)

             !Sum properties
             call mpi_allreduce(P2%M3(:,1), P2%M3(:,avg), m, mpi_double, &
                mpi_sum, mpi_comm_world, n)

             !Compute averages over n-1 processors
             P2%M3(:,1) = (P2%M3(:,avg) - P2%M3(:,1)) / dble(nproc - 1)
             P2%sgn(1) = (P2%sgn(avg) - P2%sgn(1)) / dble(nproc - 1)

             !Store physical average amongst all processors
             P2%M3(:,avg) = P2%M3(:,avg) / P2%sgn(avg) 

             !Store jackknife in the processor bin
             P2%M3(:,1)   = P2%M3(:,1) / P2%sgn(1) 

             !Compute error
             call mpi_allreduce(P2%M3(:,1)**2, P2%M3(:,err), m, mpi_double, &
                mpi_sum, mpi_comm_world, n)
             P2%M3(:,err) = P2%M3(:,err) / dble(nproc) - P2%M3(:,avg)**2 
             P2%M3(:,err) = sqrt(P2%M3(:,err) * dble(nproc-1))
#         endif

       endif
       
    end if
    
  end subroutine DQMC_Phy2_GetErr

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy2_GetIrrep_Full(P2, S)
    !
    !  Interfaced by DQMC_Phy2_GetIrrep
    !

    type(Phy2), intent(inout)   :: P2
    type(Struct), intent(inout) :: S


    ! ... Local variables ...
    integer :: i, j, ij, nb
    real(wp)  :: work(3*P2%nb),tmp(P2%nWave)

    ! ... Executable ...

    nb = P2%nb

    !make sure S%W is allocated
    if(associated(S%W))deallocate(S%W)
    allocate(S%W(nb,nb))

    !Load Full pairing matrix
    do i = 1, nb
       do j = 1, nb
          ij = S%class_b(i,j)
          S%W(i,j) = P2%M3(ij,P2%avg)
       enddo
    enddo

    !Compute waves as those diagonalizing the pairing matrix
    call lapack_dsyev('V','U', nb, S%W, nb, tmp, work, 3*nb, i)

    P2%connected=.false.

  end subroutine DQMC_Phy2_GetIrrep_Full

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy2_GetIrrep_Connected(P2, G_fun, S)
    !
    !  Interfaced by DQMC_Phy2_GetIrrep
    !

    type(phy2), intent(inout)   :: P2
    type(Struct), intent(inout) :: S
    real(wp), intent(in)        :: G_fun(S%nclass)

    integer   :: i, j, ij, n, nb
    real(wp)  :: up(S%nSite,S%nSite), dn(S%nSite,S%nSite)
    real(wp)  :: work(3*P2%nb)
    real(wp)  :: tmp(P2%nWave)


    ! ... Executable ...
    n = S%nSite
    nb = P2%nb

    !Load disconnected pairing matrix in M1.
    do i = 1, n
       do j = 1, n
          up(i,j) = S%gf_phase(i,j) * G_fun(S%D(i,j))
       enddo
    enddo
    dn = up
    call DQMC_Phy2_Pair(n, P2%M1, up, dn, S%B)

    !make sure S%W is allocated properly
    if (associated(S%W)) deallocate(S%W)
    allocate(S%W(nb,nb))

    !Expand pairing matrix to all bond pairs
    do i = 1, nb
       do j = 1, nb
          ij = S%class_b(i,j)
          S%W(i,j) = P2%M3(ij,P2%avg)
       enddo
    enddo

    !Compute connected Pairing function
    S%W = S%W - P2%M1

    !Compute waves as those diagonalizing the pairing matrix
    call lapack_dsyev('V','U', nb, S%W, nb, tmp, work, 3*n, i)

    P2%connected = .true.

  end subroutine DQMC_Phy2_GetIrrep_Connected

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy2_GetErrIrrep(P2, G_fun, S)
    use dqmc_mpi

    type(struct), intent(in)   :: S
    type(phy2), intent(inout)  :: P2
    real(wp), intent(in)       :: G_fun(S%nClass,P2%err)

    ! ... Local Variables ...
    integer :: avg, err, nbin, ibin, i, j, k, n, nproc
    real(wp)  :: factor, up(S%nSite,S%nSite), dn(S%nSite,S%nSite)
    real(wp), pointer :: MX(:,:)


    ! ... Executable ...

    nBin=P2%nBin
    avg=P2%avg
    err=P2%err
    nproc = qmc_sim%size

    factor=1.d0

    p2%M4(:,err) = 0.d0
    p2%M5(:,err) = 0.d0

    n=S%nSite

    !Fill bins with wave amplitudes using S%W from diagonalization
    do ibin = 1, avg

       !Create JK bin for Full pairing matrix (M2)
       do i = 1, P2%nb
          do j = 1, P2%nb
             k = S%class_b(i,j)
             P2%M2(i,j) = P2%M3(k,ibin)
          enddo
       enddo

       !Compute waves amplitute for bin "ibin". Full pairing function.
       call DQMC_Wave_Avg(P2%nb,P2%nWave,S%W,P2%T,P2%M2,P2%M4(1:P2%nWave,ibin),factor,P2%ldt)

       !create JK bin for G
       do i = 1, n
          do j = 1, n
             k = S%D(i,j)
             up(i,j) = S%gf_phase(i,j)*G_fun(k,ibin)
          enddo
       enddo
       dn = up

       !Construct JK bin for G*G
       call DQMC_Phy2_Pair(n, P2%M1, up, dn, S%B)
       P2%M1 = P2%M1 / n

       !Construct connected function
       P2%M1 = P2%M2 - P2%M1

       !Compute waves amplitute for bin "ibin": Connected pairing function.
       call DQMC_Wave_Avg(P2%nb,P2%nWave,S%W,P2%T,P2%M1,P2%M5(1:P2%nWave,ibin),factor,P2%ldt)

    enddo

    if (nproc > 1) then

          MX => P2%M4
          do i = 1, 2
             j = size(MX, 1)
             !Compute errorbars averaging over processors. 
#            ifdef _QMC_MPI
                call mpi_allreduce(MX(:,1)**2, MX(:,err), j, mpi_double, &
                   mpi_sum, mpi_comm_world, k)
#            endif
             MX(:,err) = MX(:,err) / dble(nproc) - MX(:,avg)**2 
             MX(:,err) = sqrt(MX(:,err) * dble(nproc-1))
             MX => P2%M5
          enddo

    else

       !Compute errors
       do i=1,P2%nWave
          !Full Pairing matrix
          P2%M4(i,err)=sqrt((nBin-1)*sum((P2%M4(i,1:nBin)-P2%M4(i,avg))**2)/nBin)
          !Connected Pairing matrix
          P2%M5(i,err)=sqrt((nBin-1)*sum((P2%M5(i,1:nBin)-P2%M5(i,avg))**2)/nBin)
       enddo

    endif


  end subroutine DQMC_Phy2_GetErrIrrep

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy2_Print(P2, wlabel, OPT)
    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !    This subroutine prints out the average and errors
    !    of pair measurements. The first three terms only.
    !
    !  Pre-assumption
    ! ===============
    !    OPT is a file handle
    !    DQMC_Phy2_GetErr was called.
    !    nWave >= 3
    !
    ! Arguments
    ! =========
    !
    type(Phy2), intent(in)    :: P2   ! Phy2
    character(*), intent(in)  :: wlabel(:)
    integer, intent(in)       :: OPT  ! Output file handle

    ! ... Local var
    integer  :: avg, err
    
    ! ... Executable ...

    if (qmc_sim%rank .ne. 0) return

    if (P2%compute) then
       avg = P2%avg
       err = P2%err
       call DQMC_Print_RealArray(0, P2%nWave, "Pair-field correlation function - accumulated:", &
            wlabel, P2%M3(:,avg:avg), P2%M3(:,err:err), OPT)
    end if
    
  end subroutine DQMC_Phy2_Print

  !--------------------------------------------------------------------!

  subroutine DQMC_phy2_WaveSymm(S, P2, Symm)

  ! INPUT
  ! W(nb,ns)     :  wave matrix to be analyzed
  ! nw           :  number of wave (number of rows of W)
  ! nb           :  number of bonds (number of column of W)
  ! wa(nw)       :  wave amplitude (expectation value out of QMC)
  ! maps(nb,ns)  :  returns how bonds are mapped under symm operations
  ! ns           :  number of symmetry operations

  ! OUTPUT
  ! wrepr(nw)  : irreducible representation of the wave
  ! wclass(nw) : group of the wave (a group is defined as the set of waves
  !              belonging to equivalent irreps)
  ! wlabel     : label of the wave (contains character of representation)
  ! nirrep     : number of irreducible representations
  ! nwclass    : number of inequivalent irreducible representations

  type(struct), intent(inout)        :: S
  type(phy2), intent(inout)          :: P2
  type(Symm_operations), intent(in)           :: Symm

  real(wp), pointer   :: W(:,:) 
  real(wp), pointer   :: wa(:) 
  integer, pointer    :: maps(:,:) 
  integer             :: iw, ir, jr, is, ib, jw, nirrep, nw, nb, ns
  integer             :: reprlabel(S%nwave),rdim(S%nwave)
  real(wp)            :: newWave(S%n_b), WaveSymm(S%nwave,S%nwave,Symm%nsymm), irrepChar(S%nwave,Symm%nsymm)
  logical             :: irrep(S%nwave)

  W=>S%W
  maps=>Symm%map_symm_p

  if(P2%connected)then
     wa=>P2%M5(:,P2%avg)
  else
     wa=>P2%M4(:,P2%avg)
  endif

  nw=S%nwave
  nb=S%n_b
  ns=Symm%nsymm

  if(.not.associated(S%wrepr)) allocate(S%wrepr(nw))

  S%wrepr=0
  nirrep=0
  rdim=0
  !Check degeneracy of eigenvalues to find the dimension of the representation
  do iw=1,nw
    if(S%wrepr(iw)==0)then
      nirrep=nirrep+1
      S%wrepr(iw)=nirrep
      rdim(nirrep)=1
      do jw=iw+1,nw
        if(S%wrepr(jw)==0.and.abs(wa(iw)-wa(jw))<1.d-6)then
          S%wrepr(jw)=nirrep
          rdim(nirrep)=rdim(nirrep)+1
        endif
      enddo
    endif
  enddo
        
  !Look at how waves transform under the action of symmetry
  WaveSymm=0.d0
  do is=1,ns
    do iw=1,nw
      !Apply symmetry "is" to wave "iw"
      do ib=1,nb
        newWave(ib)=W(maps(ib,is),iw)
      enddo
      !Dot-product between all waves (jw) and wave "iw" after transformation by "is"
      !WaveSymm(iw,jw,is)=0 means "iw", transformed by "is", has no component on "jw"
      do jw=1,nw
        WaveSymm(iw,jw,is)=sum(W(:,jw)*newWave(:))
      enddo
    enddo
  enddo

  !Check whether dimension of representation is compatible with non-zero  
  !elements in WaveSymm
  irrep=.true.
  do ir=1,nirrep
    do is=1,ns
      !Compare all pairs of waves that do not belong to the same irrep.
      !If WaveSymm is larger than 0 for any symmetry then what we found was not
      !an irrep. Accidental degeneracy.
      do jw=1,nw
        if(S%wrepr(jw)==ir)then
          do iw=1,nw
            if(S%wrepr(jw)/=ir.and.abs(WaveSymm(jw,iw,is))>1.d-10)irrep(ir)=.false.
          enddo
        endif
      enddo
    enddo
  enddo

  !Compute character of each irreducible representation
  irrepChar=0.d0
  do ir=1,nirrep
    if(irrep(ir))then
      do is=1,ns
        do iw=1,nw
          if(S%wrepr(iw)==ir) irrepChar(ir,is) = irrepChar(ir,is) + WaveSymm(iw,iw,is)
        enddo
      enddo
    endif
  enddo

  !Find which representations are equivalent (have the same character table)
  if (.not.associated(S%wclass)) allocate(S%wclass(nirrep))
  S%nwclass=0
  reprlabel=0
  do ir=1,nirrep
    if(irrep(ir))then
      do jr=1,ir-1
        if(irrep(jr))then
          if(sum((irrepChar(ir,1:ns)-irrepChar(jr,1:ns))**2)<1.d-6)exit
        endif
      enddo
      if(jr==ir)then
        S%nwclass=S%nwclass+1
        S%wclass(ir)=S%nwclass
      else
        S%wclass(ir)=S%wclass(jr)
      endif
    endif
  enddo

  !Assign to each wave a group (which is made by equivalent representations)
  do iw=1,nw
    jr=S%wrepr(iw)
    write(S%wlabel(iw),'(10(i3))') rdim(jr),(nint(irrepChar(jr,is)),is=1,ns)
  enddo

  S%nirrep=nirrep

  end subroutine

  !--------------------------------------------------------------!

  subroutine dqmc_phy2_PrintSymm(S, P2, OPT)
  use dqmc_mpi

  type(struct),   intent(inout) :: S
  type(phy2),     intent(inout) :: P2
  integer       , intent(in)    :: OPT
  integer :: iw, jw, ir, jr, ifw, avg, err, n
  real(wp)  :: Trace3(P2%err), Trace4(P2%err)

  if (qmc_sim%rank .ne. 0) return

  !aliases
  err = P2%err
  avg = P2%avg
  n   = P2%nbin  

  !Print Waves if those were found by diagonalization
  write(OPT,'(A)')" Waves coefficients (Bonds label columns and are defined above)"
  write(OPT,'(3x,100(3x,i2,2x))')(iw,iw=1,S%n_b)
  do iw=1,S%nWave
    write(OPT,'(i2,1x,20(1x,f6.3))')iw,(S%W(ir,iw),ir=1,S%n_b)
  enddo
  write(OPT,'(A)')
  write(OPT,"(76('='))")

  !Print wave amplitude ordered by symmetry
  write(OPT,'(A)')" Pair measurement: "

  !Loop over groups ...
  do ir=1,S%nwclass

    Trace3=0.d0; Trace4=0.d0
    write(OPT,'(A, i2)')' Irreducible representation',ir

    !Find a wave that represents this group ...
    do jw=1,S%nwave
      if(S%wclass(S%wrepr(jw))==ir)then
         !Set the format and write the header
         write(OPT,'(A,A)') '  Dimension and characters :',S%wlabel(jw)
         write(OPT,'(23x,A4,23x,A9)') 'Full','Connected' 
         exit
      endif
    enddo

    ! ... and print the equivalent irreps belonging to group "ir"
    Trace3=0.d0; Trace4=0.d0
    do jr=1,S%nirrep
      if(S%wclass(jr)==ir)then
        do iw=1,S%nwave 
          if(S%wrepr(iw)==jr)then
             Trace3(1:avg) = Trace3(1:avg) + P2%M4(iw,1:avg)
             Trace4(1:avg) = Trace4(1:avg) + P2%M5(iw,1:avg)
             write(OPT,'(2x,i2)',advance='no') iw
             ifw=iw
          endif
        enddo
        write(OPT,"(' :',f10.6,' +-',f10.6,6x,f10.6,' +-',f10.6)")P2%M4(ifw,avg), P2%M4(ifw,err), P2%M5(ifw,avg), P2%M5(ifw,err)
      endif
    enddo

    !Also compute and print the Trace +- err
    Trace3(err)=sqrt(sum((Trace3(1:n)-Trace3(avg))**2)*dble(n-1)/n)
    Trace4(err)=sqrt(sum((Trace4(1:n)-Trace4(avg))**2)*dble(n-1)/n)
    write(OPT,100)'Trace',Trace3(avg),Trace3(err),Trace4(avg),Trace4(err)
    write(OPT,*)

  enddo
  write(OPT,"(76('='))")

  100 format(3x,A5,3x," :",f10.6,' +-',f10.6,6x,f10.6,' +-',f10.6)

  end subroutine DQMC_phy2_PrintSymm

  !-----------------------------------------------------------------------!


end module DQMC_Phy2
