module DQMC_TDM2
#include "dqmc_include.h"

  use DQMC_UTIL
  use DQMC_STRUCT
  use DQMC_Phy2
  use DQMC_TDM1

  implicit none 
  
  !
  ! This module is to compute susceptibility.
  ! Susceptibility is stored in a 3 dimensional array. The first two 
  ! dimensions are for pair of neighbors of any two pair of particles.
  ! The third dimension is for different time slice.
  ! When aggregated, 
  !
  
  type TDM2
     real(wp), pointer :: Pair   (:,:,:)      ! working
     real(wp), pointer :: Apair  (:,:,:)      ! accumulated 
     real(wp), pointer :: Bpair  (:,:,:)      ! sus in wave form 
     real(wp), pointer :: Npair  (:,:,:)      ! no vertex sus
     
     integer  :: n
     integer  :: nClass
     integer  :: nB
     integer  :: nWave
     integer  :: L
     integer  :: itvl
     real(wp) :: dtau

     type(CCS),pointer :: Bond          ! bonds
     integer,  pointer :: D(:,:,:)      ! nonvertex info
     integer, pointer  :: SD(:,:)       ! neighbor info
     real(wp), pointer :: Wave(:,:) 
     real(wp), pointer :: T(:,:)        ! working space for getting wave
     integer  :: ldt

     character(label_len), pointer :: label(:) 
  end type TDM2
  
contains

  ! Subroutines
  ! ==================================================================
  
  subroutine DQMC_TDM2_Init(n, L, nBin, T2, S, T)
    !
    ! Purpose
    ! =======
    !    This subroutine initializes TDM2.
    !
    ! Arguments
    ! =========
    !
    type(TDM2), intent(inout) :: T2       ! time dependent measurement
    integer, intent(in)       :: n        ! No. of sites
    integer, intent(in)       :: L        ! No of time slice
    integer, intent(in)       :: nBin     ! No of Bins
    type(Struct), intent(in)  :: S
    real(wp), intent(in)      :: T(:,:)
    target  :: S, T

    ! ... local variables ...
    integer     :: i, nb1, nb2, n1, n2, r1, r2, nb, nWave

    ! ... Executable ...

    T2%n     = n
    nB       = S%n_b
    nWave    = S%nWave
    T2%nWave = S%nWave
    T2%nB    = nB
    T2%L     = L

    ! Allocate storages
    allocate(T2%pair(nB,nB,L+1))
    allocate(T2%Apair(nB,nB,L+1))
    allocate(T2%Bpair(nWave,L+1,nBin+2))
    allocate(T2%Npair(nWave,L+1,nBin+2))
    allocate(T2%D(nB,nB,n))

    T2%Wave => S%W
    T2%Bond => S%B
    T2%SD   => S%D
    T2%label=> S%wlabel
    T2%T    => T
    T2%ldt  = size(T,1)
    
    T2%pair  = ZERO
    T2%Apair = ZERO
    
    ! Build D table, which is used in computing no vertex suscetibility.
    ! D is a 3 dim array nB*nB*nClass.
    ! The meaning of entry D(i,j,k) is as follows.
    ! In a distance class k, there are two sites, say site1 and site2.  
    ! D(i,j,k) is the distance class of the ith neighbor of site1 and the
    ! jth neighbor of site2.
    ! Both i and j are indexed from 1 to 9. 

    do i = 1, n
       ! we use two sites to simulate the neighboring 
       ! structure on a reduced lattice. one is site 1, another is 
       ! the site in the reduced geometry.

       ! For all site1's neighbors and all site2's neighbors
       ! find their cooresponding 
       do nb1 = S%B%cstart(1), S%B%cstart(2)-1
          n1 = S%B%A(nb1)
          r1 = S%B%row(nb1)
          do nb2 = S%B%cstart(i), S%B%cstart(i+1)-1
             n2 = S%B%A(nb2)
             r2 = S%B%row(nb2)
             T2%D(n2, n1, i) = S%D(r2, r1)
          end do
       end do
    end do

  end subroutine DQMC_TDM2_Init

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM2_Free(T2)
    !
    ! Purpose
    ! =======
    !    This subroutine frees TDM2.
    !
    ! Arguments
    ! =========
    !
    type(TDM2), intent(inout) :: T2      ! TDM to be freed

    ! ... Executable ...

    deallocate(T2%Pair, T2%Apair, T2%Bpair, T2%Npair, T2%D)

  end subroutine DQMC_TDM2_Free

  !--------------------------------------------------------------------!
  !--------------------------------------------------------------------!

  subroutine DQMC_TDM2_Meas(T2, upt0, up0t, dnt0, dn0t, ti)
    !
    ! Purpose
    ! =======
    !    This subroutine frees TDM2.
    !
    ! Arguments
    ! =========
    !
    type(TDM2), intent(inout)    :: T2
    real(wp), intent(in)         :: up0t(:,:), upt0(:,:)
    real(wp), intent(in)         :: dnt0(:,:), dn0t(:,:)
    integer, intent(in)          :: ti

    ! ... Local scalar ...

    integer  :: n

    ! ... Executable ...

    ! Initialization
    n     = T2%n

    !  This computes the following code
    !  For all time slice ti=1,2,...L
    !     For all pair of sites (i,j)
    !        For all pair of (i's neighbors, j's neighbor) = (k,l)
    !           susnl(N_i,N_j, ti) += upt0(k,l,ti)*dnt0(i,j,ti)
    !
    !  The last time slice needs special treatment
    !      susnl(N_i,N_j, L+1) += up0t(k,l,1)*dn0t(i,j,1)
    !
    
    call DQMC_Phy2_Pair(n, T2%pair(:,:,ti), upt0(:,:), dnt0(:,:), T2%Bond)
    
    if (ti .eq. 1) then
       call DQMC_Phy2_Pair(n, T2%pair(:,:,ti), up0t(:,:), dn0t(:,:), T2%Bond)
    end if

  end subroutine DQMC_TDM2_Meas

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM2_Acc(T2, sgn)
    !
    ! Purpose
    ! =======
    !    This subroutine frees TDM2.
    !
    ! Arguments
    ! =========
    !
    type(TDM2), intent(inout)    :: T2
    real(wp), intent(in)         :: sgn

    ! ... Local scalar ...

    integer  :: n
    real(wp) :: factor

    ! ... Executable ...

    ! Initialization
    n     = (T2%L+1)*T2%nB*T2%nB

    ! Average for sus
    factor = sgn/T2%n

    ! Apair += factor*pair
    call blas_daxpy(n,factor,T2%pair,1,T2%Apair,1)

    ! accumulate sign
    T2%pair = ZERO

  end subroutine DQMC_TDM2_Acc
  
  !--------------------------------------------------------------------!

  subroutine DQMC_TDM2_Avg(T2, gnl, idx, factor)
    !
    ! Purpose
    ! =======
    !    This subroutine averges the pair measurements.
    !    which is stored in T2%sus2.
    !    The averaging process runs as follows. 
    !    1. Let T = sus2*W
    !    2. Compute b(i) = W(:,i)'*T(:,i)
    !    3. Averaging b  = b/nMeas 
    !    4. Store b into sus3(:,idx)
    !
    ! Arguments
    ! =========
    !
    type(TDM2), intent(inout) :: T2                 ! T2
    real(wp), intent(in)      :: gnl(:,:)
    real(wp), intent(in)      :: factor
    integer, intent(in)       :: idx

    ! ... BLAS function ...
    real(wp), external :: ddot

    ! ... local scalar ...
    real(wp), pointer :: susnl(:,:,:) 
    integer, pointer  :: D(:,:,:) 
    integer, pointer  :: SD(:,:) 
    
    integer  :: ti, i, j, nWave, L, nB, inbr, jnbr
 
    ! ... Executable ...
    nWave  =  T2%nWave
    L      =  T2%L 
    nB     =  T2%nB
    susnl  => T2%pair
    D      => T2%D
    SD     => T2%SD
    
    ! For suscetibility
    do j = 1, L+1
       call DQMC_Wave_Avg(nB, nWave, T2%Wave, T2%T, T2%Apair(:,:,j), &
            T2%Bpair(1:nWave,j,idx), factor, T2%ldt)
    end do
    
    ! For non vertex suscetibility
    ! we use averaged G_nl to compute the sus 
    susnl = ZERO
    do ti = 1, L+1
       do i = 1, T2%n
          do inbr = 1, nB
             do jnbr = 1, nB
                susnl(jnbr,inbr,ti) = susnl(jnbr,inbr,ti) + &
                     gnl(D(jnbr, inbr, i), ti) * gnl(SD(i,1),ti)
             end do
          end do
       end do
    end do

    ! non vertex FT 
    do j = 1, L+1
       call DQMC_Wave_Avg(nB, nWave, T2%Wave, T2%T, susnl(:,:,j), &
            T2%Npair(1:nWave,j,idx), ONE, T2%ldt)
    end do
    
    ! Sanitize 
    T2%Apair = ZERO

  end subroutine DQMC_TDM2_Avg

  !--------------------------------------------------------------------!

end module DQMC_TDM2
