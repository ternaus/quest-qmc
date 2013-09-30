module dqmc_simqmc

  use dqmc_Hubbard
  use dqmc_Gtau

  implicit none 

  ! this module will be used as a DMFT solver.
  !
  ! Data type
  ! =========
  type simqmc
     type(Hubbard)       :: Hub            ! Hubbard Model  
     integer             :: nc, na         ! n = nc*(na+1)
     type(Gtau)          :: tau
     
     ! working space for inverting M
     real(wp), pointer   :: work(:)
     integer, pointer    :: ipiv(:)
     integer             :: lwork

  end type simqmc
  
contains

  ! =============================================================== !

  subroutine DQMC_QMC_Init(qmc, nc, na, L, U, dtau, K, nwarm, npass, &
       seed, initHSF)
    !
    ! Purpose
    ! =======
    !    This subrotine initializes QMC simulation
    !
    ! Pre-assumption
    ! ==============
    !    matrix K should be nxn, where n = nc*(1+na)
    !
    ! Arguments
    ! =========
    !
    type(simqmc), intent(inout)   :: qmc           ! simulator
    integer, intent(in)           :: nc, na, L     ! scale parameters
    real(wp), intent(in)          :: U, dtau       ! Hubbard parameters
    real(wp), intent(inout)       :: K(:,:)        ! K matrix
    integer, intent(in)           :: nwarm, npass  ! simulation param
    integer, intent(in)           :: seed          ! random seed
    logical, intent(in)           :: initHSF
    target :: qmc 
    
    ! ... Local variables ...
    integer :: n, val(8), i, j
    type(Hubbard), pointer :: Hub 
    real(wp) :: temp, lambda
    integer,  parameter :: nort  = 12
    integer,  parameter :: nwrap = 12
    integer,  parameter :: nmeas= 5
    integer,  parameter :: nbin = 10
    real(wp), parameter :: diff = 0.001

    ! ... Executable ...
    
    ! set alias
    Hub => qmc%Hub

    ! scale parameters
    n     = nc*(1+na)
    Hub%n = n
    Hub%L = L
    
    ! Hubbard parameters
    Hub%n_t   = 1
    allocate(Hub%t(1))
    Hub%t(1)  = ONE
    Hub%n_mu  = 1
    allocate(Hub%mu(1))
    Hub%mu(1) = ZERO
    Hub%n_U   = 1
    allocate(Hub%U(1))
    Hub%U(1)  = U
    Hub%dtau  = dtau
    allocate(Hub%S%map(n))
    Hub%S%map= 1

    ! special case in the computing
    Hub%comp_dn = .true.
    Hub%copy_up = .false.

    ! simulation parameters
    Hub%nWarm    = nWarm
    Hub%nPass    = nPass
    Hub%nMeas    = nmeas
    if (U .gt. 6) then
       Hub%nTry  = 2
    else
       Hub%nTry  = 0
    end if
    Hub%gamma   = HALF
    Hub%nAccept = 1
    Hub%nReject = 0

    ! random seed
    Hub%idum = seed
    if (seed .eq. 0) then
       call date_and_time(VALUES=val)
       Hub%idum = val(8)*val(7)+val(6)**mod(val(5),5)
    end if
    Hub%seed = Hub%idum * (/1,2,3,4/)
    Hub%seed = mod(abs(Hub%seed), 4095)
    if (mod(Hub%seed(4),2) .eq. 0) then
       Hub%seed(4) = Hub%seed(4) + 1
    end if
    
    ! Initialize working space 
    call DQMC_WSpace_Allocate(n, 9, Hub%WS)

    ! HSF field
    if (initHSF) then
       allocate(Hub%HSF(nc,L))
       Hub%HSF = ONE
       do i = 1, Hub%L
          call ran0(nc, Hub%WS%R5(1:nc), Hub%seed)
          where(Hub%WS%R5(1:nc) .gt. HALF) Hub%HSF(:,i) = -ONE
       end do
    end if

    ! Initialize lookup table
    allocate(Hub%explook(-2:2,1))
    temp   = exp(dtau*U*HALF)    
    lambda = log(temp+sqrt(temp*temp-ONE))
    do i = -2, 2
       Hub%explook(i,1)=exp(i*lambda)
    end do

    ! Initialize V matrices
    !    The element of V(i) is either exp(nu) or exp(-nu)
    !    where nu = acosh(exp(U*dtau/2)). (see reference [1].) 
    !    The values of exp(nu) and exp(-nu) are stored in a lookup 
    !    table explook.  The decision of wheather V(i,j) is exp(nu) 
    !    or exp(-nu) is given by the list hub, which is a list 
    !    or random +1 and -1. Matrix V for spin up and down have
    !    opposite selection decision.
    ! 
    allocate(Hub%V_up(n,L))
    allocate(Hub%V_dn(n,L))
    Hub%V_up = ONE
    Hub%V_dn = ONE
    do i = 1, L
       do j = 1, nc
          Hub%V_up(j,i) = Hub%explook( int(Hub%HSF(j,i)), 1)
          Hub%V_dn(j,i) = Hub%explook(-int(Hub%HSF(j,i)), 1)
       end do
    end do

    ! Initialize Green functions
    Hub%B%n = n
    Hub%B%name = "Dense MatB"
    allocate(Hub%B%B(n,n))
    allocate(Hub%B%Bi(n,n))

    call DQMC_B_ExpInit(Hub%B, K, Hub%WS)
    
    call DQMC_SeqB_Init(n, L, nort, Hub%B, Hub%SB, Hub%WS)
    
    call DQMC_GFun_Init(n, L, Hub%G_up, Hub%V_up, Hub%WS, &
         nwrap, diff, diff, GMAT_UP)
    
    call DQMC_GFun_Init(n, L, Hub%G_dn, Hub%V_dn, Hub%WS, &
         nwrap, diff, diff, GMAT_DN)

    ! Initialize simulation range
    Hub%n_start = 1
    Hub%n_end   = nc

    Hub%S%checklist(STRUCT_INIT) = .false.

    ! Initialize Gtau
    call DQMC_Gtau_Init(n, L, TAU_T0, nOrt, nWrap, qmc%tau, Hub%B, Hub%WS)

    nullify(qmc%work)
    nullify(qmc%ipiv)

    Hub%P0%init = .false.

  end subroutine DQMC_QMC_Init

  ! =============================================================== !

  subroutine DQMC_QMC_Free(qmc)
    !
    ! Purpose
    ! =======
    !    This subrotine releases allocated memory.
    !
    ! Arguments
    ! =========
    !
    type(simqmc), intent(inout)   :: qmc           ! simulator
    target :: qmc

    ! ... Local ... 
    type(Hubbard), pointer        :: Hub

    ! ... Executable ...

    Hub => qmc%Hub
    deallocate(Hub%t, Hub%mu, Hub%U, Hub%S%map)    
    deallocate(Hub%explook, Hub%HSF, Hub%V_up, Hub%V_dn)
    call DQMC_B_Free(Hub%B)    
    call DQMC_Gfun_Free(Hub%G_up)
    call DQMC_Gfun_Free(Hub%G_dn)
    call DQMC_WSpace_Free(Hub%WS)
    call DQMC_SeqB_Free(Hub%SB)
    call DQMC_Gtau_Free(qmc%tau)

    if (associated(qmc%ipiv)) then
       deallocate(qmc%ipiv, qmc%work)
    end if

  end subroutine DQMC_QMC_Free

  ! =============================================================== !

  subroutine DQMC_measure(qmc, Aup, Adn, sgnup, sgndn)
    !
    ! Purpose
    ! =======
    !    This subrotine performs physical measurements.
    !    basically, it computes Gtau.
    !
    ! Arguments
    ! =========
    !
    type(simqmc), intent(inout)   :: qmc             ! qmc simulator
    real(wp), intent(inout)       :: Aup(:,:)        ! for working space
    real(wp), intent(inout)       :: Adn(:,:)
    real(wp), intent(out)         :: sgnup, sgndn
    target :: qmc

    ! ... Local variables ...
    integer  :: n, L, nl
    integer  :: info
    real(wp) :: lw(1)
    type(Hubbard), pointer :: Hub


    ! ... Executable ...
    
    Hub => qmc%Hub
    L   = Hub%L
    n   = Hub%n
    nl  = n*L

    ! allocate memory if necessary
    if (.not. associated(qmc%ipiv)) then
       allocate(qmc%ipiv(nL))
       call lapack_dgetri(nL, Aup, nL, qmc%ipiv, lw, -1, info)
       qmc%lwork = int(lw(1))
       allocate(qmc%work(qmc%lwork))
    end if
    
    ! construct Gup and Gdn
    sgnup = Meas_GetG(n, L, nL, Aup, Hub%B, Hub%V_up, &
         Hub%WS%R1, qmc%work, qmc%ipiv, qmc%lwork)
    sgndn = Meas_GetG(n, L, nL, Adn, Hub%B, Hub%V_dn, &
         Hub%WS%R1, qmc%work, qmc%ipiv, qmc%lwork)
    
    contains

      ! ===================================================== !

      function Meas_GetG(n, L, nL, M, B, V, W, work, ipiv, lwork) result(sgn)
        !
        ! Purpose
        ! =======
        !    This subrotine constructs matrix G,
        !        [  I            B1 ]
        !        [ -B2   I          ]
        !    M = [      -B3         ]
        !        [         ... ...  ]
        !        [          -BL  I  ] 
        !
        !    G = inv(M)
        !
        ! Arguments
        ! =========
        !
        integer, intent(in)     :: n, L, nL
        real(wp), intent(inout) :: M(:,:)       ! dim(M) = nL*nL
        type(matB), intent(in)  :: B            ! dim(B) = n*n
        real(wp), intent(in)    :: V(:,:)       ! dim(V) = n*L
        real(wp), intent(inout) :: W(:,:)       ! working space dim=n*n
        real(wp), intent(inout) :: work(:)      ! working space lwork
        integer, intent(inout)  :: ipiv(:)      ! dim(ipiv) = nL
        integer, intent(in)     :: lwork        ! dim of work
        real(wp) :: sgn

        ! ... Local variables ...
        integer :: i, j, info

        ! ... Executable ...
        
        ! construct M
        M = ZERO
        do i = 1, nL
           M(i,i) = ONE
        end do

        call DQMC_GetB(n, W, B, V(:,1), W)
        M(1:n, nL-n+1:nL) = W

        do i = 1, L-1
           call DQMC_GetB(n, W, B, V(:,i+1), W)
            M(i*n+1:(i+1)*n, (i-1)*n+1:i*n) = -W
        end do

        ! LU decomposition M
        call lapack_dgetrf(nL, nL, M, nL, ipiv, info)
        
        ! compute sign(det(M))
        sgn = ONE
        do i = 1, nL
           if (ipiv(i).ne. i) then
              sgn = -sgn
           end if
           
           if (M(i,i) .lt. ZERO) then
              sgn = -sgn
           end if
        end do
      
        ! invert M
        call lapack_dgetri(nL, M, nL, ipiv, work, lwork, info)
        
      end function Meas_GetG

      ! ===================================================== !

  end subroutine DQMC_measure
  

  ! =============================================================== !
  subroutine DQMC_Sim_Gtau(qmc, G_up, G_dn, ii, offset, sgnup, sgndn)
    !
    ! Purpose
    ! =======
    !    This subrotine compute Gup and Gdn for .
    !
    ! Pre-assumption
    ! ==============
    !    G_up and G_dn are NxN
    !    ii = 1, 2, ..., L
    !    offset = 0, 1, ..., L-1
    !
    ! Arguments
    ! =========
    !
    type(simqmc), intent(inout)   :: qmc           ! simulator
    real(wp), intent(inout)       :: G_up(:,:)
    real(wp), intent(inout)       :: G_dn(:,:)
    integer, intent(in)           :: ii, offset
    real(wp), intent(out)         :: sgnup, sgndn
    
    target  :: qmc
    pointer :: G_up
    pointer :: G_dn

    ! ... Executable ...

    ! get G(ii+offset, ii)

    call DQMC_MakeGtau(qmc%tau, qmc%Hub%G_up, qmc%Hub%G_dn, ii, offset)
    G_up => qmc%tau%upt0
    G_dn => qmc%tau%dnt0
    sgnup = qmc%Hub%G_up%sgn
    sgndn = qmc%Hub%G_dn%sgn
    
  end subroutine DQMC_Sim_Gtau

  ! =============================================================== !

  subroutine DQMC_Sim(qmc, G_avg, G_err)
    !
    ! Purpose
    ! =======
    !    This subrotine performs simulation and outputs results.
    !
    ! Pre-assumption
    ! ==============
    !    G_avg and G_err should be nc*nc*L
    !
    ! Arguments
    ! =========
    !
    type(simqmc), intent(inout)   :: qmc           ! simulator
    real(wp), intent(inout)       :: G_avg(:,:)
    real(wp), intent(inout)       :: G_err(:,:)
    target :: qmc

    ! ... Local variables ...
    integer :: i, j, k, nIter 
    type(Hubbard), pointer :: Hub 

    ! ... Executable ...

    Hub  => qmc%Hub

    ! Warmup sweep    
    do i = 1, Hub%nWarm
       ! The second parameter means no measurement should be made.
       call DQMC_Hub_Sweep(Hub, NO_MEAS0)
       call DQMC_Hub_Sweep2(Hub, Hub%nTry)
    end do
 
    ! We divide all the measurement into nBin,
    ! each having nPass/nBin pass.
    nIter  = Hub%nPass/Hub%nmeas
    do j = 1, nIter
       do k = 1, Hub%nmeas
          call DQMC_Hub_Sweep(Hub, NO_MEAS0)
          call DQMC_Hub_Sweep2(Hub, Hub%nTry)
       end do
    end do
    
    ! release memory
    call DQMC_QMC_Free(qmc)

  end subroutine DQMC_Sim
  
end module dqmc_simqmc
