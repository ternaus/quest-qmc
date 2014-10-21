program dqmc_verify

  use DQMC_2DPERL
  use DQMC_Hubbard
  use DQMC_Struct
  use DQMC_Phy0
  ! Chia-Chen: 09/06/2011
  ! MPI module is require in order to run verity in QUEST 1.0.8
  ! This is because :
  ! 1) We need to get the correct numebr of processors
  !    in DQMC_Phy0_Avg() even for serial runs.
  ! 2) The new DQMC_Phy0_Avg() can use multiprocessors
  !    to calculate observable averages and errors.
  use DQMC_MPI 
  implicit none

  !  Purpose
  ! =========
  ! This program verifies the correctness of the QUEST code 
  ! with two special cases on a 4x4 periodic lattice.
  !
  !     Case 1. One site (t=0)
  !     Case 2  No Coulomb interaction (U=0)
  !
  ! The computed values, like density and energy, are verified
  ! against theoretical values.
  !
  !  Parameters
  ! ============

  integer,      parameter  :: nx = 4, ny = 4, N = nx * ny
  real(wp),     parameter  :: t(4)  = [0.3_wp, 0.6_wp, ONE, ZERO]
  real(wp),     parameter  :: dtau  = 0.125_wp  
  real(wp),     parameter  :: U(6)  = [ONE, TWO, TWO * TWO, ZERO, -ONE, -TWO]
  real(wp),     parameter  :: mu(3) = [HALF, ZERO, -HALF]
  integer,      parameter  :: L = 12, HSF_IPT = -1, n_t = 1
  integer,      parameter  :: nWarm = 1000, nPass = 5000, nTry = 0
  integer,      parameter  :: nmeas = 12, nBin = 10, tausk = 10
  integer,      parameter  :: idum = 0, nOrth = 12, nWrap = 12
  integer,      parameter  :: ssxx = 0, fixw = 0
  real(wp),     parameter  :: errrate = 0.001_wp, difflim = 0.001_wp
  !real(wp),     parameter  :: pi = 3.141592653589793238462643_wp
  character(*), parameter  :: FMT_CMP    = &
       "(a20,f10.6,'  |',f10.6,' +-',f10.6,'  |  ',f6.2,' :',f6.2)"
  character(*), parameter  :: FMT_CONFIG = &
       "('Parameters : ',4(a,' =',f6.2,', '))"
  character(*), parameter  :: FMT_TITLE  = &
       "(20X,'Theoretical | Computed (avg +- error) |  |T-C| : error')"

  real(wp)                 :: t_up(4), t_dn(4)
  real(wp)                 :: mu_up(3), mu_dn(3)

  !
  !  Variables
  ! ===========
  
  type(Hubbard)      :: Hub 
  integer            :: i, j, k
  real(wp)           :: rho, energy_total, avg, err, one_site_occupancy
  real(wp)           :: tmp1, tmp2, tmp3, beta
  real(wp)           :: lambda(N), X(N)
  real               :: t1, t2
  character(len=30)  :: name
  integer            :: cnt(0:3) 

  ! QUEST now allows spin-dependent hopping and chemical potential.
  t_up = t
  t_dn = t
  mu_up = mu
  mu_dn = mu

  !
  !  Executable
  ! ============
  
  cnt = 0
  call cpu_time(t1)

  !Count the number of processors
  call DQMC_MPI_Init(qmc_sim,PLEVEL_1)
  
  ! Initialization
  call DQMC_Init_2DPerl(nx, ny, Hub%S, IMP_TRIANGLE)

  ! Case 1: t=0, run through mu = -0.5, 0.0, 0.5 and U = 1, 2, 4
  ! ==============================================================
  write(STDOUT, *)
  write(STDOUT,"(a,i3,a,f6.3)") "A 4x4 periodic lattice with L =", L, &
       ", dtau=", dtau
  write(STDOUT, *)
  write(STDOUT, *) "=============================="
  write(STDOUT, *) "| CASE 1. Single site (t=0)  |"
  write(STDOUT, *) "=============================="
  do i = 1, 3
     do j = 1, 6
        ! Initialize the parameter of the simulation
        call DQMC_Hub_Init(Hub, U(j:j), t_up(4:4), t_dn(4:4), mu_up(i:i), mu_dn(i:i), L, n_t, 1, 1, &
             dtau, HSF_IPT, nWarm, nPass, nMeas, nTry, nBin, tausk, idum, &
             nOrth, nWrap, fixw, errrate, difflim, HALF, 0, 0, ZERO, ZERO, ssxx, HSF_DISC)
        beta = L * dtau
     
        ! Execute
        call DQMC_Hub_Run(Hub,0)
     
        ! Check against theoretical results
        write(STDOUT, *)
        write(STDOUT, FMT_CONFIG) "t", ZERO, "mu", mu_up(i), "U", U(j), "beta", beta
        write(STDOUT, FMT_DBLINE)
        write(STDOUT, FMT_TITLE) 
        write(STDOUT, FMT_SGLINE) 
             
        ! 1. One-site density
        ! 
        !            2*exp((U/2+mu)*beta)+2*exp(2*mu*beta)
        !     rho = ---------------------------------------
        !            1+2*exp((U/2+mu)*beta)+exp(2*mu*beta)
        !
        tmp1 = exp((U(j) / TWO + mu(i)) * beta)
        tmp2 = exp(TWO * mu(i) * beta)
        tmp3 = ONE / (ONE + TWO * tmp1 + tmp2)
        rho = TWO * (tmp1 + tmp2) * tmp3
        call DQMC_Phy0_GetResult(Hub%P0, P0_DENSITY, name, avg, err)
        call Display("          Density : ", rho, avg, err)
     
        ! 2. One-site energy 
        !
        !               U*exp(2*mu*beta)                   
        !     E = ------------------------------------- - (mu + U/2)*rho
        !         1+2*exp((U/2+mu)*beta)+exp(2*mu*beta)
        !
        !    
        energy_total  = U(j) * tmp2 * tmp3 - (mu(i) + 0.5d0 * U(j)) * rho + 0.25d0 * U(j)

        call DQMC_Phy0_GetResult(Hub%P0, P0_ENERGY, name, avg, err)
        call Display("          Energy : ", energy_total, avg, err)
     
        ! 3. One-site occupancy
        !        
        !                     exp(2*mu*beta)
        !    PE = --------------------------------------
        !          1+2*exp((U/2+mu)*beta)+exp(2*mu*beta)
        !
        one_site_occupancy = tmp2 * tmp3 * Hub%U(1)
        call DQMC_Phy0_GetResult(Hub%P0, P0_NUD, name, avg, err)
        call Display(" Double occupancy : ", one_site_occupancy, avg, err)

        write(STDOUT, FMT_DBLINE) 
     end do
  end do

  ! Case 2: U=0, run through mu = -0.5, 0.0, 0.5 and t = 0.3, 0.6, 1.0
  ! ===================================================================
  ! 
  ! We compute some terms that will be used in verification.
  !
  write(STDOUT,*) 
  write(STDOUT,*) "=========================================="
  write(STDOUT,*) "| CASE 2. No Coulomb interaction (U=0)   |"
  write(STDOUT,*) "=========================================="
  write(STDOUT,*) 
  do i = 0, nx-1
     do j = 0, ny-1
        lambda(i*nx+j+1) = TWO*(cos(TWO*i*pi/nx)+cos(TWO*j*pi/ny))
     end do
  end do

  !
  ! In this case, we only need to run 1 measurement loop. Therefore, 
  ! we set nWarm = 0, nPass = 2, nMeas = 1, nBin = 1
  !     
  do i = 1, 3
     do j = 1, 3
        call DQMC_Hub_Init(Hub, U(4:4), t_up(j:j), t_dn(j:j), mu_up(i:i), mu_dn(i:i), L, n_t, 1, 1, &
             dtau, HSF_IPT, nWarm, nPass, nMeas, nTry, nBin, tausk, idum, &
             nOrth, nWrap, fixw, errrate, difflim, HALF, 0, 0, ZERO, ZERO, ssxx, HSF_DISC)
        beta = L*dtau
     
        call DQMC_Hub_Run(Hub,0)
     
        ! Check against theoretical results
        write(STDOUT, FMT_CONFIG) "t",t(j),"mu",mu(i),"U",ZERO,"beta",beta
        write(STDOUT, FMT_DBLINE)
        write(STDOUT, FMT_TITLE) 
        write(STDOUT, FMT_SGLINE) 
             
        ! 1. Density
        !                              1
        !    rho = 1.N * sum_k -----------------
        !                       1+exp(beta*x_k)
        !
        !          x_k = -t*lambda_k - mu
        !     lambda_k = 2(cos(theta_kx)+cos(theta_ky))
        !     theta_kx = 2kx*pi/N_x
        !     theta_ky = 2kx*pi/N_y
        !           kx = 0, 1, ... N_x-1
        !           ky = 0, 1, ... N_y-1
        !
        rho = ZERO
        do k = 1, N
           x(k) = ONE/(ONE+exp(beta*(-t(j)*lambda(k)-mu(i))))
           rho = rho + TWO*x(k)
        end do
        rho = rho / N
        call DQMC_Phy0_GetResult(Hub%P0, P0_DENSITY, name, avg, err)
        call Display("          Density : ", rho, avg, err)
        
     
        ! 2. One-site energy (not including chemical energy)
        !
        !                        -t*lambda_k
        !     E = 1/N * sum_k -------------------
        !                      exp(beta*x_k)+1
        !
        one_site_occupancy = -t(j) * TWO * dot_product(lambda, x) / N - mu(i) * rho
        call DQMC_Phy0_GetResult(Hub%P0, P0_ENERGY, name, avg, err)
        call Display("          Energy : ", one_site_occupancy, avg, err)
             
        write(STDOUT, FMT_DBLINE) 
     end do
  end do

  ! Statistics
  write(STDOUT,"(f6.2,a)") dble(cnt(1))/cnt(0)*100, &
       "% within 1 error bar (Expected 63.2%)."
  write(STDOUT,"(f6.2,a)") dble(cnt(1)+cnt(2))/cnt(0)*100, &
       "% within 2 error bar (Expected 86.5%)."

  call cpu_time(t2)
  write(STDOUT,*) "Running time:",  t2-t1, "(second)"

  ! clean up
  call DQMC_Hub_Free(hub)

contains 
  !===============================================================!
  !                    Supporting subroutines                     !
  !===============================================================!

  subroutine Display(name, theo, avg, err)
    !
    ! Format the print out for each case, and make some statistics
    !
    character(*), intent(in) :: name
    real(wp), intent(in)     :: theo, avg, err
    
    ! Local variabel
    integer  :: index
    real(wp) :: ratio
    
    ! Executable
    if (err /= ZERO) then
       ratio = abs(theo - avg) / err
       write(STDOUT, FMT_CMP) name, theo, avg, err, ratio, ONE
       index = ceiling(ratio)
       if (index > 3) index = 3
    else
       ratio = abs(theo - avg)
       write(STDOUT, FMT_CMP) name, theo, avg, err, ratio, ZERO
       if(ratio <= 1.0D-10) then
          index = 1
       else
          index = 3
       end if
    end if
    
    ! 
    cnt(0) = cnt(0) + 1
    cnt(index) = cnt(index) + 1

  end subroutine Display

  !===============================================================!
  !                    Supporting subroutines                     !
  !===============================================================!
  
  subroutine DQMC_Phy0_GetResult(P0, meas, name, avg, err)
    !
    ! Get results of measurements
    !
    type(Phy0), intent(in) :: P0
    integer, intent(in)    :: meas
    character(*), intent(inout) :: name
    real(wp), intent(inout):: avg, err

    name = P0_STR(meas)
    avg  = P0%meas(meas, P0%avg)
    err  = P0%meas(meas, P0%err)

  end subroutine DQMC_Phy0_GetResult

end program dqmc_verify
