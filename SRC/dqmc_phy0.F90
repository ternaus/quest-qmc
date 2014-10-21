module DQMC_Phy0
#include "dqmc_include.h"

  use DQMC_UTIL
  use DQMC_WSPACE
  use DQMC_STRUCT

  implicit none 
  
  !
  ! This module contains data structure and subroutines for some 
  ! physics measurement on Hubbard's model, including
  ! 
  !     1.  Up spin occupancy
  !     2.  Down spin occupancy
  !     3.  Potential energy
  !     4.  Kinetic energy
  !     5.  Total engery
  !     6.  Density
  !     7.  XX ferromagnetic structure factor
  !     8.  ZZ ferromagnetic structure factor
  !     9.  XX antiferromagnetic structure factor
  !    10.  ZZ antiferromagnetic structure factor
  !    11.  RMS XX AF SF
  !    12.  RMS ZZ AF SF
  !    13.  Average sign
  !    14.  Average up sign
  !    15.  Average down sign
  !    16.  Equal time Green's function
  !    17.  Density-density correlation fn: (up-up)
  !    18.  Density-density correlation fn: (up-dn)
  !    19.  XX Spin correlation function
  !    20.  ZZ Spin correlation function
  !
  ! Measurement 1-15 are scalars. 16-20 are lists of length
  ! nClass, which are the number of distinct autocorrelation terms.
  ! (see DQMC_STRUCT for more details about nClass.)
  !
  ! List of subroutines
  ! ===================
  !    1. DQMC_Phy0_Init(P0, nClass, nBin, nHist)
  !    2. DQMC_Phy0_Normalize(P0, u)
  !    3. DQMC_Phy0_Print(P0, S, OPT)
  !    4. DQMC_Phy0_GetErr(P0)
  !    5. DQMC_Phy0_Dump(P0, idx, opt)
  !    6. DQMC_PHY0_Hist(n, nBin, over, under, H, list, GetIndex)
  !    7. DQMC_Meas0(n, P0, G_up, G_dn, mu, t, sgnup, sgndn, S, up, dn)
  !
  !    *** Most subroutines are for internal use only. User program 
  !        should not need any of them.
  !
  ! Data Structure and Parameters
  ! =============================
  !    The data type Phy0 is consisted of two parts: the measurements
  !    and the histograms. 
  !    
  !    The measurements are put into bins and then analysized in the end.
  !    To make data manipulation easy, all scalar variables are put into
  !    an array S and identified by the indeces, which are defined by the
  !    parameters below. Each column of S represents a bin.
  !    List variables, like Green's function, are put into separated 
  !    arrays, in which each column is also for one bin. 
  !        
  !    There are two special bins, averge and error, which are used to
  !    store the final averaged result and error. Their indeces are 
  !    specified by 'avg' and 'err' respectively.
  !
  !    The histogram part consists of three histograms, each having
  !    nHist+2 bins. They are histogram for up occupancy (Nup), 
  !    down occupancy (Ndn), and Nup*Ndn. The additional two bins
  !    are for data exceeds the range, whose indeces are specified
  !    by 'over' and 'under' respectively.
  !

  ! Array    
  integer, parameter  :: narrays = 9

  ! Index of the array varaiables
  integer, parameter  :: IMEAS = 0
  integer, parameter  :: IGFUN = 1
  integer, parameter  :: IGFUP = 2
  integer, parameter  :: IGFDN = 3
  integer, parameter  :: ISPXX = 4
  integer, parameter  :: ISPZZ = 5
  integer, parameter  :: IAVSP = 6
  integer, parameter  :: IDEN0 = 7
  integer, parameter  :: IDEN1 = 8
  integer, parameter  :: IPAIR = 9

  ! Parameter for the index of scalar variables (IMEAS)
  integer, parameter :: P0_NUP       = 1
  integer, parameter :: P0_NDN       = 2
  integer, parameter :: P0_NUD       = 3
  integer, parameter :: P0_KE        = 4
  integer, parameter :: P0_ENERGY    = 5
  integer, parameter :: P0_DENSITY   = 6
  integer, parameter :: P0_CHIT      = 7
  integer, parameter :: P0_CV        = 8

  integer, parameter :: P0_SFERRO    = 9
  integer, parameter :: P0_SFER2     = 10
  integer, parameter :: P0_SAF       = 15
  integer, parameter :: P0_SAFSQ     = 16
  integer, parameter :: P0_SAF2      = 17
  integer, parameter :: P0_SAF2SQ    = 18

  integer, parameter :: P0_potential_energy    = 11
  integer, parameter :: P0_hopping_energy    = 12
  integer, parameter :: P0_double_occupancy    = 13
  integer, parameter :: P0_magnetization_squared    = 14

  integer, parameter :: P0_N_NO_SAF  = 14
  integer, parameter :: P0_N         = 18

  integer, parameter :: P0_SGN       = 1
  integer, parameter :: P0_SGNUP     = 2
  integer, parameter :: P0_SGNDN     = 3


  ! Name of scalar variables
  character(len=*), parameter :: P0_STR(P0_N) = (/&
       "          Up spin occupancy : ", &
       "        Down spin occupancy : ", &
       "             <U*N_up*N_dn>  : ", &
       "             Kinetic energy : ", &
       "               Total energy : ", &
       "                    Density : ", &
       "                Chi_thermal : ", &
       "              Specific heat : ", &
       "  XX Ferro structure factor : ", &
       "  ZZ Ferro structure factor : ", &
       "           Potential energy : ", &
       "             Hopping energy : ", &
       "           Double occupancy : ", &
       "      Magnetization squared : ", &
       "     XX AF structure factor : ", &
       "  Root Mean Square of XX AF : ", &
       "     ZZ AF structure factor : ", &
       "  Root Mean Square of ZZ AF : "/)

  character(len=*), parameter :: P0_SIGN_STR(3) = (/&
       "                   Avg sign : ", &
       "                Avg up sign : ", &
       "                Avg dn sign : "/)


  type Phy0
     ! Measurement part
     integer  :: nClass                     ! Number of distinct 
                                            ! autocorrelction terms
     integer  :: nBin                       ! Number of terms
     integer  :: nMeas                      ! Number of measurements
     integer  :: avg, err                   ! Index for average and error
                                            ! bins.
     integer  :: cnt                        ! Number of measurement for 
                                            ! current bin
     integer  :: idx                        ! current bin index

     integer   :: n                         ! number of sites
     real(wp)  :: beta                      ! Inverse Temperature
    
     ! Scalar array
     real(wp), pointer :: meas(:, :)        ! Scalar varaibles
     real(wp), pointer :: sign(:, :)        ! Scalar varaibles

     ! Indices
     integer :: IARR(0: narrays + 1) 
     integer :: IARRFT(1: narrays + 1)
     integer :: IARREV(1: narrays + 1)

     real(wp), pointer   :: AllProp(:, :)             ! Vector of all properties
     real(wp), pointer   :: AllPropFT(:, :)           ! Matrix of FT 
     real(wp), pointer   :: AllPropEigVal(:, :)       ! Vector of Fourier transforms
     complex*16, pointer :: AllPropEigVec(:, :, :, :)   ! Vector with the normal modes

     !Pointers to AllProp
     real(wp), pointer :: G_fun(:, :)    ! Green's function
     real(wp), pointer :: Gf_up(:, :)    ! Green's function
     real(wp), pointer :: Gf_dn(:, :)    ! Green's function
     real(wp), pointer :: SpinXX(:, :)   ! XX Spin correlation function
     real(wp), pointer :: SpinZZ(:, :)   ! ZZ Spin correlation function
     real(wp), pointer :: AveSpin(:, :)  ! Ave Spin correlation function
     real(wp), pointer :: Den0(:, :)     ! Density-density correlation 
     real(wp), pointer :: Den1(:, :)     ! up-up (0) and up-dn (1) 
     real(wp), pointer :: Pair(:, :)     ! on-site pairing

     ! working space
     real(wp), pointer :: up(:) 
     real(wp), pointer :: dn(:) 
     logical :: compSAF
     logical :: init
     logical :: initFT

  end type Phy0

contains

  ! Subroutines
  ! ==================================================================
  
  subroutine DQMC_Phy0_Init(P0, S, beta, nBin, WS)
    !
    ! Purpose
    ! =======
    !    This subroutine initializes Phy0.
    !
    !  Pre-assumption
    ! ==============
    !    nClass, nBin and nHist are positive integers.
    !
    ! Arguments
    ! =========
    !
    type(Phy0), intent(inout) :: P0      ! Phy0 to be initialized
    type(Struct), intent(in)  :: S
    integer, intent(in)       :: nBin    ! No of bins
    real(wp), intent(in)      :: beta
    type(WSpace), intent(in), target :: WS

    ! ... Local vars ...
    integer :: i, n

    ! ... Executable ...

    P0%nClass  = S%nClass
    P0%nBin    = nBin
    P0%n       = S%nsite
    P0%beta    = beta

    P0%avg     = nBin + 1
    P0%err     = nBin + 2
    P0%cnt     = 0
    P0%idx     = 1

    P0%compSAF = S%checklist(STRUCT_PHASE)
    ! count total number of data
    if (P0%compSAF) then
       P0%nMeas = P0_N
    else
       p0%nMeas = P0_N_NO_SAF
    end if

    ! Allocate storages for sign and properties
    n = P0%nmeas + P0%nClass * narrays
    allocate(P0%sign(3, P0%err))  
    allocate(P0%AllProp(n, P0%err))

    !Pointers to beginning of each array
    P0%IARR(IMEAS) = 1
    do i = 1, narrays + 1
       P0%IARR(i) = P0%nmeas + 1 + (i - 1) * P0%nClass 
    enddo
 
    P0%meas    => P0%AllProp(P0%IARR(IMEAS):P0%IARR(IMEAS + 1) - 1, :)
    P0%G_fun   => P0%AllProp(P0%IARR(IGFUN):P0%IARR(IGFUN + 1) - 1, :)
    P0%Gf_up   => P0%AllProp(P0%IARR(IGFUP):P0%IARR(IGFUP + 1) - 1, :)
    P0%Gf_dn   => P0%AllProp(P0%IARR(IGFDN):P0%IARR(IGFDN + 1) - 1, :)
    P0%SpinXX  => P0%AllProp(P0%IARR(ISPXX):P0%IARR(ISPXX + 1) - 1, :)
    P0%SpinZZ  => P0%AllProp(P0%IARR(ISPZZ):P0%IARR(ISPZZ + 1) - 1, :)
    P0%AveSpin => P0%AllProp(P0%IARR(IAVSP):P0%IARR(IAVSP + 1) - 1, :)
    P0%Den0    => P0%AllProp(P0%IARR(IDEN0):P0%IARR(IDEN0 + 1) - 1, :)
    P0%Den1    => P0%AllProp(P0%IARR(IDEN1):P0%IARR(IDEN1 + 1) - 1, :)
    P0%Pair    => P0%AllProp(P0%IARR(IPAIR):P0%IARR(IPAIR + 1) - 1, :)

    !allocate(P0%meas(P0%nMeas, nBin+2))
    !allocate(P0%sign(3, nBin+2))
    !allocate(P0%G_fun  (nClass, nBin+2))
    !allocate(P0%Gf_up  (nClass, nBin+2))
    !allocate(P0%Gf_dn  (nClass, nBin+2))
    !allocate(P0%SpinXX (nClass, nBin+2))
    !allocate(P0%SpinZZ (nClass, nBin+2))
    !allocate(P0%AveSpin(nClass, nBin+2))
    !allocate(P0%Den0   (nClass, nBin+2))
    !allocate(P0%Den1   (nClass, nBin+2))
    !allocate(P0%Pair   (nClass, nBin+2))


    ! Initialize 
    P0%meas    = ZERO
    P0%sign    = ZERO
    P0%G_fun   = ZERO
    P0%Gf_up   = ZERO
    P0%Gf_dn   = ZERO
    P0%SpinXX  = ZERO
    P0%SpinZZ  = ZERO
    P0%AveSpin = ZERO
    P0%Den0    = ZERO
    P0%Den1    = ZERO
    P0%Pair    = ZERO
    
    P0%up => WS%R5
    P0%dn => WS%R6

    P0%init = .true.

    ! 10/26/2012
    ! The following pointers will be allocated in DQMC_GetFT():
    !
    !     AllPropFT(:,:)           ! Matrix of FT 
    !     AllPropEigVal(:,:)       ! Vector of Fourier transforms
    !     AllPropEigVec(:,:,:,:)   ! Vector with the normal modes
    !
    ! However, dqmc_2dperl.F90 does not calculate momentum space observables, the status of the
    ! three pointer arrays become associated, but never allocated. This is why in /EXAMPLE/test,
    ! test.F90 generates the error message at the end of the simulation:
    !
    !    pointer being freed was not allocated
    !
    ! This happens only in /EXAMPLE/test/test.F90, not in /EXAMPLE/geom/ggeom.F90. To fix this,
    ! I added a flag that checks whether these arrays are allocated. The default value is .false.
    ! It will be changed to .true. once DQMC_GetFT() is called.
    P0%initFT = .false.

   end subroutine DQMC_Phy0_Init

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy0_Free(P0)
    !
    ! Purpose
    ! =======
    !    This subroutine frees Phy0.
    !
    ! Arguments
    ! =========
    !
    type(Phy0), intent(inout) :: P0      ! Phy0 to be freed

    ! ... Executable ...
    if (P0%init) then

       nullify(P0%meas)
       nullify(P0%G_fun, P0%Gf_up, P0%Gf_dn)
       nullify(P0%SpinXX, P0%SpinZZ, P0%AveSpin)
       nullify(P0%Den0, P0%Den1)
       nullify(P0%Pair)

       deallocate(P0%AllProp, P0%sign)

       ! P0%initFT = .true. if DQMC_GetFT() was called
       if (P0%initFT) then
          deallocate(P0%AllPropEigVal)
          deallocate(P0%AllPropEigVec)
       endif
    end if

  end subroutine DQMC_Phy0_Free

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy0_Avg(P0)
    !
    ! Purpose
    ! =======
    !    This subroutine averges the data with a bin.
    !    It also computes some measurements from others, such as
    !    density and energy.
    !
    !  Pre-assumption
    ! ==============
    !    idx is in the range of (1,nClass).
    !
    ! Arguments
    ! =========
    !
    type(Phy0), intent(inout) :: P0     ! Phy0
    
    ! ... local scalar ...
    real(wp) :: factor
    integer  :: idx, n

    ! ... Executable ...
    idx = P0%idx

    ! compute the normalization factor = 1/cnt
    if (P0%cnt == 0) then
       call DQMC_Error("Phy0 normalize: cnt = 0", 0)
    end if
    factor = ONE / P0%cnt

    ! average
    P0%meas(:, idx) = P0%meas(:, idx) * factor
    P0%sign(:, idx) = P0%sign(:, idx) * factor

    if (P0%compSAF) then
       ! The sqaure terms
       P0%meas(P0_SAFSQ, idx) = sqrt(abs(P0%meas(P0_SAFSQ, idx)))
       P0%meas(P0_SAF2SQ, idx) = sqrt(abs(P0%meas(P0_SAF2SQ, idx)))
    end if

    ! This list terms
    n = P0%nClass
    call blas_dscal(n, factor, P0%G_fun (1:n, idx), 1)
    call blas_dscal(n, factor, P0%Gf_up (1:n, idx), 1)
    call blas_dscal(n, factor, P0%Gf_dn (1:n, idx), 1)
    call blas_dscal(n, factor, P0%SpinXX(1:n, idx), 1)
    call blas_dscal(n, factor, P0%SpinZZ(1:n, idx), 1)
    call blas_dscal(n, factor, P0%Den0  (1:n, idx), 1)
    call blas_dscal(n, factor, P0%Den1  (1:n, idx), 1)
    call blas_dscal(n, factor, P0%Pair  (1:n, idx), 1)

    ! Change bin
    P0%idx = P0%idx + 1

    ! reset the counter
    p0%cnt = 0

  end subroutine DQMC_Phy0_Avg

  !--------------------------------------------------------------------!
  
  subroutine DQMC_Phy0_Print(P0, S, OPT)
    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !    This subroutine prints out the average and errors
    !    of measurements. Structure S will give labels for
    !    each autocorrelation terms.
    !
    !  Pre-assumption
    ! ===============
    !    OPT is a file handle
    !    DQMC_Phy0_GetErr was called.
    !
    ! Arguments
    ! =========
    !
    type(Phy0), intent(in)    :: P0   ! Phy0
    type(Struct), intent(in)  :: S    ! Underline lattice structure
    integer, intent(in)       :: OPT  ! Output file handle

    ! ... Local scalar ...
    integer :: nClass, avg, err

    ! ... Executable ...

    if (qmc_sim%rank /= 0) return

    nClass = P0%nClass
    avg    = P0%avg
    err    = P0%err

    ! Scalar terms
    call DQMC_Print_RealArray(0, 3, "Sign of equal time measurements:", &
         P0_SIGN_STR, P0%sign(:,avg:avg), P0%sign(:,err:err), OPT)
    
    call DQMC_Print_RealArray(0, P0%nmeas, "Equal time measurements:", &
         P0_STR, P0%meas(:,avg:avg), P0%meas(:,err:err), OPT)

    ! Function terms
    call DQMC_Print_RealArray(0, nClass, "Mean Equal time Green's function:", &
         S%clabel, P0%G_fun(:, avg:avg), P0%G_fun(:, err:err), OPT)
    
    call DQMC_Print_RealArray(0, nClass, "Up Equal time Green's function:", &
         S%clabel, P0%Gf_up(:, avg:avg), P0%Gf_up(:, err:err), OPT)
    
    call DQMC_Print_RealArray(0, nClass, "Down Equal time Green's function:", &
         S%clabel, P0%Gf_dn(:, avg:avg), P0%Gf_dn(:, err:err), OPT)
    
    call DQMC_Print_RealArray(0, nClass, &
         "Density-density correlation fn: (up-up)", &
         S%clabel, P0%Den0(:, avg:avg), P0%Den0(:, err:err), OPT)
    
    call DQMC_Print_RealArray(0, nClass, &
         "Density-density correlation fn: (up-dn)", &
         S%clabel, P0%Den1(:, avg:avg), P0%Den1(:, err:err), OPT)
    
    call DQMC_Print_RealArray(0, nClass, "XX Spin correlation function:", &
         S%clabel, P0%SpinXX(:, avg:avg), P0%SpinXX(:, err:err), OPT)
    
    call DQMC_Print_RealArray(0, nClass, "ZZ Spin correlation function:", &
         S%clabel, P0%SpinZZ(:, avg:avg), P0%SpinZZ(:, err:err), OPT)

    call DQMC_Print_RealArray(0, nClass, "Average Spin correlation function:", &
         S%clabel, P0%AveSpin(:, avg:avg), P0%AveSpin(:, err:err), OPT)
    
    call DQMC_Print_RealArray(0, nClass, "Pairing correlation function:", &
         S%clabel, P0%Pair(:, avg:avg), P0%Pair(:, err:err), OPT)
    
  end subroutine DQMC_Phy0_Print

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy0_GetErr(P0)
    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !    This subroutine computes the average and errors
    !    of measurements. 
    !
    ! Argument
    ! ========
    !
    type(Phy0), intent(inout) :: P0

    ! ... Local Scalar ...
    integer  :: i, n
    integer  :: avg, err, mpi_err
    integer  :: nproc

    ! ... Local Array
    real(wp) :: sum_sgn, sgn(P0%nBin), y(P0%nBin), data(P0%nBin)
    
    ! ... Executable ...
    n   = P0%nBin
    avg = P0%avg
    err = P0%err
    nproc = qmc_sim%size

    if (nproc == 1) then

       ! Average sign, sign up, sign down
       do i = P0_SGNUP, P0_SGNDN
          data = P0%sign(i, 1:n)
          call DQMC_JackKnife(n, P0%sign(i, avg), P0%sign(i, err), data, &
               y, sgn, sum_sgn)
       end do
       
       data = P0%sign(P0_SGN, 1:n)
       call DQMC_JackKnife(n, P0%sign(P0_SGN, avg), P0%sign(P0_SGN, err), data, &
            y, sgn, sum_sgn)
       
       ! Average other single terms
       do i = 1, P0%nMeas
          data = P0%meas(i, 1:n)
          call DQMC_SignJackKnife(n, P0%meas(i, avg), P0%meas(i, err), data, &
               y, sgn, sum_sgn)
       end do
       
       ! Average Green function
       do i = 1, P0%nClass
          data =  P0%G_fun(i, 1:n)
          call DQMC_SignJackKnife(n, P0%G_fun(i, avg), P0%G_fun(i, err), &
               data, y, sgn, sum_sgn)
       end do
       
       ! Up Green function
       do i = 1, P0%nClass
          data =  P0%Gf_up(i, 1:n)
          call DQMC_SignJackKnife(n, P0%Gf_up(i, avg), P0%Gf_up(i, err), &
               data, y, sgn, sum_sgn)
       end do
       
       ! Down Green function
       do i = 1, P0%nClass
          data =  P0%Gf_dn(i, 1:n)
          call DQMC_SignJackKnife(n, P0%Gf_dn(i, avg), P0%Gf_dn(i, err), &
               data, y, sgn, sum_sgn)
       end do
       
       ! Average correlated Density 
       do i = 1, P0%nClass
          data = P0%Den0(i, 1:n)
          call DQMC_SignJackKnife(n, P0%Den0(i, avg), P0%Den0(i, err), &
               data, y, sgn, sum_sgn)
       end do
       
       do i = 1, P0%nClass
          data = P0%Den1(i, 1:n)
          call DQMC_SignJackKnife(n, P0%Den1(i, avg), P0%Den1(i, err), &
               data, y, sgn, sum_sgn)
       end do
       
       ! Average spin
       do i = 1, P0%nClass
          data = P0%SpinXX(i, 1:n)
          call DQMC_SignJackKnife(n, P0%SpinXX(i, avg), P0%SpinXX(i, err), &
               data, y, sgn, sum_sgn)
       end do
       
       do i = 1, P0%nClass
          data = P0%SpinZZ(i, 1:n)
          call DQMC_SignJackKnife(n, P0%SpinZZ(i, avg), P0%SpinZZ(i, err), &
               data, y, sgn, sum_sgn)
       end do

       do i = 1, P0%nClass
          P0%AveSpin(i, 1:n) = ( P0%SpinZZ(i, 1:n) + 2.d0*P0%SpinXX(i, 1:n) ) / 3.d0
          data = P0%AveSpin(i, 1:n)
          call DQMC_SignJackKnife(n, P0%AveSpin(i, avg), P0%AveSpin(i, err), &
               data, y, sgn, sum_sgn)
       end do

       ! Pairing
       do i = 1, P0%nClass
          data = P0%Pair(i, 1:n)
          call DQMC_SignJackKnife(n, P0%Pair(i, avg), P0%Pair(i, err), &
               data, y, sgn, sum_sgn)
       end do

       ! Store Jackknife in bins
       do i = 1, n
          P0%sign(:,i) = (n*P0%sign(:,avg) - P0%sign(:,i)) / dble(n-1)
          P0%AllProp(:,i) = (sum_sgn*P0%AllProp(:,avg) - P0%AllProp(:,i)) / dble(n-1)
          P0%AllProp(:,i) = P0%AllProp(:,i) / P0%sign(P0_SGN,i)
       enddo

       ! Deal with error and expectation values of cv and chi_thermal properly
       P0%meas(P0_CHIT,1:avg) = P0%meas(P0_CHIT,1:avg) - P0%n * P0%beta**2 * P0%meas(P0_ENERGY,1:avg) &
          * P0%meas(P0_DENSITY,1:avg)
       P0%meas(P0_CV, 1:avg) = P0%meas(P0_CV, 1:avg) - P0%n * (P0%beta * P0%meas(P0_ENERGY, 1:avg))**2       

       P0%meas(P0_CV, err) = sqrt((n-1) * sum((P0%meas(P0_CV, 1:n) - P0%meas(P0_CV, avg))**2))
       P0%meas(P0_CHIT, err) = sqrt((n-1) * sum((P0%meas(P0_CHIT, 1:n) - P0%meas(P0_CHIT, avg))**2))

    else

       mpi_err = 0

#      ifdef _QMC_MPI
          
          n = size(P0%AllProp,1)

          !Set bin to 1 and fill the average spin
          P0%AveSpin(:,1) = ( P0%SpinZZ(:,1) + 2.d0*P0%SpinXX(:,1) ) / 3.d0

          !Average sign
          call mpi_allreduce(P0%sign(:,1), P0%sign(:,avg), 3, mpi_double, &
             mpi_sum, mpi_comm_world, mpi_err)

          !Average properties
          call mpi_allreduce(P0%AllProp(:,1), P0%AllProp(:,avg), n, mpi_double, &
             mpi_sum, mpi_comm_world, mpi_err)

          !Compute average over n-1 processors
          P0%AllProp(:,1) = (P0%AllProp(:,avg) - P0%AllProp(:,1)) / dble(nproc - 1)
          P0%sign(:,1)    = (P0%sign(:,avg) - P0%sign(:,1)) / dble(nproc - 1)

          !Store average amongst all processors
          P0%AllProp(:,avg) = P0%AllProp(:,avg) / P0%sign(P0_SGN,avg) 
          P0%sign(:,avg)    = P0%sign(:,avg) / dble(nproc)

          !Store jackknife in the processor bin
          P0%AllProp(:,1)   = P0%AllProp(:,1) / P0%sign(P0_SGN,1) 

          !Compute proper chi_thermal and C_v
          P0%meas(P0_CHIT,1:avg) = P0%meas(P0_CHIT,1:avg) - P0%n * P0%beta**2 * P0%meas(P0_ENERGY,1:avg) &
             * P0%meas(P0_DENSITY,1:avg)
          P0%meas(P0_CV,1:avg) = P0%meas(P0_CV,1:avg) - P0%n * (P0%beta * P0%meas(P0_ENERGY,1:avg))**2

          !Compute error
          call mpi_allreduce(P0%AllProp(:,1)**2, P0%AllProp(:,err), n, mpi_double, &
             mpi_sum, mpi_comm_world, mpi_err)
          P0%AllProp(:,err) = P0%AllProp(:,err) / dble(nproc) - P0%AllProp(:,avg)**2 
          P0%AllProp(:,err) = sqrt(P0%AllProp(:,err) * dble(nproc-1))

          !Compute error for sign
          call mpi_allreduce(P0%sign(:,1)**2, P0%sign(:,err), 3, mpi_double, &
             mpi_sum, mpi_comm_world, mpi_err)
          P0%sign(:,err) = P0%sign(:,err) / dble(nproc) - P0%sign(:,avg)**2 
          P0%sign(:,err) = sqrt(P0%sign(:,err) * dble(nproc-1))

          P0%meas(P0_CHIT,avg) = P0%meas(P0_CHIT,avg) - P0%n * P0%beta**2 * P0%meas(P0_ENERGY,avg) &
             * P0%meas(P0_DENSITY,avg)
      
          P0%meas(P0_CV,avg) = P0%meas(P0_CV,avg) - P0%n * (P0%beta * P0%meas(P0_ENERGY,avg))**2

#      endif

    endif
 
  end subroutine DQMC_Phy0_GetErr

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy0_Meas(n, P0, G_up, G_dn, U, mu_up, mu_dn, t_up, t_dn, sgnup, sgndn, S)
    !
    ! Purpose
    ! =======
    !    This subroutine performs some physics measurement on
    !    Hubbard model.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)          :: n            ! Number of sites
    type(Phy0), intent(inout)    :: P0           ! Phy0
    real(wp), intent(in)         :: G_up(n,n)    ! Green's function
    real(wp), intent(in)         :: G_dn(n,n)    ! for spin up and down
    real(wp), intent(in)         :: sgnup, sgndn ! Sgn for det(G_up) det(G_dn)
    real(wp), intent(in)         :: mu_up(n), mu_dn(n)  ! Chemical and Kinetic para
    real(wp), intent(in)         :: t_up(:), t_dn(:)  ! Chemical and Kinetic para
    real(wp), intent(in)         :: U(:)         ! Chemical and Kinetic para
    type(Struct), intent(in)     :: S            ! Underline structure
    target :: S

    ! ... local scalar ...

    integer  :: i, j, k, ph                      ! Loop iterator
    integer  :: tmp, idx, m                      ! Helper variables
    real(wp) :: sgn                        
    real(wp) :: var1, var2, var3          
    integer, pointer  :: start(:) 
    integer, pointer  :: r(:) 
    integer, pointer  :: A(:) 

    ! Auxiliary variable for chi_thermal and C_v
    real(wp) :: Cbar, Nbar, Tbar, un
    real(wp) :: h_up(n, n), h_dn(n, n) 
    real(wp) :: A_up(n, n), A_dn(n, n)

    ! ... executable ...

    idx = P0%idx
    tmp = P0%avg

    ! initialization
    ! Here we use avg bin as a temp variable 
    P0%meas(:,tmp)   = ZERO

    P0%G_fun(:,tmp)  = ZERO
    P0%Gf_up(:,tmp)  = ZERO
    P0%Gf_dn(:,tmp)  = ZERO
    P0%Den0(:,tmp)   = ZERO
    P0%Den1(:,tmp)   = ZERO
    P0%SpinXX(:,tmp) = ZERO
    P0%SpinZZ(:,tmp) = ZERO
    P0%Pair(:,tmp)   = ZERO

    ! Compute the site density for spin up and spin down
    do i = 1, n
       !======================================================!
       ! The density of electrons of spin up(dn) on site i    !
       ! is 1-G_up(i,i) (1-G_dn(i,i)).                        !
       ! nup (ndn) is the sum of all spin up (down) electrons.!
       !======================================================!
       P0%up(i)  = ONE - G_up(i, i)
       P0%dn(i)  = ONE - G_dn(i, i)
       P0%meas(P0_NUD, tmp) = P0%meas(P0_NUD, tmp)+ &
            P0%up(i) * P0%dn(i) * U(S%Map(i))
       !======================================================!
       ! Double occupancy P0%up(i) * P0%dn(i)
       !======================================================!
       P0%meas(P0_double_occupancy, tmp) = P0%meas(P0_double_occupancy, tmp) +&
       P0%up(i) * P0%dn(i)
       !=====================================================================!
       ! Potential energy (P0%up(i)-0.5d0) * (P0%dn(i)-0.5d0) * U(S%Map(i))
       !=====================================================================!
       P0%meas(P0_potential_energy, tmp) = P0%meas(P0_potential_energy, tmp)+ &
            (P0%up(i) - 0.5d0) * (P0%dn(i) - 0.5d0) * U(S%Map(i))
       
    end do

    P0%meas(P0_NUP, tmp) = sum(P0%up)
    P0%meas(P0_NDN, tmp) = sum(P0%dn)
    
    !=================================================================!
    ! Kinetic energy = tt*sum_{ij\sigma}(G_{ij\sigma}+G_{ji\sigma}) - 
    ! - \sum_{i\sigma} \mu_{i\sigma} (n_{i\sigma} + U_i / 2)!
    ! where site i and site j are neighbors   !
    !=================================================================!
    ! Hopping energy = tt*sum_{ij\sigma}(G_{ij\sigma}+G_{ji\sigma})

    ! set alias
    start => S%T%cstart
    r     => S%T%row
    A     => S%T%A
    
    ! loop all adj sites
    do i = 1, n  ! for each column
       do j = start(i), start(i + 1)-1 ! for each nonzero elements
          P0%meas(P0_KE, tmp) =  P0%meas(P0_KE, tmp) + &
               t_up(A(j)) * G_up(r(j), i)               + &
               t_dn(A(j)) * G_dn(r(j), i)
          P0%meas(P0_hopping_energy, tmp) = P0%meas(P0_hopping_energy, tmp) + &
               t_up(A(j)) * G_up(r(j), i)               + &
               t_dn(A(j)) * G_dn(r(j), i)
       end do
       P0%meas(P0_KE, tmp)  = P0%meas(P0_KE, tmp)            - &
            (mu_up(S%Map(i)) + 0.5d0 * U(S%map(i))) * P0%up(i) - &
            (mu_dn(S%Map(i)) + 0.5d0 * U(S%map(i))) * P0%dn(i)            
    end do

    !=================================================================!
    ! Total occupancy = nup + ndn
    !=================================================================!
    P0%meas(P0_DENSITY, tmp) = P0%meas(P0_NUP, tmp) + &                                                                  
         P0%meas(P0_NDN, tmp)      

    !=================================================================! 
    ! Magnetisation squared = 1/4 (rho - 2 double_occupancy)    
    !=================================================================! 
    P0%meas(P0_magnetization_squared, tmp) = 0.25d0 * (P0%meas(P0_density, tmp) -&
     2 * P0%meas(P0_double_occupancy, tmp))

    !=================================================================!
    ! Total energy = hopping energy + potential energy - 
    ! - sum_{i \sigma} (\mu_{up i \sigma} n_{dn i \sigma})
    !=================================================================!
    P0%meas(P0_ENERGY, tmp) = P0%meas(P0_hopping_energy, tmp) + &
         P0%meas(P0_potential_energy, tmp)
    do i = 1, n
      P0%meas(P0_ENERGY, tmp) = P0%meas(P0_ENERGY, tmp) - &
       (mu_up(S%Map(i)) * P0%up(i) + mu_dn(S%Map(i)) * P0%dn(i))
    enddo
    !=========================================!
    ! Chi_thermal 
    !=========================================!

    ! Fill h_up, h_dn with hopping matrix elements
    h_up = ZERO
    h_dn = ZERO
    do i = 1, n
       do j = start(i), start(i + 1) - 1
          h_up(r(j), i) =  -t_up(A(j))
          h_dn(r(j), i) =  -t_dn(A(j))
       end do
       h_up(i,i) =  h_up(i,i) - mu_up(S%Map(i)) - 0.5d0 * U(S%Map(i))
       h_dn(i,i) =  h_dn(i,i) - mu_dn(S%Map(i)) - 0.5d0 * U(S%Map(i))
    end do

    ! Gfun * t
    ! Fill h_up, h_dn 
    call blas_dgemm('N', 'N', n, n, n, ONE, G_up, n, h_up, n, ZERO, A_up, n)
    call blas_dgemm('N', 'N', n, n, n, ONE, G_dn, n, h_dn, n, ZERO, A_dn, n)

    ! Total number of particles
    Nbar = sum(P0%up) + sum(P0%dn)

    Tbar = 0.d0
    do i = 1, n
       Tbar = Tbar + h_up(i, i) + h_dn(i, i)
    enddo

    Cbar = 0.d0
    do i = 1, n
       Cbar = Cbar + A_up(i, i) + A_dn(i, i)
    enddo

    !< N T >
    P0%meas(P0_CHIT, tmp) = (Tbar - Cbar) * Nbar 
    do j = 1, n
       do k = 1, n
          P0%meas(P0_CHIT, tmp) = P0%meas(P0_CHIT, tmp) - G_up(j,k) * A_up(k,j) - G_dn(j,k) * A_dn(k,j)
       enddo
    enddo
    P0%meas(P0_CHIT,tmp) = P0%meas(P0_CHIT, tmp) + Cbar

    !< N U >
    P0%meas(P0_CHIT,tmp) = P0%meas(P0_CHIT, tmp) + Nbar * P0%meas(P0_NUD, tmp)
    do i = 1, n
       un = ONE
       do k = 1, n 
          un = un - G_up(i, k) * G_up(k, i)
       enddo
       P0%meas(P0_CHIT, tmp) = P0%meas(P0_CHIT, tmp) +  un*P0%dn(i)*U(S%Map(i))
       un = ONE
       do k = 1, n 
          un = un - G_dn(i, k) * G_dn(k, i)
       enddo
       P0%meas(P0_CHIT, tmp) = P0%meas(P0_CHIT, tmp) + un*P0%up(i)*U(S%Map(i))
    enddo
    P0%meas(P0_CHIT, tmp) = P0%meas(P0_CHIT, tmp) - TWO * P0%meas(P0_NUD, tmp)

    !Scale by inverse temperature
    P0%meas(P0_CHIT, tmp) = P0%beta * P0%beta * P0%meas(P0_CHIT, tmp)

    !=========================================!
    ! Specific heat
    !=========================================!

    !< T T >
    P0%meas(P0_CV, tmp) = (Tbar - Cbar)**2 
    do i = 1, n
       do j = 1, n 
          P0%meas(P0_CV, tmp) = P0%meas(P0_CV, tmp) + (h_up(i, j) - A_up(i, j)) * A_up(j, i) &
             + (h_dn(i, j) - A_dn(i, j)) * A_dn(j, i)
       enddo
    enddo
    
    !< T U >
    P0%meas(P0_CV,tmp) = P0%meas(P0_CV,tmp) + (Tbar - Cbar)*P0%meas(P0_NUD, tmp) 
    do i = 1, n 
       !un = U(S%map(i)) * P0%up(i)
       un = ZERO
       do j = 1, n
          !P0%meas(P0_CV,tmp) = P0%meas(P0_CV,tmp) + un * G_dn(j,i) * (h_dn(i,j) - A_dn(i,j))
          un = un + G_dn(j, i) * (h_dn(i, j) - A_dn(i, j))
       enddo
       P0%meas(P0_CV,tmp) = P0%meas(P0_CV,tmp) + un *  U(S%map(i)) * P0%up(i)
       !un = U(S%map(i)) * P0%dn(i)
       un = ZERO
       do j = 1, n
          !P0%meas(P0_CV,tmp) = P0%meas(P0_CV,tmp) + un * G_up(j,i) * (h_up(i,j) - A_up(i,j))
          un = un + G_up(j, i) * (h_up(i, j) - A_up(i, j))
       enddo
       P0%meas(P0_CV, tmp) = P0%meas(P0_CV, tmp) + un *  U(S%map(i)) * P0%dn(i)
    enddo

    !< U T >
    P0%meas(P0_CV, tmp) = P0%meas(P0_CV, tmp) + (Tbar - Cbar)*P0%meas(P0_NUD, tmp) 
    do i = 1, n 
       P0%meas(P0_CV, tmp) = P0%meas(P0_CV, tmp) + U(S%map(i)) * P0%up(i) * A_dn(i,i)
       !un = U(S%map(i)) * P0%up(i)
       un = ZERO
       do j = 1, n
          !P0%meas(P0_CV,tmp) = P0%meas(P0_CV,tmp) - un * A_dn(i,j) * G_dn(j,i)
          un = un + A_dn(i,j) * G_dn(j,i)
       enddo
       P0%meas(P0_CV,tmp) = P0%meas(P0_CV,tmp) - un * U(S%map(i)) * P0%up(i)

       P0%meas(P0_CV,tmp) = P0%meas(P0_CV,tmp) + U(S%map(i)) * P0%dn(i) * A_up(i,i)
       !un = U(S%map(i)) * P0%dn(i)
       un = ZERO
       do j = 1, n
          !P0%meas(P0_CV,tmp) = P0%meas(P0_CV,tmp) - un * A_up(i,j) * G_up(j,i)
          un = un + A_up(i,j) * G_up(j,i)
       enddo
       P0%meas(P0_CV,tmp) = P0%meas(P0_CV,tmp) - un * U(S%map(i)) * P0%dn(i)
    enddo

    !< U U >
    ! Redefine A_up and A_dn
    do j = 1, n 
       do i = 1, n
          un = sqrt(U(S%Map(i))*U(S%Map(j)))
          A_up(i,j) = un * (P0%up(i)*P0%up(j) - G_up(i,j)*G_up(j,i))
          A_dn(i,j) = un * (P0%dn(i)*P0%dn(j) - G_dn(i,j)*G_dn(j,i))
       enddo
       A_up(j,j) = P0%up(j) * U(S%Map(j))
       A_dn(j,j) = P0%dn(j) * U(S%Map(j))
    enddo
    ! Compute UU contribution to Cv
    do j = 1, n 
       do i = 1, n
          P0%meas(P0_CV,tmp) = P0%meas(P0_CV,tmp) + A_up(i,j)*A_dn(i,j)
       enddo
    enddo

    ! Scale by inverse T
    P0%meas(P0_CV, tmp) = P0%beta * P0%beta * P0%meas(P0_CV, tmp)

    !=====================!
    ! Autocorelation term.!
    !=====================!
    if (P0%compSAF) then
       P0%meas(P0_SAF, tmp)  = TWO*n-P0%meas(P0_NUP, tmp)-&
            P0%meas(P0_NDN, tmp)
       P0%meas(P0_SAF2,tmp)  = P0%meas(P0_SAF, tmp)
    end if

    do i = 1,n
       do j = 1,n
          var1 = G_up(i,j) * G_up(j,i) + G_dn(i,j) * G_dn(j,i)
          var2 = -TWO * G_up(i,j) * G_dn(j,i)
          var3 = G_up(i,i) * G_up(j,j) + G_dn(i,i) * G_dn(j,j) - &
                 TWO * G_up(i,i) * G_dn(j,j) - var1
          
          ! k is the index
          k = S%D(i,j)
          ph = S%gf_phase(i,j)
          P0%G_fun(k, tmp) = P0%G_fun(k, tmp) + ph*(G_up(i,j) + G_dn(i,j))
          P0%Gf_up(k, tmp) = P0%Gf_up(k, tmp) + ph*G_up(i,j) 
          P0%Gf_dn(k, tmp) = P0%Gf_dn(k, tmp) + ph*G_dn(i,j)
          P0%Den0(k, tmp)  = P0%Den0(k, tmp) + &
               P0%up(i)*P0%up(j) + P0%dn(i)*P0%dn(j) - var1
          P0%Den1(k, tmp)  = P0%Den1(k, tmp) + P0%up(i)*P0%dn(j)
          P0%SpinXX(k, tmp) = P0%SpinXX(k, tmp) + var2
          P0%SpinZZ(k, tmp) = P0%SpinZZ(k, tmp) + var3
          P0%Pair(k,tmp)  = P0%Pair(k,tmp) + G_up(i,j) * G_dn(i,j)
          
          if (P0%compSAF) then
             var1 = S%P(i)*S%P(j)
             P0%meas(P0_SAF, tmp) = P0%meas(P0_SAF, tmp) + var1 * var2
             P0%meas(P0_SAF2,tmp) = P0%meas(P0_SAF2,tmp) + var1 * var3
          end if
       end do
       ! special case for (i,i)
       k = S%D(i,i)
       var1 =  G_up(i,i) + G_dn(i,i)
       P0%Den0(k, tmp)   = P0%Den0(k, tmp)   + var1
       P0%SpinXX(k, tmp) = P0%SpinXX(k, tmp) + var1
       P0%SpinZZ(k, tmp) = P0%SpinZZ(k, tmp) + var1
    end do
    
    P0%meas(P0_SFERRO, tmp) = sum(P0%SpinXX(:,tmp))
    P0%meas(P0_SFER2,  tmp) = sum(P0%SpinZZ(:,tmp))
    
    ! Average
    P0%meas(:,tmp) = P0%meas(:,tmp) / n
    do i = 1, P0%nClass
       P0%G_fun (i, tmp) = P0%G_fun (i, tmp) / S%F(i) * HALF
       P0%Gf_up (i, tmp) = P0%Gf_up (i, tmp) / S%F(i)
       P0%Gf_dn (i, tmp) = P0%Gf_dn (i, tmp) / S%F(i)
       P0%SpinXX(i, tmp) = P0%SpinXX(i, tmp) / S%F(i)
       P0%SpinZZ(i, tmp) = P0%SpinZZ(i, tmp) / S%F(i)
       P0%Den0  (i, tmp) = P0%Den0  (i, tmp) / S%F(i) * HALF
       P0%Den1  (i, tmp) = P0%Den1  (i, tmp) / S%F(i)
       P0%Pair(i, tmp)   = P0%Pair(i, tmp) / S%F(i) * HALF
    end do
    

    if (P0%compSAF) then
       P0%meas(P0_SAFSQ, tmp) = P0%meas(P0_SAF, tmp) * P0%meas(P0_SAF, tmp)
       P0%meas(P0_SAF2SQ,tmp) = P0%meas(P0_SAF2,tmp) * P0%meas(P0_SAF2,tmp)
    end if

    ! Accumulate result to P0(:, idx)
    sgn = sgnup * sgndn
    P0%meas(:, idx) =  P0%meas(:, idx) + P0%meas(:, tmp) * sgn

    m = P0%nClass
    call blas_daxpy(m, sgn, P0%G_fun (1:m,tmp), 1, P0%G_fun (1:m,idx), 1)
    call blas_daxpy(m, sgn, P0%Gf_up (1:m,tmp), 1, P0%Gf_up (1:m,idx), 1)
    call blas_daxpy(m, sgn, P0%Gf_dn (1:m,tmp), 1, P0%Gf_dn (1:m,idx), 1)
    call blas_daxpy(m, sgn, P0%SpinXX(1:m,tmp), 1, P0%SpinXX(1:m,idx), 1)
    call blas_daxpy(m, sgn, P0%SpinZZ(1:m,tmp), 1, P0%SpinZZ(1:m,idx), 1)
    call blas_daxpy(m, sgn, P0%Den0  (1:m,tmp), 1, P0%Den0  (1:m,idx), 1)
    call blas_daxpy(m, sgn, P0%Den1  (1:m,tmp), 1, P0%Den1  (1:m,idx), 1)
    call blas_daxpy(m, sgn, P0%Pair  (1:m,tmp), 1, P0%Pair(1:m,idx), 1)

    P0%sign(P0_SGN,   idx) =  P0%sign(P0_SGN,   idx) + sgn
    P0%sign(P0_SGNUP, idx) =  P0%sign(P0_SGNUP, idx) + sgnup
    P0%sign(P0_SGNDN, idx) =  P0%sign(P0_SGNDN, idx) + sgndn
    P0%cnt = P0%cnt + 1

  end subroutine DQMC_Phy0_Meas

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy0_GetFT(P0, class, phase, ft_wgt_t, ft_wgt_g, nkt, nkg, na, nt)
    !
    ! Purpose
    ! =======
    !    This subroutine computes the FT of the correlation functions
    !    stored in phy0. It fills the bins transforming the jackknived
    !    values of the real space correlation functions. The error is
    !    estimated using the standard jackknife formula.
    !   
    ! Comments
    ! ========
    !    LAPACK could have been used in, at least, two places but this
    !    routine is not critical (only called at the end) so that
    !    having explicit matrix multiplication improves readability
    !    without degrading performance.
    !
    ! Arguments
    ! =========
    integer, intent(in)    :: na                  !number of site in primitive cell
    integer, intent(in)    :: nt                  !number of translations
    integer, intent(in)    :: nkt,nkg              !number of k-points
    integer, intent(in)    :: class(na*nt,na*nt)  !classes of pairs of sites
    integer, intent(in)    :: phase(na*nt,na*nt)  !classes of pairs of sites
    complex*16, intent(in), target :: ft_wgt_t(nt, nkt)      !fourier weights : exp(ikr)
    complex*16, intent(in), target :: ft_wgt_g(nt, nkg)    !fourier weights : exp(ikr)
    type(Phy0), intent(inout) :: P0                  !container

    ! ... Local variables ...
    real(wp) :: rwork(3*na), phcurr
    real(wp), pointer :: curr(:) 
    real(wp), pointer :: currft(:) 
    integer :: ph(na*nt,na*nt)
    integer :: nak, ip, ibin, nBin
    integer :: avg, err, ia, ja, ik, i, j, it, nk, jt
    complex*16 ::  work(2*na), U(na,na), W(na,na)
    complex*16, pointer :: ft_wgt(:,:) 
    complex*16, parameter :: ZEROZ=(0.d0, 0.d0), ONEZ=(1.d0, 0.d0)

    nBin = P0%nBin
    avg  = P0%avg
    err  = P0%err
    nk   = max(nkt, nkg)
    nak  = na * nk
 
    !Initialize pointers to Fourier transform
    do ip = 1, narrays + 1 
       P0%IARREV(ip) = (ip-1)*nak + 1
       P0%IARRFT(ip) = (ip-1)*nak*(na+1)/2 + 1
    enddo

    !Allocate storage for Fourier transform
    allocate(P0%AllPropFT(nak*(na+1)*narrays/2, P0%err))
    allocate(P0%AllPropEigVal(nak*narrays, P0%err))
    allocate(P0%AllPropEigVec(na, na, nk, narrays))
    P0%initFT = .true.

    !Loop over properties to Fourier transform
    do ip = 1, narrays

       !select k-grid with twist or not (only gfun's require twist)
       select case(ip) 
         case(IGFUN, IGFUP, IGFDN) 
            nk = nkt
            ft_wgt => ft_wgt_t
            ph = phase
         case default
            nk = nkg
            ft_wgt => ft_wgt_g
            ph = nint(ONE)
       end select

       !Fourier transform each bin and average
       do ibin = P0%avg, 1, -1

          !Pointer to property "ip" in bin "ibin"
          curr   => P0%AllProp(P0%IARR(ip):P0%IARR(ip+1)-1, ibin) 

          !Loop over inequvalent k-points
          do ik = 1, nk

             U = 0.d0
             !Compute Fourier transform matrix at "ik" for (ja,ia)
             do ia = 1, na
                do ja = ia, na
                   !sum over translations
                   do it = 1, nt
                      do jt = 1, nt
                         !Find atom which is the translation of "ja" by "it"
                         i = (it-1) * na + ia
                         j = (jt-1) * na + ja
                         !Use class that corresponds to the (i,j) pair 
                         phcurr = ph(i,j) * curr(class(i,j))
                         U(ja,ia) = U(ja,ia) + phcurr * ft_wgt(it,ik) * dconjg(ft_wgt(jt,ik))
                      enddo
                   enddo
                enddo
             enddo

             !Pointer to Fourier transform of "ip" at "ik" in "ibin"
             i = P0%IARRFT(ip) + (ik-1)*na*(na+1)/2
             do ia = 1, na
                do ja = ia, na
                   P0%AllPropFT(i, ibin) = U(ja,ia)
                   i = i + 1
                enddo
             enddo
             currft => P0%AllPropEigVal(P0%IARREV(ip)+(ik-1)*na:P0%IARREV(ip)+ik*na-1, ibin) 

             !Diagonalize Fourier transform matrix. Store eigenvalues in AllPropEigVal.
             if (ibin == P0%avg) then
                call zheev('V', 'L', na, U, na, currft, work, 2*na, rwork, it)
                !Store eigenvectors as well for the average
                P0%AllPropEigVec(:,:,ik,ip) = U
             else
                call zheev('N', 'L', na, U, na, currft, work, 2*na, rwork, it)
                call zgemm('N', 'N', na, na, na, ONEZ, U, na, P0%AllPropEigVec(:,:,ik,ip), na, ZEROZ, W, na)
                call zgemm('T', 'N', na, na, na, ONEZ, P0%AllPropEigVec(:,:,ik,ip), na, W, na, ZEROZ, U, na)
                do i = 1, na
                   currft(i) = dreal(U(i,i))
                enddo
             endif

          enddo ! Loop over k-points

       enddo ! Loop over bins

    enddo ! Loop over properties

  end subroutine DQMC_Phy0_GetFT

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy0_GetErrFT(P0)
    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !    This subroutine computes the error the Fourier Transform of all properties.
    !
    !  Pre-assumption
    ! ===============
    !    DQMC_Phy0_GetFT was called.
    !
    ! Arguments
    ! =========
    !
    type(Phy0), intent(inout) :: P0                  !container

    integer :: err, avg, n, m, nbin, i, nproc

    nbin = P0%nbin
    err  = P0%err
    avg  = P0%avg
    nproc = qmc_sim%size

    n = size(P0%AllPropEigVal, 1)
    m = size(P0%AllPropFT, 1)

    if (nproc > 1) then

       !Compute errorbars avgraging over processors
#      ifdef _QMC_MPI
          call mpi_allreduce(P0%AllPropEigVal(:,1)**2, P0%AllPropEigVal(:,err), n, mpi_double, &
             mpi_sum, mpi_comm_world, i)
          call mpi_allreduce(P0%AllPropFT(:,1)**2, P0%AllPropFT(:,err), m, mpi_double, &
             mpi_sum, mpi_comm_world, i)
#      endif
       P0%AllPropEigVal(:,err) = P0%AllPropEigVal(:,err) / dble(nproc) - P0%AllPropEigVal(:,avg)**2 
       P0%AllPropEigVal(:,err) = sqrt(P0%AllPropEigVal(:,err) * dble(nproc-1))

       P0%AllPropFT(:,err) = P0%AllPropFT(:,err) / dble(nproc) - P0%AllPropFT(:,avg)**2 
       P0%AllPropFT(:,err) = sqrt(P0%AllPropFT(:,err) * dble(nproc-1))

    else

       !Compute errorbars using bins
       do i = 1, n
          P0%AllPropEigVal(i,err) = sqrt((nbin-1)*sum((P0%AllPropEigVal(i,avg) - P0%AllPropEigVal(i,1:nbin))**2)/nbin)
       enddo
  
       do i = 1, m
          P0%AllPropFT(i,err) = sqrt((nbin-1)*sum((P0%AllPropFT(i,avg) - P0%AllPropFT(i,1:nbin))**2)/nbin)
       enddo
  
    endif

  end  subroutine DQMC_Phy0_GetErrFT

  !--------------------------------------------------------------------!
  
  subroutine DQMC_Phy0_PrintFT(P0, na, nkt, nkg, OPT)

    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !    This subroutine prints out the Fourier transform.
    !    Labels need to be produced separately.
    !
    !  Pre-assumption
    ! ===============
    !    OPT is a file handle
    !    DQMC_Phy0_GetErrFT was called.
    !
    ! Arguments
    ! =========
    !
    type(Phy0), intent(in)    :: P0        ! Phy0
    integer, intent(in)       :: OPT       ! Output file handle
    integer, intent(in)       :: na        ! number of sites in unit cell
    integer, intent(in)       :: nkt, nkg  ! number of non-equivalent k-points

    ! ... Local varaiables ...
    integer :: i, ia, ja, ii, avg, err, nakt, nakg
    real(wp), pointer :: FTptr(:,:) 
    complex*16, pointer :: Nmptr(:,:,:) 
    character(len=30), pointer :: clabel(:) 

    if (qmc_sim%rank /= 0) return

    avg    = P0%avg
    err    = P0%err
    nakt    = na*nkt
    nakg    = na*nkg

    allocate(clabel(max(nakt,nakg)*(na+1)/2))

    ! ... Executable ...
    ii = 0 
    do i = 1, nkt
       do ia = 0, na - 1
          do ja = ia, na - 1
             ii = ii + 1
             if (ia + ja == 0) then
                write(clabel(ii),'(3(i5))')i, ja, ia
             else
                write(clabel(ii),'(5x,2(i5))')ja, ia
             endif
          enddo
       enddo
    enddo
    FTptr => P0%AllPropFT(P0%IARRFT(IGFUN):P0%IARRFT(IGFUN+1)-1,:)
    call DQMC_Print_RealArray(0, nakt*(na+1)/2, "FT of Ave Equal t Green's function:", &
         clabel, FTptr(:, avg:avg), FTptr(:, err:err), OPT)

    FTptr => P0%AllPropFT(P0%IARRFT(IGFUP):P0%IARRFT(IGFUP+1)-1,:)
    call DQMC_Print_RealArray(0, nakt*(na+1)/2, "FT of Up Equal t Green's function:", &
         clabel, FTptr(:, avg:avg), FTptr(:, err:err), OPT)

    FTptr => P0%AllPropFT(P0%IARRFT(IGFDN):P0%IARRFT(IGFDN+1)-1,:)
    call DQMC_Print_RealArray(0, nakt*(na+1)/2, "FT of Dn Equal t Green's function:", &
         clabel, FTptr(:, avg:avg), FTptr(:, err:err), OPT)

    ii = 0 
    do i = 1, nkg
       do ia = 0, na-1
          do ja = ia, na-1
             ii = ii + 1
             if (ia+ja == 0) then
                write(clabel(ii),'(3(i5))')i, ja, ia
             else
                write(clabel(ii),'(5x,2(i5))')ja, ia
             endif
          enddo
       enddo
    enddo
    FTptr => P0%AllPropFT(P0%IARRFT(IDEN0):P0%IARRFT(IDEN0+1)-1,:)
    call DQMC_Print_RealArray(0, nakg*(na+1)/2, "FT of Density-density correlation fn: (up-up)", &
         clabel, FTptr(:, avg:avg), FTptr(:, err:err), OPT)

    FTptr => P0%AllPropFT(P0%IARRFT(IDEN1):P0%IARRFT(IDEN1+1)-1,:)
    call DQMC_Print_RealArray(0, nakg*(na+1)/2, "FT of Density-density correlation fn: (up-dn)", &
         clabel, FTptr(:, avg:avg), FTptr(:, err:err), OPT)

    FTptr => P0%AllPropFT(P0%IARRFT(ISPXX):P0%IARRFT(ISPXX+1)-1,:)
    call DQMC_Print_RealArray(0, nakg*(na+1)/2, "FT of XX spin correlation fn:", &
         clabel, FTptr(:, avg:avg), FTptr(:, err:err), OPT)

    FTptr => P0%AllPropFT(P0%IARRFT(ISPZZ):P0%IARRFT(ISPZZ+1)-1,:)
    call DQMC_Print_RealArray(0, nakg*(na+1)/2, "FT of ZZ spin correlation fn:", &
         clabel, FTptr(:, avg:avg), FTptr(:, err:err), OPT)

    FTptr => P0%AllPropFT(P0%IARRFT(IAVSP):P0%IARRFT(IAVSP+1)-1,:)
    call DQMC_Print_RealArray(0, nakg*(na+1)/2, "FT of Average spin correlation fn:", &
         clabel, FTptr(:, avg:avg), FTptr(:, err:err), OPT)

    FTptr => P0%AllPropFT(P0%IARRFT(IPAIR):P0%IARRFT(IPAIR+1)-1,:)
    call DQMC_Print_RealArray(0, nakg*(na+1)/2, "FT of Pairing correlation fn:", &
         clabel, FTptr(:, avg:avg), FTptr(:, err:err), OPT)

    if (na > 1) then

       ii = 0 
       do i = 1, nkt
          do ia = 1, na
             ii = ii + 1
             if (ia == 1) then
                write(clabel(ii),'(2(i5))')i, ia
             else
                write(clabel(ii),'(5x,i5)')ia
             endif
          enddo
       enddo
       FTptr => P0%AllPropEigVal(P0%IARREV(IGFUN):P0%IARREV(IGFUN+1)-1,:)
       call DQMC_Print_RealArray(0, nakt, "Eigenvalues of Ave Equal t Green's function:", &
            clabel, FTptr(:, avg:avg), FTptr(:, err:err), OPT)

       FTptr => P0%AllPropEigVal(P0%IARREV(IGFUP):P0%IARREV(IGFUP+1)-1,:)
       call DQMC_Print_RealArray(0, nakt, "Eigenvalues of Up Equal t Green's function:", &
            clabel, FTptr(:, avg:avg), FTptr(:, err:err), OPT)

       FTptr => P0%AllPropEigVal(P0%IARREV(IGFDN):P0%IARREV(IGFDN+1)-1,:)
       call DQMC_Print_RealArray(0, nakt, "Eigenvalues of Dn Equal t Green's function:", &
            clabel, FTptr(:, avg:avg), FTptr(:, err:err), OPT)

       ii = 0 
       do i = 1, nkg
          do ia = 1, na
             ii = ii + 1
             if (ia == 1) then
                write(clabel(ii),'(2(i5))')i, ia
             else
                write(clabel(ii),'(5x,i5)')ia
             endif
          enddo
       enddo
       FTptr => P0%AllPropEigVal(P0%IARREV(IDEN0):P0%IARREV(IDEN0+1)-1,:)
       call DQMC_Print_RealArray(0, nakg, "Eigenvalues of Density-density correlation fn: (up-up)", &
            clabel, FTptr(:, avg:avg), FTptr(:, err:err), OPT)

       FTptr => P0%AllPropEigVal(P0%IARREV(IDEN1):P0%IARREV(IDEN1+1)-1,:)
       call DQMC_Print_RealArray(0, nakg, "Eigenvalues of Density-density correlation fn: (up-dn)", &
            clabel, FTptr(:, avg:avg), FTptr(:, err:err), OPT)

       FTptr => P0%AllPropEigVal(P0%IARREV(ISPXX):P0%IARREV(ISPXX+1)-1,:)
       call DQMC_Print_RealArray(0, nakg, "Eigenvalues of XX spin correlation fn:", &
            clabel, FTptr(:, avg:avg), FTptr(:, err:err), OPT)

       FTptr => P0%AllPropEigVal(P0%IARREV(ISPZZ):P0%IARREV(ISPZZ+1)-1,:)
       call DQMC_Print_RealArray(0, nakg, "Eigenvalues of ZZ spin correlation fn:", &
            clabel, FTptr(:, avg:avg), FTptr(:, err:err), OPT)

       FTptr => P0%AllPropEigVal(P0%IARREV(IAVSP):P0%IARREV(IAVSP+1)-1,:)
       call DQMC_Print_RealArray(0, nakg, "Eigenvalues of Average spin correlation fn:", &
            clabel, FTptr(:, avg:avg), FTptr(:, err:err), OPT)

       FTptr => P0%AllPropEigVal(P0%IARREV(IPAIR):P0%IARREV(IPAIR+1)-1,:)
       call DQMC_Print_RealArray(0, nakg, "Eigenvalues of Pairing correlation fn:", &
            clabel, FTptr(:, avg:avg), FTptr(:, err:err), OPT)

       Nmptr => P0%AllPropEigVec(:,:,:,IGFUN)
       call DQMC_Print_EigenMode(na, nkt, "Eigenmodes of Ave Equal t Green's function (Natural orbitals):", &
             Nmptr, OPT)

       Nmptr => P0%AllPropEigVec(:,:,:,IGFUP)
       call DQMC_Print_EigenMode(na, nkt, "Eigenmodes of Up Equal t Green's function (Natural orbitals):", &
             Nmptr, OPT)

       Nmptr => P0%AllPropEigVec(:,:,:,IGFDN)
       call DQMC_Print_EigenMode(na, nkt, "Eigenmodes of Down Equal t Green's function (Natural orbitals):", &
             Nmptr, OPT)

       Nmptr => P0%AllPropEigVec(:,:,:,IDEN0)
       call DQMC_Print_EigenMode(na, nkg,"Eigenmodes of Density-density correlation fn: (up-up)", &
             Nmptr, OPT)

       Nmptr => P0%AllPropEigVec(:,:,:,IDEN1)
       call DQMC_Print_EigenMode(na, nkg, "Eigenmodes of Density-density correlation fn: (up-dn)", &
             Nmptr, OPT)

       Nmptr => P0%AllPropEigVec(:,:,:,ISPXX)
       call DQMC_Print_EigenMode(na, nkg, "Eigenmodes of XX-Spin correlation fn: ", &
             Nmptr, OPT)

       Nmptr => P0%AllPropEigVec(:,:,:,ISPZZ)
       call DQMC_Print_EigenMode(na, nkg, "Eigenmodes of ZZ-Spin correlation fn: ", &
             Nmptr, OPT)

       Nmptr => P0%AllPropEigVec(:,:,:,IAVSP)
       call DQMC_Print_EigenMode(na, nkg, "Eigenmodes of Average Spin correlation fn: ", &
             Nmptr, OPT)

       Nmptr => P0%AllPropEigVec(:,:,:,IPAIR)
       call DQMC_Print_EigenMode(na, nkg, "Eigenmodes of Pairing correlation fn: ", &
             Nmptr, OPT)

    endif

    deallocate(clabel)

  end subroutine DQMC_Phy0_PrintFT

end module DQMC_Phy0
