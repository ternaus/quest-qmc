program dqmc_test
  ! This program checks the execution time for
  ! QUEST on a 2-dimensional PEriodic Reactangular
  ! Lattice (2DPERL)

  use DQMC_2DPERL
  use DQMC_MPI
  implicit none

  integer :: t1, t2, rate

#ifdef DQMC_PROFILE
  !gfun_profile = .true.
  !matb_profile = .true.
  call profile_enable()
#endif

  call system_clock(t1)
  !count number of processors
  call DQMC_MPI_Init(qmc_sim, PLEVEL_1)

  call DQMC_Comp_2DPerl

  call system_clock(t2, rate)
  write(STDOUT,*) "Running time:",  (t2-t1) / REAL(rate), "(second)"

#ifdef DQMC_PROFILE
  !call gfun_print()
  !call matb_print()
  call profile_print()
#endif

end program dqmc_test
