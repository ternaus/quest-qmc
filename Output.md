# Input/configuration parameters. #
  1. U - (Real)
  1. t\_up - (Real)
  1. t\_down - (Real)
  1. mu\_up - (Real)
  1. mu\_down - (Real)
  1. Time slice - L - (Integer). Number of time slice.
  1. Number of sites
  1. dtau - (Real) discretization parameter.
  1. beta - inverse temperature.
  1. Number of warmup sweep - how many Monte Carlo loops need be executed for warm up.
  1. Number of measurement sweep - how many Monte Carlo loops need be executed for measurement.
  1. Frequency of measurement
  1. Random seed - (Integer) Random number seed.
  1. Frequency of recomputing G
  1. Global move number of sites
  1. Accept count
  1. Reject count
  1. Approximate accept rate
  1. gamma
  1. Type of matrix B
  1. Type of matrix HSF

# Results of the simulation. #
If a measurement is a single real number, it will be shown with three
terms:
  * name
  * average
  * error



## Sign of equal time measurements: ##
  1. Avg sign
  1. Avg up sign
  1. Avg dn sign
## Equal time measurements. ##
  1. Up spin occupancy -(Real.) in the literature is defined as `<n_up>`
  1. Down spin occupancy - (Real.) in the literature is defined as `<n_down>`
  1. <U\*N\_up\*N\_dn> -
  1. Kinetic energy -(Real.) KE = -\sum_{ij\sigma} t_{ij\sigma} (c^{\dagger}_{ij\sigma} c_{ij \sigma} + h.c) - \sum\_i (mu\_i + U\_i / 2) n\_i
. This is really weird name for this part of the energy, at least because from the elementary physics course we know that Kinetic energy should be non-negative, which is not the case here.
  1. Total Energy - This is not average of the Hamiltonian as you could guess, but! Total energy is Kinetic energy defined above and <U N\_up N\_dn>. So finally it is hopping term +` <U N_up N_dn> + mu  <n>`
  1. Density - (Real). In literature is defined as $\rho$ or `<n>`. Should be equal to `<n_up> + <n_down>.`
  1. Chi\_thermal
  1. Specific heat - specific heat.
  1. XX Ferro structure factor
  1. ZZ Ferro structure factor
  1. XX AF structure factor
  1. Root Mean Square of XX AF
  1. ZZ AF structure factor
  1. Root Mean Square of ZZ AF
  1. Mean Equal time Green's function
  1. Up Equal time Green's function
  1. Down Equal time Green's function
  1. Density-density correlation fn: (up-up)
  1. Density-density correlation fn: (up-dn)
  1. XX Spin correlation function
  1. ZZ Spin correlation function
  1. Average Spin correlation function
  1. Pairing correlation function
  1. Pair-field correlation function - accumulated:
    1. s- wave(real) - s-wave pairing correlation function. Not normalized per site.