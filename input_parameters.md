## Files. ##
  1. ofile - output file.
  1. gfile - geometry file. [Example, how to create a geometry file.](https://code.google.com/p/quest-qmc/wiki/Lieb_lattice_tutorial)

## Hubbard Model. ##
  1. mu\_up(real) - chemical potential for spins up.
  1. mu\_dn(real) - chemical potential for spins down.
  1. L (integer) - number of the time slices.
  1. dtau (real) - beta discretezation step.
  1. HSF(integer) - indicator of how Hubbard Stratonovich field(HSF) is input.
    1. -1, randomly generating HSF.
    1. 0, use HSF in memeory
    1. 1, read HSF from file
  1. bcond(real, real, real) - boundary conditions. [How to define boundary conditions.](https://code.google.com/p/quest-qmc/wiki/BoundaryConditions)

## Metropolis algorithm ##
  1. nwarm(integer) - number of warm up sweeps.
  1. npass(integer) - number of measurement sweep.
  1. ntry(integer) - number of sites to be flipped in the global sweep. This parameter controls the frequency of a type of Monte Carlo move which is needed for big U (U > 8). [Explanation, why do we need this](http://prb.aps.org/abstract/PRB/v44/i19/p10502_1)
    * A global move costs N^3
cpu time compared to
N^2 for a move of a single HS variable. So they are more expensive.
    * On the other hand, you are changing L HS variables, not just one.
  1. tausk(integer) - specifies the frequency of performing physical measurements.
  1. tdm(integer): this parameter can have values 0 or 1 an it defines if perform time dependent measurements or not.
    * 0 - we do not perform time dependent measurements.
    * 1 - we perform time dependent measurements.
  1. seed(integer) - random number seed.
## Measurement. ##
  1. nbin(integer) - determines how the computed data is divided into bins.
  1. nhist(integer) -
## Numerical. ##
  1. north(integer) - frequency of performing orthogonalization in matrix product. This is needed when we use stabilization algorithm in the calculation of the Green's function.
  1. nwrap(integer) - frequency of performing recomputing in Green's function calculation. In QUEST nwrap will be dynamically adjusted according to the errors of updating. When the difference in G is higher than "difflim" with a frequency higher than "errate", the "nwrap" is reduced. For small U's, using large nwrap is fine because matrix elements are well-behaved. At large interaction strengths, there is huge scale difference between elements of B  matrices. As a result, numerical error accumulates rather quickly if nwrap is too large.
  1. fixwrap(real) - NO IDEA. BUT! **When ASQRD is enabled and fix wrap is set to 1, north and nwrap cannot be equal.**
  1. errate(real) - tolerable error rate of recomputing.
  1. difflim(real) - tolerable difference of the matrices computed from different methods.

## Example. ##
Python template that I use to generate the input file.

```
# ==========================
# files
# ==========================
ofile = {ofile}
gfile = {gfile}
#=======================
# Hubbard Model
#=======================
mu_up = {mu}
mu_dn = {mu}
L = {nSlices}
dtau = {dtau}
HSF = -1
bcond = 0.0, 0.0, 0.0
#=======================
# Metropolis Algorithm
#=======================
nwarm = {nwarm}
npass = {npass}
ntry = {ntry}
tausk = 1
tdm = {tdm}
seed = 1234567
#=======================
# Measurements
#=======================
nbin = 5
nhist = 1
#=======================
# Numerical
#=======================
north = {north}
nwrap = 15
fixwrap = 1
errrate = 0.01
difflim = 0.000001
```