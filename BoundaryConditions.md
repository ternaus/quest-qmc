
# Introduction #
In Physics we are interested in describing properties of the large systems. Or we are trying to match our simulation with theoretical results which are normally derived for the infinite number of particles.


So, right now we work with finite lattices and than try to extrapolate our results to the infinity. There are many different boundary conditions that we can use. For example, for some models it was proved that open boundary conditions lead to the 1/N corrections, periodic boundary conditions lead to the 1 / N^2 corrections. Sometimes antiperiodic boundary conditions is better than periodic, sometimes they are not.

So for all practical purposes it makes sense to do simulation with different cases and than average over them.

# How can you change boundary conditions in QUEST? #
## Periodic and antiperiodic boundary conditions. ##
(Based on the great explanation of Simone Chiesa, nearly copy-paste from his words.).

In the input file there is a line that says something like:

```
bcond = x, y, z
```
where `x, y, z` are real numbers.

What do these numbers mean?

Those three numbers `(x,y,z)` define the phase of the hopping as you cross the cluster boundary:

  * `exp(i x pi)` as you wrap around the cluster along the 1st supercell vector
  * `exp(i y pi)` as you wrap around the cluster along the 2nd supercell vector
  * `exp(i z pi)` as you wrap around the cluster along the 3rd supercell vector

So `0.0, 0.0, 0.0` means periodic. `1.0, 1.0, 1.0` means fully antiperiodic. You can only use 1.0 or 0.0 (any combination of them) since the code still uses real variables.

## Open boundary conditions. ##
Right now only way to do this is to define big unit cell.