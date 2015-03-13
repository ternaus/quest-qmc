# Introduction #

Sometimes we need to add new measurements to the output.


# Single number. #


At the very  beginning of DQMC\_Phy0 you see a bunch of parameter
definition (P0\_NUP, P0\_DN...) which is just an emulation of C/C++ enum.

So the first step in adding an observable is adding an entry here. Note that the number of variables is either P0\_N or P0\_N\_NO\_SAF depending on the value of P0%compSAF (see DQMC\_Phy0\_Init). These constant will both need to be incremented. Something like this will do
```
  ! Parameter for the index of scalar variables (IMEAS)
  integer, parameter :: P0_NUP       = 1
  integer, parameter :: P0_NDN       = 2
  integer, parameter :: P0_NUD       = 3
  integer, parameter :: P0_KE        = 4
  integer, parameter :: P0_ENERGY    = 5
  integer, parameter :: P0_DENSITY   = 6
  integer, parameter :: P0_CHIT      = 7
  integer, parameter :: P0_CV        = 8
  integer, parameter :: P0_VLAD      = 9

  integer, parameter :: P0_SFERRO    = 10
  integer, parameter :: P0_SFER2     = 11
  integer, parameter :: P0_SAF       = 12
  integer, parameter :: P0_SAFSQ     = 13
  integer, parameter :: P0_SAF2      = 14
  integer, parameter :: P0_SAF2SQ    = 15

  integer, parameter :: P0_N_NO_SAF  = 11
  integer, parameter :: P0_N         = 15
```
Then add the corresponding "label" for printing like so
```
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
       "           Vlad new property:",&
       "  XX Ferro structure factor : ", &
       "  ZZ Ferro structure factor : ", &
       "     XX AF structure factor : ", &
       "  Root Mean Square of XX AF : ", &
       "     ZZ AF structure factor : ", &
       "  Root Mean Square of ZZ AF : "/)

```
Finally code up the actual measurement in DQMC\_Phy0\_Meas
```
P0%meas(P0_VLAD, tmp) = "some product/sum of G's"
```