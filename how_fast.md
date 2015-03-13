
# Introduction #

Easiest way to figure out how complexity of your algorithm depends on the input parameters is to do several independent simulations and plot time that each simulation uses versus different parameters.

# Results. #
**All results are only for time independent measurements!**

<img src='https://dl.dropboxusercontent.com/u/3483290/quest-qmc/t_vs_N.png' alt='time versus num_sites' height='400' width='440'>

<img src='https://dl.dropboxusercontent.com/u/3483290/quest-qmc/t_vs_L.png' alt='time versus num_slices' height='400' width='420'>

There is not enough data for accurate conclusion, but roughly, time grows as a power low in both cases, but with a different exponent.<br>
<br>
To get he power exactly we switch to the log-log scale:<br>
<br>
<img src='https://dl.dropboxusercontent.com/u/3483290/quest-qmc/log_t_vs_logN.png' alt='log (time) versus log(num_sites)' height='400' width='400'>

<img src='https://dl.dropboxusercontent.com/u/3483290/quest-qmc/logt_vs_logL.png' alt='log(time) versus log(num_slices)' height='400' width='400'>

<h1>Conclusion</h1>
<ol><li>When we increase number of sites, time grows as N^2.70.<br>
</li><li>When we increase number of time slices, time grows as L^1.43.</li></ol>

I'd say, time grows too fast. Using openMP, MPI, GPU, different type of the compiler, compiler flags will not change this power law.<br>
<br>
<h1>Set up.</h1>
<h2>Computer details.</h2>

Typical computer on the cougar cluster will have:<br>
<br>
<ol><li>Linux machines with <code>kernel 2.6.32-220.13.1.el6.x86_64</code>
</li><li>Intel Pentium 4 CPU 3.20GHz</li></ol>

<h2>Compilation details.</h2>

<ol><li>BLAS/LAPACK library - Intel mkl.<br>
</li><li>ifort compiler flags - <code>-m64 -O3 -unroll</code>
</li><li>openMP and MPI were not used.</li></ol>

<h2>Input parameters</h2>

Input file template:<br>
<br>
<pre><code># ==========================<br>
# files<br>
# ==========================<br>
ofile = {ofile}<br>
gfile = {gfile}<br>
#=======================<br>
# Hubbard Model<br>
#=======================<br>
mu_up = {mu}<br>
mu_dn = {mu}<br>
L = {nSlices}<br>
dtau = 0.1<br>
HSF = -1<br>
bcond = 0.0, 0.0, 0.0<br>
#=======================<br>
# Metropolis Algorithm<br>
#=======================<br>
nwarm = 1000<br>
npass = 50000<br>
ntry = 0<br>
tausk = 1<br>
tdm = 1<br>
#=======================<br>
# Measurements<br>
#=======================<br>
nbin = 5<br>
nhist = 1<br>
#=======================<br>
# Numerical<br>
#=======================<br>
north = 5<br>
nwrap = 15<br>
fixwrap = 1<br>
errrate = 0.01<br>
difflim = 0.000001<br>
</code></pre>

Geometry file template:<br>
<pre><code>#NDIM<br>
2<br>
#PRIM<br>
1.0 0.0<br>
0.0 1.0<br>
#SUPER<br>
{nx} 0<br>
0 {ny}<br>
#ORB<br>
s0 0.0d0 0.0d0 0.0d0 #0<br>
s1 0.5d0 0.0d0 0.0d0 #1<br>
s2 0.0d0 0.5d0 0.0d0 #2<br>
#HAMILT <br>
0 1 0.5 0.0 0.0 1.0 1.0 0.0<br>
0 1 -0.5 0.0 0.0 1.0 1.0 0.0<br>
0 2 0.0 0.5 0.0 1.0 1.0 0.0<br>
0 2 0.0 -0.5 0.0 1.0 1.0 0.0<br>
0 0 0.0 0.0 0.0 0.0 0.0 -4.0<br>
1 1 0.0 0.0 0.0 0.0 0.0 -4.0<br>
2 2 0.0 0.0 0.0 0.0 0.0 -4.0<br>
#SYMM<br>
d  0.0d0 0.0d0 0.0d0 1.0d0 0.0d0 0.0d0<br>
d  0.0d0 0.0d0 0.0d0 0.0d0 1.0d0 0.0d0<br>
c4 0.0d0 0.0d0 0.0d0 0.0d0 0.0d0 1.0d0<br>
#PHASE<br>
1 0<br>
0 1<br>
s0 0.0 0.0 0.0  1.0<br>
s1 0.5 0.0 0.0  -1.0<br>
s2 0.0 0.5 0.0  -1.0<br>
#BONDS   <br>
0 0  0.0  0.0  0.0 # 1<br>
#0 1  0.5  0.0  0.0 # 2 -2<br>
#0 2  0.0  0.5  0.0 # 3 -3<br>
#PAIR      <br>
		       1     <br>
s-wave    1.0 <br>
#END<br>
</code></pre>