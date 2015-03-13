
# Introduction #

Lieb lattice is square depleted lattice, that looks like:

<img src='https://dl.dropboxusercontent.com/u/3483290/lattice.png' alt='Smiley face' height='100' width='100'>

It is square lattice with unit cell, containing 3 lattice sites.<br>
<br>
To do simulations with it we need to write an input file that describes Lattice and Hamiltonian acting on it. In the geometry file below I use:<br>
<ol><li>Distance between sites to be equal to 0.5<br>
</li><li>Number of the unit cells in the x direction is 7.<br>
</li><li>Number of the unit cells in the y direction is 5.<br>
</li><li>Hopping for the upspin electrons is 1.<br>
</li><li>Hopping for the downspin electrons is 1.<br>
</li><li>U term on the each site site is -4.</li></ol>

<h1>Geometry file definition.</h1>
<h2>NDIM</h2>
<pre><code>#NDIM<br>
2<br>
</code></pre>

Dimension of the lattice. For our case it is 2.<br>
<br>
<h2>PRIM</h2>
<pre><code>#PRIM<br>
1.0 0.0<br>
0.0 1.0<br>
</code></pre>
Primitive lattice vectors. I've chosen distance between sites to be 0.5, and this leads to the two primitive vectors of the length 1. One pointing in the x and another pointing in the y direction.<br>
<h2>SUPER</h2>
<pre><code>#SUPER<br>
7 0<br>
0 5<br>
</code></pre>
This is supercell lattice vectors. This means that lattice can be obtained  if we use translation in the x direction up to 7 times and up to 5 times in the y direction.<br>
<br>
To create simulation for different parameters I use template that is used in the python script. Part of the template for this section looks like:<br>
<pre><code>#SUPER<br>
{nx} 0<br>
0 {ny}<br>
</code></pre>
<h2>ORB</h2>
<pre><code>#ORB<br>
s0 0.0d0 0.0d0 0.0d0 #0<br>
s1 0.5d0 0.0d0 0.0d0 #1<br>
s2 0.0d0 0.5d0 0.0d0 #2<br>
</code></pre>
This part of the file is about unit cell. 3 atoms in the unit cell, so here should be 3 lines. (1 line per atom.) The first entry is a string which serves as a label for the atom. The next three entries in the line are the position of the atom in the unit cell. QUEST automatically assigns a number to each atom (ie to each line in the #ORB section of the .geom file) beginning with zero. In the #ORB section QUEST demands three dimensional objects.<br>
<br>
<h2>HAMILT</h2>
<pre><code>#HAMILT <br>
0 1 0.5 0.0 0.0 1.0 1.0 0.0<br>
0 1 -0.5 0.0 0.0 1.0 1.0 0.0<br>
0 2 0.0 0.5 0.0 1.0 1.0 0.0<br>
0 2 0.0 -0.5 0.0 1.0 1.0 0.0<br>
0 0 0.0 0.0 0.0 0.0 0.0 -4.0<br>
1 1 0.0 0.0 0.0 0.0 0.0 -4.0<br>
2 2 0.0 0.0 0.0 0.0 0.0 -4.0<br>
</code></pre>
This block of the code defines the hoppings, on-site energies, and interaction strengths.<br>
<br>
The convention for the lines is following:<br>
Each line begins with two entries which are the (automatically assigned) atom numbers from #ORB. For problems without a basis these will just be ’0.0’. To include a hopping, the next three entries in the line should be the direction to the neighboring site. QUEST automatically makes H Hermitian, so each pair of connected sites requires only one line. The next two entries are the hopping values for the up and down electrons, which are allowed to be different in QUEST. The final entry is U , which should be set to zero for lines defining hopping. To include a U value, use ’0.0 0.0 0.0’ for<br>
the neighboring site inputs and insert a value in the last (eighth) entry.<br>
<br>
<br>
In our case we have 2 hoppings within unit cell. (first and third line) and two hoppings between different unit cells (second and forth lines)<br>
<br>
part of my template to create simulations with different parameters looks like:<br>
<pre><code>#HAMILT <br>
0 1 0.5 0.0 0.0 {t} {t} 0.0<br>
0 1 -0.5 0.0 0.0 {t} {t} 0.0<br>
0 2 0.0 0.5 0.0 {t} {t} 0.0<br>
0 2 0.0 -0.5 0.0 {t} {t} 0.0<br>
0 0 0.0 0.0 0.0 0.0 0.0 {u}<br>
1 1 0.0 0.0 0.0 0.0 0.0 {u}<br>
2 2 0.0 0.0 0.0 0.0 0.0 {u}<br>
</code></pre>

<h2>SYMM</h2>
<pre><code>#SYMM<br>
d  0.0d0 0.0d0 0.0d0 1.0d0 0.0d0 0.0d0<br>
d  0.0d0 0.0d0 0.0d0 0.0d0 1.0d0 0.0d0<br>
c4 0.0d0 0.0d0 0.0d0 0.0d0 0.0d0 1.0d0<br>
</code></pre>
This section describes point group symmetries that lattice has. This is optional, but can increase speed of the simulation. It is not used in the simulation itself, but it can be used in calculating of the measurables.<br>
<br>
QUEST supports 3 types of symmetry definitions.<br>
<ol><li>’cn’, where ’n’ is an integer, tells QUEST the lattice is symmetric under rotations by 2π/n. The six numbers following by ’cn’ specify the three Cartesian coordinates of a point belonging to the axis, following by the axis direction in the Cartesian coordinates. In this case ’c4’ is π/2 and indicates the x and y directions are equivalent. the point belonging to the axis is the origin ’0.0 0.0 0.0’ and the axis is the z direction ’0.0 0.0 1.0’.<br>
</li><li>The symmetry ’d’ is a mirror plane symmetry. It too is followed by six numbers. The first three are the Cartesian coordinates of a point in the plane, and the final three are the Cartesian coordinates of the normal to the plane.<br>
</li><li>’I’ is used for inversion symmetry. It is followed by three numbers, the Cartesian coordinates of the inversion point. In specifying #SYMM you must list all three components of the vectors even if the lattice is two dimensional.</li></ol>


Lieb lattice has same point group symmetries as square lattice, which is D4. As I understand it is enough to define only group generators(2 for D4). You can define more, but this information will be redundant.<br>
<br>
<h2>PHASE</h2>
<pre><code>1 0<br>
0 1<br>
s0 0.0 0.0 0.0  1.0<br>
s1 0.5 0.0 0.0  -1.0<br>
s2 0.0 0.5 0.0  -1.0<br>
</code></pre>
This section is optional, but important! (Based on the great explanation of Simone Chiesa, nearly copy-paste from his words.).<br>
<ol><li>If your lattice is not bipartite. This section is not applicable to the your case.<br>
</li><li>If your lattice is bipartite it can increase speed of your simulations by around 50%.</li></ol>

Idea is that QUEST takes advantage of particle-hole symmetry. We use this to tag atoms with ±1. The logic is the following: you need to first specify a cell (let me call it phase-cell) in which you can tag atoms with ±1 and that, once replicated using translations, will cover the entire lattice.<br>
<br>
The first two lines after #PHASE define the phase-cell in the same format used for SUPER: you need to insert an NDIM x NDIM matrix of integers. The phase cell will need to be contained an integer number of times in the super-cell. The code will otherwise abort. So the smaller the phase-cell the better, in some sense.<br>
<br>
You then need to specify a number of lines equal to the number of sites inside the phase cell. In the case of the Lieb lattice unit cell and phase cell are the same. Each line will contain the position of the sites inside the phase-cell followed by the tag (+/-1). Positions are always given as 3D arrays so each line must contain 4 numbers, The format of the position is the same used in ORB i.e. the position are given in Cartesian coordinates.<br>
<br>
<h2>BONDS</h2>
<pre><code>0 0  0.0  0.0  0.0  # 1<br>
1 1  0.0  0.0  0.0  # 2<br>
2 2  0.0  0.0  0.0  # 3<br>
0 1  0.5  0.0  0.0  # 4<br>
0 1  -0.5  0.0  0.0 # 5<br>
0 2  0.0  0.5  0.0  # 6<br>
0 2  0.0  -0.5  0.0 # 7<br>
</code></pre>

The first two entries in each line specify the orbital type where the up and down electrons are created. (the counter adopts a C-like convention<br>
staring from 0). The next three entries (dx, dy, dz) are the displacement<br>
indicating the actual distance  between the orbitals on which the pairs is<br>
created.<br>
<br>
What follows the "#" sign on each line is not read by the code but it is useful<br>
for specifying the following PAIR field.<br>
<font color='red'> I am not sure at all about this section.</font>
<h2>PAIR</h2>
<pre><code>#PAIR    <br>
            1    2    3    4    5    6    7<br>
s-wave     1.0  1.0  1.0  0.0  0.0  0.0  0.0<br>
</code></pre>
This specifies how the bonds in the previous section have to be combined when computing pairing amplitude. <font color='red'> I am not sure at all about this section. And I do not know how to define other pairings right now, so I would be happy if someone will do this.</font>
<h2>END</h2>
Do not forget to add:<br>
<pre><code>#END<br>
</code></pre>
at the end of the file.