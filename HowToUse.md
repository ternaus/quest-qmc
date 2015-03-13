# Introduction #

So, you installed the package. What is next? How to start seeing physics?

Right now structure of the package in confusing. So, I will try to go step by step to explain what and how.

## Download and install ##

First and most important step. You should get the package to the computer of interest and install it. Process is described in the
[Download and Install](https://code.google.com/p/quest-qmc/wiki/Download)

# Where is main program? #

Executable that we are interested in is kept in the folder:
```
quest-qmc/EXAMPLE/geom
```

it is strange and confusing, so, hopefully we will fix it is the next release.

There is a file in this folder named:

```
ggeom
```

It is, again, strange and confusing, but such is life.

## How to start the simulation? ##
First of all you go to this
```
quest-qmc/EXAMPLE/geom
```
folder.

There are two very important files:
  * in
  * square.geom

These two files specify all input parameters that you need to do the simulation.

Right now you can just start this in terminal:

```
./ggeom in
```

It will create two files:
  * test.out
  * test.tdm.out

Both are just text files with the output of the simulation. Where
  * test.out - time independent measurements.
  * test.tdm.out - time dependent measurements.

Both are described in [Output](https://code.google.com/p/quest-qmc/wiki/Output)

# Questions that may raise during this process: #

  * Question: May I use other filename and not "in", let's say you want to use "my\_best\_simulation".
  * Answer: Yes, you can do this. You can check that if you change filename from "in" to "my\_best\_simulation" and run:
```
./ggeom my_best_simulation
```
you will get files test.out and test.tdm.out and they will be exactly the same as you got with the old one "in" file.

---

  * Question: I do not like this "test.out", I want "my\_best\_simulation\_output.out". How do I do this?
  * Answer: Edit your "in". It is just text file and it has a line
```
ofile = test
```

Change it to:
```
ofile = my_best_simulation_output
```
start
```
./ggeom in
```
and you will get desired files as an output.

---

  * Question: It is square lattice by default and I want triangular  lattice. How do I proceed?
  * Answer: First you need to create a file that describes triangular lattice. You can do it yourself following [Lieb lattice tutorial](https://code.google.com/p/quest-qmc/wiki/Lieb_lattice_tutorial), but for the lattice that you are interested in. But before doing this I would recommend to check folder
```
quest-qmc/geometries/templates
```
because most commonly used geometries such as chain, square, triangular, honeycomb can be found there.

---

  * Question: Why input parameters are specified in two files and not just one?
  * Answer: It looks natural to divide input parameters into two classes. One of them is responsible for the Hamiltonian(hoppings, couplings, lattice structure, spatial symmetries) and here we name it geometry file. Second one is related to everything else(filenames, numerical parameters, etc).

---

  * Question: There are so many parameters. What do they mean?
  * Answer: Parameters with thir description can be found here: [Input parameters](https://code.google.com/p/quest-qmc/wiki/input_parameters)


---

  * Question: To summarize. What should I do I need to do to get some physical result?
  * Answer: To start actual simulation you need to define 2 input files, one related to Hamiltonian (example square.geom) another one related to the input parameters (example in) and start
```
./ggeom <filname with input parameters> 
```
and what is important! all 3 files should be in same folder.