# QUEST Example Programs

This directory contains example programs demonstrating the usage of QUEST (QUantum Electron Simulation Toolbox) for quantum Monte Carlo simulations of the Hubbard model.

## Directory Structure

- `geom/`: Main geometry-based simulation program
- `test/`: Performance testing programs
- `verify/`: Verification programs against known results

## Compilation

To compile all examples:

```bash
make
```

This will compile the programs in each subdirectory.

## Running Simulations

### Main Simulation (geom/)

The main simulation program reads a geometry file (defining lattice structure) and an input file (defining simulation parameters).

1. Navigate to the geom directory:

```bash
cd geom
```

2. Run the simulation:

```bash
./ggeom < in
```

Example files provided:

- `square.geom`: Defines a 2D square lattice
- `in`: Contains simulation parameters

You can create your own geometry and input files following these examples.

### Input File Parameters

The input file (`in`) contains several parameter sections:

- Lattice dimension
- Hubbard model parameters (U, t, μ, etc.)
- Monte Carlo algorithm settings
- Measurement settings
- Numerical parameters

### Test Programs (test/)

Contains test cases for different lattice sizes:

- `small.in`: 8x8 lattice
- `median.in`: 10x10 lattice
- `large.in`: 16x16 lattice

Run tests with:

```bash
cd test
./test < small.in
```

### Verification (verify/)

Contains programs to verify the code against known analytical results for special cases:

- Single site (t=0)
- No Coulomb interaction (U=0)

Run verification with:

```bash
cd verify
./verify
```

## Output Files

The simulation generates several output files (based on the `ofile` parameter in your input):

- `.out`: Main results
- `.geometry`: Geometry information
- `.tdm.out`: Time-dependent measurements (if enabled)
-
