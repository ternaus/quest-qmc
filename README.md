# QUEST (Quantum Electron Simulation Toolbox)

QUEST is a scientific computing package designed for quantum electron simulations using Determinant Quantum Monte Carlo (DQMC) methods.

## Prerequisites for MacOS

Install the required compilers using Homebrew:
```bash
brew install gcc
brew install gfortran
```

## Installation

1. Clone the repository:
```bash
git clone https://github.com/ternaus/quest-qmc.git
cd quest-qmc
```

2. Build the project:
```bash
# This will automatically:
# 1. Build OpenBLAS if needed
# 2. Build the DQMC library
# 3. Build the examples
make
```

## Project Structure

- `SRC/`: Source code for the DQMC library
- `EXAMPLE/`: Example applications
- `OpenBLAS/`: Contains OpenBLAS library when built

## Available Make Commands

```bash
make        # Build everything (default)
make help   # Show all available commands
make clean  # Clean all built files
```

## Running Examples

After building, you can run the examples from the EXAMPLE directory:
```bash
cd EXAMPLE
./example1  # Replace with actual example name
```

## Performance Notes

The package uses OpenBLAS for linear algebra operations, which automatically detects and utilizes available CPU cores for optimal performance.

## License

QUEST is licensed under the Dual License. See the LICENSE file for details.

## Citation

If you use this software in your research, please cite:
[Add citation information here]

## Troubleshooting

Common issues and solutions:

1. OpenMP issues on MacOS:
```bash
# If you see OpenMP-related errors, install:
brew install libomp
```

2. Compilation errors:
   - Make sure you have the latest versions of gcc and gfortran
   - Check that all paths in the Makefile are correct
