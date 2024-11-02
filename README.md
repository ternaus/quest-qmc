# QUEST (Quantum Electron Simulation Toolbox)

QUEST is a scientific computing package designed for quantum electron simulations using Determinant Quantum Monte Carlo (DQMC) methods. It supports both CPU and GPU acceleration through CUDA.

## Prerequisites

- Fortran compiler (gfortran or ifort)
- C++ compiler (g++ or icpc)
- BLAS/LAPACK libraries (OpenBLAS by default, or Intel MKL)
- Optional: CUDA toolkit for GPU acceleration
- Optional: MAGMA library for GPU-accelerated linear algebra

## Installation

1. Clone the repository:

```bash
git clone https://github.com/ternaus/quest-qmc.git

cd quest-qmc
```

2. Configure the build in the Makefile:

Key configuration options:
- `COMPILER`: Choose between `gnu` (default) or `intel` compilers
- `LAPACK`: Choose between:
  - `default` (OpenBLAS)
  - `mkl_seq` (Intel MKL sequential)
  - `mkl_par` (Intel MKL parallel)
- `MKLPATH`: Path to Intel MKL if using MKL (default: `$(MKLROOT)/lib/intel64`)
- `MAGMAPATH`: Path to MAGMA installation (optional)
- `CUDAPATH`: Path to CUDA installation (optional)

3. Build the project:

Build OpenBLAS (if using default LAPACK)

```bash
make lapack
```

```bash
make all
```


## Features

- **Multiple Backend Support**:
  - CPU-only computation
  - GPU acceleration via CUDA
  - Choice of BLAS/LAPACK implementations

- **Acceleration Features**:
  - `FLAG_CKB`: Checkboard decomposition
  - `FLAG_ASQRD`: GPU-accelerated equal-time Green's function kernel
  - `FLAG_BSOFI`: GPU-accelerated time-dependent Green's function kernel
  - `FLAG_CUDA`: Enable NVIDIA CUDA support

## Project Structure

- `SRC/`: Source code for the DQMC library
- `EXAMPLE/`: Example applications
- `OpenBLAS/`: Contains OpenBLAS library when built

## Usage

1. Configure the desired features in Makefile
2. Build the project
3. Run examples from the EXAMPLE directory

## Cleaning

To clean all built files:

```bash
make clean
```



## Important Notes

1. ASQRD method is incompatible with checkboard decomposition. Do not enable both `FLAG_ASQRD` and `FLAG_CKB` simultaneously.

2. When using CUDA support:
   - MAGMA library is required
   - Set proper paths for both CUDA and MAGMA installations

3. For Intel MKL users:
   - Ensure `MKLROOT` environment variable is set correctly
   - Choose between sequential (`mkl_seq`) or parallel (`mkl_par`) versions

## Performance Optimization

- For Intel processors, consider using the Intel compiler (`COMPILER = intel`) with MKL
- For GPU acceleration, enable CUDA support and use MAGMA
- Use parallel MKL for multi-threaded performance on CPU

## License

QUEST is licensed under the Dual License. See the LICENSE file for details.

## Citation

If you use this software in your research, please cite:
[Add citation information here]
