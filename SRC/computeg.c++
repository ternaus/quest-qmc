#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <stdexcept>

#include "blaslapack.h"
#include "profile.h"
#include "kernel_cpu.h"
#ifdef DQMC_CUDA
#include "kernel_gpu.h"
#endif

// class definition for ComputeG 
class gfun {
public:
  gfun(int n, int L, int nWrap, int fixw);
  void computeg(double *B, int nOrth, double *h, int il, 
                double *G, double *sgn, double *det);
  void stratification(double *B, int k, double *h, int il, 
                      double *tau, double *Q, double *D, double *T);
  void swapg(double *B, double *Bi, double *h, double *G);  
  void invalid_cache(int j);  
  void compute_M_cache(double *B, double *h, int &il, int k, 
                       double *M, double *work);
  ~gfun();

private:
  // parameters from the simulation
  int n;
  int L;
  int nOrth;
  int nWrap;
  int fixw;
  // memory allocation
  double *tau;
  int *ipiv;
  double *U;
  double *D;
  double *T;
  double *A;
  double *W;
  double *Db;
  int lwork;
  double *work;
#ifdef DQMC_CUDA
  int *perm;
  double *B_gpu;
  double *B_fortran;
  double *Bi_gpu;
  double *Bi_fortran;
  double *h_gpu;
  double *U_gpu;
  double *D_gpu;
  double *T_gpu;
  int nb;
  double *work_gpu;
#endif
  double **cache_M;
  bool *cache_valid;
};

// Exported functions called from Fortran
extern "C" {
  void cpp_gfun_init_(long *cpp_data, int *n, int *L, int *nWrap, int *fixw);
  void cpp_gfun_computeg_(long *cpp_data, int *il, double *sgn, double *G, double *h,
                          double *B, int *nOrth, double *det);
  void cpp_gfun_free_(long *cpp_data);
  void cpp_gfun_swapg_(long *cpp_data, double *B, double *Bi, double *h, double *G);
  void cpp_gfun_invalid_cache_(long *cpp_data, int *j);
}

void cpp_gfun_init_(long *cpp_data, int *n, int *L, int *nWrap, int *fixw)
{
  try {
    *cpp_data = (long) new gfun(*n, *L, *nWrap, *fixw);
  } catch (std::exception& e) {
    std::cerr << "Error in C++ gfun_init: " << e.what() << std::endl;
    exit(1);
  }
}

void cpp_gfun_computeg_(long *cpp_data, int *il, double *sgn, double *G, double *h,
                        double *B, int *nOrth, double *det)
{
  try {
    ((gfun *)*cpp_data)->computeg(B, *nOrth, h, *il, G, sgn, det);
  } catch (std::exception& e) {
    std::cerr << "Error in C++ gfun_computeg: " << e.what() << std::endl;
    exit(1);
  }
}

void cpp_gfun_free_(long *cpp_data)
{
  try {
    delete (gfun *)*cpp_data;
    *cpp_data = 0;
  } catch (std::exception& e) {
    std::cerr << "Error in C++ gfun_free: " << e.what() << std::endl;
    exit(1);
  }
}

void cpp_gfun_swapg_(long *cpp_data, double *B, double *Bi, double *h, double *G)
{
  try {
    ((gfun *)*cpp_data)->swapg(B, Bi, h, G);
  } catch (std::exception& e) {
    std::cerr << "Error in C++ gfun_swapg: " << e.what() << std::endl;
    exit(1);
  }
}

void cpp_gfun_invalid_cache_(long *cpp_data, int *j)
{
  try {
    ((gfun *)*cpp_data)->invalid_cache(*j);
  } catch (std::exception& e) {
    std::cerr << "Error in C++ gfun_invalid_cache: " << e.what() << std::endl;
    exit(1);
  }
}
  
// Computes log(abs(det(M))) and sign(det(M)) from a LU decomposition
void sgndet(int n, double *M, int ldm, int *ipiv, double &sgn, double &det)
{
  sgn = 1.0; det = 0.0;
  for (int i = 0; i < n; i++) {
    if (ipiv[i] != (i + 1)) sgn = -sgn;
    if (M[i * ldm + i] < 0.0) {
      sgn = -sgn;
      det += log(-M[i * ldm + i]);
    } else {
      det += log(M[i * ldm + i]); 
    }
  }  
}

// Computes M=B_{il+1}...B_1*B_L...B_{il+k}
void compute_M(int n, double *B, int L, double *h, int &il, int k, double *M, double *work)
{
  int l;
  PROFILE_BEGIN();
  il++; if (il >= L) il = 0;
  SCALEROW(n, h + il * n, B, M);
  for (l = 1; l < k; l++) {
    il++; if (il >= L) il = 0;
    DGEMM("NN", n, n, n, 1.0, B, n, M, n, 0.0, work, n);  
    SCALEROW(n, h + il * n, work, M);
  }
  PROFILE_END(profile_computem, 0);
}

void gfun::invalid_cache(int j)
{
  if (cache_valid) {
    if (j > 0) {
      if (j > L) throw(std::invalid_argument("j > n"));
      cache_valid[j - 1] = false;
    } else {
      for (int i = 0; i < L; i++)
        cache_valid[i] = false;
    }
  }
}

void gfun::compute_M_cache(double *B, double *h, int &il, int k, double *M, double *work)
{
  if (!cache_M || k == 1) {
    // Don't use the cache
    compute_M(n, B, L, h, il, k, M, work);
  } else {
    int i = il;
    bool valid = true;
    for (int l = 0; l < k && valid; l++) {
      i++; if (i >= L) i = 0;
      valid = cache_valid[i];
    }
    if (valid) {
      COPY(M, cache_M[il], n * n * sizeof(double));
      il = i;
    } else {
      i = il;
      compute_M(n, B, L, h, il, k, M, work);
      if (!cache_M[i])
#ifdef DQMC_CUDA
        CUDACHECK(cudaMalloc((void**)(cache_M + i), n * n * sizeof(double))); 
#else
        cache_M[i] = new double[n * n];
#endif
      COPY(cache_M[i], M, n * n * sizeof(double));
      for (int l = 0; l < k; l++) {
        i++; if (i >= L) i = 0;
        cache_valid[i] = true;
      }
    }
  }
}

gfun::gfun(int n, int L, int nWrap, int fixw)
{
  this->n = n;
  this->L = L;
  this->nOrth = -1; // mark as uninitialized
  this->nWrap = nWrap;
  this->fixw = fixw;
  
#ifdef DQMC_CUDA
  // Start CUBLAS library
  gpu_init();
  // Allocate buffers in the graphics memory
  CUDACHECK(cudaMalloc((void**)&B_gpu, n * n * sizeof(double)));
  CUDACHECK(cudaMalloc((void**)&Bi_gpu, n * n * sizeof(double)));
  B_fortran = Bi_fortran = NULL;
  CUDACHECK(cudaMalloc((void**)&h_gpu, n * L * sizeof(double)));
  CUDACHECK(cudaMalloc((void**)&U_gpu, n * n * sizeof(double)));
  CUDACHECK(cudaMalloc((void**)&D_gpu, n * sizeof(double)));
  CUDACHECK(cudaMalloc((void**)&T_gpu, n * n * sizeof(double)));
  CUDACHECK(cudaMalloc((void**)&A, n * n * sizeof(double)));
  CUDACHECK(cudaMalloc((void**)&W, n * n * sizeof(double)));
  CUDACHECK(cudaMalloc((void**)&perm, n * sizeof(int)));
#else
  // Allocate memory 
  A = new double[n * n];
  W = new double[n * n];
#endif

  U = new double[n * n];
  D = new double[n];
  T = new double[n * n];

  tau = new double[n];
  ipiv = new int[n];
  Db = new double[n]; 

  cache_M = NULL;
  cache_valid = NULL;

  // Allocate workspace for LAPACK
  double temp;
  lapack_dgeqp3(n, n, U, n, ipiv, tau, &temp, -1);
  lwork = (int)temp;
  lapack_dgeqrf(n, n, U, n, tau, &temp, -1);
  if (temp > lwork) lwork = (int)temp;
  lapack_dormqr("RN", n, n, n, U, n, tau, A, n, &temp, -1);  
  if (temp > lwork) lwork = (int)temp;
  lapack_dorgqr(n, n, n, U, n, tau, &temp, -1);
  if (temp > lwork) lwork = (int)temp;
  work = new double[lwork];

#ifdef DQMC_CUDA
  // Allocate workspace for MAGMA
  nb = magma_get_dgeqrf_nb(n);
  CUDACHECK(cudaMalloc((void**)&work_gpu, 
            (2*n+(n+31)/32*32)*nb * sizeof(double)));
#endif
}

// Stratification loop
void gfun::stratification(double *B, int k, double *h, int il, 
                          double *tau, double *Q, double *D, double *T)
{
  int i, l;

#ifdef DQMC_CUDA

  compute_M_cache(B, h, il, k, Q, W); // Q = B_1
  gpu_getmatrix(n, n, sizeof(double), Q, U); 
  for (i = 0; i < n; i++) ipiv[i] = 0;
  lapack_dgeqp3(n, n, U, n, ipiv, tau, work, lwork); // QRP
  for (i = 0; i < n; i++) ipiv[i]--;
  gpu_setmatrix(n, n, sizeof(double), U, Q); 
  gpu_setvector(n, sizeof(int), ipiv, perm); 

  gpu_diag(n, Q, D); // D = diag(R)

  gpu_scalerowperm(n, D, Q, perm, T); // T = D^-1*R*P

  lapack_dorgqr(n, n, n, U, n, tau, work, lwork); // build Q
  gpu_setmatrix(n, n, sizeof(double), U, Q); 

  for (l = 1; l < L / k; l++) {
    // Q = (B_l*Q)*D
    compute_M_cache(B, h, il, k, W, A); // W = B_l
    gpu_dgemm("NN", n, n, n, 1.0, W, n, Q, n, 0.0, A, n); // A = W * Q
        
    // A = A * D and compute the norm of each column
    gpu_normcol(n, A, D, work_gpu);

    // compute a permutation P that sorts the columns of A
    gpu_sort(n, work_gpu, perm);
    
    // apply the permutation
    gpu_permute(n, perm, A, Q);
        
    // standard QR
    gpu_dgeqrf(n, n, Q, n, tau, work_gpu);
    
    // built R
    COPY(W, Q, n * n * sizeof(double));
    magmablas_dswapdblk(n, nb, W, n, 1, work_gpu + n * nb, nb, 0);
    
    gpu_diag(n, W, D); // D = diag(R)

    // T = D^-1 * R * P * T
    gpu_scalerowperm(n, D, W, perm, A); // A = D^-1*R*P

    gpu_dorgqr(n, n, n, Q, n, tau, work_gpu, nb); // build Q

    DGEMM("NN", n, n, n, 1.0, A, n, T, n, 0.0, W, n); // W = A * T
    COPY(T, W, n * n * sizeof(double)); // T = W
  }

#else

  compute_M_cache(B, h, il, k, Q, W); // Q = B_1
  for (i = 0; i < n; i++) ipiv[i] = 0;
  lapack_dgeqp3(n, n, Q, n, ipiv, tau, work, lwork); // QRP
  for (i = 0; i < n; i++) ipiv[i]--;

  cpu_diag(n, Q, D); // D = diag(R)

  cpu_scalerowperm(n, D, Q, ipiv, T); // T = D^-1*R*P

  for (l = 1; l < L / k; l++) {
    // Q = (B_l*Q)*D
    compute_M_cache(B, h, il, k, A, W); // A = B_l
    lapack_dormqr("RN", n, n, n, Q, n, tau, A, n, work, lwork); // A = A * Q
        
#if 0
    // Original stratification method with QRP
    int j;
    for (j = 0; j < n; j++)
      for (i = 0; i < n; i++)
	Q[j * n + i] = A[j * n + i] * D[j];
    for (i = 0; i < n; i++) ipiv[i] = 0;
    lapack_dgeqp3(n, n, Q, n, ipiv, tau, work, lwork); // QRP
    for (i = 0; i < n; i++) ipiv[i]--;
#else
    // A = A * D and compute the norm of each column
    cpu_normcol(n, A, D, Db);

    // compute a permutation P that sorts the columns of A
    cpu_sort(n, Db, ipiv);
   
    // apply the permutation
    cpu_permute(n, ipiv, A, Q);
        
    // standard QR
    lapack_dgeqrf(n, n, Q, n, tau, work, lwork);
#endif
    cpu_diag(n, Q, D); // D = diag(R)
    
    // T = D^-1 * R * P * T
    cpu_scalerowperm(n, D, Q, ipiv, A); // A = D^-1*R*P
    DGEMM("NN", n, n, n, 1.0, A, n, T, n, 0.0, W, n); // W = A * T
    COPY(T, W, n * n * sizeof(double)); // T = W

  }
  lapack_dorgqr(n, n, n, Q, n, tau, work, lwork); // build Q

#endif
}


void gfun::computeg(double *B, int k, double *h, int il,
                    double *G, double *psgn, double *pdet)
{
  if (k != 1) {
    if (nOrth == -1) {
      // Initialization of nOrth
      nOrth = k;
      if (L % nOrth) throw(std::invalid_argument("L must be a multiple of nOrth for prepivoting"));
      if (nOrth == nWrap && fixw) {
        // Allocate space for block cache
	cache_M = new double*[L];
	cache_valid = new bool[L];
	for (int i = 0; i < L; i++) {
	  cache_M[i] = NULL;
	  cache_valid[i] = false;
	}
      }
    } else if (k != nOrth) throw(std::invalid_argument("nOrth cannot change during simulation"));
  }

  // Compute G using the ASQRD algorithm
  PROFILE_ENABLE();
  PROFILE_BEGIN();

  int i;
  
  // Build Q * D * T = B_L * B_L-1 * .. * B_1
#ifdef DQMC_CUDA
  if (B != B_fortran) {
    // Do not copy B if its pointer is the same as previous calls
    gpu_setmatrix(n, n, sizeof(double), B, B_gpu);
    B_fortran = B;
  }
  gpu_setmatrix(n, L, sizeof(double), h, h_gpu);
  stratification(B_gpu, k, h_gpu, il - 1, tau, 
                 U_gpu, D_gpu, T_gpu);
  gpu_getmatrix(n, n, sizeof(double), U_gpu, U);
  gpu_getvector(n, sizeof(double), D_gpu, D);
  gpu_getmatrix(n, n, sizeof(double), T_gpu, T);
#else
  stratification(B, k, h, il - 1, tau, U, D, T);
#endif

  // compute G = (Db^-1 * Q' + Ds * T)^-1 * Db^-1 * Q'
  
  // split D and compute det(Db)
  double sgn2 = 1.0, det2 = 0.0;
  for (i = 0; i < n; i++)
    if (fabs(D[i]) > 1) {
      Db[i] = D[i]; D[i] = 1.0;
      if (Db[i] < 0.0) {
        sgn2 = -sgn2;
        det2 += log(-Db[i]);
      } else {
        det2 += log(Db[i]); 
      }
    } else {
      Db[i] = 1.0;
    }
   
  // G = Db * Q' ; T = Db * Q' + Ds * T
  cpu_scalerowadd(n, Db, U, D, T, G);
  
  // G = T^-1 * G
  lapack_dgetrf(n, n, T, n, ipiv);
  double sgn1, det1; sgndet(n, T, n, ipiv, sgn1, det1);
  lapack_dgetrs("N", n, n, T, n, ipiv, G, n);
  
  // compute det(Q)
  double sgn3 = 1.0;
  for (i = 0; i < n - 1; i++) 
    if (tau[i] != 0) sgn3 = -sgn3;

  *psgn = sgn1 * sgn2 * sgn3; *pdet = - det1 - det2;
 
  PROFILE_END(profile_computeg, 0);
  PROFILE_DISABLE();
}

gfun::~gfun() {
  // Free allocated memory
#ifdef DQMC_CUDA  
  CUDACHECK(cudaFree(B_gpu));
  CUDACHECK(cudaFree(Bi_gpu));
  CUDACHECK(cudaFree(h_gpu));
  CUDACHECK(cudaFree(U_gpu));
  CUDACHECK(cudaFree(D_gpu));
  CUDACHECK(cudaFree(T_gpu));
  CUDACHECK(cudaFree(A));
  CUDACHECK(cudaFree(W));
  CUDACHECK(cudaFree(work_gpu));
  CUDACHECK(cudaFree(perm));
  gpu_shutdown();
#else
  delete []A;
  delete []W;
#endif
  delete []U;
  delete []D;
  delete []T;

  delete []tau;
  delete []ipiv;
  delete []Db;
  if (cache_M) {
    for (int i = 0; i < L; i++)
      if (cache_M[i])
#ifdef DQMC_CUDA  
        CUDACHECK(cudaFree(cache_M[i]));
#else
        delete []cache_M[i];
#endif
    delete []cache_M;
    delete []cache_valid;
  }
  delete []work;
}

void gfun::swapg(double *B, double *Bi, double *h, double *G) 
{
  // Computes G = h * B * G * Bi / h
  PROFILE_ENABLE();
  PROFILE_BEGIN();
#ifdef DQMC_CUDA
  /* Do not copy B or Bi if their pointers are the same as previous calls */
  if (B != B_fortran) {
    gpu_setmatrix(n, n, sizeof(double), B, B_gpu);
    B_fortran = B;
  }
  if (Bi != Bi_fortran) {
    gpu_setmatrix(n, n, sizeof(double), Bi, Bi_gpu);
    Bi_fortran = Bi;
  }
  gpu_setmatrix(n, n, sizeof(double), G, T_gpu);
  gpu_setvector(n, sizeof(double), h, h_gpu);
  // W = B * G
  DGEMM("NN", n, n, n, 1.0, B_gpu, n, T_gpu, n, 0.0, W, n);
  // G = W * Bi
  DGEMM("NN", n, n, n, 1.0, W, n, Bi_gpu, n, 0.0, T_gpu, n);
  // G = h * G / h
  SCALEROWCOL(n, h_gpu, T_gpu);
  gpu_getmatrix(n, n, sizeof(double), T_gpu, G);
#else
  // W = B * G
  DGEMM("NN", n, n, n, 1.0, B, n, G, n, 0.0, T, n);
  // G = W * Bi
  DGEMM("NN", n, n, n, 1.0, T, n, Bi, n, 0.0, G, n);
  // G = h * G / h
  SCALEROWCOL(n, h, G);
#endif
  PROFILE_END(profile_swapg, 0);
  PROFILE_DISABLE();
}
