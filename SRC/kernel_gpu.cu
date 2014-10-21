// Implementation of the different kernels in the GPU using CUDA

#include "blaslapack.h"
#include "profile.h"
#include "kernel_gpu.h"

#include <cstdio>
#include <cublas_v2.h>

// Thread block size for all kernels
#define K 64

#ifdef DQMC_PROFILE

#define GPU_PROFILE_BEGIN() if (profile_enabled) cudaDeviceSynchronize(); \
                            PROFILE_BEGIN()
#define GPU_PROFILE_END(i, n) if (profile_enabled) cudaDeviceSynchronize(); \
                              PROFILE_END(i, n)
#else

#define GPU_PROFILE_BEGIN() 
#define GPU_PROFILE_END(i, n)

#endif

cuda_exception::cuda_exception(const char *file, int line, int code) { 
  snprintf(message, sizeof(message), "CUDA error #%i at %s:%i\n",
           code, file, line);
}

const char* cuda_exception::what() const throw() { return message; }

// static stuff for the CUBLAS library

static cublasHandle_t handle = NULL;

void gpu_diag(int n, double *A, double *D)
{
  GPU_PROFILE_BEGIN();
  CUDACHECK(cublasDcopy(handle, n, A, n + 1, D, 1));
  GPU_PROFILE_END(profile_dcopy, 16.0 * n);
}

__global__ void kernel_normcol(int n, double *A, double *D, double *c)
{
  int i, tid = threadIdx.x, j = blockIdx.x * n; // working column 
  __shared__ volatile double shared[K];
  double t, w = 0.0, d = D[blockIdx.x];
  double *p = A + j + tid;
  for (i = tid; i < n; i += K) { 
    // t = A[j + i] * d;
    t = *p * d;
    // A[j + i] = t;
    *p = t;
    w += t * t;
    p += K;
  }
  shared[tid] = w;
  __syncthreads();
  if (K >= 512) { if (tid < 256) { shared[tid] = w = w + shared[tid + 256]; } __syncthreads(); } 
  if (K >= 256) { if (tid < 128) { shared[tid] = w = w + shared[tid + 128]; } __syncthreads(); } 
  if (K >= 128) { if (tid <  64) { shared[tid] = w = w + shared[tid +  64]; } __syncthreads(); }
  if (tid < 32) { 
    if (K >=  64) shared[tid] = w = w + shared[tid + 32]; 
    if (K >=  32) shared[tid] = w = w + shared[tid + 16]; 
    if (K >=  16) shared[tid] = w = w + shared[tid +  8]; 
    if (K >=   8) shared[tid] = w = w + shared[tid +  4]; 
    if (K >=   4) shared[tid] = w = w + shared[tid +  2]; 
    if (K >=   2) shared[tid] = w = w + shared[tid +  1]; 
  }   
  if (tid == 0) 
    c[blockIdx.x] = w;
}

void gpu_normcol(int n, double *A, double *D, double *c)
{
  GPU_PROFILE_BEGIN();
  kernel_normcol <<< n , K >>> (n, A, D, c);
  CUDACHECK(cudaGetLastError());
  GPU_PROFILE_END(profile_normcol, 2.0 * n * n);
}

__global__ void kernel_permute(int n, int *ipiv, double *A, double *Q)
{
  int i, j = blockIdx.x * n; // working column
  int p = ipiv[blockIdx.x] * n;
  double *pQ = Q + j + threadIdx.x;
  double *pA = A + p + threadIdx.x;
  for (i = threadIdx.x; i < n; i += K) {
    // Q[j + i] = A[p + i];
    *pQ = *pA;
    pQ += K; pA += K;
  } 
}

void gpu_permute(int n, int *ipiv, double *A, double *Q)
{
  GPU_PROFILE_BEGIN();
  kernel_permute <<< n , K >>> (n, ipiv, A, Q);
  CUDACHECK(cudaGetLastError());
  GPU_PROFILE_END(profile_permute, 16.0 * n * n);
}

__global__ void kernel_scalerowperm(int n, double *D, double *Q, int *ipiv, double *T)
{
  int i, j = blockIdx.x; // working column
  int p = ipiv[j] * n;
  double *pT = T + p + threadIdx.x;
  double *pQ = Q + j * n + threadIdx.x;
  double *pD = D + threadIdx.x;
  for (i = threadIdx.x; i <= j; i += K) {
    // T[p + i] = Q[j * n + i] / D[i]; // T = D^-1*R*P
    *pT = *pQ / *pD;
    pT += K; pQ += K; pD += K;
  }
  for (; i < n; i += K) {
    // T[p + i] = 0;
    *pT = 0.0;
    pT += K;
  }
}

void gpu_scalerowperm(int n, double *D, double *Q, int *ipiv, double *T)
{
  GPU_PROFILE_BEGIN();
  kernel_scalerowperm <<< n , K >>> (n, D, Q, ipiv, T);
  CUDACHECK(cudaGetLastError());
  GPU_PROFILE_END(profile_scalerowperm, 0.5 * n * n);
}

__global__ void kernel_scalerow(int n, double *h, double *B, double *M)
{
  int i, j = blockIdx.x * n; // working column
  double *pM = M + j + threadIdx.x;
  double *ph = h + threadIdx.x;
  double *pB = B + j + threadIdx.x;
  for (i = threadIdx.x; i < n; i += K) {
    // M[j + i] = h[i] * B[j + i];
    *pM = *ph * *pB;
    pM += K; ph += K; pB += K;
  }
}

void gpu_scalerow(int n, double *h, double *B, double *M)
{
  GPU_PROFILE_BEGIN();
  kernel_scalerow <<< n , K >>> (n, h, B, M);
  CUDACHECK(cudaGetLastError());
  GPU_PROFILE_END(profile_scalerow, 1.0 * n * n);
}

__global__ void kernel_scalerowcol(int n, double *h, double *G)
{
  // G = diag(h) * G * diag(h)^-1
  int i, j = blockIdx.x * n; // working column
  double t, f = 1.0 / h[blockIdx.x];
  double *pG = G + j + threadIdx.x;
  double *ph = h + threadIdx.x;
  for (i = threadIdx.x; i < n; i += K) {
    // G[j + i] = h[i] * G[j + i] / f;
    t = *ph * *pG;
    *pG = t * f;
    pG += K; ph += K;
  }
}

void gpu_scalerowcol(int n, double *h, double *G)
{
  GPU_PROFILE_BEGIN();
  kernel_scalerowcol <<< n , K >>> (n, h, G);
  CUDACHECK(cudaGetLastError());
  GPU_PROFILE_END(profile_scalerowcol, 2.0 * n * n);
}

void gpu_dgemm(const char *trans, int m, int n, int k, double alpha, double *a,
               int lda, double *b, int ldb, double beta, double *c, int ldc)
{
  GPU_PROFILE_BEGIN();
  cublasOperation_t transa = trans[0] == 'N' ? CUBLAS_OP_N : CUBLAS_OP_T;
  cublasOperation_t transb = trans[1] == 'N' ? CUBLAS_OP_N : CUBLAS_OP_T;
  CUDACHECK(cublasDgemm(handle, transa, transb, m, n, k, &alpha, a, lda, b, ldb, &beta, c, ldc));
  GPU_PROFILE_END(profile_dgemm, 2.0 * k * m * n);
}

void gpu_dgeqrf(int m, int n, double *dA, int ldda, double *tau, double *dT)
{
  int info;
  GPU_PROFILE_BEGIN();
  magma_dgeqrf3_gpu(m, n, dA, ldda, tau, dT, &info);
  if (info) throw(lapack_exception("GPU DGEQRF", info)); 
  GPU_PROFILE_END(profile_dgeqrf, 2.0*m*n*n - 2.0*n*n*n/3.0 + m*n + n*n + 
                  14.0*n/3.0);
}

void gpu_dorgqr(int m, int n, int k, double *da, int ldda, double *tau,
                double *dT, int nb)
{
  int info;
  GPU_PROFILE_BEGIN();
  magma_dorgqr_gpu(m, n, k, da, ldda, tau, dT, nb, &info);
  if (info) throw(lapack_exception("GPU DORGQR", info)); 
  GPU_PROFILE_END(profile_dorgqr, 4.0*m*n*k - 2.0*(m+n)*k*k + 4.0*k*k*k/3.0 + 
                                  3.0*n*k - m*k - k*k - 4.0*k/3.0);
}

void gpu_dgetrf(int m, int n, double *dA, int ldda, int *ipiv)
{
  int info;
  GPU_PROFILE_BEGIN();
  magma_dgetrf_gpu(m, n, dA, ldda, ipiv, &info);
  if (info) throw(lapack_exception("GPU DGETRF", info)); 
  GPU_PROFILE_END(profile_dgetrf, m*n*n - n*n*n/3.0 - n*n/2.0 + 5.0*n/6.0);
}

void gpu_dgetrs(const char *trans, int n, int nrhs, double *dA, int ldda, 
                int *ipiv, double *dB, int lddb)
{
  int info;
  GPU_PROFILE_BEGIN();
  magma_dgetrs_gpu(trans[0], n, nrhs, dA, ldda, ipiv, dB, lddb, &info);
  if (info) throw(lapack_exception("GPU DGETRS", info)); 
  GPU_PROFILE_END(profile_dgetrs, nrhs * (2.0*n*n - n));
}

void gpu_setvector(int n, int size, void *src, void *dst)
{
  GPU_PROFILE_BEGIN();
  CUDACHECK(cublasSetVector(n, size, src, 1, dst, 1)); 
  GPU_PROFILE_END(profile_transfer, 2.0 * n * size);
}

void gpu_getvector(int n, int size, void *src, void *dst)
{
  GPU_PROFILE_BEGIN();
  CUDACHECK(cublasGetVector(n, size, src, 1, dst, 1)); 
  GPU_PROFILE_END(profile_transfer, 2.0 * n * size);
}

void gpu_setmatrix(int m, int n, int size, void *src, void *dst)
{
  GPU_PROFILE_BEGIN();
  CUDACHECK(cublasSetMatrix(m, n, size, src, m, dst, m)); 
  GPU_PROFILE_END(profile_transfer, 2.0 * m * n * size);
}

void gpu_getmatrix(int m, int n, int size, void *src, void *dst)
{
  GPU_PROFILE_BEGIN();
  CUDACHECK(cublasGetMatrix(m, n, size, src, m, dst, m)); 
  GPU_PROFILE_END(profile_transfer, 2.0 * m * n * size);
}

void gpu_copy(void *dst, void *src, int size)
{
  GPU_PROFILE_BEGIN();
  CUDACHECK(cudaMemcpy(dst, src, size, cudaMemcpyDeviceToDevice))
  GPU_PROFILE_END(profile_dcopy, 2.0 * size);
}

__global__ void kernel_sort(int length, double *val, int *ind)
{
  const unsigned int idx = blockIdx.x;
  int i, tid = threadIdx.x;
  __shared__ volatile int shared[K];
  int w = 0;
  double myValue = val[idx];
  for (i = tid; i < length; i += K) { 
    if (myValue < val[i] || (val[i] == myValue && i > idx)) {
      w++;
    }
  }
  shared[tid] = w;
  __syncthreads();
  if (K >= 512) { if (tid < 256) { shared[tid] = w = w + shared[tid + 256]; } __syncthreads(); } 
  if (K >= 256) { if (tid < 128) { shared[tid] = w = w + shared[tid + 128]; } __syncthreads(); } 
  if (K >= 128) { if (tid <  64) { shared[tid] = w = w + shared[tid +  64]; } __syncthreads(); }
  if (tid < 32) { 
    if (K >=  64) shared[tid] = w = w + shared[tid + 32]; 
    if (K >=  32) shared[tid] = w = w + shared[tid + 16]; 
    if (K >=  16) shared[tid] = w = w + shared[tid +  8]; 
    if (K >=   8) shared[tid] = w = w + shared[tid +  4]; 
    if (K >=   4) shared[tid] = w = w + shared[tid +  2]; 
    if (K >=   2) shared[tid] = w = w + shared[tid +  1]; 
  }   
  if (tid == 0) 
    ind[w] = idx;
}

void gpu_sort(int n, double *Db, int *ipiv)
{
  // GPU version 
  int blocks = n / K;
  if (n % K) blocks++;
  GPU_PROFILE_BEGIN();
  kernel_sort <<< n , K >>> (n, Db, ipiv);
  CUDACHECK(cudaGetLastError());
  GPU_PROFILE_END(profile_sort, 8.0 * n * n);
}

void gpu_init()
{
  if (!handle) CUDACHECK(cublasCreate(&handle));
  CUDACHECK(cudaFuncSetCacheConfig(kernel_normcol, cudaFuncCachePreferL1));
  CUDACHECK(cudaFuncSetCacheConfig(kernel_permute, cudaFuncCachePreferL1));
  CUDACHECK(cudaFuncSetCacheConfig(kernel_scalerowperm, cudaFuncCachePreferL1));
  CUDACHECK(cudaFuncSetCacheConfig(kernel_scalerow, cudaFuncCachePreferL1));
  CUDACHECK(cudaFuncSetCacheConfig(kernel_scalerowcol, cudaFuncCachePreferL1));
  CUDACHECK(cudaFuncSetCacheConfig(kernel_sort, cudaFuncCachePreferL1));
}

void gpu_shutdown()
{
  if (handle) CUDACHECK(cublasDestroy(handle));
  handle = NULL;
}
