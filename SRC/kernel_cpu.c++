// Implementation of the different kernels in the CPU using OpenMP

#include "profile.h"
#include "kernel_cpu.h"

void cpu_diag(int n, double *A, double *D)
{
  int i;
  for (i = 0; i < n; i++) 
    D[i] = A[i * n + i]; // D = diag(A)
}

void cpu_normcol(int n, double *A, double *D, double *c)
{
  int i, j;
  double dot, tmp;
  PROFILE_BEGIN();
  #ifdef _OPENMP
  #pragma omp parallel for shared(n, A, D, c) private(i, j, tmp, dot) schedule(static)
  #endif
  for (j = 0; j < n; j++) {
    dot = 0.0;
    for (i = 0; i < n; i++) {
      tmp = A[j * n + i] * D[j];
      A[j * n + i] = tmp;
      dot += tmp * tmp;
    }
    c[j] = dot;
  }
  PROFILE_END(profile_normcol, 3.0 * n * n);
}

void cpu_permute(int n, int *ipiv, double *A, double *Q)
{
  int i, j, p;
  PROFILE_BEGIN();
  #ifdef _OPENMP
  #pragma omp parallel for shared(n, ipiv, Q, A) private(i, j, p) schedule(static)
  #endif
  for (j = 0; j < n; j++) {
    p = ipiv[j];
    for (i = 0; i < n; i++)
      Q[j * n + i] = A[p * n + i];
  }
  PROFILE_END(profile_permute, 0.0);
}

void cpu_scalerowperm(int n, double *D, double *Q, int *ipiv, double *T)
{
  int i, j, p;
  PROFILE_BEGIN();
  #ifdef _OPENMP
  #pragma omp parallel for shared(n, ipiv, T, Q, D) private(i, j, p) schedule(static)
  #endif
  for (j = 0; j < n; j++) {
    p = ipiv[j];
    for (i = 0; i <= j; i++)
      T[p * n + i] = Q[j * n + i] / D[i]; // T = D^-1*R*P
    for (; i < n; i++)
      T[p * n + i] = 0;
  }
  PROFILE_END(profile_scalerowperm, 0.5 * n * n);
}

void cpu_scalerow(int n, double *h, double *B, double *M)
{
  int i, j;
  PROFILE_BEGIN();
  #ifdef _OPENMP
  #pragma omp parallel for shared(n, M, B, h) private(i, j) schedule(static)
  #endif
  for (j = 0; j < n; j++)
    for (i = 0; i < n; i++)
      M[j * n + i] = B[j * n + i] * h[i]; 
  PROFILE_END(profile_scalerow, n * n);
}

void cpu_scalerowcol(int n, double *h, double *G)
{
  int i, j;
  PROFILE_BEGIN();
  #ifdef _OPENMP
  #pragma omp parallel for shared(n, G, h) private(i, j) schedule(static)
  #endif
  for (j = 0; j < n; j++) 
    for (i = 0; i < n; i++)
      G[j * n + i] = h[i] * G[j * n + i] / h[j];
  PROFILE_END(profile_scalerowcol, 2.0 * n * n);
}

void cpu_scalerowadd(int n,  double *Db, double *U, double *D, double *T, double *G)
{
  int i, j;
  double tmp;
  PROFILE_BEGIN();
  #ifdef _OPENMP
  #pragma omp parallel for shared(n, Db, U, D, T, G) private(i, j, tmp) schedule(static)
  #endif
  for (j = 0; j < n; j++)
    for (i = 0; i < n; i++) {
      tmp = U[i * n + j] / Db[i];
      G[j * n + i] = tmp;
      T[j * n + i] = tmp + D[i] * T[j * n + i];
    }
  PROFILE_END(profile_scalerowadd, 3.0 * n * n);
}

void cpu_sort(int n, double *Db, int *ipiv)
{
  int i, j, p, t;
  double tmp, max;
  PROFILE_BEGIN();

  for (j = 0; j < n; j++) {
    ipiv[j] = j;
  }
  for (j = 0; j < n - 1; j++) {
    // find column with maximum norm
    max = Db[j];
    p = j;
    for (i = j + 1; i < n; i++)
      if (Db[i] > max) {
      	max = Db[i];
      	p = i;
      }
    // swap columns
    if (p != j) {
      t = ipiv[j]; 
      ipiv[j] = ipiv[p]; 
      ipiv[p] = t;
      tmp = Db[j]; 
      Db[j] = Db[p]; 
      Db[p] = tmp;
    }
  }
  PROFILE_END(profile_sort, 0);
}
