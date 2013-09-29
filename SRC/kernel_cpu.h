
void cpu_diag(int n, double *A, double *D);
void cpu_normcol(int n, double *A, double *D, double *c);
void cpu_permute(int n, int *ipiv, double *A, double *Q);
void cpu_scalerowperm(int n, double *D, double *Q, int *ipiv, double *T);
void cpu_scalerow(int n, double *h, double *B, double *M);
void cpu_scalerowcol(int n, double *h, double *G);
void cpu_scalerowadd(int n, double *Db, double *U, double *D, double *T, double *G);
void cpu_sort(int n, double *Db, int *ipiv);

#ifndef DQMC_CUDA

#define DGEMM blas_dgemm
#define SCALEROW cpu_scalerow
#define SCALEROWCOL cpu_scalerowcol
#define COPY(d, o, s) memcpy(d, o, s)

#endif
