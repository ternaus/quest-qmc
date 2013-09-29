
#include <exception>
#include <magma.h>

class cuda_exception : public std::exception {
  char message[100];
  public:
    cuda_exception(const char *file, int line, int code);
    virtual const char* what() const throw();
};

#define CUDACHECK(a) { int code = (a); if (code) throw(cuda_exception(__FILE__, __LINE__, code)); }

void gpu_init();
void gpu_shutdown();

void gpu_diag(int n, double *A, double *D);
void gpu_normcol(int n, double *A, double *D, double *c);
void gpu_permute(int n, int *ipiv, double *A, double *Q);
void gpu_scalerowperm(int n, double *D, double *Q, int *ipiv, double *T);
void gpu_scalerow(int n, double *h, double *B, double *M);
void gpu_scalerowcol(int n, double *h, double *G);

void gpu_dgemm(const char *trans, int m, int n, int k, double alpha, double *a,
               int lda, double *b, int ldb, double beta, double *c, int ldc);
void gpu_dgeqrf(int m, int n, double *dA, int ldda, double *tau, double *dT);
void gpu_dorgqr(int m, int n, int k, double *da, int ldda, double *tau,
                double *dT, int nb);
void gpu_dgetrf(int m, int n, double *dA, int ldda, int *ipiv);
void gpu_dgetrs(const char *trans, int n, int nrhs, double *dA, int ldda,
                int *ipiv, double *dB, int lddb);

void gpu_setvector(int n, int size, void *src, void *dst);
void gpu_getvector(int n, int size, void *src, void *dst);
void gpu_setmatrix(int m, int n, int size, void *src, void *dst);
void gpu_getmatrix(int m, int n, int size, void *src, void *dst);
void gpu_copy(void *dst, void *src, int size);
void gpu_sort(int n, double *Db, int *ipiv);

#ifdef DQMC_CUDA

#define DGEMM gpu_dgemm
#define SCALEROW gpu_scalerow
#define SCALEROWCOL gpu_scalerowcol
#define COPY(d, o, s) gpu_copy(d, o, s)

#endif
