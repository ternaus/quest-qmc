#include <cstdio>

#include "blaslapack.h"
#include "profile.h"

lapack_exception::lapack_exception(const char* function, int info) { 
  snprintf(message,  sizeof(message), "LAPACK error in %s info=%i", 
  function, info); 
}

const char* lapack_exception::what() const throw() { return message; }

/* Actual interface to Fortran routines */

extern "C" {
  void dcopy_(int *n, double *dx, int *incx, double *dy, int *incy);
  void daxpy_(int *n, double *da, double *dx, int *incx, double *dy, int *incy);
  void dscal_(int *n, double *da, double *dx, int *incx);
  double ddot_(int *n, double *dx, int *incx, double *dy, int *incy);
  void dgemv_(const char *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);
  void dgemm_(const char *transa, const char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
  void dtrmm_(const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb);
  void dtrsm_(const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb);

  void dsyev_(const char *jobz, const char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info);
  void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
  void dgeqp3_(int *m, int *n, double *a, int *lda, int *jpvt, double *tau, double *work, int *lwork, int *info);
  void dorgqr_(int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *lwork, int *info);
  void dormqr_(const char *side, const char *trans, int *m, int *n, int *k, double *a, int *lda, double *tau, double *c, int *ldc, double *work, int *lwork, int *info);
  void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
  void dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
  void dgetrs_(const char *trans, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
  void dgejsv_(const char *joba, const char *jobu, const char *jobv, const char *jobr, const char *jobt, const char *jobp, int *m, int *n, double *a, int *lda, double *sva, double *u, int *ldu, double *v, int *ldv, double *work, int *lwork, int *iwork, int *info);
  void dgerfsx_(const char *trans, const char *equed, int *n, int *nrhs, double *a, int *lda, double *af, int *ldaf, int *ipiv, double *r, double *c, double *b, int *ldb, double *x, int *ldx, double *rcond, double *berr, int *n_err_bnds, double *err_bnds_norm, double *err_bnds_comp, int *nparams, double *params, double *work, int *iwork, int *info); 

  void dgeqrf_(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info);
  void dlarft_(const char *direct, const char *storev, int *n, int *k, double *v, int *ldv, double *tau, double *t, int *ldt);
}

#define N ((double)*n)
#define N2 ((double)*n * (double)*n)
#define N3 ((double)*n * (double)*n * (double)*n)
#define M ((double)*m)
#define M2 ((double)*m * (double)*m)
#define K ((double)*k)
#define K2 ((double)*k * (double)*k)
#define K3 ((double)*k * (double)*k * (double)*k)

/* routine stubs */

void blas_dcopy_(int *n, double *sx, int *incx, double *sy, int *incy) {
  PROFILE_BEGIN();
  dcopy_(n, sx, incx, sy, incy);
  PROFILE_END(profile_dcopy, 0);
}

void blas_daxpy_(int *n, double *sa, double *sx, int *incx, double *sy, int *incy) {
  PROFILE_BEGIN();
  daxpy_(n, sa, sx, incx, sy, incy);
  PROFILE_END(profile_daxpy, 2*N);
}

void blas_dscal_(int *n, double *sa, double *sx, int *incx) {
  PROFILE_BEGIN();
  dscal_(n, sa, sx, incx); 
  PROFILE_END(profile_dscal, N);
}

double blas_ddot_(int *n, double *sx, int *incx, double *sy, int *incy) {
  PROFILE_BEGIN();
  double dot = ddot_(n, sx, incx, sy, incy);
  PROFILE_END(profile_ddot, 2*N);
  return dot;
}

void blas_dgemv_(const char *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy) {
  PROFILE_BEGIN();
  dgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy); 
  PROFILE_END(profile_dgemv, 2*M*N);
}

void blas_dgemm_(const char *transa,const char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc) {
  PROFILE_BEGIN();
  dgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc); 
  PROFILE_END(profile_dgemm, 2*M*N*K);
}

void blas_dtrmm_(const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb) {
  PROFILE_BEGIN();
  dtrmm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
  if (*side == 'L' || *side == 'l') {
    PROFILE_END(profile_dtrmm, N*M2);
  } else {
    PROFILE_END(profile_dtrmm, M*N2);
  } 
}

void blas_dtrsm_(const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb) {
  PROFILE_BEGIN();
  dtrsm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
  if (*side == 'L' || *side == 'l') {
    PROFILE_END(profile_dtrsm, N*M2);
  } else {
    PROFILE_END(profile_dtrsm, M*N2);
  } 
}

// Workaround for NERSC carver buggy MKL 
/* extern "C" {
void MKL_Set_Num_Threads(int nth);
int  MKL_Get_Max_Threads(void);
} */

void lapack_dsyev_(const char *jobz, const char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info) {
  PROFILE_LWORK_BEGIN();
  // int nth = MKL_Get_Max_Threads();
  // MKL_Set_Num_Threads(1);
  dsyev_(jobz, uplo, n, a, lda, w, work, lwork, info);
  // MKL_Set_Num_Threads(nth);
  PROFILE_LWORK_END(profile_dsyev, 0);
}

void lapack_dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info) {
  PROFILE_BEGIN();
  dgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
  PROFILE_END(profile_dgesv, N * N2 - N3 / 3 - N2 / 2 + 5 * N / 6 +
                       (double)*nrhs * (2 * N2 - N));
}

void lapack_dgeqp3_(int *m, int *n, double *a, int *lda, int *jpvt, double *tau, double *work, int *lwork, int *info) {
  PROFILE_LWORK_BEGIN();
  dgeqp3_(m, n, a, lda, jpvt, tau, work, lwork, info);
  PROFILE_LWORK_END(profile_dgeqp3, 2*M*N2 - 2*N3/3 + M*N + N2 + 14*N/3);
}

void lapack_dorgqr_(int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *lwork, int *info) {
  PROFILE_LWORK_BEGIN();
  dorgqr_(m, n, k, a, lda, tau, work, lwork, info);
  PROFILE_LWORK_END(profile_dorgqr, 4*M*N*K - 2*(M+N)*K2 + 4*K3/3 + 3*N*K - M*K - K2 - 4*K/3);
}

void lapack_dormqr_(const char *side, const char *trans, int *m, int *n, int *k, double *a, int *lda, double *tau, double *c, int *ldc, double *work, int *lwork, int *info) {
  PROFILE_LWORK_BEGIN();
  dormqr_(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
  PROFILE_LWORK_END(profile_dormqr, 4*M*N*K - 2*N*K2 + 3*N*K);
}

void lapack_dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info) {
  PROFILE_BEGIN();
  dgetrf_(m, n, a, lda, ipiv, info);
  PROFILE_END(profile_dgetrf, M*N2 - N3/3 - N2/2 + 5*N/6);
}

void lapack_dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info) {
  PROFILE_LWORK_BEGIN();
  dgetri_(n, a, lda, ipiv, work, lwork, info);
  PROFILE_LWORK_END(profile_dgetri, 4*N3 - N2 + 5*N/3);
}

void lapack_dgetrs_(const char *trans, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info) {
  PROFILE_BEGIN();
  dgetrs_(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
  PROFILE_END(profile_dgetrs, (double)*nrhs * (2*N2 - N));
}

void lapack_dgejsv_(const char *joba, const char *jobu, const char *jobv, const char *jobr, const char *jobt, const char *jobp, int *m, int *n, double *a, int *lda, double *sva, double *u, int *ldu, double *v, int *ldv, double *work, int *lwork, int *iwork, int *info) {
  PROFILE_LWORK_BEGIN();
  dgejsv_(joba, jobu, jobv, jobr, jobt, jobp, m, n, a, lda, sva, u, ldu, v, ldv, work, lwork, iwork, info);
  PROFILE_LWORK_END(profile_dgejsv, 0);
}

#ifdef _SXX
void lapack_dgerfsx_(const char *trans, const char *equed, int *n, int *nrhs, double *a, int *lda, double *af, int *ldaf, int *ipiv, double *r, double *c, double *b, int *ldb, double *x, int *ldx, double *rcond, double *berr, int *n_err_bnds, double *err_bnds_norm, double *err_bnds_comp, int *nparams, double *params, double *work, int *iwork, int *info) {
  PROFILE_BEGIN();
  dgerfsx_(trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
  PROFILE_END(profile_dgerfsx, 0);
}
#endif

void lapack_dgeqrf_(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info) {
  PROFILE_LWORK_BEGIN();
  dgeqrf_(m, n, a, lda, tau, work, lwork, info);
  PROFILE_LWORK_END(profile_dgeqrf, 2*M*N2 - 2*N3/3 + M*N + N2 + 14*N/3);
}

void lapack_dlarft_(const char *direct, const char *storev, int *n, int *k, double *v, int *ldv, double *tau, double *t, int *ldt) {
  PROFILE_BEGIN();
  dlarft_(direct, storev, n, k, v, ldv, tau, t, ldt);
  PROFILE_END(profile_dlarft, 0);
}
