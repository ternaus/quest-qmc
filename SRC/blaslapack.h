#include <exception>

/* Interface to Fortran double precision BLAS and LAPACK with profiling */

extern "C" {
  void blas_dcopy_(int *n, double *sx, int *incx, double *sy, int *incy);
  void blas_daxpy_(int *n, double *sa, double *sx, int *incx, double *sy, int *incy);
  void blas_dscal_(int *n, double *sa, double *sx, int *incx);
  double blas_ddot_(int *n, double *sx, int *incx, double *sy, int *incy);
  void blas_dgemv_(const char *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);
  void blas_dgemm_(const char *transa, const char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
  void blas_dtrmm_(const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb);
  void blas_dtrsm_(const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb);

  void lapack_dsyev_(const char *jobz, const char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info);
  void lapack_dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
  void lapack_dgeqp3_(int *m, int *n, double *a, int *lda, int *jpvt, double *tau, double *work, int *lwork, int *info);
  void lapack_dorgqr_(int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *lwork, int *info);
  void lapack_dormqr_(const char *side, const char *trans, int *m, int *n, int *k, double *a, int *lda, double *tau, double *c, int *ldc, double *work, int *lwork, int *info);
  void lapack_dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
  void lapack_dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
  void lapack_dgetrs_(const char *trans, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
  void lapack_dgejsv_(const char *joba, const char *jobu, const char *jobv, const char *jobr, const char *jobt, const char *jobp, int *m, int *n, double *a, int *lda, double *sva, double *u, int *ldu, double *v, int *ldv, double *work, int *lwork, int *iwork, int *info);
  void lapack_dgerfsx_(const char *trans, const char *equed, int *n, int *nrhs, double *a, int *lda, double *af, int *ldaf, int *ipiv, double *r, double *c, double *b, int *ldb, double *x, int *ldx, double *rcond, double *berr, int *n_err_bnds, double *err_bnds_norm, double *err_bnds_comp, int *nparams, double *params, double *work, int *iwork, int *info); 

  void lapack_dgeqrf_(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info);
  void lapack_dlarft_(const char *direct, const char *storev, int *n, int *k, double *v, int *ldv, double *tau, double *t, int *ldt);

  double dlange_(const char *norm, int *m, int *n, double *a, int *lda, double *work);
}

/* Stubs for easier C++ interface */

class lapack_exception : public std::exception {
  char message[100];
  public:
    lapack_exception(const char* function, int info);
    virtual const char* what() const throw();
};

inline void blas_dcopy(int n, double *sx, int incx, double *sy, int incy) {
  blas_dcopy_(&n, sx, &incx, sy, &incy);
}

inline void blas_daxpy(int n, double sa, double *sx, int incx, double *sy, int incy) {
  blas_daxpy_(&n, &sa, sx, &incx, sy, &incy);
}

inline void blas_dscal(int n, double sa, double *sx, int incx) {
  blas_dscal_(&n, &sa, sx, &incx);
}

inline double blas_ddot(int n, double *sx, int incx, double *sy, int incy) {
  return blas_ddot_(&n, sx, &incx, sy, &incy);
}

inline void blas_dgemv(const char *trans, int m, int n, double alpha, double *a, int lda, double *x, int incx, double beta, double *y, int incy) {
  blas_dgemv_(trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); 
}

inline void blas_dgemm(const char *tab, int m, int n, int k, double alpha, double *a, int lda, double *b, int ldb, double beta, double *c, int ldc) {
  blas_dgemm_(tab, tab+1, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); 
}

inline void blas_dtrmm(const char *sutd, int m, int n, double alpha, double *a, int lda, double *b, int ldb) {
  blas_dtrmm_(sutd, sutd+1, sutd+2, sutd+3, &m, &n, &alpha, a, &lda, b, &ldb);
}

inline void blas_dtrsm(const char *sutd, int m, int n, double alpha, double *a, int lda, double *b, int ldb) {
  blas_dtrmm_(sutd, sutd+1, sutd+2, sutd+3, &m, &n, &alpha, a, &lda, b, &ldb);
}

inline void lapack_dgeqrf(int m, int n, double *a, int lda, double *tau, double *work, int lwork)
{
  int info;
  lapack_dgeqrf_(&m, &n, a, &lda, tau, work, &lwork, &info);
  if (info) throw(lapack_exception("DGEQRF", info));
}

inline void lapack_dorgqr(int m, int n, int k, double *a, int lda, double *tau, double *work, int lwork)
{
  int info;
  lapack_dorgqr_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
  if (info) throw(lapack_exception("DORGQR", info));
}

inline void lapack_dlarft(const char* ds, int n, int k, double *v, int ldv,
                   double *tau, double *t, int ldt) {
  lapack_dlarft_(ds, ds+1, &n, &k, v, &ldv, tau, t, &ldt);
}

inline void lapack_dgetrf(int m, int n, double *a, int lda, int *ipiv)
{
  int info;
  lapack_dgetrf_(&m, &n, a, &lda, ipiv, &info);
  if (info) throw(lapack_exception("DGETRF", info));
}

inline void lapack_dgetri(int n, double *a, int lda, int *ipiv, double *work, int lwork)
{
  int info;
  lapack_dgetri_(&n, a, &lda, ipiv, work, &lwork, &info);
  if (info) throw(lapack_exception("DGETRI", info));
}

inline void lapack_dgeqp3(int m, int n, double *a, int lda, int *jpvt, double *tau, double *work, int lwork)
{
  int info;
  lapack_dgeqp3_(&m, &n, a, &lda, jpvt, tau, work, &lwork, &info);
  if (info) throw(lapack_exception("DGEQP3", info));
}

inline void lapack_dormqr(const char *st, int m, int n, int k, double *a, int lda, double *tau, double *c, int ldc, double *work, int lwork)
{
  int info;
  lapack_dormqr_(st, st+1, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
  if (info) throw(lapack_exception("DORMQR", info)); 
}

inline void lapack_dgetrs(const char *trans, int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb)
{
  int info;
  lapack_dgetrs_(trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
  if (info) throw(lapack_exception("DGETRS", info)); 
}
  
inline void lapack_dgejsv(const char *job, int m, int n, double *a, int lda, double *sva, double *u,
                          int ldu, double *v, int ldv, double *work, int lwork, int *iwork)
{
  int info;
  lapack_dgejsv_(job, job+1, job+2, job+3, job+4, job+5, &m, &n, a, &lda, sva, u, &ldu, v, &ldv,
                 work, &lwork, iwork, &info);
  if (info) throw(lapack_exception("DGEJSV", info)); 
}

inline void lapack_dsyev(const char *ju, int n, double *a, int lda, double *w, 
                         double *work, int lwork)
{
  int info;
  lapack_dsyev_(ju, ju+1, &n, a, &lda, w, work, &lwork, &info);
  if (info) throw(lapack_exception("DSYEV", info)); 
}

inline double lapack_dlange(const char *norm, int m, int n, double *a, int lda, double *work)
{
  return dlange_(norm, &m, &n, a, &lda, work);
}
