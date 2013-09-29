module lapack_mod

  use blas_mod

  implicit none
  integer, parameter :: DLARNV_UNI_0_1  = 1
  integer, parameter :: DLARNV_UNI_N1_1 = 2
  integer, parameter :: DLARNV_NORMAL   = 3

#define DB kind(1.0d0)

  interface

     ! INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
     ! choose problem-dependent parameters for the local environment
     integer function ilaenv(ispec, name, opts, n1, n2, n3, n4)
       character*(*) :: name, opts
       integer   :: ispec, n1, n2, n3, n4
     end function ilaenv

     ! random number generator
     subroutine dlarnv(idist, iseed, n, x)
       integer   :: idist, n, iseed(4)
       real(DB)  :: x(n)
     end subroutine dlarnv
     
     ! compute eigenvalues/eigenvectors of symmetric A 
     subroutine dsyev(jobZ, uplo, n, A, lda, W, Work, lwork, info)
       character :: jobZ, uplo
       integer   :: n, lda, lwork, info
       real(DB)  :: A(lda,*), W(*), Work(*)
     end subroutine dsyev

     ! Linear system solver for general matrix
     subroutine dgesv(n, rhs, A, lda, pvt, B, ldb, info)
       integer   :: n, rhs, lda, ldb, info, pvt(*)
       real(DB)  :: A(lda,*), B(ldb, *)
     end subroutine dgesv

     ! QR-factorization with column pivoting
     subroutine dgeqp3(m, n, A, lda, jpvt, tau, work, lwork, info)
       integer   :: m, n, lda, lwork, info, jpvt(*)
       real(DB)  :: A(lda,*), tau(*), work(*)
     end subroutine dgeqp3

     ! generate Q-factor
     subroutine dorgqr(m, n, k, A, lda, tau, work, lwork, info)
       integer   :: m, n, k, lda, lwork, info
       real(DB)  :: A(lda,*), work(*), tau(*)
     end subroutine dorgqr
     
     ! multiply Q-factor to a matrix
     subroutine dormqr(side, trans, m, n, k, A, lda, tau, c, ldc, &
          work, lwork, info)
       character :: side, trans
       integer   :: info, k, lda, ldc, lwork, m, n
       real(DB)  :: A(lda,*), C(ldc,*), tau(*), work(*)
     end subroutine dormqr

     ! LU-decomposition
     subroutine dgetrf(m, n, A, lda, ipiv, info)
       integer   :: info, lda, m, n, ipiv(*)
       real(DB)  :: A(lda,*)
     end subroutine dgetrf

     ! Compute the inverse
     subroutine dgetri(n, A, lda, ipiv, work, lwork, info)
       integer   :: info, lda, n, ipiv(*), lwork
       real(DB)  :: A(lda,*), work(*)
     end subroutine dgetri

     ! solves a system of linear equations
     subroutine dgetrs(trans, n, nrhs, A, lda, ipiv, B, ldb, info)
       character :: trans
       integer   :: info, lda, ldb, n, nrhs, ipiv(*)
       real(DB)  :: A(lda,*), B(ldb,*)
     end subroutine dgetrs

     
  end interface

end module lapack_mod
