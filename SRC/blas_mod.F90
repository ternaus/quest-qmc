module blas_mod

  ! 
  !  This module is designed for interface checking
  !

#define DB kind(1.0d0)

#ifdef DQMC_PROFILE

#define BLAS_BEGIN if (blas_profile) call system_clock(blas_t1)
#define BLAS_END(c, t, f, n) if (blas_profile) call blas_end(c, t, f, n)

    real(DB) :: blas_dcopy_time = 0, blas_daxpy_time = 0, blas_dscal_time = 0
    real(DB) :: blas_dgemv_time = 0, blas_dgemm_time = 0, blas_dtrmm_time = 0
    real(DB) :: blas_ddot_time = 0
    integer :: blas_t1, blas_t2
    integer  :: blas_dcopy_count = 0, blas_daxpy_count = 0, blas_dscal_count = 0
    integer  :: blas_dgemv_count = 0, blas_dgemm_count = 0, blas_dtrmm_count = 0
    integer  :: blas_ddot_count = 0
    integer  :: blas_profile = 0, blas_rate
    real(DB) :: blas_dcopy_flops = 0, blas_daxpy_flops = 0, blas_dscal_flops = 0
    real(DB) :: blas_dgemv_flops = 0, blas_dgemm_flops = 0, blas_dtrmm_flops = 0
    real(DB) :: blas_ddot_flops = 0
    
  contains

    subroutine blas_end(c, t, f, n)
      integer :: c, n
      real(DB) :: t, f
      call system_clock(blas_t2, blas_rate)
      c = c + 1
      t = t + (blas_t2 - blas_t1) / REAL(blas_rate)
      f = f + n
    end subroutine blas_end

    subroutine blas_print()
       write(*,*) "DCOPY         ", blas_dcopy_count, blas_dcopy_time, blas_dcopy_flops / blas_dcopy_time / 1e6
       write(*,*) "DDOT          ", blas_ddot_count, blas_ddot_time, blas_ddot_flops / blas_ddot_time / 1e6
       write(*,*) "DAXPY         ", blas_daxpy_count, blas_daxpy_time, blas_daxpy_flops / blas_daxpy_time / 1e6
       write(*,*) "DSCAL         ", blas_dscal_count, blas_dscal_time, blas_dscal_flops / blas_dscal_time / 1e6
       write(*,*) "DGEMV         ", blas_dgemv_count, blas_dgemv_time, blas_dgemv_flops / blas_dgemv_time / 1e6
       write(*,*) "DGEMM         ", blas_dgemm_count, blas_dgemm_time, blas_dgemm_flops / blas_dgemm_time / 1e6
       write(*,*) "DTRMM         ", blas_dtrmm_count, blas_dtrmm_time, blas_dtrmm_flops / blas_dtrmm_time / 1e6
    end subroutine blas_print
    
#else

#define BLAS_BEGIN
#define BLAS_END(c, t, f, n)

  contains

#endif

     ! y = x
     subroutine blas_dcopy(n, x, incx, y, incy)
       integer   :: n, incx, incy
       real(DB)  :: x(*), y(*)
       BLAS_BEGIN
       call dcopy(n, x, incx, y, incy)
       BLAS_END(blas_dcopy_count, blas_dcopy_time, blas_dcopy_flops, n)
     end subroutine blas_dcopy

     ! y = a*x + y
     subroutine blas_daxpy(n, a, x, incx, y, incy)
       integer   :: n, incx, incy
       real(DB)  :: a, x(*), y(*)
       BLAS_BEGIN
       call daxpy(n, a, x, incx, y, incy)
       BLAS_END(blas_daxpy_count, blas_daxpy_time, blas_daxpy_flops, 2 * n)
     end subroutine blas_daxpy
          
     ! dx = da*dx
     subroutine blas_dscal(n, da, dx, incx)
       integer   :: n, incx
       real(DB)  :: da, dx(*) 
       BLAS_BEGIN
       call dscal(n, da, dx, incx)
       BLAS_END(blas_dscal_count, blas_dscal_time, blas_dscal_flops, n)
     end subroutine blas_dscal

     real(DB) function blas_ddot(n, x, incx, y, incy)
       integer   :: n, incx, incy
       real(DB)  :: x(*), y(*)
       real(DB), external :: ddot
       BLAS_BEGIN
       blas_ddot = ddot(n, x, incx, y, incy)
       BLAS_END(blas_ddot_count, blas_ddot_time, blas_ddot_flops, 2 * n)
     end function blas_ddot

     ! matrix-vector multiplication: y = alpha*A*x + beta*y
     ! DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
     subroutine blas_dgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
       character :: trans
       integer   :: m, n, lda, incx, incy
       real(DB)  :: alpha, beta, A(lda,*), x(*), y(*)
       BLAS_BEGIN
       call dgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
       BLAS_END(blas_dgemv_count, blas_dgemv_time, blas_dgemv_flops, 2 * m * n)
     end subroutine blas_dgemv


     ! matrix-matrix multiplication
     subroutine blas_dgemm(transA, transB, m, n, k, alpha, A, lda, &
          B, ldb, beta, c, ldc)
       character :: transA, transB
       integer   :: m, n, k, lda, ldb, ldc
       real(DB)  :: alpha, beta, A(lda,*), B(ldb,*), C(ldc,*)
       BLAS_BEGIN
       call dgemm(transA, transB, m, n, k, alpha, A, lda, &
          B, ldb, beta, c, ldc)
       BLAS_END(blas_dgemm_count, blas_dgemm_time, blas_dgemm_flops, 2 * k * m * n)
     end subroutine blas_dgemm

     ! triangular matrix-matrix multiplication
     subroutine blas_dtrmm (side, uplo, transA, diag, m, n, alpha, &
          A, lda, B, ldb)
       character :: side, uplo, transA, diag
       integer   :: m, n, lda, ldb
       real(DB)  :: alpha, A(lda,*), B(ldb,*)
       BLAS_BEGIN
       call dtrmm (side, uplo, transA, diag, m, n, alpha, &
          A, lda, B, ldb)
       if (UPLO .EQ. 'L') then
         BLAS_END(blas_dtrmm_count, blas_dtrmm_time, blas_dtrmm_flops, m * m * n)
       else
         BLAS_END(blas_dtrmm_count, blas_dtrmm_time, blas_dtrmm_flops, m * n * n)
       end if
     end subroutine blas_dtrmm
     
  
end module blas_mod
