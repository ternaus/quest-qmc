module DQMC_Util
  
!  use LAPACK_MOD
!  use BLAS_MOD

  implicit none 
  
  ! 
  ! This module contains basic utilities used by other DQMC codes. 
  ! It also defines basic parameters.
  ! 
  ! Subroutine List
  ! ===============
  !
  !    DQMC_MatDiff(A, B) : evaluates the difference of two matrices.
  !    DQMC_Eye(A)        : returns A as an identity matrix.
  !    DQMC_ScaleCol(n, A, D, inv) : compute A*D or A*inv(D)
  !    DQMC_ScaleRow(n, A, D, inv) : compute D*A or inv(D)*A
  !    Error(message, no) : print out an error message and stop the program.
  !    ran0(n, var, seed)     : random number generators
  !    ran2(idum) result(ran) : random number generators from old program.
  !    dumpA(A, m, n, OPT): print out the content of A
  !
  ! Parameters
  ! ==========
  !
  integer,  parameter :: WP = kind(1.0d0)  ! work precision
  real(WP), parameter :: ZERO = 0.0D0      ! constant 0
  real(WP), parameter :: ONE  = 1.0D0      ! constant 1
  real(WP), parameter :: TWO  = 2.0D0      ! constant 1
  real(WP), parameter :: HALF = 0.5D0      ! constant 1

  integer,  parameter :: STDERR = 0        ! standard error output
  integer,  parameter :: STDOUT = 6        ! standard output
  integer,  parameter :: STDIN  = 5        ! standardinput

  character(*), parameter :: FMT_STRINT  = "(a30, i12)"
  character(*), parameter :: FMT_STRDBL  = "(a30, f19.6)"
  character(*), parameter :: FMT_STR2BL  = "(a30, '(', f11.6, ',', f11.6, ')')"
  character(*), parameter :: FMT_VALERR  = "(a30, f12.6,' +- ',f12.6)"
  character(*), parameter :: FMT_INTPAR  = "(i3,i3)"
  character(*), parameter :: FMT_DBLINE  = "(76('='))"
  character(*), parameter :: FMT_SGLINE  = "(76('-'))"
  character(*), parameter :: FMT_POINT   = "('point ; dx=', i3, ' ; dy=', i3, ' :')"

  ! Preset parameters for dlarnv() call. These parameters were defined previously in lapack_mod.F90
  ! which is no longer used. 
  integer, parameter :: DLARNV_UNI_0_1  = 1
  integer, parameter :: DLARNV_UNI_N1_1 = 2
  integer, parameter :: DLARNV_NORMAL   = 3


 
  interface conjg
     module procedure conjg_real, conjg_real1, conjg_real2
  end interface conjg

  interface DQMC_JackKnife
     module procedure DQMC_JackKnife_Real, DQMC_JackKnife_Complex
  end interface DQMC_JackKnife

  interface DQMC_SignJackKnife
     module procedure DQMC_SignJackKnife_Real, DQMC_SignJackKnife_Complex
  end interface DQMC_SignJackKnife

  interface DQMC_Print_Array
     module procedure DQMC_Print_RealArray, DQMC_Print_ComplexArray
  end interface DQMC_Print_Array

contains

  !--------------------------------------------------------!
  ! Function extension for conjg
  !--------------------------------------------------------|

  function conjg_real(x) result(y)
      real*8, target, intent(in) :: x
      real*8, pointer :: y 
      y => x
  end function conjg_real

  function conjg_real1(x) result(y)
      real*8, target, intent(in) :: x(:)
      real*8, pointer :: y(:)
      y => x
  end function conjg_real1

  function conjg_real2(x) result(y)
      real*8, target, intent(in) :: x(:,:)
      real*8, pointer :: y(:,:)
      y => x
  end function conjg_real2

  !--------------------------------------------------------!
  ! Matrix computations
  !--------------------------------------------------------|

  function DQMC_MatDiff(n, A, B) result(diff)
    !
    ! Purpose
    ! =======
    !    This function computes sum(abs(A-B)).
    !
    ! Pre-assumption
    ! ==============
    !    Matrix A and B have the same dimension.
    !    On return, A = A - B, and B is untouched.
    !
    ! Arguments
    ! =========
    integer,  intent(in)    :: n         ! the order of A and B
    real(WP), intent(in)    :: A(n,n)    ! 
    real(WP), intent(in)    :: B(n,n)    !
    !
    ! ... Return value ...
    !
    real(WP) :: diff, r

    ! ... Blas function ...
    !real(WP) :: ddot
    
    ! ... Local scalar ...
    integer  :: i

    ! ... Executable ...

    !diff = ZERO
    !  
    !do i = 1, n
    !   call daxpy(n, -ONE, B(1,i), 1, A(1,i), 1)
    !   diff = diff+ddot(n, A(1,i), 1, A(1,i), 1)
    !end do

    !diff = sqrt(diff)/n/n

    !diff = maxval(abs(A-B))

    diff = 0.d0
    do i = 1, n
       r = abs(1.d0 - A(i,i)) / abs(1.d0 - B(i,i))
       diff = max( diff, abs(log(r)) ) 
    end do

  end function DQMC_MatDiff

  !--------------------------------------------------------!

  function DQMC_MatNorm(n, A) result(norm)
    !
    ! Purpose
    ! =======
    !    This function computes sum(abs(A-B)).
    !
    ! Pre-assumption
    ! ==============
    !    Matrix A and B have the same dimension.
    !    On return, A = A - B, and B is untouched.
    !
    ! Arguments
    ! =========
    integer,  intent(in)    :: n         ! the order of A and B
    real(WP), intent(in)    :: A(n,n)    ! 
    !
    ! ... Return value ...
    !
    real(WP) :: norm

    ! ... Local scalar ...
    integer  :: i, j
    real(WP) :: maxa, temp

    ! ... Executable ...

    norm = ZERO
    maxa = ZERO

    do i = 1, n
       do j = 1, n
          temp = abs(A(i,j))
          if (temp>maxa) then
             maxa = temp
          end if
       end do
    end do

    do i = 1, n
       do j = 1, n
          temp = abs(A(i,j)) / maxa
          norm = norm+ temp*temp
       end do
    end do

    norm = sqrt(norm)*maxa

  end function DQMC_MatNorm

  !--------------------------------------------------------!

  subroutine DQMC_Eye(n, A)
    ! 
    ! Purpose
    ! =======
    !    This subroutine returns A as an identity matrix.
    !
    ! Pre-assumption
    ! ==============
    !    Matrix A is square.
    !
    ! Arguments
    ! =========
    integer,  intent(in)    :: n        ! The order of A
    real(WP), intent(inout) :: A(n,n)   ! returned identity
   
    ! ... Local scalar ...
    integer  :: i
   
    ! ... Executable ...

    A = ZERO

    do i = 1,n
       A(i,i) = ONE
    end do
    
  end subroutine DQMC_Eye

  !--------------------------------------------------------!

  subroutine DQMC_Trans(n, At, A)
    !
    ! Purpose
    ! =======
    !    This subroutine returns the transpose of A
    !
    ! Argument
    ! ========
    integer, intent(in)     :: n
    real(wp), intent(inout) :: At(n,n)
    real(wp), intent(in)    :: A(n,n)
    
    ! ... local scalar 
    
    !! decide the sgn of det(Q)
    !do i = 1, n
    !   At(i,1:n) = A(1:n,i)
    !end do
    At = transpose(A)
    
  end subroutine DQMC_Trans

  !--------------------------------------------------------!
  
  subroutine DQMC_ScaleCol(n, A, D)
    ! 
    ! Purpose
    ! =======
    !    This subroutine computes A = A*D.
    !
    ! Pre-assumption
    ! ==============
    !    Matrix A is order n and D is of length n.
    !
    ! Argument
    ! ========
    integer,  intent(in)    :: n         ! The order of A
    real(WP), intent(inout) :: A(n,n)    ! 
    real(WP), intent(in)    :: D(n)      !

    ! ... Local scalar ...
    integer  :: i

    ! ... Executable ...

    ! A = A*D
    do i = 1, n
       call dscal(n, D(i), A(:, i), 1)
    end do

  end subroutine DQMC_ScaleCol
  
  !--------------------------------------------------------!
  
  subroutine DQMC_ScaleRow(n, A, D)
    ! 
    ! Purpose
    ! =======
    !    This subroutine computes A = D*A.
    !
    ! Pre-assumption
    ! ==============
    !    Matrix A is order n and D is of length n.
    !
    ! Argument
    ! ========
    integer,  intent(in)    :: n        ! The order of A
    real(WP), intent(inout) :: A(n,n)   !
    real(WP), intent(in)    :: D(n)     ! 

    ! ... Local scalar ...
    integer  :: i
    
    ! A = D*A
    do i = 1, n
       call dscal(n, D(i), A(i:,1), n)
    end do

  end subroutine DQMC_ScaleRow

  !--------------------------------------------------------!
  
  subroutine DQMC_ScaleColInv(n, A, D)
    ! 
    ! Purpose
    ! =======
    !    This subroutine computes A=A*inv(D).
    !
    ! Pre-assumption
    ! ==============
    !    Matrix A is order n and D is of length n.
    !
    ! Argument
    ! ========
    integer,  intent(in)    :: n         ! The order of A
    real(WP), intent(inout) :: A(n,n)    ! 
    real(WP), intent(in)    :: D(n)      !

    ! ... Local scalar ...
    real(WP) :: uno
    integer  :: i

    ! ... Executable ...
    uno = ONE

    ! A = A*inv(D)
    do i = 1, n
       call dscal(n, uno/D(i), A(:, i), 1)
    end do

  end subroutine DQMC_ScaleColInv

  !--------------------------------------------------------!
  
  subroutine DQMC_ScaleRowInv(n, A, D)
    ! 
    ! Purpose
    ! =======
    !    This subroutine computes A = inv(D)*A.
    !
    ! Pre-assumption
    ! ==============
    !    Matrix A is order n and D is of length n.
    !
    ! Argument
    ! ========
    integer,  intent(in)    :: n        ! The order of A
    real(WP), intent(inout) :: A(n,n)   !
    real(WP), intent(in)    :: D(n)     ! 

    ! ... Local scalar ...
    real(WP) :: uno
    integer  :: i
    
    ! ... Executable ...
    uno = ONE

    ! A = inv(D)*A
    do i = 1, n
       call dscal(n, uno/D(i), A(i:,1), n)       
    end do

  end subroutine DQMC_ScaleRowInv
 
  !--------------------------------------------------------------------!
  ! Statistics
  !--------------------------------------------------------------------!

  subroutine DQMC_SignJackKnife_Real(n, avg, err, x, y, sgn, sum_sgn)
    !
    ! Purpose
    ! =======
    !    This subroutine implements delete-1 JackKnife method for
    !    the data and sign.
    !    X = (x1, x2, ..., xn) be the input data.
    !    sgn = (sgn1, sgn2, ..., sgnn)
    !
    !    Y = (y1, y2, ..., yn) is the Jacknife resampling of X with sign.
    !    
    !    where y_i = (sum(x)-x_i)/sgn_i
    !    The JackKnife variance of X with sign is defined as 
    !
    !     n-1 
    !    ----- sqrt(sum(y_i-avg_y)^2)
    !      n
    !
    !    where avg_y = sum(y)/n
    !
    ! Arguments
    ! =========
    integer,  intent(in)    :: n
    real(wp), intent(out)   :: avg
    real(wp), intent(out)   :: err
    real(wp), intent(in)    :: x(:)
    real(wp), intent(inout) :: y(:)
    real(wp), intent(in)    :: sgn(:)
    real(wp), intent(in)    :: sum_sgn
    
    ! ... Parameter ...
    real(wp), parameter :: TOL = 1.0D-12
    
    ! ... Local variable ...
    real(wp) :: sum_x, avg_y
    
    ! ... Executable ...
    
    ! standard division
    sum_x  = sum(x)
    y      = (sum_x-x(1:n))/sgn(1:n)
    avg_y  = sum(y)/n

    y    = y - avg_y
    y    = y * y
    err  = sum(y)*(n-1)/n
    err  = sqrt(err)

    ! compute average
    avg = sum_x/sum_sgn
    
    ! If error is small enough, then regard it as 0.
    if (err .lt. TOL*abs(avg)) then
       err = ZERO
    end if
    
  end subroutine DQMC_SignJackKnife_Real
  
  !--------------------------------------------------------------------!

  subroutine DQMC_SignJackKnife_Complex(n, avg, err, x, y, sgn, sum_sgn)
    !
    ! Purpose
    ! =======
    !    This subroutine implements delete-1 JackKnife method for
    !    the data and sign.
    !    X = (x1, x2, ..., xn) be the input data.
    !    sgn = (sgn1, sgn2, ..., sgnn)
    !
    !    Y = (y1, y2, ..., yn) is the Jacknife resampling of X with sign.
    !    
    !    where y_i = (sum(x)-x_i)/sgn_i
    !    The JackKnife variance of X with sign is defined as 
    !
    !     n-1 
    !    ----- sqrt(sum(y_i-avg_y)^2)
    !      n
    !
    !    where avg_y = sum(y)/n
    !
    ! Arguments
    ! =========
    integer,  intent(in)    :: n
    complex(wp), intent(out)   :: avg
    complex(wp), intent(out)   :: err
    complex(wp), intent(in)    :: x(:)
    complex(wp), intent(inout) :: y(:)
    complex(wp), intent(in)    :: sgn(:)
    complex(wp), intent(in)    :: sum_sgn
    
    ! ... Parameter ...
    real(wp), parameter :: TOL = 1.0D-12
    
    ! ... Local variable ...
    complex(wp) :: sum_x, avg_y
    real(wp) :: rp, ip
    integer  :: i
    
    ! ... Executable ...
    
    ! standard division
    sum_x  = sum(x)
    y      = (sum_x-x(1:n))/sgn(1:n)
    avg_y  = sum(y)/n

    y    = y - avg_y
    do i = 1, n
       rp = dble(y(i))**2
       ip = aimag(y(i))**2
       err = err + dcmplx(rp, ip)
    enddo
    err  = (err*(n-1))/n
    rp   = sqrt(dble(err))
    ip   = sqrt(aimag(err))
    err  = dcmplx(rp, ip)

    ! compute average
    avg = sum_x/sum_sgn
    
    ! If error is small enough, then regard it as 0.
    if (abs(err) .lt. TOL*abs(avg)) then
       err = ZERO
    end if
    
  end subroutine DQMC_SignJackKnife_Complex

  !--------------------------------------------------------------------!
  
  subroutine DQMC_JackKnife_Real(n, avg, err, x, y, sgn, sum_sgn)
    !
    ! Purpose
    ! =======
    !    This subroutine implements delete-1 JackKnife method. Let
    !    X = (x1, x2, ..., xn) be the input data.
    !    Y = (y1, y2, ..., yn) is the Jacknife resampling of X.
    !    
    !    where y_i = (sum(x)-x_i)/(n-1)
    !    The JackKnife variance of X is defined as 
    !
    !     n-1 
    !    ----- sqrt(sum(y_i-avg_y)^2)
    !      n
    !
    !    where avg_y = sum(y)/n, which equals to avg_x
    !    
    ! Arguments
    ! =========
    integer,  intent(in)    :: n
    real(wp), intent(out)   :: err
    real(wp), intent(out)   :: avg
    real(wp), intent(in)    :: x(n)
    real(wp), intent(out)   :: y(n)
    real(wp), intent(out)   :: sgn(n)
    real(wp), intent(out)   :: sum_sgn


    ! ... Parameter ...
    real(wp), parameter :: TOL = 1.0D-12
    
    ! ... Local variable ...
    real(wp) :: sum_x
    
    ! ... Executable ...

    sum_x  = sum(x)
    ! compute average
    avg    = sum_x/ n

    ! sgn and sum sgn will be used for other analysis
    sgn    = (sum_x-x(1:n))
    sum_sgn= sum_x

    ! compute y (Jackkife sample)
    y    = sgn/(n-1)

    ! avg_y = avg_x 
    y    = y - avg
    y    = y * y
    err  = (sum(y)*(n-1))/n
    err  = sqrt(err)

    ! If error is small enough, then regard it as 0.
    if (err .lt. TOL*abs(avg)) then
       err = ZERO
    end if
    
  end subroutine DQMC_JackKnife_Real

  !--------------------------------------------------------------------!
  
  subroutine DQMC_JackKnife_Complex(n, avg, err, x, y, sgn, sum_sgn)
    !
    ! Purpose
    ! =======
    !    This subroutine implements delete-1 JackKnife method. Let
    !    X = (x1, x2, ..., xn) be the input data.
    !    Y = (y1, y2, ..., yn) is the Jacknife resampling of X.
    !    
    !    where y_i = (sum(x)-x_i)/(n-1)
    !    The JackKnife variance of X is defined as 
    !
    !     n-1 
    !    ----- sqrt(sum(y_i-avg_y)^2)
    !      n
    !
    !    where avg_y = sum(y)/n, which equals to avg_x
    !    
    ! Arguments
    ! =========
    integer,  intent(in)    :: n
    complex(wp), intent(out)   :: err
    complex(wp), intent(out)   :: avg
    complex(wp), intent(in)    :: x(n)
    complex(wp), intent(out)   :: y(n)
    complex(wp), intent(out)   :: sgn(n)
    complex(wp), intent(out)   :: sum_sgn


    ! ... Parameter ...
    real(wp), parameter :: TOL = 1.0D-12
    
    ! ... Local variable ...
    complex(wp) :: sum_x
    real(wp) :: rp, ip
    integer  :: i
    
    ! ... Executable ...

    sum_x  = sum(x)
    ! compute average
    avg    = sum_x/ n

    ! sgn and sum sgn will be used for other analysis
    sgn    = (sum_x-x(1:n))
    sum_sgn= sum_x

    ! compute y (Jackkife sample)
    y    = sgn/(n-1)

    ! avg_y = avg_x 
    y    = y - avg
    do i = 1, n
       rp = dble(y(i))**2
       ip = aimag(y(i))**2
       err = err + dcmplx(rp, ip)
    enddo
    err  = (err*(n-1))/n
    rp   = sqrt(dble(err))
    ip   = sqrt(aimag(err))
    err  = dcmplx(rp, ip)

    ! If error is small enough, then regard it as 0.
    if (abs(err) .lt. TOL*abs(avg)) then
       err = ZERO
    end if
    
  end subroutine DQMC_JackKnife_Complex
  
  !--------------------------------------------------------------------!

  subroutine DQMC_GetErr(n, err, avg, list)
    !
    ! Purpose
    ! =======
    !    This subroutine computes error of the measurements.
    !
    ! Arguments
    ! =========
    integer,  intent(in)    :: n
    real(wp), intent(out)   :: err
    real(wp), intent(in)    :: avg
    real(wp), intent(inout) :: list(n)
    
    ! ... Parameter ...
    real(wp), parameter :: TOL = 1.0D-12
    
    ! ... Local variable ...
    real(wp) :: tmp
    
    ! ... Executable ...
    
    ! compute average
    tmp  = sum(list)/n
    
    ! standard division
    list =  list - tmp
    list = list * list
    err  = (sum(list)*(n-1))/n
    err  = sqrt(err)
    
    ! If error is small enough, then regard it as 0.
    if (err .lt. TOL*abs(avg)) then
       err = ZERO
    end if
    
  end subroutine DQMC_GetErr  

  !-----------------------------------------------------------------!

  subroutine DQMC_GetErr1(n, data, avg, err)
    !
    ! Purpose
    ! =======
    !    This subroutine computes error of the measurements.
    !
    ! Arguments
    ! =========
    integer, intent(in)   :: n
    real(wp), intent(in)  :: data(n)
    real(wp), intent(out) :: avg, err
 
    ! ... Parameter ...
    real(wp), parameter :: TOL = 1.0D-12

    ! ... local vars ...
    integer :: i
    real(wp):: s
    
    ! Executable
    avg = sum(data) / n
    s = ZERO
    do i = 1, n
       s = s + (data(i)-avg)**2
    end do
    s = s / n
    err = sqrt(s)/sqrt(n-ONE)

    if (err .lt. TOL*abs(avg)) then
       err = ZERO
    end if

  end subroutine DQMC_GetErr1
   
  !--------------------------------------------------------!

  subroutine DQMC_GetErr2(n, sm, ssq, avg, err)
    !
    ! Purpose
    ! =======
    !    This subroutine computes error of the measurements.
    !
    ! Arguments
    ! =========
    integer, intent(in)   :: n
    real(wp), intent(in)  :: sm, ssq
    real(wp), intent(out) :: avg, err
 
    ! ... Parameter ...
    real(wp), parameter :: TOL = 1.0D-12

    ! ... local vars ...
    real(wp):: s
    
    ! Executable
    avg = sm / n
    s   = ssq / n - avg*avg
    err = sqrt(s) / sqrt(n-ONE)

    if (err .lt. TOL*abs(avg)) then
       err = ZERO
    end if

  end subroutine DQMC_GetErr2

  !--------------------------------------------------------!
  ! Miscellanceous
  !--------------------------------------------------------!

  subroutine DQMC_Error(message, no)
    ! 
    ! Purpose
    ! =======
    !    This subroutine prints out an error message and
    !    a message number, and then stop the program.
    !
    ! Arguments
    ! =========
    ! 
    character(*), intent(in) :: message   ! Error message
    integer, intent(in)      :: no        ! Message number
    integer                  :: dum

    ! ... Executable ...

    dum = no
    write(STDERR,*) "Error: ", message
    stop

  end subroutine DQMC_Error
  
  !--------------------------------------------------------!

  subroutine DQMC_Warning(message, no)
    ! 
    ! Purpose
    ! =======
    !    This subroutine prints out an error message and
    !    a message number.
    !
    ! Arguments
    ! =========
    ! 
    character(*), intent(in) :: message   ! Warning message
    integer, intent(in)      :: no        ! Message number
    integer                  :: dum

    ! ... Executable ...

    dum = no
    write(STDERR,*) "Warning: ", message

  end subroutine DQMC_Warning
  
  !--------------------------------------------------------!

#ifdef _QMC_MPI

  subroutine ran0(n, var, seed)
#   ifdef _QMC_MPI
#      define SIMPLE_SPRNG
#      include "sprng_f.h"
#   endif
    ! 
    ! Purpose
    ! =======
    !    Random number generator. This subroutine uses 
    !    SPRNG library.
    !
    ! Arguments
    ! =========
    ! 
    integer, intent(in)     :: n
    real(wp), intent(out)   :: var(n)     ! random number to return
    integer                 :: i
    integer, intent(in)     :: seed(4)    ! random seeds(not used)
    
    i = seed(1) !Avoid warning

    do i=1,n
       var(i) = sprng()
    enddo

  end subroutine ran0

#else

  subroutine ran0(n, var, seed)
    ! 
    ! Purpose
    ! =======
    !    Random number generator. This subroutine uses 
    !    LAPACK's random number generator dlaruv to generate
    !    a list of random numbers
    !
    ! Arguments
    ! =========
    ! 
    integer, intent(in)     :: n          ! length of the list
    real(wp), intent(out)   :: var(n)     ! random number to return
    integer, intent(inout)  :: seed(4)    ! random seeds

    ! ... local scalar ...
    integer, parameter :: max_len = 128   ! This max length is 
                                          ! defined by dlaruv
    ! ... Executable ...

    var(1:n) = ZERO
    call dlarnv(DLARNV_UNI_0_1, seed, n, var)

  end subroutine ran0

#endif

  !--------------------------------------------------------!

  subroutine ran1(n, var, seed)
    ! 
    ! Purpose
    ! =======
    !    Random number generator. This subroutine uses 
    !    LAPACK's random number generator dlaruv to generate
    !    a list of random numbers
    !
    ! Arguments
    ! =========
    ! 
    integer, intent(in)     :: n          ! length of the list
    real(wp), intent(out)   :: var(n)     ! random number to return
    integer, intent(inout)  :: seed(4)    ! random seeds

    ! ... local scalar ...
    integer, parameter :: max_len = 128   ! This max length is 
                                          ! defined by dlaruv
    ! ... Executable ...

    var(1:n) = ZERO
    call dlarnv(2, seed, n, var)

  end subroutine ran1

 !--------------------------------------------------------!

  integer function intran(L, seed)
    ! 
    ! Purpose
    ! =======
    !    Random number generator. This subroutine
    !    generates a random integer in [1:L]
    !
    ! Arguments
    ! =========
    ! 
    integer, intent(in)     :: L          ! length of the list
    integer, intent(inout)  :: seed(4)    ! random seeds

    real(wp) :: var(1)

    ! ... Executable ...
    
    call ran0(1, var, seed)
    intran = ceiling(var(1) * L)

  end function intran

 !--------------------------------------------------------!

  subroutine ranN(n, var, seed)
    ! 
    ! Purpose
    ! =======
    !    Random number generator. This subroutine uses 
    !    LAPACK's random number generator dlaruv to generate
    !    a list of random numbers
    !
    ! Arguments
    ! =========
    ! 
    integer, intent(in)     :: n          ! length of the list
    real(wp), intent(out)   :: var(n)     ! random number to return
    integer, intent(inout)  :: seed(4)    ! random seeds

    ! ... local scalar ...
    integer, parameter :: max_len = 128   ! This max length is 
                                          ! defined by dlaruv
    ! ... Executable ...

    var(1:n) = ZERO
    call dlarnv(3, seed, n, var)

  end subroutine ranN

  !-------------------------------------------------------------!

  subroutine dumpA(m, n, A, OPT)
    implicit none
    ! 
    ! Purpose
    ! =======
    !    This subroutine prints the content of matrix A.
    !    *** for internal debugging only
    !
    ! Arguments
    ! =========
    ! 
    integer, intent(in)  :: m, n          ! dimension of A
    real(wp), intent(in) :: A(1:m,1:n)    ! matrix A
    integer, intent(in)  :: OPT           ! output device

    ! ... Local variables ...
    character(20) fmt
    integer i 
    real(wp) :: temp(n)

    ! ... Executable ...

    write(fmt,"(A,I3,A)") "(",n,"F20.15)"
    
    do i=1, m
       temp = A(i,1:n)
       write (OPT, fmt) temp
    end do
    
  end subroutine dumpA

  !--------------------------------------------------------!
    
  subroutine DQMC_Print_RealArray(n, m, title, label, avg, err, OPT)
    !
    ! Purpose
    ! =======
    !    This subroutine prints out the content of avg +- err
    !   
    ! Arguments
    ! =========
    integer,  intent(in)    :: OPT, n, m
    character(*), intent(in):: title, label(:)
    real(wp), intent(in)    :: avg(:,:), err(:,:)
    
    ! ... Local variable ...
    integer  :: i, j
    
    ! ... Executable ...

    write(OPT,*) title
    if (n .gt. 0) then
       do i = 1, n
          write(OPT,*) label(i)
          do j = 1, m
             write(OPT, "(i3,e16.8,' +-',e16.8)") j-1, avg(i,j), err(i,j)
          end do
       end do
    else
       do j = 1, m
          write(OPT, "(a,e16.8,' +-',e16.8)") label(j), avg(j,1), err(j,1)
       end do
    end if
    !write(OPT,*)
    write(OPT,FMT_DBLINE)

  end subroutine DQMC_Print_RealArray

  !--------------------------------------------------------------------!
  
  subroutine DQMC_Print_ComplexArray(n, m, title, label, avg, err, OPT)
    !
    ! Purpose
    ! =======
    !    This subroutine prints out the content of avg +- err
    !   
    ! Arguments
    ! =========
    integer,  intent(in)     :: OPT, n, m
    character(*), intent(in) :: title, label(:)
    complex(wp), intent(in)  :: avg(:,:), err(:,:)
    
    ! ... Local variable ...
    integer  :: i, j
    
    ! ... Executable ...
    
    write(OPT,*) title

    if (n .gt. 0) then
       do i = 1, n
          write(OPT,*) label(i)
          do j = 1, m
             write(OPT,"(i3,'(',e15.8,' +-',e15.8,') &
                  & +i (',e15.8,' +-',e15.8,')')")&
                  j-1, dble(avg(i,j)), dble(err(i,j)), aimag(avg(i,j)), aimag(err(i,j+m))
          end do
       end do
    else
       do j = 1, m
          write(OPT,"(3x, A,'(',e15.8,' +-',e15.8,') &
               & +i (',e15.8,' +-',e15.8,')')") &
               label(j), dble(avg(j,1)), dble(err(j,1)), aimag(avg(j,1)), aimag(err(j,1))
       end do
    end if
    
    !write(OPT,*)
    write(OPT,FMT_DBLINE)

  end subroutine DQMC_Print_ComplexArray  

  !--------------------------------------------------------------------!

  subroutine DQMC_Print_EigenMode(n, m, title, value, OPT)
    !
    ! Purpose
    ! =======
    !    This subroutine prints out the eigenmodes of correlation functions
    !   
    ! Arguments
    ! =========
    integer,  intent(in)    :: OPT, n, m
    character(*), intent(in):: title
    complex*16, intent(in)  :: value(:,:,:)

    ! ... Local variable ...
    integer  :: i, ik, ja, ia

    ! ... Executable ...

    write(OPT,*) title

    i=0
    write(OPT,'(10x,100(13x,i2,12x))')(ia,ia=1,n)
    do ik = 1, m
      do ja = 1, n
        i=i+1
        if(ja==1)then
          write(OPT,'(1x,i3,2x,i3,1x)',advance="no") ik,ja
        else
          write(OPT,'(6x,i3,1x)',advance="no") ja
        endif
        do ia = 1, n
          write(OPT,"(A,f10.6,' +i',f10.6,A)",advance="no")&
           " (",real(value(ia,ja,ik)),aimag(value(ia,ja,ik)),") "
        end do
        write(OPT,'(2(i5))',advance="yes")
      end do
    end do

    write(OPT,FMT_DBLINE)

  end subroutine DQMC_Print_EigenMode

  !--------------------------------------------------------------------!

  subroutine dqmc_getFTk(value, n, nclass, class, na, nk, ft_wgt, phase, valuek)

     real(wp),    intent(in)  :: value(nclass)
     integer,     intent(in)  :: n
     integer,     intent(in)  :: nclass
     integer,     intent(in)  :: class(n,n)
     integer,     intent(in)  :: phase(n,n)
     integer,     intent(in)  :: na, nk
     complex(wp), intent(in)  :: ft_wgt(n/na, nk)
     complex(wp), intent(out) :: valuek(nk*na*(na+1)/2)

     integer     :: ik, ia, ja, it, jt, i, j, nt, naa
     complex(wp) :: U(na,na), phcurr

     nt = n / na
     naa = na*(na+1)/2

     !Loop over inequvalent k-points
     do ik = 1, nk

        U = 0.d0
        !Compute Fourier transform matrix at "ik" for (ja,ia)
        do ia = 1, na
           do ja = ia, na
              !sum over translations
              do it = 1, nt
                 do jt = 1, nt
                    !Find atom which is the translation of "ja" by "it"
                    i = (it-1) * na + ia
                    j = (jt-1) * na + ja
                    !Use class that corresponds to the (i,j) pair 
                    phcurr = phase(i,j) * value(class(i,j))
                    U(ja,ia) = U(ja,ia) + phcurr * ft_wgt(it,ik) * dconjg(ft_wgt(jt,ik))
                 enddo
              enddo
           enddo
        enddo

        !Pointer to Fourier transform of "ip" at "ik" in "ibin"
        i = (ik-1)*naa
        do ia = 1, na
           do ja = ia, na
              i = i + 1
              valuek(i) = U(ja,ia)
           enddo
        enddo

     enddo ! Loop over k-points

  end subroutine DQMC_getFTk

  !--------------------------------------------------------------------!

  subroutine DQMC_IO_open(fname, INP_UNIT, OUT_UNIT)
    !
    ! Purpose
    ! =======
    !  Find a unit for input and output file
    !   
    ! Arguments
    ! =========

     integer, intent(out) :: INP_UNIT, OUT_UNIT
     character(len=60), intent(out) :: fname

     character(len=60) :: outname

     !Open input file
     call get_command_argument(1, fname)
     call DQMC_open_file(fname, 'old', INP_UNIT)

     !Open output file
     outname=trim(adjustl(fname))//".out"
     call DQMC_open_file(outname, 'unknown', OUT_UNIT)

  end subroutine DQMC_IO_open

  !--------------------------------------------------------------------!

  subroutine DQMC_open_file(fname, fstatus, FILE_UNIT)
  
     implicit none
     character(len=*), intent(in) :: fname, fstatus
     integer, intent(out) :: FILE_UNIT
  
     logical :: unit_is_open
  
     do FILE_UNIT = 7, 99
        inquire(unit=FILE_UNIT, opened=unit_is_open)
        if (.not.unit_is_open) exit
     enddo
     open(unit=FILE_UNIT, file=fname, status=fstatus)

  end subroutine DQMC_open_file

  !--------------------------------------------------------------------!

  subroutine DQMC_count_records(n, FILE_UNIT)
  
     implicit none
     integer, intent(out):: n
     integer, intent(in) :: FILE_UNIT
  
     character(len=1) :: c
     integer :: i
  
     n = 0
     do 
        read(FILE_UNIT,'(A1)', iostat=i)c
        if (i .ne. 0) exit
        n = n + 1
     enddo
     rewind(FILE_UNIT)
  
  end subroutine DQMC_count_records

  !--------------------------------------------------------------------!

  logical function move_to_record(string,iunit)
    !
    ! Purpose
    ! =======
    !  Move to record where string is found. If not found return false.  
    !   
    ! Arguments
    ! =========
    character(len=*), intent(in) :: string
    integer, intent(in) :: iunit

    ! ... Local Variables ...
    integer  :: ios,istring
    character(len=100) :: line

    ! ... Executable ...

    rewind(iunit)

    do
       read(iunit,'(A)',iostat=ios)line
       if(ios.ne.0)then
          move_to_record=.false.
          exit
       endif
       istring=index(line,trim(adjustl(string)))
       if(istring==0)cycle
       move_to_record=.true.
       exit
    enddo
  end function
  
  !--------------------------------------------------------------------!
  
  real(wp) function get_det(a)
    !
    ! Purpose
    ! =======
    ! Directly get determinant of 3 by 3 matrix.
    !   
    ! Arguments
    ! =========
    real(wp), intent(in) :: a(3,3)

    ! ... Local variables ...
    real(wp) :: d(3)

    ! ... Executable ...
    d(1) = a(2,1)*a(3,2) - a(3,1)*a(2,2)
    d(2) = a(3,1)*a(1,2) - a(1,1)*a(3,2)
    d(3) = a(1,1)*a(2,2) - a(2,1)*a(1,2)

    get_det = d(1)*a(1,3) + d(2)*a(2,3) + d(3)*a(3,3)

  end function get_det
  
  !--------------------------------------------------------------------!
  
  subroutine get_inverse(a,inv)
    !
    ! Purpose
    ! =======
    !  Directly get inverse of 3 by 3 matrix.
    !   
    ! Arguments
    ! =========
    real(wp), intent(in) :: a(3,3)
    real(wp), intent(out) :: inv(3,3)

    ! ... Local variables ...
    real(wp) :: invdet

    ! ... Executable ...
    inv(1,1) = a(2,2)*a(3,3) - a(2,3)*a(3,2)
    inv(1,2) = a(3,2)*a(1,3) - a(3,3)*a(1,2)
    inv(1,3) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
    inv(2,1) = a(2,3)*a(3,1) - a(2,1)*a(3,3)
    inv(2,2) = a(3,3)*a(1,1) - a(3,1)*a(1,3)
    inv(2,3) = a(1,3)*a(2,1) - a(1,1)*a(2,3)
    inv(3,1) = a(2,1)*a(3,2) - a(2,2)*a(3,1)
    inv(3,2) = a(3,1)*a(1,2) - a(3,2)*a(1,1)
    inv(3,3) = a(1,1)*a(2,2) - a(1,2)*a(2,1)

    invdet = 1.d0 / (inv(3,1)*a(1,3) + inv(3,2)*a(2,3) + inv(3,3)*a(3,3))

    inv(:,:) = inv(:,:)*invdet

  end subroutine get_inverse

  !--------------------------------------------------------------------!

end module DQMC_Util
