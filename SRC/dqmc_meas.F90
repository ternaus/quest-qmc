module DQMC_Phy0
#include "dqmc_include.h"

  use DQMC_UTIL

  implicit none 
  
  !
  ! This module puts all utilities needed for physical measurements.
  !

contains

  ! Subroutines
  ! ==================================================================
  
  subroutine DQMC_AutoCorrTime(N, maxBin, data, sgn, avg, var)  
    ! 
    ! Purpose 
    ! =======
    !     This subroutine computes autocorrelation time of a measurements.
    !     From Mark's note, this will be done by calculating the 
    !     variance of the data, vari(l) as a function of bin size.
    !     the number of datum is meas*run, we will allow l, the bin
    !     size, to run from 1 to meas.
    !     
    !                 l vari(n)                 1    N             2
    !     autolth(n) = ---------     vari(n) = ---- SUM (x_m - <x>)
    !                 2 vari(1)                N-1  m=1
    !     
    !     where x_m is the bin averages.  The autocorrelation length is
    !     done this way for each Markov process and then the results are
    !     reduced to process 0.
    !     
    integer, intent(in)     :: N               ! number of data
    integer, intent(in)     :: maxBinSize          ! max number of bins
    real(wp), intent(in)    :: data(N), sgn(N)
    real(wp), intent(inout) :: avg(maxBinSize), var(maxBinSize)

    ! ... Local Variables ...
    real(wp)  :: ratio, sum1, sum2, r
    integer   :: i, j, maxBinSize, nBin, bSize

    ! ... Executable ...
    
    sum1 = sum(data)
    sum2 = sum(abs(sgn))
    avg = ZERO
    var = ZERO

    ! quick return when there is a sign problem
    if (sum2 .le. epsilon()) then
       call DQMC_Warning("Warning: sign problem", 0)
       return
    else
       ratio = sum1/sum2
    end if


    ! try different bin size
    do i = 1, maxBinSize
       nBin = N/i
       do j = 1, nBin
          sum1 = sum(data((j-1)*i+1:j*i))
          sum2 = sum(abs(sgn((j-1)*i+1:j*i)))
          if (sum2 .le. epsilon()) then
             call DQMC_Warning("Warning: sign problem", 1)
          else
             r = sum1/sum2
             var(i) = var(i) + (r - ratio)**2
             avg(i) = avg(i) + r
          end if
       end do       

       if (nBin .gt. 1) then
          var(i) = var(i)/(nBin-1)
       end if
       arv(i) = avg(i)/nBin
    end do

    if (abs(var(1)) .gt. epsilon()) then
       r = ONE/var(1)
       do i = 1, maxBinSize
          var(i) = HALF*i*var(i)*r
       end do
    else
       call DQMC_Warning("Warning: sign problem", 2)
    end if
    
  end subroutine DQMC_AutoCorrTime

end module DQMC_Phy0
