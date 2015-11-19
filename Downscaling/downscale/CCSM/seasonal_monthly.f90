!gfortran -O3 -c seasonal_monthly.f90
module seasonal_monthly

  implicit none

  private
   public sm_power,f_power1,f_power2
  public sm_meanlog,f_meanlog1,f_meanlog2

  real :: xmean(12),xsea(4)
  integer :: seas(12) = [1,1,2,2,2,3,3,3,4,4,4,1]

contains

  subroutine sm_factor(mean,xs,xm,order)

    integer :: order

    real, parameter :: tol = 1d-6

    real :: mean(12),xs(4),xm(12)

    real :: x(16),fvec(16),fjac(16,16)
    integer :: info
    integer, parameter :: lwa = (16*(16+13))/2
    real :: wa(lwa)

    if(any(xs<=0.)) then
       print*, 'ERROR: sm_factor: xs <= 0'
       print*, 'xs = ',xs
       return
    endif

    xmean = mean; xsea = xs

    x(1:2) = xs(1); x(12) = xs(1); x(3:5) = xs(2); x(6:8) = xs(3); x(9:11) = xs(4)
    x(1:12) = log(x(1:12))
    x(13:16) = 0.0

    if(order==1) then
       call hybrj1(f_factor1,16,x,fvec,fjac,16,tol,info,wa,lwa)
       if(info/=1 .and. info/=3) then
          x(13:16) = 0.01
          call hybrj1(f_factor1,16,x,fvec,fjac,16,tol,info,wa,lwa)
       endif
    elseif(order==2) then
       call hybrj1(f_factor2,16,x,fvec,fjac,16,tol,info,wa,lwa)
       if(info/=1 .and. info/=3) then
          x(13:16) = 0.01
          call hybrj1(f_factor2,16,x,fvec,fjac,16,tol,info,wa,lwa)
       endif
    endif

    !print*, x(13:16)
!    if(info/=1 .and. info/=3) then
!       if(any(abs(xs-1.0)>0.001)) then
!          print*, info
!          print*,  exp(x(1:12))
!          print*, xs
!       endif
!    endif

    xm = exp(x(1:12))

  end subroutine sm_factor

  subroutine sm_meanlog(xs,xm,order)

    integer :: order

    real, parameter :: tol = 1d-6

    real :: xs(4),xm(12)

    real :: x(16),fvec(16),fjac(16,16)
    integer :: info
    integer, parameter :: lwa = (16*(16+13))/2
    real :: wa(lwa)

    if(any(xs<=0.)) then
       print*, 'ERROR: sm_factor: xs <= 0'
       print*, 'xs = ',xs
       return
    endif

    xsea = xs

    x(1:2) = xs(1); x(12) = xs(1); x(3:5) = xs(2); x(6:8) = xs(3); x(9:11) = xs(4)
    x(1:12) = log(x(1:12))
    x(13:16) = 0.0

    if(order == 1) then
       call hybrj1(f_meanlog1,16,x,fvec,fjac,16,tol,info,wa,lwa)
       if(info/=1 .and. info/=3) then
          x(13:16) = 0.01
          call hybrj1(f_meanlog1,16,x,fvec,fjac,16,tol,info,wa,lwa)
       endif
    elseif(order==2) then
       call hybrj1(f_meanlog2,16,x,fvec,fjac,16,tol,info,wa,lwa)
       if(info/=1 .and. info/=3) then
          x(13:16) = 0.01
          call hybrj1(f_meanlog2,16,x,fvec,fjac,16,tol,info,wa,lwa)
       endif
    endif

    xm = exp(x(1:12))

  end subroutine sm_meanlog

  subroutine f_factor1(n,x,fvec,fjac,ldfjac,iflag)

    integer :: n,ldfjac,iflag
    real :: x(n),fvec(n),fjac(ldfjac,n)

    real :: lambda(4)
    integer :: j,k,l,isea

    lambda = x(13:16)

    if(iflag == 1) then
       do j=1,12
          isea = seas(j)
          l = modulo(j,12) + 1
          k = modulo(j-2,12) + 1
          fvec(j) = 2.*x(j) - x(l) - x(k) + lambda(isea)*xmean(j)*exp(x(j))
       enddo
       fvec(13) = xmean(12)*exp(x(12)) + xmean(1)*exp(x(1)) + xmean(2)*exp(x(2)) - &
            (xmean(12) + xmean(1) + xmean(2))*xsea(1)       
       fvec(14) = xmean(3)*exp(x(3)) + xmean(4)*exp(x(4)) + xmean(5)*exp(x(5)) - &
            (xmean(3) + xmean(4) + xmean(5))*xsea(2)
       fvec(15) = xmean(6)*exp(x(6)) + xmean(7)*exp(x(7)) + xmean(8)*exp(x(8)) - &
            (xmean(6) + xmean(7) + xmean(8))*xsea(3)
       fvec(16) = xmean(9)*exp(x(9)) + xmean(10)*exp(x(10)) + xmean(11)*exp(x(11)) - &
            (xmean(9) + xmean(10) + xmean(11))*xsea(4)
    elseif(iflag==2) then
       fjac = 0.0
       do j=1,12
          isea = seas(j)
          l = modulo(j,12) + 1
          k = modulo(j-2,12) + 1
          fjac(j,j) = 2. + lambda(isea)*xmean(j)*exp(x(j)) ! derivative of exp
          fjac(j,l) = -1.0
          fjac(j,k) = -1.0
          fjac(j,12+isea) = xmean(j)*exp(x(j))
       enddo
       ! derivatives of exp below:
       fjac(13,12) = xmean(12)*exp(x(12))
       fjac(13,1) = xmean(1)*exp(x(1))
       fjac(13,2) = xmean(2)*exp(x(2))
       fjac(14,3) = xmean(3)*exp(x(3))
       fjac(14,4) = xmean(4)*exp(x(4))
       fjac(14,5) = xmean(5)*exp(x(5))
       fjac(15,6) = xmean(6)*exp(x(6))
       fjac(15,7) = xmean(7)*exp(x(7))
       fjac(15,8) = xmean(8)*exp(x(8))
       fjac(16,9) = xmean(9)*exp(x(9))
       fjac(16,10) = xmean(10)*exp(x(10))
       fjac(16,11) = xmean(11)*exp(x(11))
    else
       print*, 'ERROR: iflag = ',iflag
       stop
    endif

  end subroutine f_factor1

  subroutine f_factor2(n,x,fvec,fjac,ldfjac,iflag)

    integer :: n,ldfjac,iflag
    real :: x(n),fvec(n),fjac(ldfjac,n)

    real :: lambda(4)
    integer :: j,k,l,k2,l2,isea

    lambda = x(13:16)

    if(iflag == 1) then
       do j=1,12
          isea = seas(j)
          l = modulo(j,12) + 1
          k = modulo(j-2,12) + 1
          l2 = modulo(j+1,12) + 1
          k2 = modulo(j-3,12) + 1
          fvec(j) = 6.*x(j) - 4.*x(l) - 4.*x(k) + x(l2) + x(k2) + lambda(isea)*xmean(j)*exp(x(j))
       enddo
       fvec(13) = xmean(12)*exp(x(12)) + xmean(1)*exp(x(1)) + xmean(2)*exp(x(2)) - &
            (xmean(12) + xmean(1) + xmean(2))*xsea(1)       
       fvec(14) = xmean(3)*exp(x(3)) + xmean(4)*exp(x(4)) + xmean(5)*exp(x(5)) - &
            (xmean(3) + xmean(4) + xmean(5))*xsea(2)
       fvec(15) = xmean(6)*exp(x(6)) + xmean(7)*exp(x(7)) + xmean(8)*exp(x(8)) - &
            (xmean(6) + xmean(7) + xmean(8))*xsea(3)
       fvec(16) = xmean(9)*exp(x(9)) + xmean(10)*exp(x(10)) + xmean(11)*exp(x(11)) - &
            (xmean(9) + xmean(10) + xmean(11))*xsea(4)
    elseif(iflag==2) then
       fjac = 0.0
       do j=1,12
          isea = seas(j)
          l = modulo(j,12) + 1
          k = modulo(j-2,12) + 1
          l2 = modulo(j+1,12) + 1
          k2 = modulo(j-3,12) + 1
          fjac(j,j) = 6. + lambda(isea)*xmean(j)*exp(x(j)) ! derivative of exp
          fjac(j,l) = -4.0
          fjac(j,k) = -4.0
          fjac(j,l2) = 1.0
          fjac(j,k2) = 1.0
          fjac(j,12+isea) = xmean(j)*exp(x(j))
       enddo
       ! derivatives of exp below:
       fjac(13,12) = xmean(12)*exp(x(12))
       fjac(13,1) = xmean(1)*exp(x(1))
       fjac(13,2) = xmean(2)*exp(x(2))
       fjac(14,3) = xmean(3)*exp(x(3))
       fjac(14,4) = xmean(4)*exp(x(4))
       fjac(14,5) = xmean(5)*exp(x(5))
       fjac(15,6) = xmean(6)*exp(x(6))
       fjac(15,7) = xmean(7)*exp(x(7))
       fjac(15,8) = xmean(8)*exp(x(8))
       fjac(16,9) = xmean(9)*exp(x(9))
       fjac(16,10) = xmean(10)*exp(x(10))
       fjac(16,11) = xmean(11)*exp(x(11))
    else
       print*, 'ERROR: iflag = ',iflag
       stop
    endif

  end subroutine f_factor2

  subroutine f_meanlog1(n,x,fvec,fjac,ldfjac,iflag)

    integer :: n,ldfjac,iflag
    real :: x(n),fvec(n),fjac(ldfjac,n)

    real :: lambda(4)
    integer :: j,k,l,isea

    lambda = x(13:16)

    if(iflag == 1) then
       do j=1,12
          isea = seas(j)
          l = modulo(j,12) + 1
          k = modulo(j-2,12) + 1
          fvec(j) = 2.*x(j) - x(l) - x(k) + lambda(isea)*exp(x(j))
       enddo
       fvec(13) = exp(x(12)) + exp(x(1)) + exp(x(2)) - 3.0*xsea(1)       
       fvec(14) = exp(x(3)) + exp(x(4)) + exp(x(5)) - 3.0*xsea(2)
       fvec(15) = exp(x(6)) + exp(x(7)) + exp(x(8)) - 3.0*xsea(3)
       fvec(16) = exp(x(9)) + exp(x(10)) + exp(x(11)) - 3.0*xsea(4)
    elseif(iflag==2) then
       fjac = 0.0
       do j=1,12
          isea = seas(j)
          l = modulo(j,12) + 1
          k = modulo(j-2,12) + 1
          fjac(j,j) = 2. + lambda(isea)*exp(x(j)) ! derivative of exp
          fjac(j,l) = -1.0
          fjac(j,k) = -1.0
          fjac(j,12+isea) = exp(x(j))
       enddo
       ! derivatives of exp below:
       fjac(13,12) = exp(x(12))
       fjac(13,1) = exp(x(1))
       fjac(13,2) = exp(x(2))
       fjac(14,3) = exp(x(3))
       fjac(14,4) = exp(x(4))
       fjac(14,5) = exp(x(5))
       fjac(15,6) = exp(x(6))
       fjac(15,7) = exp(x(7))
       fjac(15,8) = exp(x(8))
       fjac(16,9) = exp(x(9))
       fjac(16,10) = exp(x(10))
       fjac(16,11) = exp(x(11))
    else
       print*, 'ERROR: iflag = ',iflag
       stop
    endif

  end subroutine f_meanlog1

  subroutine f_meanlog2(n,x,fvec,fjac,ldfjac,iflag)

    integer :: n,ldfjac,iflag
    real :: x(n),fvec(n),fjac(ldfjac,n)

    real :: lambda(4)
    integer :: j,k,l,k2,l2,isea

    lambda = x(13:16)

    if(iflag == 1) then
       do j=1,12
          isea = seas(j)
          l = modulo(j,12) + 1
          k = modulo(j-2,12) + 1
          l2 = modulo(j+1,12) + 1
          k2 = modulo(j-3,12) + 1
          fvec(j) = 6.*x(j) - 4.*x(l) - 4.*x(k) + x(l2) + x(k2) + lambda(isea)*exp(x(j))
       enddo
       fvec(13) = exp(x(12)) + exp(x(1)) + exp(x(2)) - 3.0*xsea(1)       
       fvec(14) = exp(x(3)) + exp(x(4)) + exp(x(5)) - 3.0*xsea(2)
       fvec(15) = exp(x(6)) + exp(x(7)) + exp(x(8)) - 3.0*xsea(3)
       fvec(16) = exp(x(9)) + exp(x(10)) + exp(x(11)) - 3.0*xsea(4)
    elseif(iflag==2) then
       fjac = 0.0
       do j=1,12
          isea = seas(j)
          l = modulo(j,12) + 1
          k = modulo(j-2,12) + 1
          l2 = modulo(j+1,12) + 1
          k2 = modulo(j-3,12) + 1
          fjac(j,j) = 6. + lambda(isea)*exp(x(j)) ! derivative of exp
          fjac(j,l) = -4.0
          fjac(j,k) = -4.0
          fjac(j,l2) = 1.0
          fjac(j,k2) = 1.0
          fjac(j,12+isea) = exp(x(j))
       enddo
       ! derivatives of exp below:
       fjac(13,12) = exp(x(12))
       fjac(13,1) = exp(x(1))
       fjac(13,2) = exp(x(2))
       fjac(14,3) = exp(x(3))
       fjac(14,4) = exp(x(4))
       fjac(14,5) = exp(x(5))
       fjac(15,6) = exp(x(6))
       fjac(15,7) = exp(x(7))
       fjac(15,8) = exp(x(8))
       fjac(16,9) = exp(x(9))
       fjac(16,10) = exp(x(10))
       fjac(16,11) = exp(x(11))
    else
       print*, 'ERROR: iflag = ',iflag
       stop
    endif

  end subroutine f_meanlog2

end module seasonal_monthly
