!test schem for filling in seasonal cycle when winter solar radiation is zero
!gfortran -O3 -o test_fill test_fill.f90 -I/Users/djlorenz/Research/netcdf/include -I/Users/djlorenz/lib -L/Users/djlorenz/Research/netcdf/lib -L/Users/djlorenz/lib -lnetcdff -lnetcdf  -llapack -lblas -lmylib
program test_fill

  use netcdf
  use constants
  use write_data
  implicit none

  integer, parameter :: n = 40
  real :: x(n,3)

  integer :: j

  x = missing
  do j=1,n
     x(j,1) = sin(2.*pi*j/20.)
!     x(j,1) = (j-20)*(j-37)*(j-5)
!     if(j<38 .or. j>40) then
     if(j<32 .and. j>0) then
        x(j,3) = x(j,1)
     endif


  enddo
  
  !x = cshift(x,dim=1,shift=-3)

  call fillin(n,x(:,3),x(:,2))


  call line(sy=x)









end program test_fill

! assume periodic
! assume only one region that's missing
subroutine fillin(n,x,xf)

  use constants
  implicit none

  integer :: n
  real :: x(n),xf(n)
  
  integer :: nf,j1,j2,j
  real :: x1(2),x2(2)
  real :: y(n),a(n,n)
  integer :: ipiv(n),info

  xf = x
  a = 0.0

  if(any(x==missing)) then
     if(x(1)==missing) then
        ! find j2 bound
        j = 2
        do 
           if(x(j)/=missing) exit
           j = j + 1
        enddo
        j2 = j
        ! find j1 bound
        j = n
        do 
           if(x(j)/=missing) exit
           j = j - 1
        enddo
        j1 = j
        nf = j2 - 1 + n - j1
        x1(2) = x(j1-1)
        x1(1) = x(j1)
        x2(1) = x(j2)
        x2(2) = x(j2 + 1)
     else
        ! find j1 bound
        j = 2
        do 
           if(x(j)==missing) exit
           j = j + 1
        enddo
        j1 = j - 1
        ! find j2 bound
        do j = j1+1,n
           if(x(j)/=missing) exit
        enddo
        if(j<=n) then
           j2 = j
           nf = j2 - j1 - 1
        else
           j2 = 1
           nf = n - j1
        endif
        if(j1 > 1) then
           x1(2) = x(j1-1)
        else
           x1(2) = x(n)
        endif
        x1(1) = x(j1)
        x2(1) = x(j2)
        if(j2 < n) then
           x2(2) = x(j2 + 1)
        else
           x2(2) = x(1)
        endif
     endif
  endif

  if(nf == 1) then
     y(1) = (4.0*x1(1) - 1.0*x1(2) + 4.0*x2(1) - 1.0*x2(2))/6.0
  else
     y = 0.0
     do j=1,nf
        a(j,j) = 6.0
        if(j+1<=nf) then
           a(j,j+1) = -4.0
        else
           y(j) = y(j) + 4.0*x2(1)
        endif
        if(j-1>=1) then
           a(j,j-1) = -4.0
        else
           y(j) = y(j) + 4.0*x1(1)
        endif
        if(j+2<=nf) then
           a(j,j+2) = 1.0
        else
           if(j==nf) then
              y(j) = y(j) - 1.0*x2(2)
           else
              y(j) = y(j) - 1.0*x2(1)
           endif
        endif
        if(j>=3) then
           a(j,j-2) = 1.0
        else
           if(j==1) then
              y(j) = y(j) - 1.0*x1(2)
           else
              y(j) = y(j) - 1.0*x1(1)
           endif
        endif
     enddo
     call sgesv(nf,1,a,n,ipiv,y,n,info)
  endif

  if(j1<j2) then
     do j=j1+1,j2-1
        xf(j) = y(j-j1)
     enddo
  else
     do j=j1+1,n
        xf(j) = y(j-j1)
     enddo
     do j=1,j2-1
        xf(j) = y(n-j1+j)
     enddo
  endif


end subroutine fillin
