!gfortran -O3 -o gdd gdd.f90 -I/Users/djlorenz/Research/netcdf/include -I/Users/djlorenz/lib -L/Users/djlorenz/Research/netcdf/lib -L/Users/djlorenz/lib -lnetcdff -lnetcdf -llapack -lblas -lmylib
program gdd0

  use netcdf
  use constants
  use read_data
  use nc_write
  implicit none

  integer, parameter :: year1 = 1901, year2 = 2011
  integer, parameter :: nvar=6

  real, parameter :: minlon = -173.0, maxlon = -48.0
  real, parameter :: minlat = 10.0, maxlat = 80.0

  real, allocatable :: x(:,:,:,:),y(:,:),xmon(:,:,:,:)
  real :: ym(nvar,12),ymean(12),ystd(12),yastd(12),tstd(12)
  real :: gdd

  real, parameter :: days(12) = &
       [31.0,28.25,31.0,30.0,31.0,30.0,31.0,31.0,30.0,31.0,30.0,31.0]

  integer :: ilon,ilat,iyear,imonth,im,ip,i0
  real :: Lm,L0,Lp,am,ap

  ! get prcp, Tmax Tmin, vap:
  file_expr = '/Users/djlorenz/Research/jack/data/tmax.nc '// &
       '/Users/djlorenz/Research/jack/data/tmin.nc '
  vars = 'tmax tmin '
  call read_nc(minlon,maxlon,.false.,minlat,maxlat,missing,missing, &
       [all,all,all,year1],[all,all,all,year2],1,'month')
  allocate(x(nlon,nlat,12,nyear),y(12,nyear),xmon(nlon,nlat,nvar,12))
  x = merge(0.5*(x1(:,:,1,:,:) + x1(:,:,2,:,:)),missing, &
       x1(:,:,1,:,:)/=missing .and. x1(:,:,2,:,:)/=missing)
  deallocate(x1)

  do ilat=1,nlat
!     print*, ilat,nlat
     do ilon=1,nlon
        if(x(ilon,ilat,1,1)/=missing) then
           y = x(ilon,ilat,:,:)
           !
           ymean = sum(y,dim=2)/nyear
           do iyear=1,nyear
              y(:,iyear) = y(:,iyear) - ymean
           enddo
           ystd = sqrt(sum(y**2,dim=2)/(nyear - 1))
           ! one way:
!           call tstd_annual(ymean,yastd,2)
!           tstd = sqrt(0.176175791289*days*ystd**2 + yastd**2)
           ! another way:
           do imonth=1,12
              i0 = imonth
              im = modulo(i0-2,12)+1
              ip = modulo(i0,12)+1
              L0 = days(i0); Lm = days(im); Lp = days(ip)
              am = 2.0*(ymean(i0) - ymean(im))/(L0 + Lm)
              ap = 2.0*(ymean(ip) - ymean(i0))/(Lp + L0)
              yastd(imonth) = L0*sqrt((5.*ap**2 + 6.*ap*am + 5.*am**2)/192.)
           enddo
           tstd = sqrt(0.178211466801*days*ystd**2 + yastd**2)
           ! end ways
           do imonth=1,12
              ym(1,imonth) = days(imonth)*gdd(ymean(imonth),tstd(imonth),0.0)
              ym(2,imonth) = days(imonth)*gdd(ymean(imonth),tstd(imonth),5.0)
              ym(3,imonth) = days(imonth)*max(ymean(imonth),0.0)
              ym(4,imonth) = days(imonth)*max(ymean(imonth)-5.0,0.0)
           enddo
           ym(5,:) = ystd
           ym(6,:) = yastd
           !
           xmon(ilon,ilat,:,:) = ym
        else
           xmon(ilon,ilat,:,:) = missing
        endif
     enddo
  enddo

  call ncopen(1,'gdd_obs.nc','gdd0 gdd5 gdd0m gdd5m tmstd tastd',lons=lons,lats=lats, &
       tunit='M',title='growing degree days estimated from monthly values')
  do imonth=1,12
     call ncwrite(1,r1=xmon(:,:,:,imonth))
  enddo
  call ncclose(1)

end program gdd0

function gdd(Tmean,sigma,T0)

  use constants
  implicit none

  real :: gdd
  real :: Tmean,sigma,T0

!  gdd = sigma/sqrt(2.*pi)*exp(-0.5*((Tmean-T0)/sigma)**2) + &
!       0.5*(Tmean - T0)*(1.0 + erf((Tmean - T0)/(sqrt(2.)*sigma)))
  gdd = sigma/sqrt(2.*pi)*exp(-0.5*((Tmean-T0)/sigma)**2) + &
       0.5*(Tmean - T0)*erfc((T0 - Tmean)/(sqrt(2.)*sigma))

end function gdd

!temperature st. dev. from annual cycle:
subroutine tstd_annual(xm,tstd,order)

  implicit none

  real :: xm(12),tstd(12)
  integer :: order
  real :: xd(365)

  integer, parameter :: days(12) = [31,28,31,30,31,30,31,31,30,31,30,31]
  integer, parameter :: sd(12) = [1,32,60,91,121,152,182,213,244,274,305,335]
  integer, parameter :: ed(12) = [31,59,90,120,151,181,212,243,273,304,334,365]

  integer :: iday,imonth

  if(order==1) then
     call mon1daily(xm,xd)  
  elseif(order==2) then
     call mon2daily(xm,xd)  
  else
     print*, 'ERROR: tstd_annual: order must be 1 or 2'
     stop
  endif


  tstd = 0.0
  do imonth=1,12
     do iday=sd(imonth),ed(imonth)
        tstd(imonth) = tstd(imonth) + (xd(iday) - xm(imonth))**2
     enddo
     tstd(imonth) = sqrt(tstd(imonth)/(days(imonth)-1))
  enddo

end subroutine tstd_annual

! minimize second derivative:
subroutine mon2daily(xm,xd)

  implicit none

  real :: xm(12),xd(365)

  real :: a(377,377),x(377)

  integer :: j,k,imonth
  integer :: ipiv(377),info

  integer, parameter :: days(12) = [31,28,31,30,31,30,31,31,30,31,30,31]
  integer, parameter :: sd(12) = [1,32,60,91,121,152,182,213,244,274,305,335]
  integer, parameter :: ed(12) = [31,59,90,120,151,181,212,243,273,304,334,365]

  a = 0.0; x = 0.0
  do j=1,365
     a(j,j) = 6.0
     k = modulo(j,365) + 1
     a(j,k) = -4.0
     k = modulo(j-2,365) + 1
     a(j,k) = -4.0
     k = modulo(j+1,365) + 1
     a(j,k) = 1.0
     k = modulo(j-3,365) + 1
     a(j,k) = 1.0
  enddo
  do imonth=1,12
     a(sd(imonth):ed(imonth),365+imonth) = 1.0
     a(365+imonth,sd(imonth):ed(imonth)) = 1.0
  enddo

  x(366:377) = days*xm

  call sgesv(377,1,a,377,ipiv,x,377,info)

  xd = x(1:365)

end subroutine mon2daily

! minimize second derivative:
subroutine mon1daily(xm,xd)

  implicit none

  real :: xm(12),xd(365)

  real :: a(377,377),x(377)

  integer :: j,k,imonth
  integer :: ipiv(377),info

  integer, parameter :: days(12) = [31,28,31,30,31,30,31,31,30,31,30,31]
  integer, parameter :: sd(12) = [1,32,60,91,121,152,182,213,244,274,305,335]
  integer, parameter :: ed(12) = [31,59,90,120,151,181,212,243,273,304,334,365]

  a = 0.0; x = 0.0
  do j=1,365
     a(j,j) = 2.0
     k = modulo(j,365) + 1
     a(j,k) = -1.0
     k = modulo(j-2,365) + 1
     a(j,k) = -1.0
  enddo
  do imonth=1,12
     a(sd(imonth):ed(imonth),365+imonth) = 1.0
     a(365+imonth,sd(imonth):ed(imonth)) = 1.0
  enddo

  x(366:377) = days*xm

  call sgesv(377,1,a,377,ipiv,x,377,info)

  xd = x(1:365)

end subroutine mon1daily
