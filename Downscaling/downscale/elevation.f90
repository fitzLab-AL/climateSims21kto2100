!gfortran -O3 -o elevation elevation.f90 -I/Users/djlorenz/Research/netcdf/include -I/Users/djlorenz/lib -L/Users/djlorenz/Research/netcdf/lib -L/Users/djlorenz/lib -lnetcdff -lnetcdf  -llapack -lblas -lmylib
program elevation

  use netcdf
  use constants
  use read_data
  use nc_write
  implicit none

  real, parameter :: minlon = -173.0, maxlon = -48.0
  real, parameter :: minlat = 10.0, maxlat = 80.0

  integer, parameter :: nlone=21600, nlate=10800
  real, parameter :: dl = 1./120.
  integer*2 :: ix(nlone)
  real :: lonse(nlone),latse(nlate),xin(nlone,nlate)
  real, allocatable :: x(:,:,:)

  integer :: ilon,ilat

  do ilon=1,nlone
     lonse(ilon) = -180. + dl/2. + dl*(ilon-1)
  enddo
  do ilat=1,nlate
     latse(ilat) = dl/2. + dl*(ilat-1)
  enddo

  open(1,file='/Users/djlorenz/Research/downscaling2/elevation/elevation.dat', &
       access='stream',form='unformatted')
  do ilat=1,nlate
     read(1) ix
     xin(:,ilat) = merge(real(ix),0.0,ix/=-500)
  enddo
  close(1)

! get final grid:
  file_expr = '/Users/djlorenz/Research/jack/data/wind.nc '
  vars = 'wind '
  call read_nc(minlon,maxlon,.false.,minlat,maxlat,missing,missing, &
       [all,all,all,2001],[all,all,all,2001],1,'month')
  allocate(x(nlon,nlat,2))

  call box(lonse,latse,nlone,nlate,lons,lats,nlon,nlat,xin,x(:,:,1),1,0)

  x(:,:,2) = 101.3*(1.0 - 0.0065*x(:,:,1)/293.)**5.26

  call ncopen(1,'elevation.nc','e ps',vunits='m kPa', &
       vnames='elevation; surface pressure;',lons=lons,lats=lats)
  call ncwrite(1,r1=x)
  call ncclose(1)



end program elevation
