!course monthly anomalies to high resolution full variables (i.e. not anomalies) 
!gfortran -O3 -o ps_final ps_final.f90 -I/Users/djlorenz/Research/netcdf/include -I/Users/djlorenz/lib -L/Users/djlorenz/Research/netcdf/lib -L/Users/djlorenz/lib -lnetcdff -lnetcdf  -llapack -lblas -lmylib
program ps0

  use netcdf
  use constants
  use read_data
  use nc_write
  implicit none

  ! 3 == 1981-1990, 0 == 1951-1960
  integer, parameter :: time1 = -21, time2 = 0, ntime = time2-time1+1
  
  real, parameter :: minlon = -173.0, maxlon = -48.0
  real, parameter :: minlat = 10.0, maxlat = 80.0

  integer :: nlonc,nlatc
  real, allocatable :: lonsc(:),latsc(:)

  real, allocatable :: e(:,:),xmiss(:,:),dx(:,:,:)
  real, allocatable :: ice(:,:,:),ice2(:,:)
  integer*2, allocatable :: ix(:,:,:)
  real, parameter :: sf=0.01, ao=0.0
  real :: temp

  integer :: itime,ilat,ilon,imonth

  ! get ice topo:
  file_expr = '/Users/djlorenz/Research/jack/downscale/ECBilt-CLIO/ice_topo/ice_topo.nc '
  vars = 't '
  call read_nc(minlon,maxlon,.false.,minlat,maxlat,missing,missing, &
       [all,all,1,2001],[all,all,1,2001],1,'month')
  allocate(ice(nlon,nlat,nlev))
  ice = xm(:,:,:,1,1,1)
  deallocate(xm)
  nlonc = nlon; nlatc = nlat
  allocate(lonsc(nlonc),latsc(nlatc))
  lonsc = lons; latsc = lats

  do itime=1,nlev
     ice(:,:,itime) = ice(:,:,itime) - ice(:,:,nlev)
  enddo


  ! get CRU wind for landsea mask:
  file_expr = '/Users/djlorenz/Research/jack/data/wind.nc '
  vars = 'wind '
  call read_nc(minlon,maxlon,.false.,minlat,maxlat,missing,missing, &
       [all,all,all,2001],[all,all,all,2001],1,'month')
  allocate(xmiss(nlon,nlat))
  xmiss = x1(:,:,1,1,1)
  deallocate(x1)

  ! get observed hi-res climatology:
  file_expr = '/Users/djlorenz/Research/jack/downscale/elevation.nc '
  vars = 'e '
  call read_nc(minlon,maxlon,.false.,minlat,maxlat,missing,missing, &
       [all,all,1,2001],[all,all,1,2001],1,'month')
  allocate(e(nlon,nlat),ice2(nlon,nlat),ix(nlon,nlat,ntime))
  e(:,:) = x1(:,:,1,1,1)
  deallocate(x1)

  call ncopen(1,'/Users/djlorenz/Research/jack/final/ECBilt/ps.nc','ps ', &
       vunits='kPa',vnames='surface pressure;', &
       lons=lons,lats=lats,nlev=ntime,minlev=real(time1), &
       title='Downscaled ECBilt runs. Surface pressure every thousand years.', &
       deflate=1,sf=[sf],ao=[ao])


  do itime=1,ntime
     print*, itime
     call bilinear(lonsc,latsc,nlonc,nlatc,lons,lats,nlon,nlat, &
          ice(:,:,itime),ice2,1,0)
     ice2 = merge(101.3*(1.0 - 0.0065*(ice2+e)/293.)**5.26,missing, &
          ice2/=missing .and. xmiss/=missing)
     do ilat=1,nlat
        do ilon=1,nlon
           if(ice2(ilon,ilat)/=missing) then
              temp = (ice2(ilon,ilat)-ao)/sf
              if(abs(temp)>32767.) then
                 print*, 'ERROR: conversion to 2 byte: out of range. ',temp
                 stop
              endif
              ix(ilon,ilat,itime) = nint(temp)
           else
              ix(ilon,ilat,itime) = missing_short
           endif
        enddo
     enddo
  enddo
  call ncwrite(1,im=ix)
  call ncclose(1)

end program ps0
