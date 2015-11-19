!course monthly anomalies to high resolution full variables (i.e. not anomalies) 
!gfortran -O3 -o prcp_final prcp_final.f90 read_CCSM.o write_final.o -I/Users/djlorenz/Research/netcdf/include -I/Users/djlorenz/lib -L/Users/djlorenz/Research/netcdf/lib -L/Users/djlorenz/lib -lnetcdff -lnetcdf  -llapack -lblas -lmylib
program prcp0

  use netcdf
  use constants
  use read_data
  use read_CCSM
  use write_final
  implicit none

  integer, parameter :: year1 = 1901, year2 = 2011

  ! 3 == 1981-1990, 0 == 1951-1960
  integer, parameter :: decade1 = -2200, decade2 = 3, ndecade = decade2-decade1+1
  
  real, parameter :: minlon = -173.0, maxlon = -48.0
  real, parameter :: minlat = 10.0, maxlat = 80.0

  real, allocatable :: x(:,:,:),xmiss(:,:,:),dx(:,:,:)
  real, allocatable :: xmon(:,:,:,:)
  integer*2, allocatable :: ix(:,:,:)
  real, parameter :: sf=0.1, ao=0.0
  real :: temp

  integer :: idecade,ilat,ilon,imonth

  ! get grid
  call get_CCSM(minlon,maxlon,minlat,maxlat,0,0,'PRECT')

  ! get CRU wind for landsea mask:
  file_expr = '/Users/djlorenz/Research/jack/data/wind.nc '
  vars = 'wind '
  call read_nc(minlon,maxlon,.false.,minlat,maxlat,missing,missing, &
       [all,all,all,2001],[all,all,all,2001],1,'month')
  allocate(xmiss(nlon,nlat,12))
  xmiss = x1(:,:,1,:,1)
  deallocate(x1)

  ! get observed hi-res climatology:
  ! get prcp, Tmax Tmin, vap:
  file_expr = '/Users/djlorenz/Research/jack/data/prcp.nc '
  vars = 'prcp '
  call read_nc(minlon,maxlon,.false.,minlat,maxlat,missing,missing, &
       [all,all,all,year1],[all,all,all,year2],1,'month')
  allocate(x(nlon,nlat,12))
  x = sum(x1(:,:,1,:,:),dim=4)/nyear
  x = merge(x,missing,abs(x)<=1e30)
  deallocate(x1)

  allocate(xmon(nlonc,nlatc,ndecade,12),dx(nlon,nlat,12),ix(nlon,nlat,12))

  open(1,file='prcpL_monthly.dat',access='stream',form='unformatted')
  read(1) xmon
  close(1)

  call ncopen(1,'/Users/djlorenz/Research/jack/final/CCSM/prcp.nc','prcp ', &
       vunits='mm',vnames='mean monthly precipitation accumulation;', &
       lons=lons,lats=lats,title='Downscaled TraCE CCSM runs. Monthly means by decade.', &
       deflate=1,sf=[sf],ao=[ao])

  do idecade=1,ndecade
     print*, idecade
     call bilinear(lonsc,latsc,nlonc,nlatc,lons,lats,nlon,nlat, &
          xmon(:,:,idecade,:),dx,12,0)
     dx = merge(x*dx,missing,x/=missing .and. dx/=missing .and. xmiss/=missing)
     do imonth=1,12
        do ilat=1,nlat
           do ilon=1,nlon
              if(dx(ilon,ilat,imonth)/=missing) then
                 temp = (dx(ilon,ilat,imonth)-ao)/sf
                 if(abs(temp)>32767.) then
                    print*, 'ERROR: conversion to 2 byte: out of range. ',temp
                    stop
                 endif
                 ix(ilon,ilat,imonth) = nint(temp)
              else
                 ix(ilon,ilat,imonth) = missing_short
              endif
           enddo
        enddo
     enddo
     call ncwrite(1,im=ix)
  enddo
  call ncclose(1)

end program prcp0
