!gfortran -O3 -o gdd gdd.f90 -I/Users/djlorenz/Research/netcdf/include -I/Users/djlorenz/lib -L/Users/djlorenz/Research/netcdf/lib -L/Users/djlorenz/lib -lnetcdff -lnetcdf  -llapack -lblas -lmylib
program temp0

  use netcdf
  use constants
  use read_data
  use nc_write
  implicit none

  integer, parameter :: year1 = 1950, year2 = 2005

  ! future scenarios:
  integer, parameter :: year11 = 2006, year22 = 2100
  character (len = 40) :: scens(2) = [character (len = 40) :: 'rcp45','rcp85']

  real, parameter :: minlon = -173.0, maxlon = -48.0
  real, parameter :: minlat = 10.0, maxlat = 80.0

  integer, parameter :: nmodel = 12
  character (len = 60) :: models(nmodel) = [character (len=60) :: &
       'ACCESS1-3','CanESM2','CESM1-CAM5','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-CM3', &
       'GISS-E2-R','HadGEM2-ES','inmcm4','IPSL-CM5A-MR','MIROC5','MRI-CGCM3']
  
  real, parameter :: days(12) = [31.,28.25,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.]
  
  real, allocatable :: tmstd(:,:,:),T(:,:,:,:)
  real :: temp,gdd,g,Tmean(0:13),tastd(12),tstd(12)

  integer*2, allocatable :: ix(:,:,:,:)
  real, parameter :: sf=0.1, ao=0.0

  character (len=400) :: fname
  integer :: imodel,iscen,ilat,ilon,imonth,iyear,imp,imm,ivar,im,ip,i0
  real :: Lm,L0,Lp,am,ap

  
  file_expr = '/Users/djlorenz/Research/jack/downscale/CCSM/gdd_obs.nc '
  vars = 'tmstd '
  call read_nc(missing,missing,.false.,missing,missing,missing,missing, &
       [all,all,all,2001],[all,all,all,2001],1,'month')
  allocate(tmstd(nlon,nlat,12),ix(nlon,nlat,2,12))
  tmstd = x1(:,:,1,:,1)
  deallocate(x1)

  ! do models:
  do imodel = 1,nmodel
     print*, imodel,trim(models(imodel))
     
     file_expr = '/Users/djlorenz/Research/jack/final/cmip5/historical/'// &
          trim(models(imodel))//'/temp.nc'
     vars = 'tmax tmin'
     call read_nc(minlon,maxlon,.false.,minlat,maxlat,missing,missing, &
          [all,all,all,year1],[all,all,all,year2],1,'month')
     allocate(T(nlon,nlat,12,nyear))
     T = merge(0.5*sum(x1,dim=3),missing,x1(:,:,1,:,:)/=missing)
     deallocate(x1)

     fname = '/Users/djlorenz/Research/jack/final/cmip5/historical/'//trim(models(imodel))// &
          '/gdd.nc'
     call ncopen(1,trim(fname),'gdd0 gdd5', &
          vunits='C',vnames='growing degree days base 0C; growing degree days base 5C;', &
          lons=lons,lats=lats,title='Downscaled CMIP5 models, monthly mean.', &
          deflate=1,sf=[sf],ao=[ao],syear=year1,tunit='M')

     do iyear=1,nyear
        do ilat=1,nlat
           do ilon=1,nlon
              if(T(ilon,ilat,1,1)/=missing) then
                 Tmean(1:12) = T(ilon,ilat,:,iyear)
                 Tmean(0) = T(ilon,ilat,12,max(iyear-1,1))
                 Tmean(13) = T(ilon,ilat,1,min(iyear+1,nyear))
                 do imonth=1,12
                    i0 = imonth
                    im = modulo(i0-2,12)+1
                    ip = modulo(i0,12)+1
                    L0 = days(i0); Lm = days(im); Lp = days(ip)
                    am = 2.0*(Tmean(i0) - Tmean(i0-1))/(L0 + Lm)
                    ap = 2.0*(Tmean(i0+1) - Tmean(i0))/(Lp + L0)
                    tastd(imonth) = L0*sqrt((5.*ap**2 + 6.*ap*am + 5.*am**2)/192.)
                 enddo
                 tstd = sqrt(0.178211466801*days*tmstd(ilon,ilat,:)**2 + tastd**2)
                 do imonth=1,12
                    g = days(imonth)*gdd(Tmean(imonth),tstd(imonth),0.0)
                    temp = (g-ao)/sf
                    if(abs(temp)>32767.) then
                       print*, 'ERROR: conversion to 2 byte: out of range. ',temp
                       stop
                    endif
                    ix(ilon,ilat,1,imonth) = nint(temp)
                    g = days(imonth)*gdd(Tmean(imonth),tstd(imonth),5.0)
                    temp = (g-ao)/sf
                    if(abs(temp)>32767.) then
                       print*, 'ERROR: conversion to 2 byte: out of range. ',temp
                       stop
                    endif
                    ix(ilon,ilat,2,imonth) = nint(temp)
                 enddo
              else
                 ix(ilon,ilat,:,:) = missing_short
              endif
           enddo
        enddo

        do imonth=1,12
           call ncwrite(1,i1=ix(:,:,:,imonth))
        enddo

     enddo

     call ncclose(1)
     deallocate(T)

     ! future scenarios:
     do iscen = 1,2
        
        file_expr = '/Users/djlorenz/Research/jack/final/cmip5/'//trim(scens(iscen))//'/'// &
             trim(models(imodel))//'/temp.nc'
        vars = 'tmax tmin'
        call read_nc(minlon,maxlon,.false.,minlat,maxlat,missing,missing, &
             [all,all,all,year11],[all,all,all,year22],1,'month')
        allocate(T(nlon,nlat,12,nyear))
        T = merge(0.5*sum(x1,dim=3),missing,x1(:,:,1,:,:)/=missing)
        deallocate(x1)
        
        fname = '/Users/djlorenz/Research/jack/final/cmip5/'//trim(scens(iscen))//'/'// &
             trim(models(imodel))//'/gdd.nc'
        call ncopen(1,trim(fname),'gdd0 gdd5', &
             vunits='C',vnames='growing degree days base 0C; growing degree days base 5C;', &
             lons=lons,lats=lats,title='Downscaled CMIP5 models, monthly mean.', &
             deflate=1,sf=[sf],ao=[ao],syear=year11,tunit='M')
        
        do iyear=1,nyear
           do ilat=1,nlat
              do ilon=1,nlon
                 if(T(ilon,ilat,1,1)/=missing) then
                    Tmean(1:12) = T(ilon,ilat,:,iyear)
                    Tmean(0) = T(ilon,ilat,12,max(iyear-1,1))
                    Tmean(13) = T(ilon,ilat,1,min(iyear+1,nyear))
                    do imonth=1,12
                       i0 = imonth
                       im = modulo(i0-2,12)+1
                       ip = modulo(i0,12)+1
                       L0 = days(i0); Lm = days(im); Lp = days(ip)
                       am = 2.0*(Tmean(i0) - Tmean(i0-1))/(L0 + Lm)
                       ap = 2.0*(Tmean(i0+1) - Tmean(i0))/(Lp + L0)
                       tastd(imonth) = L0*sqrt((5.*ap**2 + 6.*ap*am + 5.*am**2)/192.)
                    enddo
                    tstd = sqrt(0.178211466801*days*tmstd(ilon,ilat,:)**2 + tastd**2)
                    do imonth=1,12
                       g = days(imonth)*gdd(Tmean(imonth),tstd(imonth),0.0)
                       temp = (g-ao)/sf
                       if(abs(temp)>32767.) then
                          print*, 'ERROR: conversion to 2 byte: out of range. ',temp
                          stop
                       endif
                       ix(ilon,ilat,1,imonth) = nint(temp)
                       g = days(imonth)*gdd(Tmean(imonth),tstd(imonth),5.0)
                       temp = (g-ao)/sf
                       if(abs(temp)>32767.) then
                          print*, 'ERROR: conversion to 2 byte: out of range. ',temp
                          stop
                       endif
                       ix(ilon,ilat,2,imonth) = nint(temp)
                    enddo
                 else
                    ix(ilon,ilat,:,:) = missing_short
                 endif
              enddo
           enddo
           
           do imonth=1,12
              call ncwrite(1,i1=ix(:,:,:,imonth))
           enddo
           
        enddo
        
        call ncclose(1)
        deallocate(T)
        
     enddo
     
  enddo




end program temp0

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
