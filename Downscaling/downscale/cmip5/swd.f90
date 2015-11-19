!gfortran -O3 -o swd swd.f90 -I/Users/djlorenz/Research/netcdf/include -I/Users/djlorenz/lib -L/Users/djlorenz/Research/netcdf/lib -L/Users/djlorenz/lib -lnetcdff -lnetcdf  -llapack -lblas -lmylib
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

  integer :: nlonc,nlatc,nyearc
  real, allocatable :: lonsc(:),latsc(:),xave(:,:,:)

  integer :: nlonf,nlatf
  real, allocatable :: lonsf(:),latsf(:),xfave(:,:,:)
  real, allocatable :: xf(:,:,:),xmiss(:,:,:)
  real :: temp
  integer*2, allocatable :: ixf(:,:)

  real, parameter :: sf=0.1, ao=0.0
  
  character (len=400) :: fname
  integer :: imodel,iscen,iyear,imonth,ivar,ilon,ilat

  ! get CRU wind for landsea mask:
  file_expr = '/Users/djlorenz/Research/jack/data/wind.nc '
  vars = 'wind '
  call read_nc(minlon,maxlon,.false.,minlat,maxlat,missing,missing, &
       [all,all,all,2001],[all,all,all,2001],1,'month')
  allocate(xmiss(nlon,nlat,12))
  xmiss(:,:,:) = x1(:,:,1,:,1)
  deallocate(x1)

  ! get observed hi-res climatology:
  ! get prcp, Tmax Tmin, vap:
  file_expr = '/Users/djlorenz/Research/jack/data/srbmean.nc '
  vars = 'swd '
  call read_nc(minlon,maxlon,.false.,minlat,maxlat,missing,missing, &
       [all,all,all,2001],[all,all,all,2001],1,'month')
  nlonf = nlon; nlatf = nlat
  allocate(lonsf(nlonf),latsf(nlatf))
  lonsf = lons; latsf = lats
  allocate(xfave(nlonf,nlatf,12))
  xfave = x1(:,:,1,:,1)
  deallocate(x1)
  allocate(xf(nlonf,nlatf,12),ixf(nlonf,nlatf))
  
  ! do models:
  do imodel = 1,nmodel
     print*, imodel,trim(models(imodel))
     
     file_expr = '/Users/djlorenz/Research/jack/cmip5jack/historical/'// &
          trim(models(imodel))//'/rsds*.nc'
     vars = 'rsds'
     call read_nc(minlon,maxlon,.false.,minlat,maxlat,missing,missing, &
          [all,all,all,year1],[all,all,all,year2],1,'month')
     nlonc = nlon; nlatc = nlat; nyearc = nyear
     allocate(lonsc(nlonc),latsc(nlatc))
     lonsc = lons; latsc = lats
     allocate(xave(nlonc,nlatc,12))

     ! check sign:
     temp = sum(x1(:,:,1,:,1),x1(:,:,1,:,1)/=missing)
     if(temp < 0.) then
        print*, 'FLIP SIGN 1'
        stop
     endif

     xave = sum(x1(:,:,1,:,:),mask=x1(:,:,1,:,:)/=missing,dim=4) / count(x1(:,:,1,:,:)/=missing,dim=4)
     
     if(models(imodel) == 'HadGEM2-ES') then
        if(any(x1(:,:,1,12,nyearc)==missing)) then
           ! dec 2005 is missing
           x1(:,:,1,12,nyearc) = merge( (x1(:,:,1,11,nyearc) &
                / xave(:,:,11)) * xave(:,:,12), 0.0, xave(:,:,11)/=0.0)
        endif
     endif
     if(any(x1(:,:,1,:,:)==missing)) then
        print*, 'ERROR: missing values'
        !do iyear=1,nyear
        !   do imonth=1,12
        !      do ivar=1,2
        !         do ilat=1,nlat
        !            do ilon=1,nlon
        !               if(x1(ilon,ilat,imonth,iyear)==missing) then
        !                  print*, ilon,ilat,imonth,iyear+year1 - 1
        !               endif
        !            enddo
        !         enddo
        !      enddo
        !   enddo
        !enddo
        stop
     endif
     
     fname = '/Users/djlorenz/Research/jack/final/cmip5/historical/'//trim(models(imodel))
     call system('mkdir -p '//trim(fname))
     fname = '/Users/djlorenz/Research/jack/final/cmip5/historical/'//trim(models(imodel))// &
          '/swd.nc'

     call ncopen(1,trim(fname),'swd ', &
          vunits='Wm-2',vnames='mean monthly downward shortwave radiation at the surface;', &
          lons=lonsf,lats=latsf,title='Downscaled CMIP5 models, monthly mean.', &
          deflate=1,sf=[sf],ao=[ao],syear=year1,tunit='M')

     do iyear = year1,year2

        do imonth = 1,12
           x1(:,:,1,imonth,iyear-year1+1) = merge(x1(:,:,1,imonth,iyear-year1+1) / &
                xave(:,:,imonth), 0.0, xave(:,:,imonth) /= 0.0)
        enddo

        call bilinear(lonsc,latsc,nlonc,nlatc,lonsf,latsf,nlonf,nlatf, &
             x1(:,:,1,:,iyear-year1+1),xf,12,0)
        
        xf = merge(xf * xfave, missing, xfave /= missing .and. xmiss /= missing)
        
        do imonth = 1,12
           do ilat = 1,nlatf
              do ilon = 1,nlonf
                 if( xf(ilon,ilat,imonth) /= missing ) then
                    temp = (xf(ilon,ilat,imonth) - ao)/sf
                    if(abs(temp)>32767.) then
                       print*, 'ERROR: conversion to 2 byte: out of range. ',temp
                       stop
                    endif
                    ixf(ilon,ilat) = nint(temp)
                 else
                    ixf(ilon,ilat) = missing_short
                 endif
              enddo
           enddo
           call ncwrite(1,i1=ixf)
        enddo
        
     enddo
     call ncclose(1)
     
     ! future scenarios:
     do iscen = 1,2
        file_expr = '/Users/djlorenz/Research/jack/cmip5jack/'//trim(scens(iscen))//'/'// &
             trim(models(imodel))//'/rsds*.nc'
        vars = 'rsds'
        call read_nc(minlon,maxlon,.false.,minlat,maxlat,missing,missing, &
             [all,all,all,year11],[all,all,all,year22],1,'month')

        ! check sign:
        temp = sum(x1(:,:,1,:,1),x1(:,:,1,:,1)/=missing)
        if(temp < 0.) then
           print*, 'FLIP SIGN 1',trim(scens(iscen))
           stop
        endif

        !if(iscen==1 .and. models(imodel) == 'ACCESS1-3') then
        !   ! july 2050 tasmin is missing
        !   x1(:,:,2,7,45) = 0.5*(x1(:,:,2,6,45) + x1(:,:,2,8,45) &
        !        - xave(:,:,2,6) - xave(:,:,2,8)) + xave(:,:,2,7)
        !endif
        if(any(x1(:,:,1,:,:)==missing)) then
           print*, 'ERROR: missing values, '//trim(scens(iscen))
           !do iyear=1,nyear
           !   do imonth=1,12
           !      do ivar=1,2
           !         do ilat=1,nlat
           !            do ilon=1,nlon
           !               if(x1(ilon,ilat,imonth,iyear)==missing) then
           !                  print*, ilon,ilat,imonth,iyear+year11 - 1
           !               endif
           !            enddo
           !         enddo
           !      enddo
           !   enddo
           !enddo
           stop
        endif
        
        fname = '/Users/djlorenz/Research/jack/final/cmip5/'//trim(scens(iscen))//'/'// &
             trim(models(imodel))
        call system('mkdir -p '//trim(fname))
        fname = '/Users/djlorenz/Research/jack/final/cmip5/'//trim(scens(iscen))//'/'// &
             trim(models(imodel))//'/swd.nc'
        
        call ncopen(1,trim(fname),'swd ', &
             vunits='Wm-2',vnames='mean monthly downward shortwave radiation at the surface;', &
             lons=lonsf,lats=latsf,title='Downscaled CMIP5 models, monthly mean.', &
             deflate=1,sf=[sf],ao=[ao],syear=year11,tunit='M')
        
        do iyear = year11,year22
           
           do imonth = 1,12
              x1(:,:,1,imonth,iyear-year11+1) = merge(x1(:,:,1,imonth,iyear-year11+1) / &
                   xave(:,:,imonth), 0.0, xave(:,:,imonth) /= 0.0)
           enddo
           
           call bilinear(lonsc,latsc,nlonc,nlatc,lonsf,latsf,nlonf,nlatf, &
                x1(:,:,1,:,iyear-year11+1),xf,12,0)
           
           xf = merge(xf * xfave, missing, xfave /= missing .and. xmiss /= missing)
           
           do imonth = 1,12
              do ilat = 1,nlatf
                 do ilon = 1,nlonf
                    if( xf(ilon,ilat,imonth) /= missing ) then
                       temp = (xf(ilon,ilat,imonth) - ao)/sf
                       if(abs(temp)>32767.) then
                          print*, 'ERROR: conversion to 2 byte: out of range. ',temp
                          stop
                       endif
                       ixf(ilon,ilat) = nint(temp)
                    else
                       ixf(ilon,ilat) = missing_short
                    endif
                 enddo
              enddo
              call ncwrite(1,i1=ixf)
           enddo
           
        enddo
        call ncclose(1)
        
     enddo

     deallocate(xave,lonsc,latsc)
     
  enddo

end program temp0
  
