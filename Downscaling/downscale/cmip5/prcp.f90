!gfortran -O3 -o prcp prcp.f90 -I/Users/djlorenz/Research/netcdf/include -I/Users/djlorenz/lib -L/Users/djlorenz/Research/netcdf/lib -L/Users/djlorenz/lib -lnetcdff -lnetcdf  -llapack -lblas -lmylib
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

  integer, parameter :: nquantile = year2 - year1 + 1
  real :: qs(2,nquantile)

  integer, parameter :: nmodel = 12
  character (len = 60) :: models(nmodel) = [character (len=60) :: &
       'ACCESS1-3','CanESM2','CESM1-CAM5','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-CM3', &
       'GISS-E2-R','HadGEM2-ES','inmcm4','IPSL-CM5A-MR','MIROC5','MRI-CGCM3']

  integer :: nlonc,nlatc,nyearc
  real, allocatable :: lonsc(:),latsc(:),xave(:,:,:)

  integer :: nlonf,nlatf,nyearf
  real, allocatable :: lonsf(:),latsf(:),xfave(:,:,:),xprcp(:,:,:,:)
  real, allocatable :: xbox(:,:,:,:),quantile(:,:,:,:,:)
  real, allocatable :: xf(:,:,:),xmiss(:,:,:)
  real :: temp,junk(year22-year11+1)
  integer*2, allocatable :: ixf(:,:)
  real, allocatable :: x1max(:,:,:,:)

  real, parameter :: sf=0.1, ao=0.0
  real, parameter :: maxprcp = 2100.

  character (len=400) :: fname
  integer :: imodel,iscen,iyear,imonth,ivar,ilon,ilat,j

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
  file_expr = '/Users/djlorenz/Research/jack/data/prcp.nc '
  vars = 'prcp '
  call read_nc(minlon,maxlon,.false.,minlat,maxlat,missing,missing, &
       [all,all,all,year1],[all,all,all,year2],1,'month')
  nlonf = nlon; nlatf = nlat; nyearf = nyear
  allocate(lonsf(nlonf),latsf(nlatf))
  lonsf = lons; latsf = lats
  allocate(xfave(nlonf,nlatf,12))
  xfave = merge(sum(x1(:,:,1,:,:),dim=4)/nyear, missing, .not.any(x1(:,:,1,:,:)==missing,dim=4))
  allocate(xprcp(nlonf,nlatf,12,nyearf))
  xprcp = x1(:,:,1,:,:)
  deallocate(x1)
  allocate(xf(nlonf,nlatf,12),ixf(nlonf,nlatf))
  
  ! do models:
  do imodel = 1,nmodel
     print*, imodel,trim(models(imodel))
     
     file_expr = '/Users/djlorenz/Research/jack/cmip5jack/historical/'// &
          trim(models(imodel))//'/pr*.nc'
     vars = 'pr'
     call read_nc(minlon,maxlon,.false.,minlat,maxlat,missing,missing, &
          [all,all,all,year1],[all,all,all,year2],1,'month')
     nlonc = nlon; nlatc = nlat; nyearc = nyear
     allocate(lonsc(nlonc),latsc(nlatc))
     lonsc = lons; latsc = lats
     allocate(xave(nlonc,nlatc,12))

     allocate(x1max(nlonc,nlatc,4,12))

     ! temporary xave in this case (for missing values):
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

     !x1max(:,:,1,:) = maxval(x1(:,:,1,:,:),dim=4)

     ! box average the observations:
     allocate(xbox(nlonc,nlatc,12,nyearf))
     call box(lonsf,latsf,nlonf,nlatf,lonsc,latsc,nlonc,nlatc,xprcp,xbox,12*nyearf,0)
!     call box_smooth2(lonsf,latsf,nlonf,nlatf,lonsc,latsc,nlonc,nlatc,xprcp,xbox, &
!          12*nyearf,1.0,0)
!     call bilinear(lonsf,latsf,nlonf,nlatf,lonsc,latsc,nlonc,nlatc,xprcp,xbox,12*nyearf,0)
     ! quantile adjust and save quantile info for future scenarios:
     allocate(quantile(2,nlonc,nlatc,12,nquantile))
     do imonth=1,12
        do ilat = 1,nlatc
           do ilon = 1,nlonc
              if(.not.any(xbox(ilon,ilat,imonth,:)==missing)) then
                 call quantile_adj(nyearc,x1(ilon,ilat,1,imonth,:), &
                      xbox(ilon,ilat,imonth,:),nquantile, &
                      quantile(:,ilon,ilat,imonth,:),1)
              else
                 x1(ilon,ilat,1,imonth,:) = missing
              endif
           enddo
        enddo
       ! print*, imonth,maxval(xbox(:,:,imonth,:))
     enddo

     !x1max(:,:,2,:) = maxval(x1(:,:,1,:,:),dim=4)

     !open(99,file='quantile.dat',access='stream',form='unformatted')
     !do imonth=1,12
     !   do ivar=1,2
     !      do j=1,nquantile
     !         write(99) quantile(ivar,:,:,imonth,j)
     !      enddo
     !   enddo
     !enddo
     !close(99)
     !call ctlwrite('quantile.ctl','q0 q1',lons=lonsc,lats=latsc,nlev=nquantile, &
     !     nt=12,tunit='M')
     !stop
     
     xave = merge( sum(xbox,dim=4,mask=xbox/=missing)/count(xbox/=missing,dim=4), missing, &
          count(xbox/=missing,dim=4) > 0)

     !do imonth=1,12
     !   print*, imonth,maxval(xave(:,:,imonth))
     !enddo
     
     
     fname = '/Users/djlorenz/Research/jack/final/cmip5/historical/'//trim(models(imodel))
     call system('mkdir -p '//trim(fname))
     fname = '/Users/djlorenz/Research/jack/final/cmip5/historical/'//trim(models(imodel))// &
          '/prcp.nc'

     call ncopen(1,trim(fname),'prcp ', &
          vunits='mm',vnames='mean monthly precipitation accumulation;', &
          lons=lonsf,lats=latsf,title='Downscaled CMIP5 models, monthly mean.', &
          deflate=1,sf=[sf],ao=[ao],syear=year1,tunit='M')

     do iyear = year1,year2

        do imonth = 1,12
!           x1(:,:,1,imonth,iyear-year1+1) = merge(x1(:,:,1,imonth,iyear-year1+1) / &
!                xave(:,:,imonth), 0.0, xave(:,:,imonth) /= 0.0 .and. &
!                x1(:,:,1,imonth,iyear-year1+1) /= missing)
           x1(:,:,1,imonth,iyear-year1+1) = merge(merge(x1(:,:,1,imonth,iyear-year1+1) / &
                xave(:,:,imonth), missing, x1(:,:,1,imonth,iyear-year1+1) /= missing), &
                0.0, xave(:,:,imonth) /= 0.0)
        enddo

        call bilinear(lonsc,latsc,nlonc,nlatc,lonsf,latsf,nlonf,nlatf, &
             x1(:,:,1,:,iyear-year1+1),xf,12,0)

        !print*, iyear,maxval(xf),minval(xf,xf/=missing)

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
             trim(models(imodel))//'/pr*.nc'
        vars = 'pr'
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
        
        !x1max(:,:,3,:) = maxval(x1(:,:,1,:,:),dim=4)

        ! quantile adjust based on quantiles found previously:
        do imonth=1,12
           do ilat = 1,nlatc
              do ilon = 1,nlonc
                 if(.not.any(xbox(ilon,ilat,imonth,:)==missing)) then
                    qs = quantile(:,ilon,ilat,imonth,:)
                    call quantile_adj(nyear,x1(ilon,ilat,1,imonth,:), &
                         junk,nquantile, &
                         qs,2)
                 else
                    x1(ilon,ilat,1,imonth,:) = missing
                 endif
              enddo
           enddo
        enddo

        !x1max(:,:,4,:) = maxval(x1(:,:,1,:,:),dim=4)
       
     !open(99,file='pmax.dat',access='stream',form='unformatted')
     !write(99) x1max
     !close(99)
     !call ctlwrite('pmax.ctl','p0 pa0 p1 pa1',lons=lonsc,lats=latsc, &
     !     nt=12,tunit='M')
     !stop
 
        fname = '/Users/djlorenz/Research/jack/final/cmip5/'//trim(scens(iscen))//'/'// &
             trim(models(imodel))
        call system('mkdir -p '//trim(fname))
        fname = '/Users/djlorenz/Research/jack/final/cmip5/'//trim(scens(iscen))//'/'// &
             trim(models(imodel))//'/prcp.nc'
        
        call ncopen(1,trim(fname),'prcp ', &
             vunits='mm',vnames='mean monthly precipitation accumulation;', &
             lons=lonsf,lats=latsf,title='Downscaled CMIP5 models, monthly mean.', &
             deflate=1,sf=[sf],ao=[ao],syear=year11,tunit='M')
        
        do iyear = year11,year22
           
           do imonth = 1,12
              !x1(:,:,1,imonth,iyear-year11+1) = merge(x1(:,:,1,imonth,iyear-year11+1) / &
              !     xave(:,:,imonth), 0.0, xave(:,:,imonth) /= 0.0 .and. &
              !     x1(:,:,1,imonth,iyear-year11+1) /= missing)
              x1(:,:,1,imonth,iyear-year11+1) = merge(merge(x1(:,:,1,imonth,iyear-year11+1) / &
                   xave(:,:,imonth), missing, x1(:,:,1,imonth,iyear-year11+1) /= missing), &
                   0.0, xave(:,:,imonth) /= 0.0)
           enddo
           
           call bilinear(lonsc,latsc,nlonc,nlatc,lonsf,latsf,nlonf,nlatf, &
                x1(:,:,1,:,iyear-year11+1),xf,12,0)

           !print*, iyear,maxval(xf),minval(xf,xf/=missing)
           
           xf = merge(xf * xfave, missing, xfave /= missing .and. xmiss /= missing)

           do imonth = 1,12
              do ilat = 1,nlatf
                 do ilon = 1,nlonf
                    if( xf(ilon,ilat,imonth) /= missing ) then
                       if(xf(ilon,ilat,imonth) > maxprcp) then
                          xf(ilon,ilat,imonth) = maxprcp
                          !print*, 'max'
                       endif
                       temp = (xf(ilon,ilat,imonth) - ao)/sf
                       if(abs(temp)>32767.) then
                          print*, 'ERROR: conversion to 2 byte: out of range. ',temp
                          print*, ilon,ilat,imonth
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

     deallocate(xave,lonsc,latsc,quantile,xbox)
     
     deallocate(x1max)

  enddo

end program temp0
  
! assumes no missing values:
subroutine quantile_adj(nt,x,y,nq,quantile,ical)

  implicit none

  integer :: nt,nq,ical
  real :: x(nt),y(nt),quantile(2,nq)

  ! work
  real :: xsort(nq),ysort(nq),xa,ya,slope,q
  integer :: it,j,k,klo,khi

  if(ical == 1) then
     if(nq/=nt) then
        print*, 'ERROR: quantile_adj, ',nt,nq
        stop
     endif
     xsort = x
     call shell0(nt,xsort)
     quantile(1,:) = xsort
     ysort = y
     call shell0(nt,ysort)
     quantile(2,:) = ysort
  else
     xsort = quantile(1,:)
     ysort = quantile(2,:)
  end if

  xa = sum(xsort)/nq
  xsort = xsort - xa
  ya = sum(ysort)/nq
  ysort = ysort - ya
  slope = sum(xsort*ysort)/sum(xsort**2)

  do it=1,nt
     q = x(it)
     if(q <= 0.) then
        x(it) = 0.
     elseif(q < quantile(1,1)) then
        x(it) = quantile(2,1)*x(it)/quantile(1,1)
     elseif(q > quantile(1,nq)) then
        x(it) = quantile(2,nq) + slope*(x(it) - quantile(1,nq))
     else
        j = 3
        klo=1
        khi=nq
1       if(khi-klo>1) then
           k=(khi+klo)/2  
           if(quantile(1,k)>q)then  
              khi=k  
           else  
              klo=k  
           endif
           goto 1  
        endif
        x(it) = &
             quantile(2,klo) + &
             (quantile(2,khi)-quantile(2,klo))/ &
             (quantile(1,khi)-quantile(1,klo))* &
             (q-quantile(1,klo))
     endif
  end do


end subroutine quantile_adj
