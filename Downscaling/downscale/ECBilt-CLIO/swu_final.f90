!course monthly anomalies to high resolution full variables (i.e. not anomalies) 
!gfortran -O3 -o swu_final swu_final.f90 read_ECBilt.o write_final.o -I/Users/djlorenz/Research/netcdf/include -I/Users/djlorenz/lib -L/Users/djlorenz/Research/netcdf/lib -L/Users/djlorenz/lib -lnetcdff -lnetcdf  -llapack -lblas -lmylib
program swu0

  use netcdf
  use constants
  use read_data
  use read_ECBilt
  use write_final
  use nc_write, only : ctlwrite
  implicit none

  integer, parameter :: year1 = 1901, year2 = 2011

  ! 3 == 1981-1990, 0 == 1951-1960
  integer, parameter :: decade1 = -2100, decade2 = -1, ndecade = decade2-decade1+1
  
  real, parameter :: minlon = -173.0, maxlon = -48.0
  real, parameter :: minlat = 10.0, maxlat = 80.0

  real, allocatable :: x(:,:,:),xrel(:,:,:),xmiss(:,:,:),dx(:,:,:)
  real, allocatable :: xmon(:,:,:,:)
  integer*2, allocatable :: ix(:,:,:)
  real, parameter :: sf=0.1, ao=0.0
  real :: temp,xrel2(12)

  integer, parameter :: nlatT=360
  integer :: lat0
  real :: latsT(nlatT),insol(nlatT,ndecade,12)

  integer :: idecade,ilat,ilon,imonth

  ! get insolation: 
  do ilat=1,nlatT
     latsT(ilat) = 180./nlatT*(ilat - nlatT/2 - 0.5) 
  enddo
  do ilat=1,nlatT
     if(latsT(ilat) >= minlat) exit
     lat0 = ilat
  enddo

  open(1,file='../../berger1978/insol_high.dat',access='stream',form='unformatted')
  read(1) insol
  close(1)

  ! get grid
  call get_ECBilt(minlon,maxlon,minlat,maxlat,-1,-1,'prcp')

  ! get CRU wind for landsea mask:
  file_expr = '/Users/djlorenz/Research/jack/data/wind.nc '
  vars = 'wind '
  call read_nc(minlon,maxlon,.false.,minlat,maxlat,missing,missing, &
       [all,all,all,2001],[all,all,all,2001],1,'month')
  allocate(xmiss(nlon,nlat,12))
  xmiss = x1(:,:,1,:,1)
  deallocate(x1)

  ! get observed hi-res climatology:
  file_expr = '/Users/djlorenz/Research/jack/data/srbmean.nc '
  vars = 'swd swn'
  call read_nc(minlon,maxlon,.false.,minlat,maxlat,missing,missing, &
       [all,all,all,2001],[all,all,all,2001],1,'month')
  allocate(x(nlon,nlat,12),xrel(nlon,nlat,12))
  x = x1(:,:,1,:,1) - x1(:,:,2,:,1)
  deallocate(x1)

  do imonth=1,12
     do ilat=1,nlat
        do ilon=1,nlon
           if(insol(lat0+ilat,ndecade,imonth)/=0.) then
              xrel(ilon,ilat,imonth) = x(ilon,ilat,imonth)/ &
                   insol(lat0+ilat,ndecade,imonth)
           else
!              xrel(ilon,ilat,imonth) = missing
              xrel(ilon,ilat,imonth) = 0.0
           endif
        enddo
     enddo
  enddo
!  do ilat=1,nlat
!     do ilon=1,nlon
!        xrel2 = xrel(ilon,ilat,:)
!        call fillin(12,xrel2,xrel(ilon,ilat,:))
!     enddo
!  enddo

!  open(1,file='swd_rel.dat',access='stream',form='unformatted')
!  write(1) xrel
!  close(1)
!  call ctlwrite('swd_rel.ctl','r',lons=lons,lats=lats,nt=12,tunit='M')

  allocate(xmon(nlonc,nlatc,ndecade,12),dx(nlon,nlat,12),ix(nlon,nlat,12))

  open(1,file='swu_monthly.dat',access='stream',form='unformatted')
  read(1) xmon
  close(1)

  call ncopen(1,'/Users/djlorenz/Research/jack/final/ECBilt/swu.nc','swu ', &
       vunits='Wm-2',vnames='mean monthly upward shortwave radiation at the surface;', &
       lons=lons,lats=lats,title='Downscaled TraCE ECBilt runs. Monthly means by decade.', &
       deflate=1,sf=[sf],ao=[ao])

  do idecade=1,ndecade
     print*, idecade
     call bilinear(lonsc,latsc,nlonc,nlatc,lons,lats,nlon,nlat, &
          xmon(:,:,idecade,:),dx,12,0)
     dx = merge(xrel**dx,missing,x/=missing .and. dx/=missing .and. xmiss/=missing)
     do imonth=1,12
        do ilat=1,nlat
           do ilon=1,nlon
              if(dx(ilon,ilat,imonth)/=missing) then
                 dx(ilon,ilat,imonth) = dx(ilon,ilat,imonth)* &
                      insol(lat0+ilat,idecade,imonth)
              endif
           enddo
        enddo
     enddo
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

end program swu0

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
  
  if(any(x==missing)) then
     a = 0.0
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
     
  endif
  
  
end subroutine fillin
