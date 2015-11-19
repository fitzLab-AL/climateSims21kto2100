!seasonal precipitation => monthly precipitation
!gfortran -O3 -o prcp_monthly prcp_monthly.f90 read_ECBilt.o hybrj.o seasonal_monthly.o -I/Users/djlorenz/Research/netcdf/include -I/Users/djlorenz/lib -L/Users/djlorenz/Research/netcdf/lib -L/Users/djlorenz/lib -lnetcdff -lnetcdf  -llapack -lblas -lmylib
program prcp

  use netcdf
  use constants
  use read_data
  use nc_write
  use read_ECBilt
  use write_data
  use seasonal_monthly
  implicit none

  integer, parameter :: year1 = 1901, year2 = 2011

  ! 3 == 1981-1990, 0 == 1951-1960
  integer, parameter :: decade1 = -2100, decade2 = -1, ndecade = decade2-decade1+1
  
  real, parameter :: minlon = -173.0, maxlon = -48.0
  real, parameter :: minlat = 10.0, maxlat = 80.0

  real, allocatable :: x(:,:,:),xbox(:,:,:)
  real, allocatable :: xsea(:,:,:,:),xmon(:,:,:,:)
  real :: xtemp(12,5),mean(12)

  integer :: idecade,ilat,ilon

  ! get grid
  call get_ECBilt(minlon,maxlon,minlat,maxlat,-1,-1,'prcp')

  ! get observed climatology:
  allocate(xbox(nlonc,nlatc,12))
  ! get prcp, Tmax Tmin, vap:
  file_expr = '/Users/djlorenz/Research/jack/data/prcp.nc '
  vars = 'prcp '
  call read_nc(minlon,maxlon,.false.,minlat,maxlat,missing,missing, &
       [all,all,all,year1],[all,all,all,year2],1,'month')
  allocate(x(nlon,nlat,12))
  x = sum(x1(:,:,1,:,:),dim=4)/nyear
  x = merge(x,missing,abs(x)<=1e30)
  call box(lons,lats,nlon,nlat,lonsc,latsc,nlonc,nlatc,x,xbox,12,0)
  deallocate(x,x1)

  allocate(xsea(nlonc,nlatc,ndecade,4),xmon(nlonc,nlatc,ndecade,12))

  open(1,file='prcpL.dat',access='stream',form='unformatted')
  read(1) xsea
  close(1)

  do idecade=1,ndecade
     print*, idecade
     do ilat=1,nlatc
        do ilon=1,nlonc
           if(xbox(ilon,ilat,1)/=missing .and. xsea(ilon,ilat,idecade,1)/=missing) then
              call sm_factor(xbox(ilon,ilat,:),xsea(ilon,ilat,idecade,:),xmon(ilon,ilat,idecade,:),2)
!              call sea1month(xsea(ilon,ilat,idecade,:),xmon(ilon,ilat,idecade,:))
           else
              xmon(ilon,ilat,idecade,:) = missing
           endif
        enddo
     enddo
  enddo

  print*, any(isnan(xmon))

  call ctlwrite('prcpL_monthly.ctl','prcp',lons=lonsc,lats=latsc, &
       nlev=ndecade,minlev=0.01*decade1,dlev=0.01,nens=12)
  open(1,file='prcpL_monthly.dat',access='stream',form='unformatted')
  write(1) xmon
  close(1)

end program prcp

! minimize first derivative:
subroutine sea1month(xs,xm)

  implicit none

  real :: xs(4),xm(12)

  real :: a(16,16),x(16)

  integer :: j,k
  integer :: ipiv(16),info

  a = 0.0; x = 0.0
  do j=1,12
     a(j,j) = 4.0
     k = modulo(j,12) + 1
     a(j,k) = -2.0
     k = modulo(j-2,12) + 1
     a(j,k) = -2.0
  enddo
  a(1:2,13) = 1.0
  a(12,13) = 1.0
  a(3:5,14) = 1.0
  a(6:8,15) = 1.0
  a(9:11,16) = 1.0
  a(13,1:2) = 1.0
  a(13,12) = 1.0
  a(14,3:5) = 1.0
  a(15,6:8) = 1.0
  a(16,9:11) = 1.0

  ! 3.0 because of average
  x(13:16) = 3.0*xs

  call sgesv(16,1,a,16,ipiv,x,16,info)

  xm = x(1:12)

end subroutine sea1month

subroutine sea1month_trans(xs,xm)

  implicit none

  real :: xs(4),xm(12)

  real :: a(16,16),x(16)

  integer :: j,k
  integer :: ipiv(16),info

  a = 0.0; x = 0.0
  do j=1,12
     a(j,j) = 4.0
     k = modulo(j,12) + 1
     a(j,k) = -2.0
     k = modulo(j-2,12) + 1
     a(j,k) = -2.0
  enddo
  a(1:2,13) = 1.0
  a(12,13) = 1.0
  a(3:5,14) = 1.0
  a(6:8,15) = 1.0
  a(9:11,16) = 1.0
  a(13,1:2) = 1.0
  a(13,12) = 1.0
  a(14,3:5) = 1.0
  a(15,6:8) = 1.0
  a(16,9:11) = 1.0

  ! 3.0 because of average
  x(13:16) = 3.0*log(xs)

  call sgesv(16,1,a,16,ipiv,x,16,info)

  xm = exp(x(1:12))

end subroutine sea1month_trans

